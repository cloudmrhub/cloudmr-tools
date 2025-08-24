Here’s a clear, implementation-level recipe for computing a GRAPPA g-factor map, distilled from the two MATLAB files you shared.

# Inputs (what you must have)

* **Undersampled k-space** `data`: shape `(Nc, Mx, My, Mz)` (Nc coils).
* **Calibration region** `calib`: fully sampled ACS k-space, shape `(Nc, Nx, Ny, Nz)`.
* **Noise covariance** `C` (a.k.a. `noise`): shape `(Nc, Nc)`, estimated from noise pre-scan or corners of k-space.
* **Acceleration** `R = [Rx, Ry, Rz]` with the assumption **Rx = 1** (undersampling only in phase encodes).
* **Kernel size** `kernel = [kx, ky, kz]` (odd sizes recommended).
* **Optional SVD cutoff** `tol` for regularized pseudoinverse.&#x20;

---

# Output (what you compute)

* **Reconstructed image** (coil-combined).
* **g-factor map** `g`: shape `(Mx, My, Mz)`, the voxelwise noise-amplification factor (normalized by the net acceleration).&#x20;

---

# Step-by-step algorithm

## 1) Normalize dimensions and defaults

* If `R` or `kernel` are given as 2D, append a trailing `1` to make them 3D.
* Choose the pseudoinverse:

  * If `tol` is provided: `pinv_reg(A) = pinv(A, tol * ||A||₂)`
  * Else: `pinv_reg = @pinv`.
* Assert **Rx == 1** (x must be fully sampled).&#x20;

## 2) Prepare a container for **k-space kernels**

* Allocate `W_k` to hold all per-coil GRAPPA kernels in k-space.
  In the reference code this is initialized as `W = zeros([Nc, Nx*Ny*Nz, Nc])` and later reshaped to `(Nc_out, Nx, Ny, Nz, Nc_in)`.&#x20;

## 3) Enumerate kernel “types” (the missing-sample patterns)

* For Cartesian GRAPPA with undersampling only along y/z, there are `prod(R(2:3)) - 1` distinct hole positions relative to the acquired lines.
* Loop over `type = 1 … (prod(R(2:3)) - 1)`. Each `type` corresponds to a different relative (yy, zz) offset of a missing point within an `R(2) × R(3)` cell.&#x20;

## 4) From ACS, collect training **source** and **target** samples for this type

Call the indexing helper once on a fully “true” mask of size `(Nc, Nx, Ny, Nz)` to get:

* `trg`: the linear indices of ACS target points matching this type (one target per coil).
* `src`: the linear indices of all ACS source neighborhoods (all coils × kernel support) associated with each target.
  This is exactly what `grappa_get_indices(kernel, samp, pad, R, type)` does, where:
* `pad = floor(R .* kernel / 2)` protects borders so neighborhoods stay inside the ACS.
* Internally it builds the relative kernel source grid, chooses the target offset given `type`, finds all valid targets (respecting undersampling periodicity), and tiles across coils.&#x20;

**Intuition:** `src` stacks (coil, kx, ky, kz) neighborhood samples for each target location; `trg` points to the coil-specific target value to predict.

## 5) Train the GRAPPA weights for this type

Solve a linear least-squares mapping from sources to targets using ACS data:

$$
\textbf{W}_\text{type} \;=\; \underbrace{\text{calib}(\text{trg})}_{[Nc \times N_\text{trg}]}\;
\cdot\; \text{pinv\_reg}\!\Big(\underbrace{\text{calib}(\text{src})}_{[Nc\cdot kx\cdot ky\cdot kz \times N_\text{trg}]}\Big).
$$

This yields, for each **output coil**, a set of complex weights that linearly combine the neighborhood samples from **all coils** to predict the missing point of this type.&#x20;

## 6) Place the trained weights into a global **k-space kernel image**

* Determine the hole offset `(yy, zz)` for `type`.
* Call the indexer again with a **single-point mask** (a delta at k-space center) shifted by that `(yy, zz)` to obtain where in a kernel image each neighborhood sample lands.
* Reshape `weights` to `(Nc_out, Nc_in, kx·ky·kz)` and **write** them into the appropriate positions of `W_k` for all output coils and all input coils.
  This step builds, in k-space, the **impulse response** of the GRAPPA operator for every **input→output coil pair** at the right relative offsets.&#x20;

## 7) Convert k-space kernels to **image-space kernels**

* Reshape to `(Nc_out, Nx, Ny, Nz, Nc_in)`.
* Flip and circularly shift along x, y, z to convert the relative indexing convention into a convolution kernel centered at the origin (the code applies `flip` plus `circshift` by one voxel).
* **Pad** in k-space so kernel arrays match the final image grid `(Mx, My, Mz)`.
* **Insert identity** at the k-space center: set the voxel at DC to `eye(Nc)` so acquired samples pass through unchanged.
* **IFFT** along spatial dims to get **image-space kernels** `W_im`, scaling by `sqrt(Mx·My·Mz)` (unitary FFT convention).
  After this, `W_im` has shape `(Nc_in, Nc_out, Mx, My, Mz)` (or equivalently `(Nc_out, Nc_in, …)` depending on the final permute). Each voxel now contains an **Nc\_out×Nc\_in matrix** that linearly maps input coil images to an output (combined) image.&#x20;

## 8) Form coil images and apply the image-space kernels

* Compute image-space coil images from undersampled data: `img = IFFT(data)` along the spatial axes.
* Apply the kernels voxelwise and **sum over input coils**:

$$
x(\mathbf{r}) \;=\; \sum_{c=1}^{Nc} \big[W_{\text{im}}(\mathbf{r})\big]_{:,c}\; \cdot\; img_c(\mathbf{r}).
$$

In the reference, this is implemented as an elementwise multiply `W_im .* img` followed by a sum over the coil dimension.&#x20;

## 9) Recover the **effective linear combination vector** per voxel

For g-factor you need, for each voxel, the **linear map** that turns the vector of noisy coil data into the final scalar image value. The code reduces the Nc\_out×Nc\_in kernel and the coil image content to an **Nc×1 vector** `a( r )`:

* Reshape `W_im` to `(Nc_out, Nc_in, Nvox)`.
* Let `y(r)` be the coil image vector at voxel `r`.
* Compute $a(r)$ as the **effective combination weights** by collapsing output channels with the reconstructed image (the code uses `conj(data)` as a per-voxel factor, then sums over output coils against `W_im` to obtain an Nc×1 vector). The result captures the net linear operator $x(r)=a(r)^{H} y(r)$.&#x20;

## 10) Compute **noise variance** of the reconstructed voxel

Given the noise covariance $C$ between coils, the reconstructed noise variance at voxel $r$ is:

$$
\sigma_\text{acc}^2(r) \;=\; a(r)^{H}\, C \, a(r).
$$

This is exactly what the code computes as `sum(W .* (conj(noise * W)), 1)` after building $a(r)$ into `W` (naming overlap).&#x20;

## 11) Compute the **reference (unaccelerated) variance**

You also need the variance you would have with **fully sampled data** and an optimal coil combination at the same voxel. The code computes

$$
\sigma_\text{ref}^2(r) \;=\; s(r)^{H} \, C \, s(r),
$$

where $s(r)$ is the (complex) coil image vector (serving as a practical proxy for sensitivity-weighted optimal combination in this implementation). This appears in the code as `p = sum(squeeze(data) .* conj(noise * squeeze(data)), 1)`.&#x20;

## 12) Assemble the **g-factor map**

Finally,

$$
g(r) \;=\; \frac{\sqrt{\sigma_\text{acc}^2(r)\,/\,\sigma_\text{ref}^2(r)}}{\prod R}.
$$

Reshape `g` to `(Mx, My, Mz)`. By definition, $g \ge 1$, and values >>1 indicate strong noise amplification due to the reconstruction geometry and kernel.&#x20;

---

# Practical notes and checks

* **Indexing logic matters.** The helper `grappa_get_indices` assumes undersampling only in y/z and tiles indices across coils; if Rx≠1 it throws an error. Use `pad = floor(R .* kernel / 2)` to keep neighborhoods valid inside ACS.&#x20;
* **Kernel centering.** The flip+circshift and the identity insertion at DC align the discrete convolution properly and preserve acquired samples.&#x20;
* **Calibration quality.** Ensure ACS is large enough relative to `kernel` and `R`; if not, regularization (`tol`) becomes critical.&#x20;
* **Noise covariance.** Use a well-conditioned $C$. If you only have per-coil noise SDs, start with a diagonal $C$; if you have prewhitened data, set $C \approx I$. The g-map depends directly on $C$.&#x20;
* **2D vs 3D.** The same flow extends to 3D with `Rz>1` and `kz>1`; otherwise set `Rz=kz=1`.&#x20;

If you’d like, I can translate these steps into a compact, commented Python or MATLAB skeleton you can drop into your pipeline.
