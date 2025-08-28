import numpy as np
from scipy.io import savemat
def print_sub(M, IND, *, base=0, order='F'):
    """
    Mimics the MATLAB helper:
        function O = print_sub(M, IND)
    Prints (i,j,k)->M[i,j,k] for each linear index in IND, then a summary.
    Returns O = array of coordinates (N x 3).

    Parameters
    ----------
    M : np.ndarray, shape (A,B,C)
        3D array to index.
    IND : sequence of ints
        Linear indices into M. If base=1, these should be 1-based (MATLAB style).
        If base=0, these are 0-based (Python style).
    base : int, {0,1}, default 0
        Index base of IND.
    order : {'F','C'}, default 'F'
        Linearization order used to interpret IND.
        'F' = Fortran/MATLAB (first axis fastest), 'C' = NumPy default.

    Returns
    -------
    O : np.ndarray, shape (N,3)
        Coordinates (i,j,k) for each IND.
    """
    if M.ndim != 3:
        raise ValueError("M must be 3D (A,B,C).")
    A, B, C = M.shape
    IND = np.asarray(IND, dtype=np.int64).ravel()
    if base not in (0, 1):
        raise ValueError("base must be 0 or 1.")

    O = np.empty((IND.size, 3), dtype=int)

    for n, idx in enumerate(IND):
        L = idx - (1 if base == 1 else 0)  # convert to 0-based
        if L < 0:
            raise IndexError("Index < 1 encountered with base=1.")
        if order == 'F':
            # Fortran order: L = i + A*j + A*B*k
            i = L % A
            r = L // A
            j = r % B
            k = r // B
        elif order == 'C':
            # C order: L = ((i*B) + j)*C + k
            i = L // (B * C)
            r = L % (B * C)
            j = r // C
            k = r % C
        else:
            raise ValueError("order must be 'F' or 'C'.")

        if not (0 <= i < A and 0 <= j < B and 0 <= k < C):
            raise IndexError("Linear index out of bounds for M with given order.")

        O[n] = (i, j, k)
        print(f"{i},{j},{k} -> {M[i, j, k]}")

    # Summary like MATLAB’s display
    mu = O.mean(axis=0)
    sd = O.std(axis=0, ddof=1) if O.shape[0] > 1 else np.full(3, np.nan)
    mu_str = " ".join(f"{x:g}" for x in mu)
    sd_str = " ".join(f"{x:g}" for x in sd)
    print(f"a total of {IND.size} mean: {mu_str}, std: {sd_str}")

    return O

def pinv_matlab(A, tol=None):
    """
    MATLAB-like pseudoinverse using SVD, with MATLAB's default tolerance rule.
    Works for real/complex, single/double. Returns same dtype as A.

    tol:
      - If None: tol = max(size(A)) * eps(class(A)) * max(singular_values)
                  (same as MATLAB default)
      - If a number x: treat as *absolute* tolerance (MATLAB's pinv(A, x))

    Notes:
      - Keeps computations and outputs in A's dtype (e.g., complex64),
        to better match MATLAB single-precision results.
    """
    A = np.asarray(A)
    # Ensure SVD runs in the same precision as A to match MATLAB 'single' better
    dtype = A.dtype

    # SVD (full_matrices=False matches MATLAB's economy SVD for pinv)
    U, S, Vh = np.linalg.svd(A, full_matrices=False)
    # S is always real (singular values). Use EPS of the *real* part type.
    eps = np.finfo(A.real.dtype).eps

    if tol is None:
        # MATLAB default: max(m,n) * eps(norm(A)) ~ max(m,n) * eps * max(S)
        tol_eff = max(A.shape) * eps * (S[0] if S.size else 0.0)
    else:
        # MATLAB pinv(A, TOL) treats TOL as absolute cutoff
        tol_eff = float(tol)

    # Invert singular values above tolerance
    Sinv = np.where(S > tol_eff, 1.0 / S, 0.0).astype(A.real.dtype, copy=False)

    # X = V * diag(Sinv) * U^H (conjugate transpose)
    X = (Vh.conj().T * Sinv) @ U.conj().T

    # Cast back to original dtype (e.g., complex64) to match MATLAB single behavior
    return X.astype(dtype, copy=False)


def make_pinv_reg(tol=None):
    """
    Factory to match your 'pinv_reg' usage.

    - If tol is None: use MATLAB default tolerance inside pinv_matlab.
    - If tol is provided: you wanted MATLAB-like *relative* tol times ||A||_2.
      MATLAB's pinv(A,TOL) expects absolute TOL, so we translate:
          absolute_tol = tol * ||A||_2
    """
    if tol is None:
        return lambda A: pinv_matlab(A, None)
    else:
        return lambda A: pinv_matlab(A, tol * np.linalg.norm(A, 2))


import numpy as np

import numpy as np

def grappa_get_indices(kernel, samp, pad, R, type_, offset=0):
    """
    Python version with 0-based 'type_'.
    Returns 0-based Fortran-linear indices suitable for arr.ravel(order='F').

    src: (Nc*kx*ky*kz, N_targets)
    trg: (Nc,               N_targets)
    """
    kx, ky, kz = map(int, kernel)
    px, py, pz = map(int, pad)
    Rx, Ry, Rz = map(int, R)
    if Rx != 1:
        raise ValueError("x-direction must be fully sampled (Rx == 1).")

    samp = (samp != 0)
    Nc, dx, dy, dz = samp.shape

    # -------- Python indexing for type_ (0..Ry*Rz-2), skip acquired (0,0)
    missing_types = Ry * Rz - 1
    if not (0 <= type_ < missing_types):
        raise ValueError("Type parameter is inconsistent with R (expected 0..Ry*Rz-2).")

    # Enumerate missing offsets in row-major (y fast inside z), skipping (0,0)
    # [(dy,dz) for dz in 0..Rz-1, for dy in 0..Ry-1, if not (dy==0 and dz==0)]
    idx = 0
    dy0 = dz0 = None
    for zz in range(Rz):
        for yy in range(Ry):
            if yy == 0 and zz == 0:
                continue
            if idx == type_:
                dy0, dz0 = yy, zz
                break
            idx += 1
        if dy0 is not None:
            break
    # dy0, dz0 are in 0..Ry-1 / 0..Rz-1 and never (0,0)

    # -------- limits for targets (respect padding)
    x_rng = np.arange(px, dx - px, dtype=int)
    y_rng = np.arange(py, dy - py, dtype=int)
    z_rng = np.arange(pz, dz - pz, dtype=int)

    # -------- kernel source mask on a single-coil (dx,dy,dz) grid
    src_mask = np.zeros((dx, dy, dz), dtype=bool)
    xs = np.arange(0, Rx * kx, Rx, dtype=int)
    ys = np.arange(0, Ry * ky, Ry, dtype=int)
    zs = np.arange(0, Rz * kz, Rz, dtype=int)
    xs = xs[xs < dx]; ys = ys[ys < dy]; zs = zs[zs < dz]
    src_mask[np.ix_(xs, ys, zs)] = True
    k_idx = np.flatnonzero(src_mask.flatten(order="F")).astype(np.int64)  # len = kx*ky*kz

    # -------- target index inside the kernel block for this type (0-based)
    x_t = Rx * ((kx + 1) // 2) - 1
    y_t = Ry * ((ky + 1) // 2 - 1) + dy0
    z_t = Rz * ((kz + 1) // 2 - 1) + dz0
    trg_single = np.ravel_multi_index((x_t, y_t, z_t), (dx, dy, dz), order="F").astype(np.int64)

    # relative offsets for all kernel points
    k_rel = (k_idx - trg_single).astype(np.int64)  # (kx*ky*kz,)

    # -------- all valid target positions (shift samp by dy0,dz0 inside padded box)
    trg_mask = np.zeros((dx, dy, dz), dtype=bool)
    sub = samp[0][np.ix_(x_rng, y_rng, z_rng)]
    sub = np.roll(sub, shift=dy0, axis=1)  # y
    sub = np.roll(sub, shift=dz0, axis=2)  # z
    trg_3d = np.flatnonzero(sub.flatten(order="F")).astype(np.int64)

    # -------- map trg_3d back to full grid coordinates in Fortran-linear space
    # the sub-block starts at (px,py,pz); convert local linear -> global linear
    # Build full-grid linear indices of those target sites:
    # compute their multi-indices first (in the sub-block), then shift by pad and re-linearize
    Ax, Ay, Az = len(x_rng), len(y_rng), len(z_rng)
    # unravel in the sub-block (dx->Ax etc.), Fortran order
    i = trg_3d % Ax
    r = trg_3d // Ax
    j = r % Ay
    k = r // Ay
    xg = x_rng[i]
    yg = y_rng[j]
    zg = z_rng[k]
    trg_full = np.ravel_multi_index((xg, yg, zg), (dx, dy, dz), order="F").astype(np.int64)

    # -------- per-target sources on single-coil grid
    src_one = (k_rel[:, None] + trg_full[None, :]).astype(np.int64)  # (kpts, N_targets)

    # -------- replicate over coils to (Nc,dx,dy,dz) Fortran-linear
    coil_idx = np.arange(Nc, dtype=np.int64)[:, None]  # (Nc,1)
    trg = (trg_full[None, :] * Nc + coil_idx).astype(np.int64)       # (Nc, N_targets)
    src = (src_one[None, :, :] * Nc + coil_idx[:, None, :]).astype(np.int64)  # (Nc,kpts,N_targets)
    src = src.reshape(Nc * k_idx.size, trg.shape[1], order="C")

    if offset:
        trg += offset
        src += offset
    return src, trg
import numpy as np

def grappa_get_indices_matlab_order(kernel, samp, pad, R, type_, offset=0):
    """
    Returns src (Nc*kpts, N_targets) and trg (Nc, N_targets) with
    EXACT MATLAB row ordering:
      - kernel points: x fastest, then y, then z (from find(...))
      - rows: coil varies fastest inside each kernel tap
    All indices 0-based, to be used with arr.ravel(order='F').
    """
    kx, ky, kz = map(int, kernel)
    px, py, pz = map(int, pad)
    Rx, Ry, Rz = map(int, R)
    if Rx != 1:
        raise ValueError("x-direction must be fully sampled (Rx == 1).")

    samp = (samp != 0)
    Nc, dx, dy, dz = samp.shape

    # --- MATLAB's [yy,zz] = ind2sub([Ry,Rz], type+1), skipping (1,1)
    missing = Ry * Rz - 1
    if not (0 <= type_ < missing):
        raise ValueError("type_ must be 0..(Ry*Rz-2) (Python 0-based).")
    type1 = type_ + 1  # MATLAB 1-based index into the [Ry,Rz] grid, skipping (1,1)
    yy = (type1 % Ry) + 1 if (type1 % Ry) != 0 else Ry
    zz = (type1 - 1) // Ry + 1
    if yy == 1 and zz == 1:
        raise ValueError("Internal: acquired position picked; expected missing.")
    dy0, dz0 = yy - 1, zz - 1  # back to 0-based offsets

    # ---- limits for targets (respect padding)
    x_rng = np.arange(px, dx - px, dtype=int)
    y_rng = np.arange(py, dy - py, dtype=int)
    z_rng = np.arange(pz, dz - pz, dtype=int)

    # ---- single-coil kernel mask exactly like MATLAB
    src_mask = np.zeros((dx, dy, dz), dtype=bool)
    xs = np.arange(0, Rx * kx, Rx, dtype=int)  # 0..Rx*kx-1 step Rx
    ys = np.arange(0, Ry * ky, Ry, dtype=int)
    zs = np.arange(0, Rz * kz, Rz, dtype=int)
    src_mask[np.ix_(xs, ys, zs)] = True

    # k_idx: find(mask) in MATLAB → Fortran order linear indices (0-based here)
    k_idx = np.flatnonzero(src_mask.flatten(order="F")).astype(np.int64)  # length kpts
    kpts = k_idx.size

    # target inside kernel block (0-based)
    x_t = Rx * ((kx + 1) // 2) - 1
    y_t = Ry * ((ky + 1) // 2 - 1) + dy0
    z_t = Rz * ((kz + 1) // 2 - 1) + dz0
    trg_single = np.ravel_multi_index((x_t, y_t, z_t), (dx, dy, dz), order="F").astype(np.int64)

    # relative offsets
    k_rel = (k_idx - trg_single).astype(np.int64)  # (kpts,)

    # ---- all valid target positions (shift samp as MATLAB: circshift by [0 0 yy-1 zz-1])
    sub = samp[0][np.ix_(x_rng, y_rng, z_rng)]
    sub = np.roll(sub, shift=dy0, axis=1)
    sub = np.roll(sub, shift=dz0, axis=2)

    # Fortran-linear indices in the sub-block
    trg_sub = np.flatnonzero(sub.flatten(order="F")).astype(np.int64)  # (N_targets,)

    # Convert sub-block linear → (x,y,z) → full-grid linear (Fortran)
    Ax, Ay, Az = len(x_rng), len(y_rng), len(z_rng)
    i = trg_sub % Ax
    r = trg_sub // Ax
    j = r % Ay
    k = r // Ay
    xg = x_rng[i]; yg = y_rng[j]; zg = z_rng[k]
    trg_full = np.ravel_multi_index((xg, yg, zg), (dx, dy, dz), order="F").astype(np.int64)

    # ---- src_one = (kpts, N_targets) single-coil indices (MATLAB src = bsxfun(@plus, ...))
    src_one = (k_rel[:, None] + trg_full[None, :]).astype(np.int64)    # (kpts, N_targets)

    # ---- MATLAB replication EXACTLY:
    # src = bsxfun(@plus, (src(:)'-1)*nc+1, (0:nc-1)'); src = reshape(src, [], N_targets);
    vec = src_one.flatten(order='F')                 # src(:) in MATLAB
    base = (vec * Nc).astype(np.int64)               # (src(:)-1)*nc + 1  →  0-based: src(:)*Nc
    coils = np.arange(Nc, dtype=np.int64)[:, None]   # (Nc,1)
    src_mat = base[None, :] + coils                  # (Nc, kpts*N_targets)
    src = np.reshape(src_mat, (Nc * kpts, -1), order='F')  # (Nc*kpts, N_targets)

    # ---- trg replication (MATLAB: trg = bsxfun(@plus, (trg-1)*nc+1, (0:nc-1)'))
    trg_base = (trg_full * Nc).astype(np.int64)      # (N_targets,)
    trg = trg_base[None, :] + coils                  # (Nc, N_targets)

    if offset:
        src += offset
        trg += offset
    return src, trg





def circshift4(arr, shifts):
    # shifts is a 4-tuple for axes (0,1,2,3); positive = shift toward higher indices
    out = arr
    for ax, s in enumerate(shifts):
        if s != 0:
            out = np.roll(out, shift=int(s), axis=ax)
    return out

def grappa_gfactor(data, calib, noise, R, kernel, tol=None, debug=False):
    """
    rewrite of grappa_gfactor.m in pyhton

    Inputs:
      data   : (Nc, Mx, My, Mz) undersampled k-space (complex)
      calib  : (Nc, Nx, Ny, Nz) ACS k-space (complex)
      noise  : (Nc, Nc) coil noise covariance
      R      : (Rx, Ry[, Rz])
      kernel : (kx, ky[, kz])
      tol    : relative SVD cutoff (like MATLAB pinv(A, tol*norm(A,2)))
      debug  : print debugging information

    Returns:
      recon_img_coils : (Nc, Mx, My, Mz) image-domain unaliased reconstruction
      g               : (Mx, My, Mz) SOS-combined g-factor
    """
    
    Nc, Mx, My, Mz = data.shape
    _,  Nx, Ny, Nz = calib.shape

    # Promote to 3D
    if len(R) == 2:       R = (R[0], R[1], 1)
    if len(kernel) == 2:  kernel = (kernel[0], kernel[1], 1)

    if debug:
        print(f"[DBG] Nc={Nc}  M=({Mx},{My},{Mz})  N=({Nx},{Ny},{Nz})  R={R}  kernel={kernel}")

    # Pseudoinverse choice
    pinv_reg = make_pinv_reg(tol=tol)

    
    # MATLAB: W = zeros([Nc, Nx*Ny*Nz, Nc]);
    W = np.zeros((Nc, Nx * Ny * Nz, Nc), dtype=np.complex128)
    
    # MATLAB: for type = 1:prod(R(2:end))-1
    # Python: type_ = 0 .. prod(R(2:))-2
    n_types = int(np.prod(R[1:])) - 1 #or Cartesian GRAPPA with undersampling only along y/z, there are prod(R(2:3)) - 1 distinct hole positions relative to the acquired lines.
    for type_ in range(n_types):
        #Loop over type = 1 … (prod(R(2:3)) - 1). Each type corresponds to a different relative (yy, zz) offset of a missing point within an R(2) × R(3) cell.
        pad     = tuple(np.floor(np.array(R) * np.array(kernel) / 2).astype(int))  # floor(R.*kernel/2)
        # pad = (2, 5, 0)

        # Full-true sampling mask for the call (same as true([Nc Nx Ny Nz]) in MATLAB):
        samp = np.ones((Nc, Nx, Ny, Nz), dtype=bool)
        src, trg = grappa_get_indices_matlab_order(kernel=kernel, samp=samp, pad=pad, R=R, type_=type_)
        calib_vec = calib.ravel(order='F')
        # 3) Gather VALUES exactly like MATLAB linear indexing
        calib_trg = calib_vec[trg]                 # shape: (Nc, N_targets)
        calib_src = calib_vec[src]                 # shape: (Nc*kx*ky*kz, N_targets)
        A = calib_src
        B = calib_trg
        weights = B @ pinv_reg(A)     # shape: (Nc, N_targets, Nc)
        
        
        kpts = calib_src.shape[0] // Nc
        weights = np.reshape(weights, (Nc, Nc, kpts), order='F')

        # ---- (2) KERNEL-PLACEMENT INDICES (circshifted center) ----
        # [yy, zz] = ind2sub(R(2:3), type+1)  (MATLAB 1-based)
        iy, iz = np.unravel_index(type_ + 1, (R[1], R[2]), order='F')  # 0-based
        yy, zz = iy + 1, iz + 1
        mask = np.zeros((Nc, Nx, Ny, Nz), dtype=bool)
        cx = Nx//2
        cy = Ny//2
        cz = Nz//2
        mask[:, cx, cy, cz] = True
        samp_kernel = circshift4(mask, (0, 0, 1 - yy, 1 - zz))
        src_kernel, _ = grappa_get_indices_matlab_order(kernel=kernel, samp=samp_kernel, pad=pad, R=R, type_=type_)

        # ---- (3) WRITE INTO W USING DIRECT (i,j,c) MAPPING (no flatten views) ----
        
        s = np.asarray(src_kernel, dtype=np.int64)     # already 0-based per your note

        i =  s % Nc                             # coil index
        j = (s // Nc)                           # spatial (Nx*Ny*Nz) index in Fortran flattening

        # One target-coil slice at a time
        for c in range(Nc):
            # W[i, j, c] = weights[c, :, :].ravel(order='F')   # length must equal len(src)
            W[i,j,c]= np.expand_dims(weights[c, :, :].ravel(order='F'),axis=-1)
        # from scipy.io import savemat
        # savemat(f'/g/src{type_:02d}.mat', {'W': W})  # Save W for debugging
    # ---- (4) APPLY WEIGHTS TO DATA ----
    
        
    W = np.stack(
    [W[:, :, c].reshape(Nc, Nx, Ny, Nz, order='F') for c in range(Nc)],
    axis=-1
    )  # -> (Nc, Nx, Ny, Nz, Nc)
        
        
        # --- Flip, pad, permute, center DC, and inverse FFT (MATLAB port) ---

    # 1) Flip each spatial axis and circshift by +1 (dims 2..4 in MATLAB → axes 1..3 here)
    W = np.roll(np.flip(W, axis=1), shift=1, axis=1)  # flip dim-2 (Nx), then shift
    W = np.roll(np.flip(W, axis=2), shift=1, axis=2)  # flip dim-3 (Ny), then shift
    W = np.roll(np.flip(W, axis=3), shift=1, axis=3)  # flip dim-4 (Nz), then shift

    # 2) Pad to (Mx, My, Mz) centered (pre = ceil, post = floor)
    dx, dy, dz = Mx - Nx, My - Ny, Mz - Nz
    pad_width = [
        (0, 0),                                      # Nc (leading)        → no pad
        (int(np.ceil(dx/2)), int(np.floor(dx/2))),   # Nx → Mx (pre, post)
        (int(np.ceil(dy/2)), int(np.floor(dy/2))),   # Ny → My
        (int(np.ceil(dz/2)), int(np.floor(dz/2))),   # Nz → Mz
        (0, 0),                                      # Nc (trailing)       → no pad
    ]
    W = np.pad(W, pad_width, mode="constant", constant_values=0)

    # 3) Permute: [5,1,2,3,4] → bring last Nc to front → (Nc_out, Nc_in, Mx, My, Mz)
    W = np.transpose(W, (4, 0, 1, 2, 3))

    # 4) Put identity at k-space center voxel (MATLAB: floor(end/2+1) → Python: L//2)
    W[:, :, Mx // 2, My // 2, Mz // 2] = np.eye(Nc, dtype=W.dtype)

    # 5) ifft along spatial dims 3:5 with sqrt(Mx*My*Mz) scaling (MATLAB ifftdim)
    sp_axes = (2, 3, 4)

    # Helper that mirrors MATLAB-style centered IFFT along given axes
    def ifftdim(X, axes):
        return np.fft.fftshift(
            np.fft.ifftn(np.fft.ifftshift(X, axes=axes), axes=axes, norm="ortho"),
            axes=axes
        )

    Nsp = np.prod([W.shape[a] for a in sp_axes])
    W = ifftdim(W, sp_axes) * np.sqrt(Nsp) 
    
    print(f"[DBG] W after IFFT shape: {W.shape}") if debug else None

    data = ifftdim(data, axes=(1,2,3))

    data = np.expand_dims(data, axis=0)   # shape (1, Nc, Mx, My, Mz)
    
    res = np.sum(W * data, axis=1)   # sum over Nc_in
    data = np.squeeze(res)
    

    Nvox = Mx * My * Mz

    # 1) Flatten spatial dims (Fortran order to match MATLAB)
    Wv = np.reshape(W, (Nc, Nc, Nvox), order='F')                   # (Nc, Nc, Nvox)
    d  = np.reshape(np.conj(data), (Nc, 1, Nvox), order='F')        # (Nc, 1,  Nvox)

    # 2) W = squeeze(sum(d .* Wv, 1));  % sum over Nc_in (first dim in MATLAB code)
    #    Here we sum over the first axis (Nc) because d and Wv share that axis.
    Wv = np.sum(d * Wv, axis=0)                                     # (Nc, Nvox)

    # 3) W = sum(W .* conj(noise * W), 1);
    #    noise @ Wv -> (Nc, Nvox); elementwise multiply and sum over coils
    Wv = np.sum(Wv * np.conj(noise @ Wv), axis=0)                   # (Nvox,)

    # 4) p = sum(squeeze(d) .* conj(noise * squeeze(d)), 1);
    d2 = np.squeeze(d, axis=1)                                      # (Nc, Nvox)
    p  = np.sum(d2 * np.conj(noise @ d2), axis=0)                   # (Nvox,)

    # 5) g = reshape(sqrt(abs(W)./abs(p)), [Mx, My, Mz]) / prod(R);
    G = np.reshape(np.sqrt(np.abs(Wv) / np.abs(p)), (Mx, My, Mz), order='F') / np.prod(R)

    # 6) data = reshape(conj(data), Nc, Mx, My, Mz);   % restore original shape
    # data = np.reshape(np.conj(d2), (Nc, Mx, My, Mz), order='F')
    
    
    
    return G