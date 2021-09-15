import demag_nonequidistant as dn

if False:
    nx, ny = 5, 7
    dx, dy =  1.5e-9, 2.2e-9
    z_spacing = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9, 1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9]
else:
    nx, ny = 2, 3
    dx, dy =  1.e-9, 2e-9
    z_spacing = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9]

N_demag = dn.setup_demagtensor(nx = nx, ny = ny, dx = dx, dy = dy,
        dz_list = z_spacing)

# print(N_demag)
print(N_demag.shape)

def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True

if module_exists('arrayfire'):
    import arrayfire as af
    af_N_demag = af.from_ndarray(N_demag)
    # print(1e30 * af_N_demag)
    N_fft = af.fft2_r2c(af_N_demag)
    # print(1e30 * N_fft)
    print("N_fft.shape  ", N_fft.shape)

# ==== magnumaf ====

if module_exists('magnumaf'):
    import magnumaf as maf
    mesh = maf.NonequiMesh(nx, ny, dx, dy, z_spacing)
    demag = maf.NonequiDemagField(mesh, verbose = False, caching = False, nthreads = 6)
    mafNfft = demag.get_Nfft()
    # print(1e30 * mafNfft)
    print("mafNfft.shape", mafNfft.shape)
    maxdiff = af.max(N_fft - mafNfft)
    print("maxdiff=", maxdiff)

    if maxdiff != 0.0:
        raise ValueError("N_fft values differ!, maxdiff is", maxdiff, "but should be 0.0")
