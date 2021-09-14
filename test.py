import demag_nonequidistant as dn

# N_demag = dn.setup_demagtensor(nx = 2, ny = 3, dx = 1e-9, dy = 2e-9,
#         dz_list = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9])

nx, ny = 5, 7
dx, dy =  1.5e-9, 2.2e-9
z_spacing = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9, 1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9]

# nx, ny = 2, 3
# dx, dy =  1.e-9, 2e-9
# z_spacing = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9]

N_demag = dn.setup_demagtensor(nx = nx, ny = ny, dx = dx, dy = dy,
        dz_list = z_spacing)

# print(N_demag)
print(N_demag.shape)

import arrayfire as af
af_N_demag = af.from_ndarray(N_demag)
# print(1e30 * af_N_demag)
N_fft = af.fft2_r2c(af_N_demag)
# print(1e30 * N_fft)
print(N_fft.shape)

# ==== magnumaf ====

def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True

if module_exists('magnumaf'):
    import magnumaf as maf
    mesh = maf.NonequiMesh(nx, ny, dx, dy, z_spacing)
    demag = maf.NonequiDemagField(mesh, verbose = True, caching = False, nthreads = 6)
    mafNfft = demag.get_Nfft()
    # print(1e30 * mafNfft)
    print(mafNfft.shape)
    
    print(af.max(N_fft - mafNfft))
