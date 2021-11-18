import unittest
import demag_nonequi_layers as dnl
import numpy as np

def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True

if module_exists('arrayfire'):
    import arrayfire as af
if module_exists('magnumaf'):
    import magnumaf as maf

class TestNonequiDemag(unittest.TestCase):
    nx, ny = 4, 8
    dx, dy =  1.e-9, 2e-9
    z_spacing = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9]

    if module_exists('arrayfire') and module_exists('magnumaf'):
        maf_mesh = maf.NonequiMesh(nx, ny, dx, dy, z_spacing)
        maf_demag = maf.NonequiDemagField(maf_mesh, verbose = False, caching = False, nthreads = 6)
        N_fft2d_maf = maf_demag.get_Nfft()

        def test_N_arrayfire_fft2d_eq_to_magnumaf(self):
            N_numpy = dnl.N_demag(nx = self.nx, ny = self.ny, dx = self.dx, dy = self.dy, dz_list = self.z_spacing)
            N_numpy_af_fft2d = af.fft2_r2c(af.from_ndarray(N_numpy))
            max_abs_diff = af.max(af.abs(N_numpy_af_fft2d - self.N_fft2d_maf))
            self.assertEqual(max_abs_diff, 0.0)

        def test_N_numpy_fft2d_almost_eq_to_magnumaf(self):
            Ndemag_fft2d_numpy = dnl.N_demag_fft2d(nx = self.nx, ny = self.ny, dx = self.dx, dy = self.dy, dz_list = self.z_spacing)
            max_abs_diff = np.max(np.abs(Ndemag_fft2d_numpy - self.N_fft2d_maf.to_ndarray()))
            self.assertAlmostEqual(max_abs_diff, 0.0, places = 40)

        def test_Heff_in_Apm_almost_eq_to_magnumaf(self):
            Nfft2d_numpy = dnl.N_demag_fft2d(nx = self.nx, ny = self.ny, dx = self.dx, dy = self.dy, dz_list = self.z_spacing)
            m = np.zeros((self.nx, self.ny, len(self.z_spacing), 3))
            m[:, :, :, 0] = 1 # set homogen in x dir

            # magnum.af:
            Ms = 1.0 * maf.Constants.mu0
            m0 = af.from_ndarray(m)
            state = maf.State(maf.Mesh(0,0,0,0,0,0), Ms, m0)
            Heff_af = self.maf_demag.H_in_Apm(state)

            # TODO enable when impl
            # Heff_numpy = dnl.Heff_in_Apm(Nfft2d_numpy, m)
            # max_abs_diff_Heff = np.max(np.abs(Heff_numpy - Heff_af.to_ndarray()))
            # self.assertAlmostEqual(max_abs_diff_Heff, 0.0, places = 40)

    else:
        print("Info: arrayfire or magnum.af not found, skipping depending tests.")


if __name__ == '__main__':
    unittest.main()
