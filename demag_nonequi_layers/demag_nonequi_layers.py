import numpy as np
import demag_nonequi_layers.demag_nonequi_layers_cpp as dnlc
from typing import List

# or reexport: from demag_nonequi_layers_cpp import setup_demagtensor as N_demag
def N_demag(nx: int, ny: int, dx: float, dy: float, dz_list: List[float]):
    return dnlc.setup_demagtensor(nx = nx, ny = ny, dx = dx, dy = dy,
        dz_list = dz_list)

def rfft2_along_xy(N_demag: np.ndarray):
    """Applies a 2D FFT of the demag tensor along the x and y axes."""
    # Note: arrayfire uses equivalent to axes = (1, 0);
    # Note: axis = (0, 1) would reduce second dim, not first dim as af
    return np.fft.rfft2(N_demag, axes = (1, 0))

def N_demag_fft2d(nx: int, ny: int, dx: float, dy: float, dz_list: List[float]):
    return rfft2_along_xy(N_demag(nx = nx, ny = ny, dx = dx, dy = dy,
        dz_list = dz_list))

def Heff_in_Apm(Nfft2d: np.ndarray, m: np.ndarray):
    """Calculate the effective demagnetization field."""
    raise RuntimeError("Not yet imlemented.")
