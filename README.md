Demagnetization tensor of non-equidistant magnetic layers
===

A small standalone project calculating the demagnetization tensor from [1] in multi-threaded C++ and exposing it as a numpy array.

### Installation:
Either use the [Dockerfile](Dockerfile) or for local installation:

`$ sudo pip install .`

### Usage:
```python
import demag_nonequi_layers as dnl
N = dnl.N_demag(nx = 128, ny = 64, dx = 1.0e-9, dy = 2.5e-9, dz_list = [1.0e-9, 2.5e-9, 3.7e-9, 4.8e-9])
```

### References:
This method is also available in [magnum.af](https://github.com/magnum-af/magnum.af), with an additional GPU-accelerated effective field method.

[1] P. Heistracher, F. Bruckner, C. Abert, C. Vogler, and D. Suess, “Hybrid FFT algorithm for fast demagnetization field calculations on non-equidistant magnetic layers,” *Journal of Magnetism and Magnetic Materials*, vol. 503, p. 166592, June 2020.
