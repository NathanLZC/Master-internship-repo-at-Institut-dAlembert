import tamaas as tm
import tamaas.utils
import matplotlib.pyplot as plt
import numpy as np
import time

# parallel computation
tm.initialize(16)

# Initialize threads and fftw
tm.set_log_level(tm.LogLevel.info)  # Show progression of solver

# discretization
n = 256

# Surface generator
sg = tm.SurfaceGeneratorFilter2D([n, n])
sg.random_seed = 1

# Spectrum
sg.spectrum = tm.Isopowerlaw2D()

# Parameters
sg.spectrum.q0 = 4
sg.spectrum.q1 = 4
sg.spectrum.q2 = 128
sg.spectrum.hurst = 0.8


# Generating surface
surface = sg.buildSurface()
surface /= tm.Statistics2D.computeSpectralRMSSlope(surface)

np.save("more_rough_surface.npy", surface)

plt.imshow(surface)
