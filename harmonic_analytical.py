import numpy as np
# import pandas as pd
# from matplotlib.ticker import FormatStrFormatter
import scipy.constants

def zk(bhw, k, d):
    # return 1 / (np.exp(beta * hbar * omega / 2) - np.exp(- beta * hbar * omega / 2))
    return pow((np.exp(0.5 * k * bhw)) / (np.exp(k * bhw) - 1), d)

def dzk(bhw, k, d):
    return -0.5 * k * d * zk(bhw, k, d) * \
                (1 + np.exp(- k * bhw)) / (1 - np.exp(- k * bhw))

def boson_energy(bhw, N, d):
    zs = np.zeros(N + 1)
    dzs = np.zeros(N + 1)
    zs[0] = 1
    dzs[0] = 0
    for m in range(1, N + 1):
        sig_z = 0.0
        sig_dz = 0.0
        for k in range(1, m + 1):
            sig_dz += dzk(bhw, k, d) * zs[m - k] + zk(bhw, k, d) * dzs[m - k]
            sig_z += zk(bhw, k, d) * zs[m - k]
        zs[m] = (1 / m) * sig_z
        dzs[m] = (1 / m) * sig_dz
    return (-1 / zs[N] * dzs[N])


def analytical_energy_bhw(num_bosons, bhw):
    N = num_bosons
    dim = 3

    return boson_energy(bhw, N, dim)

def analytical_energy(num_bosons, temperature_kelvin):

    spring_constant = 1.21647924E-8
    mass = 1.0
    omega = np.sqrt(spring_constant / mass)

    # Boltzmann = scipy.constants.Boltzmann
    Boltzmann = 1
    # hbar = scipy.constants.hbar
    hbar = 1

    temperature = temperature_kelvin * 3.1668152e-06 # see units.py of i-pi for conversion to atomic units

    beta = (1 / (Boltzmann * temperature))

    bhw = beta * hbar * omega

    return analytical_energy_bhw(num_bosons, bhw) * hbar * omega
