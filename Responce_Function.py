import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sp
import scipy.fft as fft
from tqdm import tqdm


def refractive_index(w, w_TO, w_LO, gamm, epsilInf):
    n_w = np.sqrt(1 + ((w_LO ** 2 - w_TO ** 2) / (w_TO ** 2 - (sp.hbar * w) ** 2 - 1j * sp.hbar * w * gamm)) * epsilInf)
    return n_w


def complex_function(n_w, c, thickness, w, n_g):
    g_w = (2 / (n_w + 1)) * (c * (np.exp(-1j * 2 * sp.pi * w * thickness * (n_g - n_w) / c) - 1) / (
            -1j * 2 * sp.pi * w * thickness * (n_g - n_w)))
    return g_w


def electrooptic_coeff(r_e, w, w_TO, gamm, Faust_Henry_coeff):
    r41_w = r_e * (1 + Faust_Henry_coeff * (1 / (1 - (((sp.hbar * w) ** 2 - 1j * sp.hbar * w * gamm) / (w_TO ** 2)))))
    return r41_w


# __ Documentation ___
# w - frequency (Question with dimension)
# n_w - variable for refractive index
# thickness - thickness of the crystal (d in formula). Question with dimension
# sp.hbar - Plank constant from scipy lib
# sp.pi - pi constant from scipy lib


# __ Constants ___
w_TO = 177  # transverse optical phonon frequency
w_LO = 206  # longitudinal optical phonon frequency
gamm = 3.01  # lattice damping
epsilInf = 6.7  # high-frequency dielectric constant
Faust_Henry_coeff = 0.02  # represents the ratio between the ionic and the electronic part of the electro-optic effect
n_g = 2.8528  # group refractive index for 800 nm
r_e = 1  # ??? electronic nonlinearity is assumed to be constant at the mid- and far-infrared frequencies
c = 1  # ??? speed of light ???  Question with dimension
max_thz = 200 # Change here interested range of THz from 0

plt.rcParams.update({'figure.max_open_warning': 0}) # Error Catching (opened plots)

# Building the function
print("Loading...")
for i in range(1, max_thz+1):
    thickness = i
    w = np.linspace(1, max_thz, max_thz)
    n_w = refractive_index(w, w_TO, w_LO, gamm, epsilInf)
    ComplexResponse_function = complex_function(n_w, c, thickness, w, n_g) * electrooptic_coeff(r_e, w, w_TO, gamm,
                                                                                                Faust_Henry_coeff)
#  Plotting the function

    fig = plt.figure()
    ax = fig.add_subplot()
    fig.suptitle('Response function R(w) for different thicknesses of ZnTe crystal\n \n' + str(i) + ' micrometers',
                 fontname="Times New Roman", fontweight="bold")

    plt.ylim(-0.025, 0.015) # Set the limit for y axis

    plt.xlabel('THz', fontname="Times New Roman", fontweight="bold")
    plt.ylabel('R(w)', fontname="Times New Roman", fontweight="bold")

    plt.yticks(fontsize=10, fontname="Times New Roman")
    plt.xticks(fontsize=10, fontname="Times New Roman")

    # plt.plot(w, np.real(fft.fft(ComplexResponse_function)))  # Real space
    plt.plot(w, np.imag(fft.fft(ComplexResponse_function)))  # Imaginary space

    name = str(i) + '.png'
    plt.savefig(name)

# Making a gif file

images = []
for i in range(1, max_thz+1):
    images.append(imageio.imread(str(i) + '.png'))
imageio.mimsave('movie.gif', images)

print("Finished. Check your code filepath.")