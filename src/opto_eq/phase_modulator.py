import numpy as np

def phase_modulator(E_field, V, wavelength = 1.5e-6, crystal_cut='X'):
    # Since our simulation is Transverse Ex and Ey, our phase-modulator is automatically Z-Propagating as Ez = 0
    # We must choose between X-Cut or Y-Cut. X-Cut means Y surface, so apply phase to Ey.
    # Y-Cut means X surface, so apply phase to Ex.
    # X-Cut has r13 and Y-Cut has r22 config. LiNbO3 crystal is well-defined for all such configs. Pockels is the
    # linear electro-optic effect of materials in presence of electric-field.

    # Kerr effect is ignored. It is negligible compared to Pockels, and is highly non-linear electro-optic effect
    # of materials in presence of electric-field.
    # Phase modulator design is single symmetry, single-waveguide, therefore it is polarization insensitive,
    # and polarization will apply to either Ex or Ey, depending on Cut.

    # Sources:
    # Site: https://www.fiberoptics4sale.com/blogs/wave-optics/guided-wave-electro-optic-modulators
    # Site: https://www.fiberoptics4sale.com/blogs/wave-optics/pockels-effect
    # Paper: Titanium-Diffused Lithium Niobate Waveguide Devices by R.C. Alferness
    # Paper:  Lithium niobate thin film electro-optic modulator by Jikun Liu
    # Paper: Recent Progress of LiNbO3 Based Electro-optic Modulators with NRZ Coding in High Speed Photonic Networks
    # Paper:  Broadband integrated optical modulators:achievements and prospects by V M Petrov

    # Parameters:
    n_o = 2.2               # Ordinary refractive index of LiNbO3
    n_e = 2.14              # Extraordinary refractive index of LiNbO3, only used for non Z-propagating.
    d = 24.133e-6           # Separation of the Electrodes, also called Se sometimes
    L = 5.588e-2            # Length of LiNbO3 film
    Gamma = 0.8             # Physical/Spatial Overlap factor between electrode and waveguide
    r13 = 10.12e-12         # For X-Cut Config
    r22 = 3.4e-12           # For Y-Cut Config

    # Mathematical physical formulation. Kerr effect would be added here if necessary.
    if crystal_cut == "X":
        Vpi = (wavelength * d) / (2 * n_o ** 3 * r13 * Gamma * L)
        # print(f"Vπ = {Vpi:.3f} V")
        # Phase-shift formulation.
        phi = (np.pi * V) / Vpi

        # Apply via jones matrix
        pm = np.array([
            [1, 0],
            [0, np.exp(1j * phi)]
        ])
        return np.transpose(pm @ np.transpose(E_field))
        # return E

    elif crystal_cut == "Y":
        Vpi = (wavelength * d) / (2 * n_o ** 3 * r22 * Gamma * L)
        # print(f"Vπ = {Vpi} V")
        # Phase-shift formulation.
        phi = (np.pi * V) / Vpi

        # Apply via jones matrix
        pm = np.array([
            [np.exp(1j * phi), 0],
            [0, 1]
        ])

        return np.transpose(pm @ np.transpose(E_field))

    else:
        raise Exception("Wrong or improper crystal orientation!")
