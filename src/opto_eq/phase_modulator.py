import numpy as np


class PhaseModulator:
    def __init__(self, crystal_cut='X', modulation="DC", params=None):
        """
        Initialize a Phase Modulator using the Pockels effect in LiNbO3.

        Parameters:
            crystal_cut (str): 'X' or 'Y'. Determines modulation axis.
            modulation (str): 'DC' or 'RF'. Only 'DC' is implemented currently.
            params (dict, optional): Dictionary of modulator parameters.
                Accepts keys: 'wavelength', 'n_o', 'n_e', 'd', 'L', 'Gamma', 'r13', 'r22'.

        Raises:
            RuntimeError: If unknown crystal_cut, modulation, or parameter key is provided.

        Important Notes:
            - RF modulation is not yet implemented.
            - Based on Z-propagation with Ex, Ey transverse field simulation.
            - Based on well-known Pockels effect in LiNbO3.
        Explanation:
            - Since our simulation is Transverse Ex and Ey, our phase-modulator is automatically Z-Propagating as Ez = 0
            - We must choose between X-Cut or Y-Cut. X-Cut means Y surface, so apply phase to Ey.
            - Y-Cut means X surface, so apply phase to Ex.
            - X-Cut has r13 and Y-Cut has r22 config. LiNbO3 crystal is well-defined for all such configs. Pockels is the
            - linear electro-optic effect of materials in presence of electric-field.

            - Kerr effect is ignored. It is negligible compared to Pockels, and is highly non-linear electro-optic effect
            - of materials in presence of electric-field.
            - Phase modulator design is single symmetry, single-waveguide, therefore it is polarization insensitive,
            - and polarization will apply to either Ex or Ey, depending on Cut.

        Sources:
            - Site: https://www.fiberoptics4sale.com/blogs/wave-optics/guided-wave-electro-optic-modulators
            - Site: https://www.fiberoptics4sale.com/blogs/wave-optics/pockels-effect
            - Paper: Titanium-Diffused Lithium Niobate Waveguide Devices by R.C. Alferness
            - Paper:  Lithium niobate thin film electro-optic modulator by Jikun Liu
            - Paper: Recent Progress of LiNbO3 Based Electro-optic Modulators with NRZ Coding in High Speed Photonic Networks
            - Paper:  Broadband integrated optical modulators:achievements and prospects by V M Petrov
        """
        self.default_params = {
            "wavelength": 1550e-9,   # Wavelength (in nm)
            "n_o": 2.2,     # Ordinary refractive index of LiNbO3
            "n_e": 2.14,     # Extraordinary refractive index of LiNbO3, only used for non Z-propagating.
            "d": 24.133e-6,       # Separation of the Electrodes, also called Se sometimes
            "L": 5.588e-2,       # Length of LiNbO3 film
            "Gamma": 0.8,   # Physical/Spatial Overlap factor between electrode and waveguide
            "r13": 10.12e-12,     # For X-Cut Config
            "r22": 3.4e-12     # For Y-Cut Config
        }
        if params is None:
            params = self.default_params

        for key in params:
            if key not in self.default_params:
                raise RuntimeError(f"Unidentified Parameter: {key}")

        if crystal_cut not in ["X", "Y"]:
            raise RuntimeError(f"Unidentified Crystal Cut: {crystal_cut}")
        else:
            self.crystal_cut = crystal_cut

        if modulation not in ["DC", "RF"]:
            raise RuntimeError(f"Unidentified Modulation Type: {modulation}")
        elif modulation == "RF":
            raise NotImplementedError("RF modulation not implemented yet.")
        else:
            self.modulation = modulation

        for key, value in params.items():
            setattr(self, key, value)

        self.__Vpi = self.get_vpi()

    def __repr__(self):
        return f"PhaseModulator(cut={self.crystal_cut}, modulation={self.modulation}, Vpi={self.Vpi:.2f} V)"

    def get_vpi(self):
        """
        Gets half-angle voltage of the modulator
        :return: Voltage Vπ
        """
        if self.crystal_cut == "X":
            bot = (2 * self.n_o ** 3 * self.r13 * self.Gamma * self.L)
            if bot == 0:
                raise ZeroDivisionError("Invalid Parameters: Denominator of Vπ calculation is Zero!")
            return (self.wavelength * self.d) / bot
        else:
            bot = (2 * self.n_o ** 3 * self.r22 * self.Gamma * self.L)
            if bot == 0:
                raise ZeroDivisionError("Invalid Parameters: Denominator of Vπ calculation is Zero!")
            return (self.wavelength * self.d) / bot

    def get_phi(self, V):
        """
        Gets the phi value of the modulator
        :param V: Modulation voltage (in Volts)
        :return: Phase shift φ (in radians)
        """
        if np.size(V) == 1:
            return (np.pi * V) / self.__Vpi

    @property
    def Vpi(self):
        return self.__Vpi

    def modulate(self, E_field, V):
        """
        Modulate the E_field depending on Modulation Type
        :param E_field:
        :param V: Modulation DC Voltage
        :return: Modulated E_Field as numpy array of (N, 2)
        """
        if E_field.shape[-1] != 2:
            raise ValueError("E_field must be a 2D array with shape (N, 2)")

        if self.modulation == "DC" and np.size(V) == 1:
            phi = (np.pi * V) / self.__Vpi
            if self.crystal_cut == "X":
                pm = np.array([
                    [1, 0],
                    [0, np.exp(1j * phi)]
                ])
                return np.transpose(pm @ np.transpose(E_field))
            else:
                pm = np.array([
                    [np.exp(1j * phi), 0],
                    [0, 1]
                ])
                return np.transpose(pm @ np.transpose(E_field))
        elif self.modulation == "RF" and np.size(V) > 1:
            # Is not implemented yet, therefore return E_field
            return E_field
