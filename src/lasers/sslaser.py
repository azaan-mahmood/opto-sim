import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=14199
# E is a column matrix, use column stack E = np.column_stack((Ex, Ey))
# Modeling and simulation of pulsed Er:YAG laser system
# Spectroscopic properties of Er:Yb by Lincong Rao


class SolidStateLaser:
    def __init__(self, wavelength, polarization_azimuth, polarization_ellipticity=None, frequency=None, power_dbm=1.0):
        self.tau1 = 1.071e-3 # Lifetime of lower energy level (Handbook Series on Semiconductor Parameters), also called radiative life time
        self.tau2 = 6.381e-3  # Lifetime of upper energy level
        self.tau_c = 14.4e-6 # Photon Decay or also called as Fluorescence Lifetime or Photon Lifetime
        self.N0 = 8.59e20 # Total number of atoms of Er+3:Yb+3
        self.sigma12 = 0.98e-20  # Absorption cross-section
        self.sigma21 = 1.01e-20  # Emission cross-section
        self.g2byg1 = 1 # Degeneracy Rate.
        self.Rp = self.g2byg1*(1/self.tau2)  # Pump rate (Minimum required for population inversion)
        self.alpha = 10  # Gain/Loss coefficient to damp I growth
        self.Gamma = 1.0  # Confinement factor, confinement of photons in the laser
        self.c = 3e8  # Speed of light in vacuum (m/s)
        self.h = 6.626e-34  # Planck's constant (J·s)

        self.wavelength = wavelength
        self.polarization_azimuth = polarization_azimuth
        if frequency is None:
            self.frequency = self.c / self.wavelength
        else:
            self.frequency = frequency
        if polarization_ellipticity is None:
            self.polarization_ellipticity = 0
        else:
            self.polarization_ellipticity = polarization_ellipticity
        self.power_dbm = power_dbm
        self.power_mw = 10**(power_dbm/10)
        self.I_0 = self.power_mw / (self.h * self.frequency)
        self.power_out, self.final_photons = self.out_pow()

    def out_pow(self):
        N1_0 = 2e23 # Initial population of lower energy level
        N2_0 = 0  # Initial population of upper energy level
        y0 = [N1_0, N2_0, self.I_0]
        t_span = [0, 1e-6]
        t_eval = np.linspace(0, 1e-6, 1000)
        # sol = solve_ivp(self.rate, t_span, y0, t_eval=t_eval, method='RK45')
        sol = solve_ivp(
            self.rate,
            t_span,
            y0,
            t_eval=t_eval,
            method='BDF',
            rtol=1e-6,
            atol=1e-9
        )
        # Integrate photon density over time
        integrated_I = abs(np.trapz(sol.y[2], sol.t))
        average_I = integrated_I / (t_span[1] - t_span[0])
        # Average power output
        power_out = average_I * self.h * self.frequency
        return power_out, abs(sol.y[2][-1])

    def rate(self, t, y):
        N1, N2, I = y
        dN1_dt = (N2 - self.g2byg1 * N1) * self.c * I * self.sigma12 + N2/self.tau2 - N1/self.tau1
        dN2_dt = -N2*self.sigma12*I*self.c - N2/self.tau2 + self.Rp*(self.N0 - N2)
        dI_dt = self.c*I*self.sigma12*(N2-self.g2byg1*N1) - self.alpha*I/self.tau_c
        return [dN1_dt, dN2_dt, dI_dt]

    def get_electric_field(self, t=0, over_period=False, normalize = True):
        """
        Get the electric field vector at time t or over one period of the oscillation.

        Parameters:
        t (float): Time in seconds. Ignored if over_period  is True.
        over_period (bool): If True, return the electric field vectors over one period.
        normalize (bool): If True, returns electric field as a norm, otherwise as product of
                          E = E*sqrt(power)

        Returns:
        numpy array: If over_period is False, return the electric field vector [Ex, Ey, Ez] at time t.
                     If over_period is True, return a tuple (t, Ex, Ey) where t is an array of time points,
                     Ex and Ey are arrays of the electric field components over one period.
                     Returns a normalized E Matrix which only contains direction and not amplitude
        """
        E0 = np.sqrt(2 * self.power_out)  # amplitude of the electric field
        omega = 2 * np.pi * self.frequency
        phi = self.polarization_azimuth
        chi = self.polarization_ellipticity

        if over_period:
            t = np.linspace(0, 2 * np.pi / self.frequency, 1000)
            E = np.array([self._calculate_electric_field(ti, E0, omega, phi, chi) for ti in t])
            if normalize:
                E_normalized = E / np.linalg.norm(E)
                return E_normalized
            elif not normalize:
                return E
        else:
            return self._calculate_electric_field(t, E0, omega, phi, chi)

    def _calculate_electric_field(self, t, E0, omega, phi, chi):
        """
        Calculate the electric field vector at time t.

        Parameters:
        t (float): Time in seconds.
        E0 (float): Amplitude of the electric field.
        omega (float): Angular frequency of the electric field.
        phi (float): Polarization azimuth angle.
        chi (float): Polarization ellipticity angle.

        Returns:
        numpy array: Electric field vector [Ex, Ey].
        """
        # From Gerd Keiser Book (Optical Fiber Communications)
        # Electric Field Ex (Generally called slow axis or ordinary ray).
        # Electric Field Ey (Generally called fast axis or extraordinary ray).
        # Assume Ez = 0 because of its orthogonal propagation
        # Not using exp(), was causing logical errors!
        # Might figure out later #ADD TO_DO!
        Ex = E0 * (np.cos(omega * t + phi) + 1j*np.sin(omega*t + phi))
        Ey = E0 * (np.cos(omega * t + chi) + 1j*np.sin(omega*t + chi))
        return np.array([Ex, Ey])

    def plot_photon_density(self, t, I):
        """Plot the photon density over time."""
        plt.figure(figsize=(8, 6))
        plt.plot(t, I, label='Photon Density (I)')
        plt.xlabel('Time (s)')
        plt.ylabel('Photon Density')
        plt.title('Photon Density vs Time')
        plt.grid(True)
        plt.legend()
        plt.show()

    def __str__(self):
        return (f"Polarized Light Source: λ={self.wavelength:.2e} m,"
                f" φ={self.polarization_azimuth:.4f} rad, f={self.frequency:.2e} Hz, "
                f"Pdbm={10 * np.log10(self.power_out):.2f} dBm, "
                f"Pout={self.power_out:.3f} mW, "
                f"Photon Density={self.I_0}")