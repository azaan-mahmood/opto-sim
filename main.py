from src.lasers import sslaser as laser
from src.viewers.stokes import compute_stokes_parameters, poincare
from src.viewers import polarimeter
from src.opto_eq import polarizer
import numpy as np

source = laser.SolidStateLaser(
    wavelength=1550e-9,  # Laser wavelength
    polarization_azimuth=np.pi,  # 45° polarization
    polarization_ellipticity=np.pi/4,
    power_dbm=-5,  # arbitrary power unit
    frequency=5e6
)

E = source.get_electric_field(normalize=False, over_period=True)
S0, S1, S2, S3 = compute_stokes_parameters(E)
psi = 0.5 * np.arctan2(S2, S1)
chi = 0.5 * np.arcsin(S3 / S0)
print(f"------------------Stokes of Source---------------------")
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
print(f"Psi (polarization azimuth) = {np.rad2deg(psi):.3f}°")
print(f"Chi (polarization ellipticity) = {np.rad2deg(chi):.3f}°")
psi = 0.5 * np.arctan2(S2, S1)
chi = 0.5 * np.arcsin(S3 / S0)
print(f"------------------Stokes after 45---------------------")
E = polarizer(E, polarization="45")
S0, S1, S2, S3 = compute_stokes_parameters(E)
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
print(f"Psi (polarization azimuth) = {np.rad2deg(psi):.3f}°")
print(f"Chi (polarization ellipticity) = {np.rad2deg(chi):.3f}°")
polarimeter(E)
poincare(S1, S2, S3)