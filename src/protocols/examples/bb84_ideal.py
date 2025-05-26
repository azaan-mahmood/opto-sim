import numpy as np
from src.lasers.sslaser import SolidStateLaser
from src.opto_eq import cable, optics, phase_modulator
from src.viewers import fields, polarimeter
from src.detectors import apd
import random


Vpi = 3.757  # Phase modulator Vpi as defined


def simulate_bb84(num_bits, show_pol = False):
    """
    Simulation for an ideal BB84, that does not contain Eve.
    # Example usage
    # num_bits = 1000
    # params = simulate_bb84(num_bits, show_pol = False)
    # print(f"Quantum Bit Error Rate (QBER): {params[7] * 100:.2f}%")
    :param num_bits: Number of bits for the simulation
    :param show_pol: Show polarization for each bit.
    :return: list containing [alice_bits, alice_bases, bob_bits, bob_bases, sifted_indices, sifted_alice, sifted_bob, qber]


    """
    # Initialize lists to store Alice's and Bob's data
    alice_bits, alice_bases = [], []
    bob_bits, bob_bases = [], []

    detector = apd(
        wavelength=1550e-9, quantum_efficiency=0.9, gain=10, excess_noise_factor=10,
        load_resistance=50, temperature=25, frequency=40
    )
    alice_laser = SolidStateLaser(
        wavelength=1550e-9,
        polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2,
        power_dbm=-5,
        frequency=5e6
    )

    eve_laser = SolidStateLaser(
        wavelength=1550e-9,
        polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2,
        power_dbm=-5,
        frequency=5e6
    )

    for _ in range(num_bits):
        # Alice's random choices
        alice_basis = random.choice(['C', 'X'])
        alice_bit = random.randint(0, 1)
        # phase_alice = alice_phase_shift(alice_basis, alice_bit)
        if alice_basis == 'C':
            if alice_bit == 0:
                phase_alice = Vpi / 2
            else:
                phase_alice = 3 * Vpi / 2
        else:  # 'X' basis
            if alice_bit == 0:
                phase_alice = 0
            else:
                phase_alice = Vpi

        # Generate photon with Alice's encoding
        E = alice_laser.get_electric_field(normalize=False, over_period=True)
        pout = alice_laser.power_out

        # Apply Alice Phase Modulation
        E = optics.polarizer(E, '45')  # Initial polarization
        E = phase_modulator(E, phase_alice, 'X')  # Alice's phase shift

        if show_pol:
            polarimeter(E, title=f"Bit Number {num_bits},Alice Bit/Basis = {alice_bit} / {alice_basis}")

        # Channel transmission (QC). Dispersion = False = Ideal QC
        E, _ = cable(
            fiber_length=100, E=E, dispersion=False, pin=pout,
            attenuation_factor=0.182, temperature=25, num_bends=10
        )

        # Bob's random basis choice
        bob_basis = random.choice(['C', 'X'])
        if bob_basis == 'C':
            phase_bob = 0
        else:
            phase_bob = Vpi / 2

        # Apply Bob Phase Modulation
        E = phase_modulator(E, phase_bob, 'X')  # Bob's phase shift
        # Measurement
        Ex, Ey = optics.pbs(E)
        photons_x = detector.output(E=Ex, area=1, bandwidth=1e6)
        photons_y = detector.output(E=Ey, area=1, bandwidth=1e6)
        ratio = abs((photons_x - photons_y) / (photons_x + photons_y))
        # Determine Bob's bit
        if ratio > 0.001:
            if photons_x > photons_y:
                bob_bit = 0
            elif photons_y > photons_x:
                bob_bit = 1
        else:
            # If photon counts are equal, this might suggest a +45 or -45 polarization state
            # randomly assign a bit
            bob_bit = random.randint(0, 1)

        if show_pol:
            polarimeter(E, title=f"Bit Number {num_bits},Bob Bit/Basis = {bob_bit} / {bob_basis}")
        # Store results
        alice_bits.append(alice_bit)
        alice_bases.append(alice_basis)
        bob_bits.append(bob_bit)
        bob_bases.append(bob_basis)

    # Sift keys to retain matching bases
    sifted_indices = [i for i in range(num_bits) if alice_bases[i] == bob_bases[i]]
    sifted_alice = [alice_bits[i] for i in sifted_indices]
    sifted_bob = [bob_bits[i] for i in sifted_indices]

    # Calculate QBER using a portion of the sifted key
    qber = 0.0
    if len(sifted_alice) > 0:
        sample_size = min(len(sifted_alice) // 2, 100)  # Check up to 100 bits
        if sample_size > 0:
            error_count = 0
            for i in random.sample(range(len(sifted_alice)), sample_size):
                if sifted_alice[i] != sifted_bob[i]:
                    error_count += 1
            qber = error_count / sample_size

    return [alice_bits, alice_bases, bob_bits, bob_bases, sifted_indices, sifted_alice, sifted_bob, qber]
