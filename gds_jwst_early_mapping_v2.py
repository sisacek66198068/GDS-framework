import numpy as np

print("=== GDS JWST EARLY MAPPING v2 ===")
print()

# toy structure window fraction
f_start = 0.40247
f_end   = 0.425685

# Early-universe calibration window:
# assume daughter-region birth maps to first 1 Gyr,
# not full 13.8 Gyr.
T_early = 1.0  # Gyr

age_start = f_start * T_early
age_end   = f_end * T_early

# matter-era approximate inversion:
# t(z) ≈ 17.3 Gyr / (1+z)^(3/2)
# rough high-z LCDM-like scaling
C = 17.3

z_start = (C / age_start)**(2.0/3.0) - 1.0
z_end   = (C / age_end)**(2.0/3.0) - 1.0

print(f"Early-calibrated age start [Gyr] = {age_start:.6f}")
print(f"Early-calibrated age end [Gyr]   = {age_end:.6f}")
print()
print(f"Approx redshift start            = {z_start:.6f}")
print(f"Approx redshift end              = {z_end:.6f}")
print()

if 8 <= z_end <= 20 or 8 <= z_start <= 20:
    print("RESULT:")
    print("GDS early structure window maps into JWST-relevant high-redshift range.")
else:
    print("RESULT:")
    print("Mapping still misses JWST range; needs calibrated CLASS/GDS background.")

print()
print("NOTE:")
print("This is an early-universe calibration estimate, not a final cosmological fit.")
