import numpy as np

print("=== GDS JWST REDSHIFT MAPPING TEST ===")
print()

# Toy normalized structure window
f_start = 0.40247
f_end   = 0.425685

# Approximate age of universe in Gyr
t0 = 13.8

# Map toy fraction to cosmic age
age_start = f_start * t0
age_end   = f_end * t0

print(f"Mapped cosmic age start [Gyr]   = {age_start:.6f}")
print(f"Mapped cosmic age end [Gyr]     = {age_end:.6f}")
print()

# Very crude matter-era inversion:
# t(z) ≈ t0 / (1+z)^(3/2)
# so 1+z ≈ (t0/t)^(2/3)

z_start = (t0 / age_start)**(2.0/3.0) - 1.0
z_end   = (t0 / age_end)**(2.0/3.0) - 1.0

print(f"Approx mapped redshift start    = {z_start:.6f}")
print(f"Approx mapped redshift end      = {z_end:.6f}")
print()

print("NOTE:")
print("This is a crude normalization test, not yet a calibrated cosmological fit.")
print("For JWST-level comparison, the next step is calibrated mapping using CLASS/GDS background evolution.")
