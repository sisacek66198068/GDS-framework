import numpy as np

# toy-to-physical scaling estimate

# reference stellar BH merger
# choose ~30 solar masses effective remnant

M_solar = 1.98847e30
G_SI = 6.67430e-11
c = 299792458.0

M_bh = 30.0 * M_solar

# Schwarzschild radius
R_s = 2.0 * G_SI * M_bh / c**2

# use same ratio as toy model
ratio = 25.0
R_core = R_s / ratio

# echo delay estimate
delta_t = 2.0 * R_s / c * np.log(ratio)

print("=== GDS LIGO PHYSICAL ECHO ESTIMATE ===")
print()

print(f"Black hole mass [kg]            = {M_bh:.12e}")
print(f"Schwarzschild radius [m]        = {R_s:.12e}")
print(f"Finite core radius [m]          = {R_core:.12e}")
print(f"Echo delay [seconds]            = {delta_t:.12e}")
print()

if delta_t < 1.0:
    print("RESULT:")
    print("Echo should appear shortly after merger ringdown.")
    print("This is directly testable in LIGO residuals.")
else:
    print("RESULT:")
    print("Echo delay is long and requires extended residual search.")
