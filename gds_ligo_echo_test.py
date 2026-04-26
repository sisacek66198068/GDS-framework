import numpy as np

# ============================================================
# GDS LIGO ECHO TEST
# Finite-core black hole vs classical singularity
#
# Goal:
# estimate gravitational-wave echo delay caused by
# reflection from finite saturated core instead of singularity
# ============================================================

# toy units first
G = 1.0
c = 1.0
M = 1.0

# Schwarzschild-like horizon
R_s = 2.0 * G * M / c**2

# GDS finite core from our model
R_core = 0.08

print("=== GDS LIGO ECHO TEST ===")
print()

print(f"Schwarzschild-like radius R_s   = {R_s:.12e}")
print(f"Finite core radius R_core       = {R_core:.12e}")
print()

# ------------------------------------------------
# Echo delay estimate
#
# Toy approximation:
# delta_t ~ 2 * integral(dr / (1 - R_s/r))
#
# simplified logarithmic estimate:
#
# delta_t_echo ~ 2 * R_s * ln(R_s / R_core)
#
# (order-of-magnitude diagnostic)
# ------------------------------------------------

delta_t_echo = 2.0 * R_s * np.log(R_s / R_core)

print(f"Estimated echo delay (toy)      = {delta_t_echo:.12e}")
print()

# dimensionless comparison
ratio = R_s / R_core
print(f"Horizon/Core ratio              = {ratio:.12e}")
print()

if R_core < R_s:
    print("RESULT:")
    print("Finite core lies inside horizon.")
    print("Reflection from finite core can generate delayed GW echoes.")
    print("This is absent for a true singularity.")
else:
    print("RESULT:")
    print("Core is outside horizon.")
    print("No true black-hole echo interpretation.")

print()
print("NEXT STEP:")
print("Compare predicted delay with LIGO ringdown residual structure.")
