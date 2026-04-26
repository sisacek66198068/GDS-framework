import numpy as np

# ============================================================
# GDS Hawking Evaporation Compatibility Test
#
# Goal:
# compare standard evaporation tendency with
# finite-core remnant interpretation
# ============================================================

# toy natural units
G = 1.0
c = 1.0
hbar = 1.0
kB = 1.0

# initial BH mass
M_initial = 1.0

# finite core radius from GDS
R_core = 0.08

print("=== GDS HAWKING EVAPORATION TEST ===")
print()

# Schwarzschild radius
def Rs(M):
    return 2.0 * G * M / c**2

# Hawking temperature (toy)
def T_H(M):
    # proportional form
    return 1.0 / (8.0 * np.pi * G * M)

# estimate mass corresponding to finite core
# if Rs = R_core  -> remnant threshold

M_remnant = R_core * c**2 / (2.0 * G)

T_initial = T_H(M_initial)
T_remnant = T_H(M_remnant)

print(f"Initial BH mass M_initial       = {M_initial:.12e}")
print(f"Initial Schwarzschild radius    = {Rs(M_initial):.12e}")
print()

print(f"Finite core radius R_core       = {R_core:.12e}")
print(f"Remnant threshold mass          = {M_remnant:.12e}")
print()

print(f"Initial Hawking temperature     = {T_initial:.12e}")
print(f"Remnant-scale temperature       = {T_remnant:.12e}")
print()

ratio = T_remnant / T_initial
print(f"Temperature increase factor     = {ratio:.12e}")
print()

print("INTERPRETATION:")
print("Standard Hawking evaporation drives temperature upward as mass decreases.")

if M_remnant > 0:
    print("In GDS, evaporation does not require collapse to zero mass.")
    print("A finite remnant threshold naturally appears when horizon approaches the finite core scale.")
    print("This replaces singular end-state by a finite-core endpoint.")
else:
    print("No remnant threshold detected.")
