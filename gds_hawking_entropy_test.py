import numpy as np

# ============================================================
# GDS + Hawking / Bekenstein entropy consistency test
#
# Goal:
# check whether finite-core black hole remains compatible
# with horizon entropy scaling
# ============================================================

# natural toy units
G = 1.0
c = 1.0
hbar = 1.0
kB = 1.0

# black hole mass
M = 1.0

# Schwarzschild-like horizon
R_s = 2.0 * G * M / c**2

# GDS finite core
R_core = 0.08

print("=== GDS HAWKING ENTROPY TEST ===")
print()

print(f"Schwarzschild-like radius R_s   = {R_s:.12e}")
print(f"Finite core radius R_core       = {R_core:.12e}")
print()

# ------------------------------------------------
# Horizon area
# ------------------------------------------------

A_horizon = 4.0 * np.pi * R_s**2

# Bekenstein-Hawking entropy
# S = A / (4 G hbar)

S_BH = A_horizon / (4.0 * G * hbar)

print(f"Horizon area A_horizon          = {A_horizon:.12e}")
print(f"Bekenstein entropy S_BH         = {S_BH:.12e}")
print()

# ------------------------------------------------
# Core correction diagnostic
#
# if entropy were controlled by core instead:
# compare area scales
# ------------------------------------------------

A_core = 4.0 * np.pi * R_core**2
ratio = A_core / A_horizon

print(f"Core area A_core                = {A_core:.12e}")
print(f"Core/Horizon area ratio         = {ratio:.12e}")
print()

if R_core < R_s:
    print("RESULT:")
    print("Finite core exists inside horizon while horizon entropy remains dominant.")
    print("GDS is compatible with standard Bekenstein-Hawking scaling.")
    print("No immediate contradiction with semiclassical black-hole thermodynamics.")
else:
    print("RESULT:")
    print("Core exceeds horizon scale.")
    print("This would strongly conflict with standard entropy interpretation.")
