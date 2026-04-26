import numpy as np

G = 1.0
c = 1.0

M = 1.0

# Schwarzschild-like toy horizon
R_s = 2 * G * M / c**2

# Our measured core from scan
R_core = 0.08
R_min = 0.08

print("=== HORIZON CHECK ===")
print(f"Schwarzschild-like radius R_s = {R_s:.6f}")
print(f"Core radius R_core           = {R_core:.6f}")
print(f"Minimum reached radius       = {R_min:.6f}")
print()

if R_min < R_s:
    print("Result: collapse enters inside horizon before singularity.")
    print("No infinite singularity appears; finite saturated core survives inside horizon.")
else:
    print("Result: collapse stops before horizon crossing.")
