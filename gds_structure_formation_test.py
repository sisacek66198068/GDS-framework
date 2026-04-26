import numpy as np

G = 1.0

# values from daughter-region run
rho_B = 1.744907941871e-04
H_B   = 1.177389014658e-01

t_ff = 1.0 / np.sqrt(G * rho_B)
t_H  = 1.0 / H_B

ratio = t_ff / t_H

print("=== GDS STRUCTURE FORMATION TEST ===")
print()
print(f"Late daughter density rho_B   = {rho_B:.12e}")
print(f"Late Hubble-like H_B          = {H_B:.12e}")
print()
print(f"Free-fall time t_ff           = {t_ff:.12e}")
print(f"Expansion time t_H            = {t_H:.12e}")
print(f"Ratio t_ff / t_H              = {ratio:.12e}")
print()

if t_ff < t_H:
    print("RESULT:")
    print("Gravity can overcome expansion locally.")
    print("Structure formation is dynamically allowed.")
    print("Proto-galaxy / daughter-universe formation is supported.")
else:
    print("RESULT:")
    print("Expansion dominates over collapse.")
    print("Large-scale structure formation is weak in this setup.")
