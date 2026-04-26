import numpy as np

M_solar = 1.98847e30
G = 6.67430e-11
c = 299792458.0

ratio = 25.0

masses = [10, 20, 30, 50, 70, 100, 150]

print("=== GDS LIGO ECHO MASS TABLE ===")
print("M_solar   R_s_m             R_core_m         echo_delay_s       echo_delay_ms")

for m in masses:
    M = m * M_solar
    R_s = 2.0 * G * M / c**2
    R_core = R_s / ratio
    dt = 2.0 * R_s / c * np.log(ratio)
    print(f"{m:7.1f}   {R_s: .12e}   {R_core: .12e}   {dt: .12e}   {dt*1000: .6f}")
