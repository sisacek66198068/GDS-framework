import numpy as np
import matplotlib.pyplot as plt

M_solar = 1.98847e30
G = 6.67430e-11
c = 299792458.0

ratio = 25.0

masses = np.linspace(5, 200, 300)
echo_ms = []

for m in masses:
    M = m * M_solar
    R_s = 2.0 * G * M / c**2
    dt = 2.0 * R_s / c * np.log(ratio)
    echo_ms.append(dt * 1000.0)

echo_ms = np.array(echo_ms)

plt.figure(figsize=(9,5))
plt.plot(masses, echo_ms)
plt.xlabel("Black hole remnant mass [solar masses]")
plt.ylabel("Predicted echo delay [ms]")
plt.title("GDS Finite-Core Prediction: Ringdown Echo Delay")
plt.grid(True)
plt.tight_layout()
plt.savefig("gds_ligo_echo_delay_prediction.png", dpi=200)

print("=== GDS LIGO ECHO PREDICTION PLOT ===")
print("SAVED: gds_ligo_echo_delay_prediction.png")
print()
print("Key prediction:")
print("Echo delay scales approximately linearly with black-hole remnant mass.")
