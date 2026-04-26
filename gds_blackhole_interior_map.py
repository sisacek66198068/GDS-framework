import numpy as np
import matplotlib.pyplot as plt

dt = 1e-4
T = 20.0
t = np.arange(0, T, dt)

G = 1.0
c = 1.0
M0 = 1.0
R_s = 2 * G * M0 / c**2

rho_c = 50.0
n = 3.0
A = 4.0

R_core = 0.08
k_bounce = 0.35
transfer_strength = 0.08

R_A, V_A = 3.0, -0.15
R_B, V_B = 0.05, 0.0
M_A, M_B = M0, 0.0

probe_r = 0.20
probe_v = 0.0
returned = False

hist = []

def rho_eff(M, R):
    return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

def Gamma(rho):
    return 1.0 + A / (1.0 + (rho / rho_c)**n)

for ti in t:
    rhoA = rho_eff(M_A, R_A)
    gamA = Gamma(rhoA)

    grav_A = -G * M_A / (R_A**2 + R_core**2)
    sat = (rhoA / rho_c)**n / (1.0 + (rhoA / rho_c)**n)
    bounce = k_bounce * sat / (R_A + R_core)
    gate_R = np.exp(-(R_A / (3 * R_core))**2)

    acc_A = grav_A / gamA + bounce

    dM = transfer_strength * sat * gate_R * M_A * dt
    dM = min(dM, M_A)

    M_A -= dM
    M_B += dM

    acc_B = 0.0
    if M_B > 0:
        acc_B += G * M_B / (R_B**2 + R_core**2) * (sat * gate_R)
        acc_B += 0.03 * M_B / (R_B + R_core)

    V_A += acc_A * dt
    R_A += V_A * dt

    if R_A < R_core:
        R_A = R_core
        V_A = abs(V_A) * 0.15

    V_B += acc_B * dt
    R_B += V_B * dt

    if M_B > 0:
        probe_acc = G * M_B / (probe_r**2 + R_core**2) + 0.02 * probe_r
        probe_v += probe_acc * dt
        probe_r += probe_v * dt
        if probe_r <= R_core:
            returned = True

    H_B = V_B / R_B if R_B > 0 else 0.0
    flux = dM / dt

    hist.append([
        ti, R_A, V_A, M_A, rhoA, gamA, sat, gate_R,
        flux, M_B, R_B, V_B, H_B, probe_r, probe_v
    ])

hist = np.array(hist)

ti, RA, VA, MA, rhoA, gamA, sat, gate, flux, MB, RB, VB, HB, pr, pv = hist.T

def first_time(mask):
    idx = np.where(mask)[0]
    return float(ti[idx[0]]) if len(idx) else None

t_horizon = first_time(RA < R_s)
t_sat_10 = first_time(sat > 0.10)
t_sat_50 = first_time(sat > 0.50)
t_gate = first_time(gate > 0.10)
t_transfer = first_time(MB > 1e-6)
t_core = first_time(RA <= R_core * 1.000001)
t_B_expand = first_time(RB > 0.051)

print("=== GDS BLACK-HOLE INTERIOR MECHANISM MAP ===")
print(f"Schwarzschild-like horizon R_s       = {R_s:.6e}")
print(f"Core radius R_core                   = {R_core:.6e}")
print()
print("EVENT TIMES")
print(f"horizon crossing time                = {t_horizon}")
print(f"saturation 10% time                  = {t_sat_10}")
print(f"saturation 50% time                  = {t_sat_50}")
print(f"interlink gate activation time       = {t_gate}")
print(f"mass transfer start time             = {t_transfer}")
print(f"core contact time                    = {t_core}")
print(f"B expansion visible time             = {t_B_expand}")
print()
print("FINAL / EXTREME METRICS")
print(f"R_A_min                              = {RA.min():.12e}")
print(f"rho_A_max                            = {rhoA.max():.12e}")
print(f"Gamma_A_min                          = {gamA.min():.12e}")
print(f"Gamma_A_max                          = {gamA.max():.12e}")
print(f"saturation_max                       = {sat.max():.12e}")
print(f"gate_max                             = {gate.max():.12e}")
print(f"flux_max                             = {flux.max():.12e}")
print(f"M_A_final                            = {MA[-1]:.12e}")
print(f"M_B_final                            = {MB[-1]:.12e}")
print(f"mass_total_final                     = {(MA[-1]+MB[-1]):.12e}")
print(f"R_B_final                            = {RB[-1]:.12e}")
print(f"V_B_final                            = {VB[-1]:.12e}")
print(f"H_B_late_mean                        = {HB[int(0.8*len(HB)):].mean():.12e}")
print(f"probe_final_radius                   = {pr[-1]:.12e}")
print(f"probe_final_velocity                 = {pv[-1]:.12e}")
print(f"probe_returned_to_core               = {returned}")
print()
print("INTERPRETATION")
if (RA.min() < R_s and RA.min() > 0 and MB[-1] > 1e-6 and RB[-1] > RB[0] and not returned):
    print("Confirmed mechanism:")
    print("horizon crossing -> finite saturated core -> interlink gate -> conserved mass transfer -> FRW-like daughter expansion -> no geodesic return")
else:
    print("Mechanism not fully confirmed in this run.")

np.savetxt(
    "gds_blackhole_interior_map.csv",
    hist,
    delimiter=",",
    header="t,RA,VA,MA,rhoA,GammaA,saturation,gate,flux,MB,RB,VB,HB,probe_r,probe_v",
    comments=""
)

# 1 radius map
plt.figure(figsize=(9,5))
plt.plot(ti, RA, label="A radius")
plt.plot(ti, RB, label="B radius")
plt.axhline(R_s, linestyle="--", label="Schwarzschild-like horizon")
plt.axhline(R_core, linestyle=":", label="finite core")
plt.xlabel("time [toy units]")
plt.ylabel("radius")
plt.title("Black-Hole Interior Map: Horizon, Core, Daughter Expansion")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("bh_map_01_radius.png", dpi=200)

# 2 density/saturation
plt.figure(figsize=(9,5))
plt.plot(ti, rhoA, label="effective density A")
plt.axhline(rho_c, linestyle="--", label="critical density")
plt.yscale("log")
plt.xlabel("time [toy units]")
plt.ylabel("density")
plt.title("Density Growth and Saturation Threshold")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("bh_map_02_density.png", dpi=200)

# 3 saturation/gate/flux
plt.figure(figsize=(9,5))
plt.plot(ti, sat, label="saturation")
plt.plot(ti, gate, label="interlink gate")
plt.plot(ti, flux / max(flux.max(), 1e-30), label="normalized mass flux")
plt.xlabel("time [toy units]")
plt.ylabel("activation")
plt.title("Interlink Activation: Saturation, Gate, Mass Flux")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("bh_map_03_activation_flux.png", dpi=200)

# 4 mass conservation
plt.figure(figsize=(9,5))
plt.plot(ti, MA, label="mass A")
plt.plot(ti, MB, label="mass B")
plt.plot(ti, MA+MB, "--", label="total mass")
plt.xlabel("time [toy units]")
plt.ylabel("mass")
plt.title("Conserved Transfer From Parent Interior to Daughter Region")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("bh_map_04_mass_transfer.png", dpi=200)

# 5 Hubble-like expansion
plt.figure(figsize=(9,5))
plt.plot(ti, HB, label="H_B = V_B / R_B")
plt.xlabel("time [toy units]")
plt.ylabel("H_B")
plt.title("FRW-like Expansion Rate of Daughter Region")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("bh_map_05_HB.png", dpi=200)

# 6 geodesic non-return
plt.figure(figsize=(9,5))
plt.plot(ti, pr, label="probe radius in B")
plt.axhline(R_core, linestyle=":", label="return/core boundary")
plt.xlabel("time [toy units]")
plt.ylabel("probe radius")
plt.title("Geodesic Non-Return Test")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("bh_map_06_geodesic_nonreturn.png", dpi=200)

print()
print("SAVED:")
print("  gds_blackhole_interior_map.csv")
print("  bh_map_01_radius.png")
print("  bh_map_02_density.png")
print("  bh_map_03_activation_flux.png")
print("  bh_map_04_mass_transfer.png")
print("  bh_map_05_HB.png")
print("  bh_map_06_geodesic_nonreturn.png")
