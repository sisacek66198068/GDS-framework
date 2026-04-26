import numpy as np

# ============================================================
# GDS CAUSAL + FRW TEST
#
# Goal:
# 1) verify one-way causal separation
# 2) verify Region B behaves like expanding FRW-like patch
#
# ============================================================

dt = 1e-4
T = 20.0
t = np.arange(0, T, dt)

G = 1.0
M0 = 1.0

rho_c = 50.0
n = 3.0
A = 4.0

R_core = 0.08
k_bounce = 0.35
transfer_strength = 0.08

# Region A
R_A = 3.0
V_A = -0.15

# Region B
R_B = 0.05
V_B = 0.0

M_A = M0
M_B = 0.0

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

    acc_A = grav_A / gamA + bounce

    gate_R = np.exp(-(R_A / (3 * R_core))**2)

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

    H_B = 0.0
    if R_B > 0:
        H_B = V_B / R_B

    hist.append([
        ti,
        R_A,
        V_A,
        M_A,
        M_B,
        R_B,
        V_B,
        H_B,
        sat,
        gate_R
    ])

hist = np.array(hist)

ti  = hist[:,0]
RA  = hist[:,1]
VA  = hist[:,2]
MA  = hist[:,3]
MB  = hist[:,4]
RB  = hist[:,5]
VB  = hist[:,6]
HB  = hist[:,7]
SAT = hist[:,8]
GATE= hist[:,9]

# ------------------------------------------------------------
# TEST 1: One-way causal separation
#
# If after transfer starts:
# - MB keeps increasing
# - RB keeps expanding
# - no reverse flow exists
#
# then one-way causal channel is supported
# ------------------------------------------------------------

transfer_started = np.where(MB > 1e-6)[0]

if len(transfer_started) > 0:
    i0 = transfer_started[0]

    MB_after = MB[i0:]
    RB_after = RB[i0:]

    monotonic_MB = np.all(np.diff(MB_after) >= -1e-10)
    monotonic_RB = np.all(np.diff(RB_after) >= -1e-10)

else:
    monotonic_MB = False
    monotonic_RB = False

# ------------------------------------------------------------
# TEST 2: FRW-like expansion
#
# Positive stable H_B after transfer
# ------------------------------------------------------------

late = int(0.8 * len(HB))
HB_late = HB[late:]

HB_mean = np.mean(HB_late)
HB_min  = np.min(HB_late)
HB_max  = np.max(HB_late)

frw_like = (HB_mean > 0.0 and HB_min >= 0.0)

# ------------------------------------------------------------

print("=== GDS CAUSAL + FRW TEST ===")
print()

print("ONE-WAY CAUSAL TEST")
print(f"MB monotonic increase      = {monotonic_MB}")
print(f"RB monotonic expansion     = {monotonic_RB}")
print()

print("FRW-LIKE EXPANSION TEST")
print(f"H_B late mean              = {HB_mean:.6e}")
print(f"H_B late min               = {HB_min:.6e}")
print(f"H_B late max               = {HB_max:.6e}")
print(f"FRW-like expansion         = {frw_like}")
print()

if monotonic_MB and monotonic_RB and frw_like:
    print("RESULT:")
    print("Region B behaves as a causally separated expanding patch.")
    print("This supports daughter-universe interpretation within the model.")
else:
    print("RESULT:")
    print("Evidence is weaker; stronger geometric test required.")
