import numpy as np

dt = 1e-4
T = 20.0
t = np.arange(0, T, dt)

G = 1.0
rho_c = 50.0
n = 3.0
A = 4.0
R_core = 0.08
k_bounce = 0.35
transfer_strength = 0.08

def rho_eff(M, R):
    return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

def Gamma(rho):
    return 1.0 + A / (1.0 + (rho / rho_c)**n)

def evolve_source(M, R, V):
    rho = rho_eff(M, R)
    gam = Gamma(rho)
    sat = (rho / rho_c)**n / (1.0 + (rho / rho_c)**n)
    gate = np.exp(-(R / (3 * R_core))**2)

    grav = -G * M / (R**2 + R_core**2)
    bounce = k_bounce * sat / (R + R_core)
    acc = grav / gam + bounce

    dM = transfer_strength * sat * gate * M * dt
    dM = min(dM, M)

    M -= dM
    V += acc * dt
    R += V * dt

    if R < R_core:
        R = R_core
        V = abs(V) * 0.15

    return M, R, V, dM, sat, gate

def run_independent():
    M1, R1, V1 = 1.0, 3.0, -0.15
    M2, R2, V2 = 2.5, 4.5, -0.12

    MB1, RB1, VB1 = 0.0, 0.05, 0.0
    MB2, RB2, VB2 = 0.0, 0.05, 0.0

    H1, H2 = [], []

    for _ in t:
        M1, R1, V1, dM1, sat1, gate1 = evolve_source(M1, R1, V1)
        M2, R2, V2, dM2, sat2, gate2 = evolve_source(M2, R2, V2)

        MB1 += dM1
        MB2 += dM2

        accB1 = 0.0
        if MB1 > 0:
            accB1 += G * MB1 / (RB1**2 + R_core**2) * (sat1 * gate1)
            accB1 += 0.03 * MB1 / (RB1 + R_core)

        accB2 = 0.0
        if MB2 > 0:
            accB2 += G * MB2 / (RB2**2 + R_core**2) * (sat2 * gate2)
            accB2 += 0.03 * MB2 / (RB2 + R_core)

        VB1 += accB1 * dt
        RB1 += VB1 * dt
        VB2 += accB2 * dt
        RB2 += VB2 * dt

        H1.append(VB1 / RB1 if RB1 > 0 else 0.0)
        H2.append(VB2 / RB2 if RB2 > 0 else 0.0)

    return {
        "M1": M1, "M2": M2,
        "MB1": MB1, "MB2": MB2,
        "RB1": RB1, "RB2": RB2,
        "VB1": VB1, "VB2": VB2,
        "H1": np.array(H1),
        "H2": np.array(H2),
        "total": M1 + M2 + MB1 + MB2,
    }

def run_shared():
    M1, R1, V1 = 1.0, 3.0, -0.15
    M2, R2, V2 = 2.5, 4.5, -0.12

    MD, RD, VD = 0.0, 0.05, 0.0
    H = []

    for _ in t:
        M1, R1, V1, dM1, sat1, gate1 = evolve_source(M1, R1, V1)
        M2, R2, V2, dM2, sat2, gate2 = evolve_source(M2, R2, V2)

        MD += dM1 + dM2

        accD = 0.0
        if MD > 0:
            trigger = sat1 * gate1 + sat2 * gate2
            accD += G * MD / (RD**2 + R_core**2) * trigger
            accD += 0.03 * MD / (RD + R_core)

        VD += accD * dt
        RD += VD * dt

        H.append(VD / RD if RD > 0 else 0.0)

    return {
        "M1": M1, "M2": M2,
        "MD": MD,
        "RD": RD,
        "VD": VD,
        "H": np.array(H),
        "total": M1 + M2 + MD,
    }

ind = run_independent()
sh = run_shared()

late = int(0.8 * len(t))

H1m = ind["H1"][late:].mean()
H2m = ind["H2"][late:].mean()
H1s = ind["H1"][late:].std()
H2s = ind["H2"][late:].std()

Hsm = sh["H"][late:].mean()
Hss = sh["H"][late:].std()

print("=== DUAL BLACK HOLE: SHARED VS INDEPENDENT ===")
print()

print("INDEPENDENT MODE")
print(f"B1 mass final                 = {ind['MB1']:.12e}")
print(f"B2 mass final                 = {ind['MB2']:.12e}")
print(f"B1 radius final               = {ind['RB1']:.12e}")
print(f"B2 radius final               = {ind['RB2']:.12e}")
print(f"H_B1 late mean                = {H1m:.12e}")
print(f"H_B2 late mean                = {H2m:.12e}")
print(f"H_B1 late std                 = {H1s:.12e}")
print(f"H_B2 late std                 = {H2s:.12e}")
print(f"Delta independent H           = {abs(H1m-H2m):.12e}")
print(f"Total mass independent        = {ind['total']:.12e}")
print()

print("SHARED MODE")
print(f"Shared mass final             = {sh['MD']:.12e}")
print(f"Shared radius final           = {sh['RD']:.12e}")
print(f"H_shared late mean            = {Hsm:.12e}")
print(f"H_shared late std             = {Hss:.12e}")
print(f"Total mass shared             = {sh['total']:.12e}")
print()

print("COMPARISON")
print(f"Independent total daughter mass = {(ind['MB1']+ind['MB2']):.12e}")
print(f"Shared daughter mass            = {sh['MD']:.12e}")
print(f"Mass difference                 = {abs((ind['MB1']+ind['MB2'])-sh['MD']):.12e}")
print()

if abs((ind['MB1']+ind['MB2'])-sh['MD']) < 1e-10 and Hsm > 0:
    print("RESULT:")
    print("Both descriptions conserve identical transferred mass.")
    print("Shared mode yields one common Hubble-like expansion number.")
    print("Independent mode yields separate H values.")
    print("This supports using shared H as the diagnostic for one daughter universe.")
else:
    print("RESULT:")
    print("Shared-vs-independent comparison is inconclusive.")
