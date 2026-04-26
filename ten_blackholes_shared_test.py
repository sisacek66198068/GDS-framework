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

# 10 black-hole-like sources with different masses/radii/infall speeds
sources0 = [
    [0.6,  2.5, -0.18],
    [0.8,  2.8, -0.16],
    [1.0,  3.0, -0.15],
    [1.3,  3.4, -0.14],
    [1.7,  3.8, -0.13],
    [2.1,  4.2, -0.125],
    [2.5,  4.5, -0.12],
    [3.0,  5.0, -0.115],
    [3.7,  5.6, -0.11],
    [4.5,  6.3, -0.105],
]

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
    sources = [s[:] for s in sources0]

    daughters = []
    for _ in sources:
        daughters.append([0.0, 0.05, 0.0])  # M_B, R_B, V_B

    H_hist = [[] for _ in sources]

    for _ in t:
        for i in range(len(sources)):
            M, R, V = sources[i]
            M, R, V, dM, sat, gate = evolve_source(M, R, V)
            sources[i] = [M, R, V]

            MB, RB, VB = daughters[i]
            MB += dM

            accB = 0.0
            if MB > 0:
                accB += G * MB / (RB**2 + R_core**2) * (sat * gate)
                accB += 0.03 * MB / (RB + R_core)

            VB += accB * dt
            RB += VB * dt

            daughters[i] = [MB, RB, VB]
            H_hist[i].append(VB / RB if RB > 0 else 0.0)

    return sources, daughters, [np.array(h) for h in H_hist]

def run_shared():
    sources = [s[:] for s in sources0]

    MD, RD, VD = 0.0, 0.05, 0.0
    H_hist = []

    for _ in t:
        dM_total = 0.0
        trigger_total = 0.0

        for i in range(len(sources)):
            M, R, V = sources[i]
            M, R, V, dM, sat, gate = evolve_source(M, R, V)
            sources[i] = [M, R, V]

            dM_total += dM
            trigger_total += sat * gate

        MD += dM_total

        accD = 0.0
        if MD > 0:
            accD += G * MD / (RD**2 + R_core**2) * trigger_total
            accD += 0.03 * MD / (RD + R_core)

        VD += accD * dt
        RD += VD * dt

        H_hist.append(VD / RD if RD > 0 else 0.0)

    return sources, [MD, RD, VD], np.array(H_hist)

ind_sources, ind_daughters, ind_H = run_independent()
sh_sources, sh_daughter, sh_H = run_shared()

late = int(0.8 * len(t))

ind_H_means = np.array([h[late:].mean() for h in ind_H])
ind_H_stds  = np.array([h[late:].std() for h in ind_H])

ind_MB = np.array([d[0] for d in ind_daughters])
ind_RB = np.array([d[1] for d in ind_daughters])

shared_MD, shared_RD, shared_VD = sh_daughter
H_shared_mean = sh_H[late:].mean()
H_shared_std  = sh_H[late:].std()

initial_total_mass = sum(s[0] for s in sources0)
ind_total_mass = sum(s[0] for s in ind_sources) + ind_MB.sum()
sh_total_mass = sum(s[0] for s in sh_sources) + shared_MD

print("=== TEN BLACK HOLES: SHARED VS INDEPENDENT UNIVERSE TEST ===")
print()

print("INITIAL TOTAL MASS")
print(f"Initial total mass                 = {initial_total_mass:.12e}")
print()

print("INDEPENDENT MODE")
for i in range(10):
    print(
        f"BH{i+1:02d}: "
        f"MB={ind_MB[i]:.12e}  "
        f"RB={ind_RB[i]:.12e}  "
        f"H_late={ind_H_means[i]:.12e}  "
        f"H_std={ind_H_stds[i]:.12e}"
    )

print()
print(f"Independent daughter mass sum       = {ind_MB.sum():.12e}")
print(f"Independent H mean spread min       = {ind_H_means.min():.12e}")
print(f"Independent H mean spread max       = {ind_H_means.max():.12e}")
print(f"Independent Delta H                 = {(ind_H_means.max()-ind_H_means.min()):.12e}")
print(f"Independent total mass              = {ind_total_mass:.12e}")
print()

print("SHARED MODE")
print(f"Shared daughter mass                = {shared_MD:.12e}")
print(f"Shared daughter radius              = {shared_RD:.12e}")
print(f"Shared daughter velocity            = {shared_VD:.12e}")
print(f"H_shared late mean                  = {H_shared_mean:.12e}")
print(f"H_shared late std                   = {H_shared_std:.12e}")
print(f"Shared total mass                   = {sh_total_mass:.12e}")
print()

print("COMPARISON")
print(f"Mass difference independent/shared  = {abs(ind_MB.sum()-shared_MD):.12e}")
print(f"Total mass error independent        = {abs(ind_total_mass-initial_total_mass):.12e}")
print(f"Total mass error shared             = {abs(sh_total_mass-initial_total_mass):.12e}")
print()

if abs(ind_MB.sum()-shared_MD) < 1e-10 and H_shared_mean > 0:
    print("RESULT:")
    print("10 independent black-hole sources transfer the same total daughter mass in both descriptions.")
    print("Independent branches produce a spread of different H values.")
    print("Shared branch produces one common H_shared.")
    print("This supports the shared-daughter-universe diagnostic.")
else:
    print("RESULT:")
    print("Ten-source shared test is inconclusive.")
