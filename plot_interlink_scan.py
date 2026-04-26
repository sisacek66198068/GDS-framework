import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("gds_interlink_scan_results.csv")

print("=== SUMMARY ===")
print(df.groupby(["rho_c","n","A"])["success"].mean().reset_index().to_string(index=False))

plt.figure(figsize=(8,5))
plt.hist(df["MB_final"], bins=40)
plt.xlabel("Final transferred mass M_B")
plt.ylabel("count")
plt.title("Distribution of Transferred Mass in GDS Interlink Scan")
plt.grid(True)
plt.tight_layout()
plt.savefig("scan_MB_hist.png", dpi=200)

plt.figure(figsize=(8,5))
plt.scatter(df["MB_final"], df["RB_final"], s=12)
plt.xlabel("Final transferred mass M_B")
plt.ylabel("Final expansion radius R_B")
plt.title("Expansion Strength vs Transferred Mass")
plt.grid(True)
plt.tight_layout()
plt.savefig("scan_MB_vs_RB.png", dpi=200)

plt.figure(figsize=(8,5))
plt.scatter(df["rhoA_max"], df["MB_final"], s=12)
plt.xscale("log")
plt.xlabel("Maximum effective density in A")
plt.ylabel("Final transferred mass M_B")
plt.title("Density Trigger vs Interlink Transfer")
plt.grid(True)
plt.tight_layout()
plt.savefig("scan_density_vs_MB.png", dpi=200)

print("SAVED:")
print("  scan_MB_hist.png")
print("  scan_MB_vs_RB.png")
print("  scan_density_vs_MB.png")
