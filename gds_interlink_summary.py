import pandas as pd

df = pd.read_csv("gds_interlink_scan_results.csv")

print("=== GDS INTERLINK FINAL SUMMARY ===")
print(f"Total scan cases              : {len(df)}")
print(f"Successful interlink cases    : {int(df['success'].sum())}")
print(f"Success rate                  : {df['success'].mean():.6f}")
print(f"Minimum RA_min                : {df['RA_min'].min():.6e}")
print(f"Maximum rhoA_max              : {df['rhoA_max'].max():.6e}")
print(f"Maximum transferred mass MB   : {df['MB_final'].max():.6e}")
print(f"Mean transferred mass MB      : {df['MB_final'].mean():.6e}")
print(f"Maximum expansion radius RB   : {df['RB_final'].max():.6e}")
print(f"Mean expansion radius RB      : {df['RB_final'].mean():.6e}")
print(f"Mass conservation min         : {df['mass_total'].min():.12e}")
print(f"Mass conservation max         : {df['mass_total'].max():.12e}")

R_s = 2.0
R_min = df["RA_min"].min()
print()
print("=== HORIZON DIAGNOSTIC ===")
print(f"Schwarzschild-like R_s        : {R_s:.6e}")
print(f"Minimum reached radius        : {R_min:.6e}")
print(f"Inside horizon?               : {R_min < R_s}")
