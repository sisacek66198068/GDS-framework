import numpy as np

# ============================================================
# GDS JWST EARLY STRUCTURE TEST
#
# Goal:
# compare daughter-universe early structure window
# with observed surprisingly early massive galaxies
# ============================================================

print("=== GDS JWST EARLY STRUCTURE TEST ===")
print()

# from our previous structure window result
t_start = 8.0494
t_end   = 8.5137

# normalized toy total evolution time
T_total = 20.0

# fraction of evolution when structure becomes possible
f_start = t_start / T_total
f_end   = t_end / T_total

print(f"Structure window start          = {t_start:.6f}")
print(f"Structure window end            = {t_end:.6f}")
print(f"Total modeled evolution         = {T_total:.6f}")
print()

print(f"Normalized start fraction       = {f_start:.12e}")
print(f"Normalized end fraction         = {f_end:.12e}")
print()

# comparison concept:
# LCDM expects later efficient structure growth
# GDS allows earlier localized collapse

if f_start < 0.5:
    print("RESULT:")
    print("Structure formation begins before half of total evolution time.")
    print("This supports accelerated early structure formation.")
    print("This is qualitatively consistent with JWST early massive-galaxy tension.")
else:
    print("RESULT:")
    print("Structure formation begins too late for strong JWST tension relief.")

print()
print("NEXT STEP:")
print("Map this normalized timing to observational redshift z ~ 8–15.")
