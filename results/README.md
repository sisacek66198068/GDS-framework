# GDS numerical evidence — publication-extracted CSV / Číselné výsledky z publikací

These four CSV files transcribe selected **reported values from public PDF publications**. They are **not** original simulation output files, rerun calculations, independent replications or raw observational datasets.

| CSV | Contents | Primary publication and location |
|---|---|---|
| [gds_gigaall_published.csv](gds_gigaall_published.csv) | GIGAALL best -log L for LCDM, GDS/GSC and STRING-LIKE | *GDS v3 — A Unified Dynamical Field Framework*, section 6.3, Table 6 |
| [gds_collapse_rk4_published.csv](gds_collapse_rk4_published.csv) | RK4 step-size vs minimum radius, kept separately by publication | *From Singularity to Stability*, section 3; *Regenerative Metastability*, section 4.1.2, Table 2 |
| [gds_finite_core_mass_published.csv](gds_finite_core_mass_published.csv) | Total and saturated-core masses, radii and density for **scaled modeled cases** labeled Sgr A* and M82 | *From Singularity to Stability*, sections 4.1 and 4.2 |
| [regenerative_metastability_n3000_published.csv](regenerative_metastability_n3000_published.csv) | Robust N3000 p999 and maximum comparisons | *Emergent Metastability Corridors*, section 5.1, Table 1 |

The original PDF papers are linked from [GDS publication evidence map](../publications/GDS_EVIDENCE_MAP.md).

## Scientific caveats / Vědecké výhrady

- The **Sgr A*** and **M82** rows describe numerical *scaled model cases*, **not independently measured astrophysical core masses**. Do not represent these as direct mass measurements of real black holes.
- The two published RK4 convergence sequences have different rounding and sample points; they are **kept as separate source-labeled series** rather than silently reconciled. Negative coarse-step minimum radii in one publication are preserved as reported and require numerical interpretation.
- The source PDFs do not consistently specify the physical units of the RK4 minimum-radius entries; the CSV marks units as `unspecified`.
- The best-fit GIGAALL difference does not alone establish model selection or statistical significance without parameter-count/complexity assessment and independent reproduction.
- For scientific reuse, consult the complete paper, methodology, figures, and limitations. Original terminal CSVs should replace or supplement these extracted tables once their provenance is verified.

**Status:** public-paper numerical extracts, not raw benchmark datasets. / **Stav:** číselné výpisy z veřejných publikací, nikoli surová data z výpočtů.
