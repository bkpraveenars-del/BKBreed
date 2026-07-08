# Example data provenance

The four CSVs in `inst/extdata/` are **simulated** demonstration datasets,
generated with a fixed random seed so results are reproducible. They mimic
arid-zone crop trials (pearl millet / wheat context) with realistic means,
variances and built-in trait correlations and G×E structure. They are for
teaching and package demonstration only — not real experimental records.

| File | Design | Structure |
|---|---|---|
| `pearlmillet_rbd.csv` | RBD | 12 genotypes × 3 reps, 6 traits |
| `factorial_rbd.csv` | Factorial RBD | 3 nitrogen × 4 varieties × 3 reps |
| `augmented_alpha.csv` | Augmented alpha-lattice | 4 checks + 20 test entries, 2 reps × 4 blocks |
| `mlt_stability.csv` | Multi-location trial | 10 genotypes × 5 environments × 3 reps |

Regenerate with the seed-based generator used during development
(seed = 2026). Built-in signal: yield positively linked to tillers, panicle
length and test weight; genotypes given differing environmental sensitivities
(bᵢ) so stability analysis is meaningful.
