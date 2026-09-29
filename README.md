# SPD: standardized prognostic dependence

R code and data for the manuscript "A Standardized Diagnostic for Prognosis-Driven Treatment Deviation in Per-Protocol Pharmacoepidemiologic Studies", submitted to *Pharmacoepidemiology and Drug Safety*.

SPD is the cause-specific Cox coefficient for treatment deviation per standard deviation of a baseline prognostic score. All data in this repository are simulated; there are no patient data.

## Requirements

R with the packages survival, data.table, future and future.apply (tested with R 4.3.3 and survival 3.5-8):

```r
install.packages(c("survival", "data.table", "future", "future.apply"))
```

`run_all_checks.R` needs only base R and survival.

## Quick check

```bash
Rscript run_all_checks.R
```

This takes a few seconds. It recomputes the reported values that come from stored data and prints each one next to the value in the manuscript:

- Table 3 and eTable S5, from the stored replicate-level results of the diagnostic-complementarity simulation
- Table 4 and eTables S6–S7, from the stored evaluation cohort, predicted probabilities and 2,000 bootstrap replicates
- eTable S8 and eFigure S5, from the stored event-ordering cohort
- Figures 1 and 2

## Simulations

```bash
Rscript run_all_simulations.R
```

This takes about 20 minutes on two cores (the number of parallel workers is set with `N_CORES`, default 3). It runs:

- the estimator-performance simulation for Table 2 and eTables S1–S3 (`scripts/07`, three settings)
- the model-based versus robust standard-error comparison for eTable S4 (`scripts/09`)
- eFigures S1–S4 from these results (`scripts/08`)
- the diagnostic-complementarity simulation for Table 3 and eTable S5 (`scripts/11`)

Results and figures are written to `output/`, which is not tracked by git.

## Where each result comes from

| Result | Script | Input |
|---|---|---|
| Table 2, eTables S1–S3 | 07 | simulated |
| eTable S4 | 09 | simulated |
| eFigures S1–S4 | 08 | output of 07 and 09 |
| Table 3, eTable S5 | 16 (stored replicates); 11 (new run) | `results/complementarity` |
| Figure 1 | 14 | none |
| Figure 2 | 14 | `results/complementarity` |
| Table 4, eTables S6–S7 | 12 | `results/worked_example` |
| eTable S8, eFigure S5 | 13 | `results/additional_checks` |

## Agreement with the manuscript

- Scripts 12, 13 and 16 reproduce the stored results to within 1e-7, and therefore every value printed in the manuscript.
- Script 09 reproduces eTable S4 exactly.
- Scripts 07 and 11 regenerate the simulated data. The random draws differ from those behind the reported tables, so the results agree with them to within Monte Carlo error but not digit for digit. The reported eTables S1–S4 are stored in `results/original_tables`.
- The predicted probabilities in `results/worked_example` come from the main-effects logistic regression and histogram gradient-boosting models described in Supplementary Methods S1.
- In the stored worked-example files, the known-DGM reference-risk score is labelled `Oracle reference-risk logit` (`oracle_probability`).
- The worked example and the complementarity simulation use Breslow ties; the estimator simulation uses Efron ties. Simulated event times are continuous, so the choice does not change the estimates.

## Earlier version

The original submission used a worked example based on the MIMIC-IV Clinical Database Demo. That code is kept in the repository history under the tag `v1-original-submission`.

## License

MIT; see `LICENSE`.
