# Prob-toxicity-values

## Step I: Deriving Bayesian benchmark doses (BMDs) and WHO/IPCS BMDs
1) Bayesian BMD
   - Input files:
     - All data files are in in repository: BBMD input
   - Script: Benchmark dose.Rmd
   - Output files:
     - bmd samples CAS by row.csv
     - BBMD quantile.csv
     - hed samples CAS by row.csv
2) WHO/IPCS BMD
   - Script: WHO default benchmark dose.Rmd
   - Output files:
     - WHO bmd samples CAS by row.csv
     - WHO BMD quantile.csv
     - WHO hed samples CAS by row.csv
3) Generate Bayesian model averaging benchmark dose plots
   - Input files:
     - All data files are in in repository: BBMD input
     - Script: BBMD dose-response data.Rmd
   - Output files:
     - All plots in repository: dose-response plots (Supplementary Materials File S3)

## Step II: Human toxicokinetic (TK) variability
  ### Generate steady-state concentrations (Css) and TK parameters for urine excretion from httk R package.
  ### Reference: Ring et al. 2017 (PMID: 28628784)
- Input file: 19 overlapping chemical list.csv
- Script: Human TK variability.R
- Output file:
  - Css samples CAS by row.csv
  - GFRxFub samples CAS by row.csv

## Step III: Human toxicodynamic (TD) variability
  ### Use the population-based data from lymphoblastoid cell line assays.
  ### Reference: Ford et al. 2022 (PMID: 36006120)
- Input files:
  - 19 overlapping chemical list.csv
  - tdvf_ford et al 2022.csv
- Script: Human TD variability.R
- Output file: td gsd distribution CAS by row.csv

## Step IV: Deriving toxicity values (HDMI and biomonitoring equivalents in blood and urine)
1) Chemical-specific toxicity values
   - Input files:
     - hed samples CAS by row.csv
     - Css samples CAS by row.csv
     - td gsd distribution CAS by row.csv
     - GFRxFub samples CAS by row.csv
     - regulatory reference dose.csv
   - Script: Chemical-specific prob toxicity values.R
   - Output files:
     - BEMI blood.csv
     - BEMI urine.csv
     - HDMI.csv
     - BEMI blood quantile.csv
     - BEMI urine quantile.csv
     - HDMI quantile.csv
2) WHO/IPCS approximate toxicity values
   - Input files:
     - WHO hed samples CAS by row.csv
     - regulatory reference dose.csv
   - Script: WHO prob toxicity values.R
   - Output files:
     - WHO approximate HDMI.csv
     - WHO approximate HDMI quantile.csv
## Analyses
1) Compare chemical-specific approach and WHO/IPCS approach in HDMI
   - Input files:
     - regulatory reference dose.csv
     - bmd samples CAS by row.csv
     - WHO bmd samples CAS by row.csv
     - HDMI.csv
     - WHO approximate HDMI.csv
   - Script:
   - Output files:
     - the most sensitive endpoint BBMD HDMI plot.pdf (Manuscript Figure 2)
     - all endpoints normalized BBMD HDMI ggridges plot.pdf (Supplementary Fig. S1)
2) Compare differences in HDMI, BMDs, and human variability values between two approaches
   - Input files:
     - regulatory reference dose.csv
     - BBMD quantile.csv
     - WHO approximate HDMI quantile.csv
     - HDMI quantile.csv
     - BMD quantile.csv
     - WHO BMD quantile.csv
   - Script: Differences in HDMI BMD UFh.R
   - Output file: Differences in HDMI BMD UFh bar plot.pdf (Manuscript Figure 3)
3) Degree of uncertainty between two approaches
   - Input files:
      - BBMD quantile.csv
      - WHO BMD quantile.csv
      - WHO approximate HDMI quantile.csv
      - HDMI quantile.csv
   - Script: Degree of uncertainty.R
   - Output file: Degree of uncertainty uncertainty violin plot.pdf (Manuscript Figure 4)
4) Compare chemical-specific toxicity values with exposure data
   - Input files:
      - bimonitoring data.csv
      - BEMI blood.csv
      - BEMI urine.csv
      - regulatory reference dose.csv
      - HDMI.csv
      - Expocast distribution by row.csv
   - Script: Toxicity values with exposure data ggridges.R
   - Output files:
      - the most sensitive HDMI and BEMI in blood and urine ggridges plot.pdf (Manuscript Figure 5)
      - all endpoints HDMI and BEMI in blood and urine ggridges plot.pdf (Supplementary Fig. S2)
5) Estimate Margin of Exposure (MOE)
   - Input files:
      - expocast upper bound.csv
      - WHO approximate HDMI quantile.csv
      - HDMI quantile.csv
      - bimonitoring data.csv
      - BEMI blood quantile.csv
   - Script: Margin of exposure.R
   - Output file: MOE plot.pdf (Manuscript Figure 6)
