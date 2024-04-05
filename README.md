# Prob-toxicity-values

## Step I: Deriving Bayesian benchmark doses (BMDs) and WHO/IPCS BMDs
1) Bayesian BMD
   - Input files:
     - All data files are in in repository: BBMD input
   - **Script: Benchmark dose.Rmd**
   - Output files:
     - bmd samples CAS by row.csv
     - BBMD quantile.csv
     - hed samples CAS by row.csv
2) WHO/IPCS BMD
   - **Script: WHO default benchmark dose.Rmd**
   - Output files:
     - WHO bmd samples CAS by row.csv
     - WHO BMD quantile.csv
     - WHO hed samples CAS by row.csv
3) Generate Bayesian model averaging benchmark dose plots
   - Input files:
     - All data files are in in repository: BBMD input
   - **Script: BBMD dose-response data.Rmd**
   - Output files:
     - All plots in repository: dose-response plots (Supplementary Materials File S3)

## Step II: Human toxicokinetic (TK) variability
  ### Generate steady-state concentrations (C<sub>ss</sub>) and TK parameters for urine excretion from the population-based _httk_ R package.
  ### Reference: Ring et al. 2017 (PMID: 28628784)
- Input file: 19 overlapping chemical list.csv
- **Script: Human TK variability.R**
- Output files:
  - Css samples CAS by row.csv
  - GFRxFub samples CAS by row.csv

## Step III: Human toxicodynamic (TD) variability
  ### Use the population-based data from lymphoblastoid cell line assays.
  ### Reference: Ford et al. 2022 (PMID: 36006120)
- Input files:
  - 19 overlapping chemical list.csv
  - tdvf_ford et al 2022.csv
- **Script: Human TD variability.R**
- Output file: td gsd distribution CAS by row.csv

## Step IV: Deriving toxicity values 
  ### HD<sub>M</sub><sup>I</sup> and biomonitoring equivalents (BE<sub>M</sub><sup>I</sup>) in blood and urine
1) Chemical-specific toxicity values
   - Input files:
     - hed samples CAS by row.csv
     - Css samples CAS by row.csv
     - td gsd distribution CAS by row.csv
     - GFRxFub samples CAS by row.csv
     - regulatory reference dose.csv
   - **Script: Chemical-specific prob toxicity values.R**
   - Output files:
     - BEMI blood.csv
     - BEMI urine.csv
     - HDMI.csv
     - BEMI blood quantile.csv
     - BEMI urine quantile.csv
     - HDMI quantile.csv
2) WHO/IPCS approximate HD<sub>M</sub><sup>I</sup>
   - Input files:
     - WHO hed samples CAS by row.csv
     - regulatory reference dose.csv
   - **Script: WHO prob toxicity values.R**
   - Output files:
     - WHO approximate HDMI.csv
     - WHO approximate HDMI quantile.csv
## Step V: Analyses
1) Compare BMD uncertainty distributions between Bayesian model averaging approach and WHO/IPCS approximate approach
   - **Compare BMD uncertainty distributions, differences of the median BMD, degree of uncertainty**
     - Input files:
       - regulatory reference dose.csv
       - bmd samples CAS by row.csv
       - WHO bmd samples CAS by row.csv
       - regulatory reference dose.csv
       - BBMD quantile.csv
       - WHO BMD quantile.csv
     - **Script: BMD plot.R**
     - Output files:
       - BMD plot.pdf **(Manuscript Figure 2)**
2) Compare human variability uncertainty distributions (AF<sub>intra</sub>) between chemical-specific approach and WHO/IPCS approximate approach
   - **Compare AF<sub>intra</sub> uncertainty distributions, differences of the median AF<sub>intra</sub>, degree of uncertainty**
     - Input files:
       - regulatory reference dose.csv
       - TKTDVF01.csv
       - WHO approximate intraspecies factor.csv
       - TKTDVF01 quantile.csv
       - WHO approximate intraspecies factor quantile.csv
     - **Script: AFintra plot.R**
     - Output files:
       - Intraspecies uncertainty plot.pdf **(Manuscript Figure 3)**
3) Compare HD<sub>M</sub><sup>I</sup> uncertainty distributions between chemical-specific approach and WHO/IPCS approximate approach
   - **Compare HD<sub>M</sub><sup>I</sup> uncertainty distributions, differences of the median HD<sub>M</sub><sup>I</sup>, degree of uncertainty**
     - Input files:
       - regulatory reference dose.csv
       - HDMI.csv
       - WHO approximate HDMI.csv
       - HDMI quantile.csv
       - WHO approximate HDMI quantile.csv
     - **Script: HDMI plot.R**
     - Output files:
       - HDMI plot.pdf **(Manuscript Figure 4)**
4) Compare chemical-specific toxicity values (HD<sub>M</sub><sup>I</sup> and BE<sub>M</sub><sup>I</sup> in blood and urine) with exposure data
   - **Exposure data from the U.S. EPA ExpoCast prediction, NHANES, Canadian Health Measure Surveys**
     - Input files:
       - bimonitoring data.csv
       - BEMI blood.csv
       - BEMI urine.csv
       - regulatory reference dose.csv
       - HDMI.csv
       - Expocast distribution by row.csv
     - **Script: Toxicity values with exposure data.R**
     - Output files:
       - Chemical-specific toxicity values with exposure data boxplot.pdf **(Supplementary Fig. S1)**
5) Estimate Margin of Exposure (MOE)
   - Input files:
      - expocast upper bound.csv
      - WHO approximate HDMI quantile.csv
      - HDMI quantile.csv
      - bimonitoring data.csv
      - BEMI blood quantile.csv
      - BEMI urine quantile.csv
   - **Script: MOE.R**
   - Output file: MOE plot (Expocast and biomonitoring).pdf **(Manuscript Figure 5)**
6) Compare biomonitoring data and C<sub>ss</sub>/U<sub>ss</sub>
   - Intput files:
     - expocast median and upper bound.csv
     - Css samples CAS by row.csv
     - GFRxFub samples CAS by row.csv
     - bimonitoring data.csv
   - **Script: Expocast to Css and Uss.R**
   - Output file: Expocast Css/Uss vs biomonitoring data plot.pdf **(Supplementary Fig. S2)**
