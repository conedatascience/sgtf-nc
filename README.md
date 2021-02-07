# Rapid Impact Analysis of B 1.1.7 Variant on the Spread of SARS-CoV-2 in North Carolina

This analysis examines the potential impact of the B1.1.7 variant on the tranmission of SARS-CoV-2 in North Carolina.

## Running Analysis

Ensure that the following package is installed (non-CRAN):

```r
install.packages("remotes")
remotes::install_github("conedatascience/nccovid")
```

### Data

The data reviewed to run this analysis are available in the **data-raw** directory.

### Analysis

In order to reproduce the analysis, run:

* **R/generate_variant.R** to fit the UK data for the new variant

* **R/simulate_outcomes.R** to run the stochastic, discrete ODEs 

