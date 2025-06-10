# wdz.cut: An R Function for Threshold Effect Analysis

This repository contains an R function `wdz.cut()` for performing threshold effect analysis, also known as segmented regression. It automatically finds the optimal inflection point (cutpoint) in the relationship between a continuous variable and an outcome, and uses a likelihood ratio test to check the significance of the threshold effect.

This function was developed by Dongzhe Wu.

## Features

-   Automatically detects the optimal cutpoint in a dataset.
-   Supports `glm` (binomial, gaussian), `cph` (Cox), and `lrm` (logistic) model objects.
-   Performs a likelihood ratio test to assess the statistical significance of the threshold effect.
-   Outputs a clean, formatted table summarizing the results, perfect for academic papers.

## How to Use

No installation is required. You can load this function directly from GitHub into your R session using the `source()` command.

Run the following command in your R console:

```R
source("[https://raw.githubusercontent.com/2nixisst/WDZ/main/2025Threshold.Effect.Analysis.R](https://raw.githubusercontent.com/2nixisst/WDZ/main/2025Threshold.Effect.Analysis.R)")
