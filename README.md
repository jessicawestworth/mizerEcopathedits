# Balanced Harvesting in the Celtic Sea: Using Mizer Modelling to Determine Optimal Fishing Regimes

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This repository contains the dynamic multi-species size-spectrum modeling framework used to evaluate the ecological, conservation, and food provisioning outcomes of **Species-and-Size-level Balanced Harvesting (ssBH)** compared to status quo fishing practices within the Celtic Sea ecosystem.

The framework for setting up and calibrating all multi-species size spectrum models used the [mizer](https://sizespectrum.org/mizer/) and [mizerEcopath](https://gustavdelius.github.io/mizerEcopath/) packages.

## Project Overview
This repository includes the code for data processing, model creation, model projection, results evalution and sensivity analysis for the dissertation report "BALANCED HARVESTING IN THE CELTIC SEA: USING MIZER MODELLING TO DETERMINE OPTIMAL FISHING REGIMES". Traditional fisheries management often relies on Total Allowable Catches (TACs) derived from single-species MSY models coupled with strict size-selective limits. While designed to optimize yield-per-recruit, this approach heavily biases mortality toward large-bodied, high-trophic-level predators, destabilizing ecosystem structure and removing highly fecund individuals that buffer recruitment variability.

**Balanced Harvesting (BH)** offers an alternative by distributing fishing mortality across the widest possible range of species and sizes in proportion to their natural productivity or production. This repository focuses specifically on **species-and-size-level Balanced Harvesting (ssBH)**—where fishing mortality scales directly with biological *production* (biomass produced per unit time), acting as a density-dependent safeguard against stock collapse.

### Research Questions Addressed:
1. What are the ecological and provisioning effects of different fishing intensities within an entirely ssBH regime?
2. How does a hybrid or entirely ssBH regime perform compared to the status quo regime?
3. How does a transition to a hybrid or entirely ssBH regime redistribute yields across the food web?
4. What are the underlying ecological processes driving these shifts when transitioning from the status quo?

## Methodology & Model Workflow
The framework modifies the standard `mizer` implementation (Scott et al., 2014) to enhance biological realism. The modeling process consists of:

1. **Status Quo Baseline (2012–2024):** Construction of a non-interacting, steady-state allometric model calibrated using extensive empirical data from the Celtic Sea.
2. **Dynamic Calibration:** Interspecific interactions (predation, growth, biomass flux) are enabled, and the model is simulated forward to a steady state.
3. **Scenario Projections:** Execution of hundreds of multi-species simulations exploring fully ssBH and "hybrid" fishing regimes (ssBH implemented alongside status quo practices) across a gradient of fishing intensities.
4. **Robustness Testing:** Due to exhaustive data utilization preventing independent dataset validation, model stability and parameter uncertainty are evaluated using **Morris and Regional Sensitivity Analyses**.

## Data Sources
The model is highly parameterized and tuned using regional empirical data from the Celtic Sea ecosystem:
* **Biomass & Landings:** Reference Ecopath model (Lauria, 2012; 2016) and fisheries landings data (European Commission, 2024).
* **Surveys & Life History:** Length distributions, maturity schedules, and age-length data (ICES 2025b; Silva et al., 2024).
* **Trophic Interactions:** Diet-matrix (Lauria, 2016) and stomach contents empirical data (Delius, 2026).

## Files and Accompanying Report Sections
The R code files and the accompanying report section are described below:

**Data Processing**:
* 7.1.2.1 ICES Rectangle Area Mapping: inst/ICES Statistical Rectangles to ICES Areas.qmd
* 7.1.2.2 Processing Length Distribution: inst/Processing DATRAS Fishing Survey.qmd
* 7.2.2.3 Processing Stratified Age at Length data: 
* 7.1.2.4 Processing Boarfish Data:inst/Length Stratify Boarfish Data.qmd
* 7.1.2.5 Calculating w_mat and w_mat25
* 7.1.3. Processing Landings: Processing Landings.qmd
* 7.1.3.1 Discard Rates and Total Yield Calculation: inst/Calculate Biomass Discard Rate.qmd
* 7.1.3.2 Size-specific Mortaltiy Distributions:  inst/Processing Landings.qmd
* 7.1.3.3 Gear Categorization:  inst/Processing Landings.qmd
* 7.1.5 Calculating the interaction matrix: R/reduceEcopathDiet.R
* 7.1.6 Reformulating the diet matrix: inst/Processing Ecopath Diet Matrix.qmd

**Data Matching**:
* Age Matching: R/plotAge.R & simulateAge.R & AgeDensity.R
* Length Distribution Matching: R/plot_catch.R
* Visual Matching Shiny App: R/ecopath_tune.R

**Model Celtic Sea Creation Workflow**:
* 3.1-3.8: vignettes/Celtic Sea Status Quo Model.qmd

**Balanced Harvest**:
* 3.9 Calculating ssBH: R/flux.R
* 3.10.2 Gradual Implementation: R/make_blended_ssBH_FMort.R

**Model Projection and Results**:
* 3.10 & 3.11 Evaluation Metrics: inst/create_alpha_c_sims.R
* 3.12 Sensitivity Analysis: inst/Sensitivity Analysis.R

## Model Creation Data Files
1. diets.rda
2. stomach_data_fit.rda
3. survey_length_distribution.rda
4. life_history_fishbase.rda
5. fishing_deaths.rda
6. cs_age_size.rda

## Non-interacting Celtic Sea Status Quo model
inst/final7.rds

## Installation
You can install the required version of mizerEcopath and mizer from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("sizespectrum/mizer", ref="6353b27")
remotes::install_github("jessicawestworth/mizerEcopathedits"", ref="dissertation")
```
