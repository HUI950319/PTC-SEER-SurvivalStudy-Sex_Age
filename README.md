 # Sex Disparities in Papillary Thyroid Cancer Survival: Divergent Patterns of Relative and Absolute Effects Across the Age Spectrum

## Overview

This repository contains the complete analytical code and processed data for a comprehensive survival analysis study examining sex disparities in papillary thyroid cancer (PTC) outcomes. The study investigates how relative and absolute effects of sex on cancer-specific survival vary across different age groups using advanced statistical methods.

## Study Background

**Objective**: To re-evaluate the prognostic impact of sex on cancer-specific survival (CSS) in papillary thyroid cancer patients and determine whether age modifies the effect of sex on CSS.

**Data Source**: Surveillance, Epidemiology, and End Results (SEER) database (2004-2015)

**Sample Size**: 77,349 PTC patients
- Female: 61,197 (79.1%)
- Male: 16,152 (20.9%)

**Primary Outcome**: Cancer-specific survival (CSS)

## Repository Contents

### Files Description

| File | Size | Description |
|------|------|-------------|
| `Rcode.R` | 119KB (3,248 lines) | Complete statistical analysis code with comprehensive English comments |
| `Processed_SEER_data.rds` | 631KB | Preprocessed SEER cancer registry data ready for analysis |

### Data Description

**`Processed_SEER_data.rds`**
- **Format**: R data frame in RDS format
- **Patients**: 77,349 papillary thyroid cancer cases
- **Time Period**: Diagnosed 2004-2015, followed through 2020
- **Variables**: Demographics, tumor characteristics, treatment information, survival outcomes
- **Key Variables**:
  - `Sex`: Male/Female
  - `Age`: Age at diagnosis (continuous and categorical)
  - `time`: Follow-up time in months
  - `DSS`: Disease-specific survival indicator
  - `event`: Competing risks event indicator
  - Comprehensive set of confounding variables (race, ethnicity, marital status, income, tumor grade, size, extension, etc.)

## Statistical Methods

### Core Analytical Approaches

1. **Cox Proportional Hazards Models**
   - Relative effects (hazard ratios)
   - Multiple adjustment strategies

2. **Fine-Gray Competing Risks Models**
   - Absolute effects (cumulative incidence functions)
   - Accounts for competing causes of death

3. **Restricted Cubic Spline (RCS) Analysis**
   - Non-linear age-sex interaction modeling
   - Bootstrap confidence intervals (n=1,000)

4. **Formal Interaction Testing**
   - Multiplicative scale (Cox models)
   - Additive scale (RERI, AP, SI indices)

### Causal Inference Methods

- **Unadjusted Analysis**: Kaplan-Meier and Aalen-Johansen estimators
- **Direct Standardization**: Outcome regression adjustment


## Code Structure

### Main Analysis Sections

```r
# 1. Package Loading and Function Definitions
#    - Comprehensive survival analysis packages
#    - Custom functions for causal inference methods

# 2. Cox Proportional Hazards Analysis
#    - Survival curves comparison
#    - RCS analysis for non-linear age effects
#    - Subgroup analysis by age categories
#    - Formal interaction testing

# 3. Competing Risks Analysis  
#    - Cumulative incidence functions
#    - Fine-Gray regression models
#    - Sensitivity analysis for robustness

# 4. Visualization and Results Export
#    - Forest plots generation
#    - RCS effect plots
#    - Comprehensive figure creation
```

### Key Functions

- **`.get.obj()`**: Main wrapper for different adjustment methods
- **`.get.sub.dif.times()`**: Calculate survival differences across subgroups
- **`.get.UM.haz()`**: Univariate and multivariate hazard analysis
- **`rcs.cox1_boot()`**: RCS analysis with Cox models and bootstrap CI
- **`rcs.poi_boot()`**: RCS analysis with Poisson models for absolute rates
- **`.get.interTable()`**: Formal interaction analysis table generation

## Key Findings

### Patient Characteristics
- **Age Differences**: Males diagnosed at older age (median 51 vs 47 years)
- **Tumor Aggressiveness**: Males showed more aggressive clinicopathological features
  - Larger tumor size
  - Higher rates of extrathyroidal extension
  - More multifocal disease
  - Greater lymph node involvement
  - More advanced TNM staging

### Survival Analysis Results

#### Relative Effects (Hazard Ratios)
- **Overall Effect**: Males had 46% higher cancer-specific mortality risk (HR=1.46, 95% CI: 1.26-1.69)
- **Age-Dependent Patterns**: Relative effects decreased significantly with advancing age

#### Absolute Effects (Survival Differences)
- **10-year Survival Difference**: 0.46% (95% CI: 0.25%-0.67%)
- **15-year Survival Difference**: 0.77% (95% CI: 0.31%-1.22%)
- **Age-Dependent Patterns**: Absolute effects increased with advancing age

#### Interaction Analysis
- **Multiplicative Scale**: Significant interaction (p<0.001)
- **Additive Scale**: Evidence of positive interaction
- **Clinical Implication**: Divergent patterns of relative vs absolute effects across age spectrum

## Requirements

### R Version and Packages

```r
# Required R version: ≥ 4.0.0
# Key packages:
- survival          # Cox regression and Kaplan-Meier
- tidycmprsk        # Competing risks analysis
- adjustedCurves    # Causal inference methods
- riskRegression    # Advanced survival modeling
- tidyverse         # Data manipulation and visualization
- gtsummary         # Statistical tables
- forestploter      # Forest plot generation
- patchwork         # Plot composition
- interactionR      # Interaction analysis
```

### System Requirements
- **Memory**: ≥8GB RAM recommended for bootstrap procedures
- **Processors**: Multi-core CPU recommended for parallel bootstrap (R=1,000)
- **Storage**: ~2GB for complete analysis outputs and figures


## Citation

If you use this code or data in your research, please cite:

```
[Manuscript Under Review]
Sex Disparities in Papillary Thyroid Cancer Survival: Divergent Patterns of 
Relative and Absolute Effects Across the Age Spectrum
```


## Contact Information

For questions about the analysis methods or code implementation, please contact:
- **Primary Investigator**: [Hui Ouyang]
- **Institution**: [Xiangya Hospital of central South University]
- **Email**: [ouyanghui95@163.com]

## Acknowledgments

- **Data Source**: Surveillance, Epidemiology, and End Results (SEER) Program
- **Statistical Methods**: Advanced survival analysis and causal inference techniques
- **R Community**: Comprehensive packages enabling reproducible research


**Note**: This repository contains research code for academic purposes. The processed SEER data is de-identified and publicly available. All analyses follow appropriate ethical guidelines for cancer registry data usage. 
