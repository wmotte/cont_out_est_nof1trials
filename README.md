# Continuous Outcome Estimation in N-of-1 Trials

This repository contains the code and data for the research paper:

**Defelippe V, Bartoš F, Wagenmakers EJ, Braun KPJ, Jansen FE, Otte WM.** *Continuous outcome estimation in N-of-1 trials for accelerated decision-making.* **Epilepsia.** 2026 Feb 5. doi: [10.1002/epi.70119](https://doi.org/10.1002/epi.70119)

## Overview

N-of-1 trials are single-patient crossover studies that aim to determine treatment efficacy for individual patients through multiple alternating periods of active treatment and placebo/control. This study demonstrates how continuous Bayesian outcome estimation can accelerate clinical decision-making in N-of-1 trials by allowing for earlier trial termination when sufficient evidence is obtained.

### Key Findings

- **Trial duration reduction**: Original trial duration could be reduced by 9.5% to 35% in analyzed examples
- **Dual analytical approach**: Combining Bayesian hypothesis testing (Bayes factors) with minimal clinically important difference (MCID) thresholds provides complementary evidence for decision-making
- **Adaptive methodology**: Continuous updating of treatment effect estimates as new data becomes available enables more efficient trials while maintaining statistical rigor

## Background

Traditional N-of-1 trials evaluate treatment efficacy only after predetermined cycles are completed. However, in chronic neurological conditions like epilepsy where recurrent events are measured continuously, treatment responders (or non-responders) can potentially be identified earlier. This research explores whether continuous Bayesian analysis can safely reduce trial duration and patient burden.

## Methodology

The study employs:

1. **Bayesian generalized linear models** with Poisson distributions for count data (seizures)
2. **Bayes factors (BF)** to quantify relative evidence for treatment effect vs. no effect
3. **Posterior probability distributions** for minimal clinically important differences (MCID)
4. **Sequential updating** of models with each new data point

## Repository Structure

### Main Analysis Scripts

The numbered R scripts represent the sequential analysis pipeline:

#### Data Visualization
- `00__plot_data__nicotine.R` - Visualize nicotine patch trial data (Willoughby et al.)
- `01__plot_data__headache.R` - Visualize headache treatment trial data (Williams et al.)

#### Simulation Study
- `02__simulate_refractory.R` - Simulate refractory epilepsy N-of-1 trial data
- `03__model_refractory.R` - Bayesian modeling of simulated data
- `04__plot_refractory.R` - Visualization of simulation results

#### Continuous Updating Analyses
- `05__updating__nicotine.R` - Sequential analysis of nicotine trial
- `06__updating__hippocampal_stim.R` - Sequential analysis of hippocampal stimulation trial
- `07__updating__headache.R` - Sequential analysis of headache trial
- `08__updating__simulate_refractory.R` - Sequential analysis of simulated trial

#### Bayes Factor Analyses
- `09__BF__simulate_refractory.R` - Bayes factors for simulated data
- `09a__BF__simulate_refractory_prior.R` - Prior sensitivity for simulated data
- `10__BF__headache.R` - Bayes factors for headache trial
- `10a__BF__headache_prior.R` - Prior sensitivity for headache trial
- `11__BF__nicotine.R` - Bayes factors for nicotine trial
- `11a__BF__nicotine_prior.R` - Prior sensitivity for nicotine trial
- `12__BF__hippocampal_stim.R` - Bayes factors for hippocampal stimulation
- `12a__BF__hippocampal_stim_prior.R` - Prior sensitivity for hippocampal stimulation

#### Sensitivity Analyses
- `13__sensitivity_MCID_thresholds.R` - Sensitivity to different MCID thresholds
- `14__sensitivity_prior_specifications.R` - Sensitivity to prior distributions
- `15__sensitivity_treatment_efficacy.R` - Sensitivity to treatment effect sizes

### Data

The `datasets/` folder contains reconstructed data from published N-of-1 trials:

- **Nicotine trial** (Willoughby et al., 2003) - Transdermal nicotine patches for refractory epilepsy
- **Headache trial** (Williams et al., 1986) - Somatostatin analog (SMS) for prolactinoma-related headaches
- **Hippocampal stimulation** (Tellez-Zenteno et al., 2006) - Electrical stimulation for mesial temporal lobe epilepsy
- **Hankey dataset** - Additional clinical data

Data available in multiple formats: `.csv`, `.xlsx`, and visualizations (`.png`)

### Output

Each analysis script generates outputs in corresponding directories:
- `out.00.plot.nicotine/` through `out.15.sensitivity.efficacy/`

These contain:
- Statistical model results
- Figures and visualizations
- Posterior probability distributions
- Bayes factor calculations

## Installation

### Requirements

- R (version ≥ 4.0 recommended)
- Required R packages:
  - `brms` - Bayesian regression models
  - `bridgesampling` - Bayes factor calculation via bridge sampling
  - `ggplot2` - Visualization
  - `readr` - Data import
  - Additional packages as specified in individual scripts

### Setup

```bash
# Clone the repository
git clone https://github.com/wmotte/cont_out_est_nof1trials.git
cd cont_out_est_nof1trials

# Install R dependencies (example for key packages)
# In R:
install.packages(c("brms", "ggplot2", "readr", "bridgesampling"))
```

## Usage

### Running the Full Analysis Pipeline

Scripts are numbered sequentially and can be run in order:

```bash
# Visualize original trial data
Rscript 00__plot_data__nicotine.R
Rscript 01__plot_data__headache.R

# Run simulation study
Rscript 02__simulate_refractory.R
Rscript 03__model_refractory.R
Rscript 04__plot_refractory.R

# Perform continuous updating analyses
Rscript 05__updating__nicotine.R
Rscript 06__updating__hippocampal_stim.R
Rscript 07__updating__headache.R
Rscript 08__updating__simulate_refractory.R

# Calculate Bayes factors
Rscript 09__BF__simulate_refractory.R
# ... etc.
```

### Running Individual Analyses

Each script can be executed independently:

```bash
# Example: Run nicotine trial analysis
Rscript 05__updating__nicotine.R
```

## Key Concepts

### Bayes Factor (BF)

The Bayes factor quantifies the relative evidence for the alternative hypothesis (H1: treatment effect exists) vs. the null hypothesis (H0: no treatment effect):

- **BF > 1**: Evidence favoring treatment effect
- **BF < 1**: Evidence favoring no effect
- **Log BF > 3.4**: Strong evidence for treatment effect
- **Log BF > 4.6**: Extreme evidence for treatment effect

### Minimal Clinically Important Difference (MCID)

The MCID represents the smallest treatment effect that is clinically meaningful to patients. In epilepsy trials:
- **50% seizure reduction**: Traditional MCID
- **30% seizure reduction**: Accepted for refractory epilepsy by regulatory agencies

The Bayesian approach calculates the posterior probability that the treatment effect exceeds the MCID threshold.

### Continuous Updating

Rather than waiting until trial completion, the Bayesian model is updated each time new data becomes available (e.g., daily), allowing:
- Real-time monitoring of treatment effects
- Early stopping for efficacy when strong evidence emerges
- Early stopping for futility when treatment clearly doesn't work
- Reduced patient burden and exposure to ineffective treatments

## Study Examples

### 1. Simulated Refractory Epilepsy Trial
- Patient with 5 seizures/day baseline
- 5 treatment cycles (placebo vs. verum)
- Net seizure reduction: 50%
- **Result**: Trial could stop at Day 152 instead of Day 168 (9.5% reduction)

### 2. Nicotine Patch Trial (Willoughby et al., 2003)
- 33-year-old woman with refractory frontal lobe epilepsy
- Transdermal nicotine vs. placebo patches
- Complete seizure suppression on nicotine
- **Result**: Trial could stop on Day 52 instead of Day 84 (38% reduction, ~1 month earlier)

### 3. Headache Trial (Williams et al., 1986)
- 38-year-old woman with prolactinoma-related headaches
- Somatostatin analog (SMS) vs. placebo
- 4-day trial with multiple daily measurements
- **Result**: Evidence for treatment effect emerged on Day 3-4

### 4. Hippocampal Stimulation (Tellez-Zenteno et al., 2006)
- 4 adults with refractory mesial temporal lobe epilepsy
- Electrical stimulation on vs. off
- Long-term follow-up data (48 months)
- **Result**: Demonstrates how evidence evolves with changing interventions

## Citation

If you use this code or methodology in your research, please cite:

```
Defelippe V, Bartoš F, Wagenmakers EJ, Braun KPJ, Jansen FE, Otte WM.
Continuous outcome estimation in N-of-1 trials for accelerated decision-making.
Epilepsia. 2026 Feb 5. doi: 10.1002/epi.70119.
```

## Authors

- **Victoria Defelippe** - Department of Pediatric Neurology, UMC Utrecht Brain Center, University Medical Center Utrecht
- **František Bartoš** - Department of Psychology, University of Amsterdam
- **Eric-Jan Wagenmakers** - Department of Psychology, University of Amsterdam
- **Kees P. J. Braun** - Department of Pediatric Neurology, UMC Utrecht Brain Center
- **Floor E. Jansen** - Department of Pediatric Neurology, UMC Utrecht Brain Center
- **Willem M. Otte** - Department of Pediatric Neurology, UMC Utrecht Brain Center (Corresponding author: w.m.otte@umcutrecht.nl)

## Funding

This project was funded by:
- MING Foundation
- Health-Holland (Grant Number: LSHM21053-SGF)

## License

This project is licensed under the terms specified in the LICENSE file.

## Related Resources

- **JASP Software**: User-friendly Bayesian statistics software - [jasp-stats.org](https://jasp-stats.org/)
- **brms Package**: Bayesian regression models in R - [paul-buerkner.github.io/brms/](https://paul-buerkner.github.io/brms/)
- **Original Paper**: [https://doi.org/10.1002/epi.70119](https://doi.org/10.1002/epi.70119)

## Clinical Implications

This methodology can help:
1. **Reduce patient burden** by shortening trial duration when evidence is clear
2. **Minimize placebo exposure** for responders
3. **Identify non-responders earlier** to avoid prolonged ineffective treatment
4. **Improve trial efficiency** while maintaining statistical rigor
5. **Enable personalized medicine** through adaptive trial designs

## Future Directions

Potential extensions of this work include:
- Incorporation of seizure frequency fluctuations and cyclical patterns
- Adjustment for carryover effects in specific interventions
- Integration with digital health platforms for real-time data collection
- Application to other chronic neurological conditions
- Development of patient-specific MCID thresholds through shared decision-making

## Contact

For questions or collaboration inquiries:
- Willem M. Otte: w.m.otte@umcutrecht.nl
- Victoria Defelippe: v.m.defelippe-2@umcutrecht.nl

## Acknowledgements

We thank the original investigators who published the N-of-1 trial data that enabled these re-analyses:
- Willoughby JO, Pope KJ, Eaton V (2003) - Nicotine trial
- Williams G, Ball J, Bloom S, Joplin GF (1986) - Headache trial
- Tellez-Zenteno JF, McLachlan RS, Parrent A, Kubu CS, Wiebe S (2006) - Hippocampal stimulation trial

---

**Keywords**: N-of-1 trials, Bayesian statistics, epilepsy, adaptive clinical trials, precision medicine, seizures, personalized treatment, statistical modeling
