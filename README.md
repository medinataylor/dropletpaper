## Scripts and input files for: Medina et al., *Inactivation of surrogate respiratory bacteria is driven by osmotic pressure change in drying saliva microdroplets*

### ResAM Model

Author: Beiping Luo  
Contact: beiping.luo@env.ethz.ch  
OrciD: https://orcid.org/0000-0003-1629-881X  

- `ReSAM_50.f` — ResAM model written in Fortran 77

### Input files for artificial saliva at 50%

- `droplet.dat` - droplet parameters
- `image_size.dat` - droplet size from filming
- `parameters.dat` - model parameters
- `ML.dat` - chemical species concentration

### Output ResAM data

- `resam_output_dropletpaper.zip`
  - `as30_ecoli` — artificial saliva droplet at 30% RH
  - `as40_3010` — artificial saliva droplet at 40% RH
  - `as50_0.1x_2909` — 0.1x artificial saliva droplet at 50% RH
  - `as50_ecoli` — artificial saliva droplet at 50% RH
  - `as60` — artificial saliva droplet at 60% RH
  - `as70_default` — artificial saliva droplet at 70% RH
  - `as80_2910` — artificial saliva droplet at 80% RH
  - `as5010x_2609` — 10x artificial saliva droplet at 50% RH
  - `pbs40_0303` — PBS droplet at 40% RH
  - `pbs50_1212` — PBS droplet at 60% RH
  - `rh_diff_3_new` — fast change RH experiment
  - `rh_slow_0702` — slow change RH experiment

### R Script to process ResAM output data to obtain osmotic pressure and bacteria predictions

Author: Taylor Medina  
Contact: taylor.medina@epfl.ch  
OrciD: https://orcid.org/0000-0002-3272-0139  

- `osmoticpressure.R` — written in R
