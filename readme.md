# PhenologyBigBud 🌿

## Overview
This repository has data and code used in analyses for a manuscript entitled: "Broad-scale geographic structure of phenological sensitivity emerges from cue–forcing interactions" that is  
a collaborative effort and that leverages PhenoVision leaf annotatons generated from a recent paper "Grady, E. L., E. G. Denny, C. E. Seltzer, J. Deck4, D. Li, R. Dinnage* and R. P. Guralnick* (co-senior authors).  
[Accepted]. Turning a new leaf: PhenoVision provides leaf phenology data at the global scale. Applications in Plant Sciences. We store intermediate data outputs here and point  
those looking for the full leaf annotation data to a Zenodo repository where those data reside.  

## Key Questions
This work addresses three key questions:
1) Where is temperature a dominant driver of phenological timing, and where is its influence constrained by other cues? 
2) How does seasonal climate context, especially relative temperature seasonality, interact with local cue environments to shape sensitivity to warming?
3) Do spatial gradients that organize phenological timing across geography reflect the same sensitivities expressed through interannual climate variability?

## Data Sources
Our data sources include:
- Leaf annotations (filtered to breaking leaf buds)
- 2017-2025 Daily climate at the global scale from the ERA5-Land renanalysis dataset (9km resolution)
- Global 100m Terrestrial Human Footprint (HFP-100) v1.2
- A relative seasonality raster layer that has removed latitudinal and elevational signal derived from Worldclim Bio4 at 1km resolution
- Ensemble Digital Terrain Model (EDTM) of the world at 250m resolution

## Methods
- This repository contains code and data for a global analysis of leaf budbreak phenology using AI-generated annotations from iNaturalist observations (Phenovision; Grady et al., 2026). We model departures from Hopkins' Bioclimatic Law — a geographic baseline for phenological timing based on latitude and elevation — to isolate how temperature forcing, photoperiod, temperature seasonality, and human footprint shape spring leafing across the Northern Hemisphere (21–68°N, 2017–2025). Using mixed-effects models fit to hundreds of thousands of observation-level records, we test how cue environments constrain or amplify phenological responses to warming, and project these relationships globally to map temperature sensitivity and photoperiod gating strength across temperate regions.

## Repository Structure
PhenologyBigBud/
├── data/
├── scripts/
└── outputs/

## Getting Started
- All scripts used are in R and follow a logical sequence to run in order to derive analysis outputs
- The key script that generates the main datafile used, annotating photoperiod, GDD measurements and relative temperature seasonality is the first script
- Dependencies and needed inputs set in the main first script

## Contributors
- Rob Guralnick, Lindsay Campell, Tzeresa Crimmins, Daijiang Li and Michael Beltiz

## Citation
- Will post a citation to manuscript when we preprint this.

## License
- e.g., MIT, CC-BY
