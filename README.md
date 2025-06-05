This repository documents the code associated with the manuscript "Exposure of western United States bird communities to predicted high severity fire". It is written primarily 

# System Requirements
### OS Requirements
This software has only been tested on macOS: Sonoma 14.3 but all package depndencies are available for both operating systems and scripts should be fully executable for other operating systems.

### R Dependencies
This software is primarily implemented in R (tested on v4.3.3) and package dependencies are primarily statistical and spatial data packages:

```
    betareg
    broom
    cowplot
    dplyr
    emmeans
    furrr
    ggplot2
    JuliaCall
    patchwork
    purrr
    rnaturalearth
    sf
    terra
    tidyr
    tidyterra
```

# Installation
#### Github
Download software from github:
```
git clone https://github.com/karinorman/ebirdFire
```

#### Install
Install dependencies using R package `devtools`, this should only take a few minutes if install is done before the data step below.
```
devtools::install()
```

#### Data
Download the zip file of data products from Zenodo [10.5281/zenodo.15414728](https://zenodo.org/records/15414728) and place them in folder `data/`.

# Analysis
Analysis scripts can be found in `analysis/`. Scripts `00` and `01` create data pieces already included in the data folder. 
Subsequent scripts step through the analysis in order of their prefix number. Interim data products to execute any script in the analysis folder are all included in the data download above. 
As all inputs and outputs of analysis scripts are included in the data download, they can be easily used to check reproducibility.
Individual scripts should only take a few minutes to run on a standard desktop, with the exception of `00` and `01` which may take a couple hours.
