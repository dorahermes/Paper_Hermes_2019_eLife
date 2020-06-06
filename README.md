[![Contributors](https://img.shields.io/github/contributors/badges/shields)(https://github.com/badges/shields/graphs/contributors)
[![@Dora](http://img.shields.io/twitter/follow/dora_hermes.svg?style=social)](https://twitter.com/dora_hermes?lang=en)

# Paper Hermes 2019 eLife
This repository contains functions and scripts that were used in developing the OV model and generate the results and figures from [Hermes et al, 2019](https://elifesciences.org/articles/47035).
If you use this code as a part of any publications, please please cite this work.

Code (c) Dora Hermes, Kendrick Kay and Jonathan Winawer

Please direct any comments about the code via the issues page associated with this GitHub repository, or via email to dorahermes@gmail.com

## Download data
The data can be downloaded from https://osf.io/q4yad/.
Unzip it and you should have all the data in the folder /data/.
Move this directory to /root/directory/of/Paper_Hermes_2019_eLife/.
The data can also be downloaded by running the code ```gammaModelDownloadData.m``` in /root/directory/of/Paper_Hermes_2019_eLife/.

## Software dependencies
### Make sure that the following toolboxes are downloaded and in the path:
1. MATLAB utility functions written by Kendrick Kay (https://github.com/kendrickkay/knkutils)

   Put the root directory and all sub-directories in MATLAB path.
1. FileTrip (https://github.com/fieldtrip/fieldtrip) 

   It is **not** recommended to put all sub-directories into MATLAB path, but rather put the /root/of/fieldtrip/ in the path, then run ```ft_defaults```.  Check [this](http://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path) out.
1. VISTASOFT of VISTA lab at Stanford University (https://github.com/vistalab/vistasoft)

   Put the root directory and all sub-directories in MATLAB path.
1. JSONio (https://github.com/gllmflndn/JSONio)

   Put the root directory and all sub-directories in MATLAB path.
1. Chronux (http://chronux.org/)

   Put the root directory and all sub-directories in MATLAB path.

### Matlab toolboxes:
1. curve_fitting_toolbox
1. optimization_toolbox
1. signal_toolbox
1. statistics_toolbox

## Reproducing figures
The analyses data can be loaded and the figures can be reproduced with the functions in:
/make_figures/


## Reproducing analyses
Functions to analyses the data and fit the models can be found here: /processing/.
Please note that the raw iEEG data are not included. Interested reader should refer to the reference [Hermes, 2019](https://elifesciences.org/articles/47035) or script ```preprocess_ECoGDataBeforeShare.m``` for the detailed steps of pre-processing the raw data.


## Acknowledgements
This code was made available with the support of the Netherlands Organization for Scientific Research www.nwo.nl under award number 016.VENI.178.048 to Dora Hermes and the National Institute Of Mental Health of the National Institutes of Health under Award Number R01MH111417 to Natalia Petridou and Jonathan Winawer. The ECoG data collection was facilitated by the Stanford Human Intracranial Cognitive Electrophysiology Program (SHICEP) and the Epilepsy Team at the UMC Utrecht Brain Center.

## Reference
Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035. DOI: https://doi.org/10.7554/eLife.47035
