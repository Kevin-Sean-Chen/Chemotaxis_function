# Chemotaxis_function

This is a repo for behavioral modeling for worm odor-navigation reported in the paper: https://arxiv.org/abs/2311.07117 (Chen et al. 2023). 
The repo analyzes behavioral trajectories that are output from the high-throughput worm tracking pipeline in the Leifer lab (https://github.com/Kevin-Sean-Chen/leifer-Behavior-Triggered-Averaging-Tracker-new)
The corresponding data for figure making and analysis can be found in figshare: [10.6084/m9.figshare.24764403](https://doi.org/10.6084/m9.figshare.24764403) 

## Making figures

Please download the data set from figure, adjust the corresponding downloaded directory, and run scripts in the `figures` folder. 
Figures should be self-contained in this folder. To generate figures from the original data, please check the `fig_note.txt` document in figshare folder describes specific scripts and data to use for analysis from scratch.

## Model demonstration

To validate the generative model and inference procedure, we use code in the `demo` folder to simulate behavioral output given random odor concentration and infer ground-truth parameters.
Run `mGLM.m` for the dPAW model proposed in our paper. A simplified turning model is shown in `bGLM.m`.

## Model fitting

We use scripts in `inference` folder and optimize objectives in the `functions` folder to fit models. Specifically, we use `Chemotaxis_inference_pop.m` and `classify_learn3.m` to fit MLE of parameters and to run cross-validation,
use `mutant_fit.m` to fit to different worm strains, and use `opto_mGLM.m` to fit to perturbation data.

## Track analysis

To analyze tracks and pre-process data for model fitting, we use files in the `scripts` folder. Specifically, `BRW_WV_CI.m` is used to compute different chemotaxis indices, 
`OdorFlow_3.m` is used to visualize behavioral trajectories given parameters. We call `track2data.m` functions in the `function` folder to save as Data structure that contains fields for analysis.
To load from raw data, please include https://github.com/Kevin-Sean-Chen/leifer-Behavior-Triggered-Averaging-Tracker-new repo to use additional functions.

## Citation

If you use the design files and/or codes provided in this repository, please cite:
> Kevin S. Chen, Jonathan W. Pillow, Andrew M. Leifer. (2023). Olfactory learning alters navigation strategies and behavioral variability in C. elegans arXiv:2311.07117

## Lisence
Design files and codes provided in this repository are free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
