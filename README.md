# Code underlying: Weak transcription factor clustering at binding sites can facilitate information transfer from molecular signals

This repository contains the code and data accompanying the paper:
Mijatović, T., Kok, A. R., Zwanikken, J. W., & Bauer, M. (2025). 'Weak transcription factor clustering at binding sites can facilitate information transfer from molecular signals'. Available as arXiv preprint: [arXiv:2505.07641](https://arxiv.org/abs/2505.07641).

Work performed by Zwanikken Group and Bauer Group at TU Delft.

## Abstract 
Transcription factor concentrations provide signals to cells that allow them to regulate gene expression to make correct cell fate decisions. Calculations for noise bounds in gene regulation suggest that clustering or cooperative binding of transcription factors decreases signal-to-noise ratios at binding sites. However, clustering of transcription factor molecules around binding sites is frequently observed. We develop two complementary models for clustering transcription factors at binding site sensors that allow us to study information transfer from a signal to a variable relevant to development, namely future cell fates, on the example of the signal morphogen Bicoid in early fly development. We find that weak cooperativity or clustering can allow for maximal information transfer, especially about the relevant variable. Why clustering is beneficial for transfer depends on how long the binding site sensor measures the signal before affecting downstream gene expression: for short measurements, clustering allows for the implementation of a switch, while for long measurements, weak clustering allows the sensor to access maximal developmental information provided in a nonlinear signal. Finally, we find that clustering can facilitate binding site sensors to achieve an optimal bound to molecular sensing, the information bottleneck (IB) bound. This IB optimality appears consistent with a biologically intuitive optimization.

## Overview
In the paper we discuss two different models for clustering transcription factors at a single binding site sensor: 
1. A stochastic model based on the Hill function 
2. A mechanistic model based on an Ising-like model 
This repository contains the code for generating the data from the models used for the figures in the paper. There two main folders; one for each model. Within the two directories, `HillModel\` and `IsingModel\`, all information can be found about the code and data in the READMEmd for the respective models. 

For the Bicoid mean and variance, we use data kindly provided by Thomas Gregor, from the publication: [T. Gregor, D. W. Tank, E. F. Wieschaus, and W. Bialek, *Cell* 130, 153 (2007)](https://doi.org/10.1016/j.cell.2007.05.025). 

## Repository Structure:
```
├── HillModel/ 	
│   ├── Code/
│	│    ├── Figures.ipynb
|	│    └── DataGeneration.py
│   ├── Data/
│	│    ├── ExperimentData/..
|	│    └── ModelData/..
│   ├── Figures/..
│   ├── HillModel.yml
│   └── README.md 
├── IsingModel/ 
│   ├── Analysis/..
│   ├── AnalyticalModel/..
│   ├── Data/
│	│    ├── Measurement
|	│    └── TimeSweep
│   ├── Figures/..
│   ├── Parameters/
│	│    ├── Bcd_data/..
|	│    └── ..
│   ├── Scripts/
│	│    ├── BicoidDataAnalysisScripts.ipynb
│	│    ├── MetadataCreationScripts.ipynb
│	│    ├── AnalyticalCalculationScripts.ipynb
|	│    └── FigureCreationScripts.ipynb
│   ├── StochasticModel/..
|	└── README.md 
├── LICENSE.md
└── README.md 
```
## License
Distributed under the X License. See `LICENSE.md for more information.

## Contact
For any questions or inquiries, please contact Marianne Bauer (m.s.bauer@tudelft.nl)