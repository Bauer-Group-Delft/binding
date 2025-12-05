# Ising Model
This repository contains the code for stochastic Ising-like model used in the paper:

 Mijatović, T., Kok, A. R., Zwanikken, J. W., & Bauer, M. (2025). 'Weak transcription factor clustering at binding sites can facilitate information transfer from molecular signals'. Available as arXiv preprint: [arXiv:2505.07641](https://arxiv.org/abs/2505.07641).

## Repository Structure:
```
...
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

## Code Overview 
- Parameters: Code for calculating parameters and writing metadata files.
	- The data for the mean and variance of Bicoid concentrations along the embryo in the folder Bcd_data was kindly provided to us by Thomas Gregor, from the publication: Gregor T, Tank DW, Wieschaus EF, Bialek W. Probing the limits to positional information. Cell. 2007 Jul 13;130(1):153-64. doi: 10.1016/j.cell.2007.05.025. PMID: 17632062; PMCID: PMC2253670.
- AnalytcialModel: Code for calculations on the analytical model.
- StochasticModel: Code for simulation of the stochastic model.
- Data: Metadata and data files storing results of the simulations and calculations. 
	- Measurement folder cotains data for an averaging time of 10 minutes used for calculating the mutual information.
	- TimeSweep folder contains data for different averaging times used for calculating the standard deviation as a funciton of averaging time.
	- (Meta)datafiles for the simulations with higher concentration sampling were too big too upload to github, but can be easily recreated using the supplied scripts (see workflow below, only analytical calculations required).
- Analysis: Code for analysing siulation and claculation results and calculating mutual inforations.
- Scripts: Scripts that use the code and datafiles to write (meta)data files and analyse the results.
- Figures: The figures created using Scripts/FigureCreationScripts.ipynb.

Note that some of the in-code doccumentation might be outdated.

 **Parameter Glossary**
`k_off`: off rate for Bicoid
`k_on`: on rate for Bicoid
`h_tf`: chemical potential of the bulk, referred to as mu in the paper
`e_b`: binding energy of Bicoid to a binding site
`J` : interaction energy between two Bicoid molecules bound to adjacent binding sites

## Usage
1. Create metadata files using Scripts/MetadataCreationScripts.ipynb.
2. Calculate analytical results using Scripts/AnalytcialCalculationScripts.ipynb.
3. Run the stochastic simulations in parallel using StochasticModel/StochasticSimulationRun. Note that when running StochasticModel/StochasticSimulationRun the name of the metatdatafile must be supplied as a command line argument.
4. Recreate figures from article using Scripts/FigureCreationScripts.ipynb.

Stochastic simulations were run in paralel on the DelftBlue super computer.

Most of the data required for generating the figures is already present in the Data folder, so steps 1-3 are optional when attempting to regenerate the figures from the Manuscript. The generated figures can also be found in the figures folder.

## Packages used
- abc
- copy
- datetime
- IPython.display
- itertools
- math
- matplotlib (matplotlib.pyplot, matplotlib.lines)
- mpi4py
- numpy
- os
- pickle
- random
- re
- scipy (scipy.optimize, scipy.sparse, scipy.sparse.linalg, scipy.ndimage, scipy.signal, scipy.io)
- seaborn
- sympy
- sys

See "requirements.txt" for the packages required for running the FigureCreationScripts.ipynb file in the Scripts folder. Note that extra packages from the list above may need to be installed for running the rest of the code.

## Contact
For any questions or inquiries, please contact Marianne Bauer (m.s.bauer@tudelft.nl)
