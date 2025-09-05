# Hill Model
This repository contains the code for stochastic Hill function model used in the paper:

 Mijatović, T., Kok, A. R., Zwanikken, J. W., & Bauer, M. (2025). 'Weak transcription factor clustering at binding sites can facilitate information transfer from molecular signals'. Available as arXiv preprint: [arXiv:2505.07641](https://arxiv.org/abs/2505.07641).

## Code Overview 
**Main Scripts**: 
- `DataGeneration.py`: Generates the data from the Hill model and Information bottleneck section in the paper. More information can be found in the code itself. 
- `Figures.ipynb`: Generates the figures from our paper related to the experiment data and Hill model based on the model data generated from `DataGeneration.py`. More information can be found in the notebook itself.

**Dependencies**: Stated in `environment.yml` or `requirements.txt`

## Data Overview 
**Experimental Data:** 
- `TG_normBcd.mat`: Experiment data of the Bicoid gradient, used and published with permission from Thomas Gregor from the publication [T. Gregor, D. W. Tank, E. F. Wieschaus, and W. Bialek, *Cell* 130, 153 (2007)](https://doi.org/10.1016/j.cell.2007.05.025)

 **Model Data**
- After running `DataGeneration.py` for a given set of parameters, a new directory is generated containing the following data: 
	- `parameters.json`: Contains all relevant parameters used for the model. The directory is also named after the (range) of parameters used. 
	- `MIs_all.csv`: Contains all $I(C;s)$ and $I(C;x)$ for different $h$,$k$ combinations for each $\tau$
		- Columns: `h`,`k`, `tau`, `I_CX`, `I_CS`
	- `MIs_optimal.csv`: The values of of $h$,$k$ that give the highest $I(C;x)$ for a given $I(C;s)$ for each $\tau$. 
		- Columns: `tau`, `optimal_h`,`optimal_k`, `I_CS_value`, `optimal_I_CX`
	- `MIs_optimalh.csv`: Contains the values of the optimal I(C;x) and I(C;s) per $h$ and corresponding $k$ for each $\tau$
		- Columns: `tau`, `h`,`optimalh_k_ICS`, `optimalh_I_CS`, `optimalh_k_ICX`, `optimalh_I_CX`
	- `I_CS_IB.csv` & `I_CX_IB.csv`: optimal information bottleneck. 
- Which data is used for which figure is given in `Figures.ipynb`. 

 **Parameter Glossary**
`x`: normalized position along embryo (output)
`s`: normalized Bicoid concentration (input)
`C`: normalized hunchback expression (compressed input)

`koff`: off rate for Bicoid \[s$^{-1}$\]
`kon`: on rate for Bicoid \[s$^{-1}$\]
`k`: dissociation equilibrium constant of the binding devices (koff/kon) \[-\]
`h`: cooperativity parameter / Hill coefficient of the binding sensors \[-\]
`tau`: averaging (or measurement) time \[s\]
`S_bins`: number of bins for discretization of the concentrations
`C_bins`: number of bins discretization of the occupation levels
`C_cov`: minimal \<`C`\> coverage for the optimal h,k combinations

`hillmode`: sets the type of calculation for the Hill model we want to use
- 'koff': sets koff as constant, with a variable kon defined by k
- 'kon': sets kon as constant, with a variable koff defined by k
- 'non-coop': a check for a non-cooperative binding devices (with koff as constant)

## Repository Structure:
```
├── HillModel/ 	
│   ├── Code/
│	     ├── Figures.ipynb
|	     └── DataGeneration.py
│   ├── Data/
│	     ├── ExperimentData/..
|	     └── ModelData/..
│   ├── Figures/..
│   └── HillModel.yml	
...
```

## Getting Started 
1. Clone the repository: (put in the correct reference)
   `git clone [https://github.com/yourusername/your-repo.git](https://github.com/yourusername/your-repo.git)`
2. Install dependencies if needed :
   `python conda env create -f environment.yml`
    or
    `pip install -r requirements.txt`

## Usage
1. (optional) Model data can be generated using `DataGeneration.py`; model parameters need to be changed in the python file before running. 
2. Run `Figures.ipynb` to regenerate the paper figures. Optionally new figures can be created by changing the data directory for each type of figure. Figures are stored in `../Figures/`. 

## Contact
For any questions or inquiries, please contact Marianne Bauer (m.s.bauer@tudelft.nl)