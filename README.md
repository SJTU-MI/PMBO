# PMBO
Polymer Multiobjective Bayesian Optimization (PMBO) assiciated by physical descriptors
## Description
Polymer Multiobjective Bayesian Optimization (PMBO)  framework for exploring multifunctional polymers, including four components: 1) Polymer Library, which stores all polymer candidates; 2) Feature extraction, which based on physical feature engineering or graph descriptors generation; 3) Polymer properties simulator, which calculates polymer properties at different scales using DFT, MD and ML. 4) Multi-objective Bayesian algorithm, which predicts polymer properties and recommends potentially multifunctional polymers.Please refer to our work "Discovery of Multifunctional Polymers Via Physical Descriptors-Guided Multi-objective Bayesian Optimization" for additional details.
![Framework](https://github.com/SJTU-MI/PMBO/blob/main/Framework.png)
## Installation
### Files loading and environment setup:

To download, clone this repository:<br>
````
git clone https://github.com/SJTU-MI/PMBO.git
````

To run most code in this repository, the relevant anaconda environment can be installed from environment.yml. To build this environment, run：<br>
````
cd ./PMBO
conda env create -f environment.yml
conda activate pmbo
````
Additionally, polymer physical feature engineering and properties calculations can be accessed at other GitHub repositories of [APFEforPI](https://github.com/SJTU-MI/APFEforPI) and [RadonPy](https://github.com/RadonPy/RadonPy).
## Try the desired parts of the project:
### Code in the pmbo folder
**gp.py**: Gaussian process regression model <br>
**Acq_fun.py**: Acquisition functions such as EI(expected improvement) and UCB(upper confidence bound) <br>
**hypervolume.py**: Calculation of hypervolume <br>
**utility.py**: Utility functions such as Pareto front allocation and data pre/post processing <br>
**optimize.py**: Core of multi-objective Bayesian optimization <br>
**log.py**: PMBO Logo <br>
### Tutorial
**MBO_tutorial.ipynb**: A case of multi-objective optimization for multifunctional polymers discovery <br>
**example.csv**:Benchmark dataset for testing (input file) <br>
**cal_data.csv**: MBO recommended polymers and their observed properties (output file) <br>
**HV.csv**: Optimized convergence curve evaluated by hypervolume (output file) <br>
## Authors

| **AUTHORS** |Xiang Huang, Shenghong Ju            |
|:-------------:|--------------------------------------------------|
| **VERSION** | V1.0 / July,2023                               |
| **EMAILS**  | shenghong.ju@sjtu.edu.cn                         |
| **GROUP HOME**  | https://ju.sjtu.edu.cn/en/                         |

## Attribution
This work is under BSD-2-Clause License. Please, acknowledge use of this work with the appropiate citation to the repository and research article.
