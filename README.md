# PMBO
Polymer Multiobjective Bayesian Optimization (PMBO) assiciated by physical descriptors
## Description
Polymer Multiobjective Bayesian Optimization (PMBO)  framework for exploring multifunctional polymers, including four components: 1) Polymer Library, which stores all polymer candidates; 2) Feature extraction, which based on physical feature engineering or graph descriptors generation; 3) Polymer properties simulator, which calculates polymer properties at different scales using DFT, MD and ML. 4) Multi-objective Bayesian algorithm, which predicts polymer properties and recommends potentially multifunctional polymers.Please refer to our work "Discovery of Multifunctional Polymers in Constrained Chemical Space Via Physical Descriptors-Guided Multi-objective Bayesian Optimization" for additional details.
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
## Authors

| **AUTHORS** |Xiang Huang, Shenghong Ju            |
|-------------|--------------------------------------------------|
| **VERSION** | V1.0 / July,2023                               |
| **EMAILS**  | shenghong.ju@sjtu.edu.cn                         |

## Attribution
This work is under BSD-2-Clause License. Please, acknowledge use of this work with the appropiate citation to the repository and research article.
