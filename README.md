# Code documentation for the Fish_performance_LCA repository

## Overview of the repository

This repository contains the code used to reproduce the results of the manuscript: *Jouannais.P,Gibertoni.P, Bartoli.M, Pizzol.M, LCA to evaluate the environmental opportunity costs of biological performances in fish farming (under review)* 

**Cite this repository:**
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7759944.svg)](https://doi.org/10.5281/zenodo.7759944)


 
### Overview of folders and files:

**Data**

+ **_elemental_contents.csv_** includes average elemental compositions (N, P, C) of biochemical classes  (Lipid, Phospholipids, proteins, Carbohydrates) 


+ **_Ingredient_composition_withNP.csv_** includes the fish feed composition (wheat, oil etc.) and the calculated composition in term of biochemical classes (Lipid, proteins, carbohydrates, ash, water)+ NP

+ **_Technosphere_Fish_17_03.csv_** is the foreground technosphere matrix for the fish farms and used for the calculations.

+ **_AH_combi_1.json_**  foreground database importable to Brightway2 (1).


+ **_Micro_for_combi.json_**  foreground database importable to Brightway2 (2).



**Environment**

+ **env_bw_windows.yml** File needed to create the virtual environment on WINDOWS.
+ **env_bw_ubuntu_full.yml** File needed to create the virtual environment on UBUNTU.


**Outputs_Fish**

All outputs of csv and xlsx types are saved in this folder.

**Background_mc_fish**

When running the simulations, the MonteCarlo iterations for the background are saved as pickle objects in this folder.

**Scripts Fish **

+ seven **.py** files: python scripts including the model itself and needed to run the simulations. 

Files, scripts, and their functions'interconnections are mapped below.  
<br>  

<img src="Code map_Fishmodel.jpg"
     alt="Markdown Monster icon"
     style="float: left; margin-right: 10px;" />  
<br>  




**Functions_for_physical_and_biological_calculations_3nd**

Contains functions to calculate values related to anaerobic digestion and valorization of dead fish and sludge


**Technosphere_matrix_modifications_fish_3rd**

Contains functions which modify foreground technospheres  with new mortalities, FCR, excretion etc. 

**Main_functions_fish_scenarios_FCRs**

Script containing the main functions perfoming the simulations 
of the fish farm LCAs


**Simulate_Fish_LCA**

Script which calls the functions to run simulations and produce the final results


 See *Reproducing results from the article.*



**Prepare_project_3nd** 

Creates the Brightway2 project and loads the foreground databases in it. Imports your local version of ecoinvent 3.8 consequential in the new project and loads bioshpere3.

See *Reproducing results from the article.*

**Contrib_plot** 

Contribution analyis


+ one **.R** file: **R_plot** to plot results and get statistics.



<br>

### Reproducing results from the article

*Requirements*

+ Miniconda or Anaconda
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

+ A local version of the ecoinvent 3.8 consequential database

+ A python interface (e.g., Spyder) and a R interface (e.g., Rstudio)


*Step by step procedure:*

1. **Download or clone a local copy of the repository. Keep the folders structure.**

2. **Prepare a conda environment with all needed packages**

+ From terminal or Anaconda/Miniconda terminal access the "Environment" folder. The path depends on where you saved the repository:

```
 cd <yourpathtothefolder/Environment>
```

+ Create the conda environment with all necessary packages using the .yml file corresponding to your OS.

**For Windows:**


```
conda env create --file env_bw_windows.yml
```

+ Activate the newly created environment:

```
conda activate env_bw_windows
```

+ Install the ray package from pip:
```
pip install ray
```

**For Ubuntu:**


```
conda env create --file env_bw_ubuntu.yml
```

+ Activate the newly created environment:

```
conda activate env_bw_ubuntu
```

+ Install the ray package from pip:
```
pip install ray
```

For MACOS, you should try to install the environment from the ubuntu file and, in case of issues, complete the environment by installing the problematic packages manually. 




3. **Set up the Brigtway2 project**

+ In the Scripts directory, open the file **Prepare_project_3nd.py** in a text editor or python interface and change the value of the variable ```ei38dir``` by specifying the directory where the ecoinvent files are on your drive: ```ei38dir = <yourpathtoecoinventfiles>```. 

+ From the python interface or from command, execute the whole script to prepare the Brightway2 project (```python Prepare_project_3nd.py```) .

4. **Run the simulations using the model** 

Mono-dimensional and multi-dimensional samplings are done separately with two different scripts.

+ In the Scripts directory, open the file **Simulate_Fish_LCA** 

+ Change the number of simulations if you want and run the simulation by executing the whole script from the python interface or from command line (```python Simulate.py```). 

+ Wait for all the simulations to be finished  The script will export excel and csv files in the folder  "Outputs_Fish" and "Background_mc_fish"

5. **Plot the figures based on the generated excel files **
  
+ Use **R_plot** and **Contrib_plot** after changing the paths of the files to read in the code editors.


<br>  

