These are Python scripts and data which will allow one to perform various simulations to explore cooperativity between viruses

# Set up
Must have Python 3.12 or later version installed on your computer.
Must have the following Python libraries installed:
 - numpy
 - matplotlib
 - pandas
 - openpyxl
 - lmfit

You can install these using two methods:
1. pip
 - Navigate to the direction in which this repository is located. Should be something like /users/my_name/.../HCMV_Dose_Infection
 - run the following commands: "pip install numpy", "pip install matplotlib", "pip install pandas", "pip install openpyxl", "pip install lmfit"
2. Anaconda (conda)
   - You must first have Anaconda, a Python environment mamanger, installed on your device.
   - Navigate to the folder in which this repository is located. Should be something like /users/my_name/.../HCMV_Dose_Infection
   - Run the following command: *conda env create -f hcmv_env.yaml -n hcmv_env*
   - After it finishes, run: *conda activate hcmv_env*. Now you can execute the Python code.
   - **You must navigate to HCMV_Dose_Infection and activate the environment every time you want to run the code**

# Running a simulation

A config file, called "config.yml" is provided which supplies key parameters to the Python scripts that control the simulation. **DO NOT MOVE THIS FILE**

This config file can be opened and modified with most text editors.

To choose which simulation you want to run, simply change the entry next to simul_name, the second line of config.yml, to one of: 
    - 'clump', 'comp', 'acc_dam', 'null', 'clump_comp' or 'clump_acc_dam'. 
Make sure you include the single quotation marks so that Python recognizes it as a string, which it should be. Save the config file **DO NOT CHANGE THE FILE NAME OR MOVE IT**.

To choose which experimental dataset you want to compare the simulation against, change the value of sheet to either 0,1,2,3,4,5, or 6. No single quotations this time. The names of the datasets are in the config file.

To modify the simulation parameters themselves, such as vMax, gamma, kappa, beta, b, etc., navigate to the appropriate simulation, appropriate sheet, and variable you want to change, change it, then save the config file. **DO NOT CHANGE THE FILE NAME OR MOVE IT**. Example 1:

> You are running a *null* simulation against the *2020_05_29 TR_GFP_fibroblast* dataset and want to change the value of b from 0.1 to 0.15.
> - Navigate to the *NULL_PARAMETERS* section of the config file. 
> - Navigate to the *b* subsection under *NULL_PARAMETERS*. 
> - Find the entry of *b* corresponding to *# sheet=2 2020_05_29 TR_GFP_fibroblast*. Delete 0.1 and change it to 0.15.
> - Save the config file. **DO NOT CHANGE THE FILE NAME OR MOVE IT**

**- - - - Running a simulation - - - -**
- If running from the terminal/command prompt, ensure you are in the correct folder this repository is located in. Should be something like: /users/my_name/.../HCMV_Dose_Infection.
- Ensure that you have changed config.yml to what you want.
- If you installed using the Anaconda (conda) method, ensure that you have activated the conda environment.
- Type *python run_simulation.py* and press enter.

# Plotting a simulation
There are two scripts for plotting simulation results, *make_simple_simul_plot* and *plot_clump_simul*. Both make a plot comparing the simulation results against an experimental dataset. But *plot_clump_simul* shows extra information about the clump distribution.

**Setting up *make_simple_simul_plot.py*.**
- Open this file in a Python IDE or a text editor.
- In line 18, where it says "file", type the path to the simulation results file you want to plot. These paths will always have the formula:
> simulation_results/simulation_type/filename.csv
> Example: *"simulation_results/null/NullSimul_2022_11_02_TB_GFP_fib_s=50_vMax=800.0_b=0.1_n=5.csv"*
- Make sure that the path is surrounded in quotations " ", save the file. **DO NOT CHANGE THE FILE NAME OR MOVE IT**
- If you only want to see the figure and not save it, comment out the *fig.savefig* line at the bottom using a #

**Setting up *plot_clump_simul.py*.**
- Open this file in a Python IDE or a text editor.
- In line 18, where it says "file_simul", type the path to the simulation results file you want to plot. These paths will always have the formula:
> simulation_results/clump_simulation_type/filename.csv
> Example: *"simulation_results/clump/ClumpSimul_2022_11_02_TB_GFP_fib_s=100_vMax=1000.0_poly_norm_n=3.csv"*
- In line 19, where it says "file_simul_size", type the path to the .json file of clump sizes which corresponds to the simulation. You can identify the matching files by parameters listed in the file name (ignore "n" and "run"). These .json files are contained in a subfolder called "clump_information".
> Example: *"simulation_results/clump/clump_information/ClumpSimul_2022_11_02_TB_GFP_fib_s=100_vMax=1000.0_poly_norm_run=1_CLUMP.json"*
- Make sure that the path is surrounded in quotations " ", save the file. **DO NOT CHANGE THE FILE NAME OR MOVE IT**
- If you only want to see the figure and not save it, comment out the *fig.savefig* line at the bottom using a #

**- - - - Running the plotting scripts - - - -**
- If running from the terminal/command prompt, ensure you are in the correct folder this repository is located in. Should be something like: /users/my_name/.../HCMV_Dose_Infection.
- Ensure that you have changed config.yml to what you want.
- If you installed using the Anaconda (conda) method, ensure that you have activated the conda environment.
- Type either *python make_simple_simul_plot.py* or *python plot_clump_simul.py* depending on what you want