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
   - Run the following command "conda env create -f hcmv_env.yaml -n hcmv_env"
   - After it finishes, run "conda activate hcmv_env". Now you can execute the Python code.
   - **You must navigate to HCMV_Dose_Infection and activate the environment every time you want to run the code**

# Running a simulation

A config file, called "config.yml" is provided which supplies key parameters to the Python scripts that control the simulation. **DO NOT MOVE THIS FILE**

This config file can be opened and modified with most text editors.

To choose which simulation you want to run, simply change the entry next to simul_name, the second line of config.yml, to one of: 
    - 'clump', 'comp', 'acc_dam', 'null', 'clump_comp' or 'clump_acc_dam'. 
Make sure you include the single quotation marks so that Python recognizes it as a string, which it should be. Save the config file **DO NOT CHANGE THE FILE NAME OR MOVE IT**.

To choose which experimental dataset you want to compare the simulation against, change the value of sheet to either 0,1,2,3,4,5, or 6. No single quotations this time. The names of the datasets are in the config file.

To modify the simulation parameters themselves, such as vMax, gamma, kappa, beta, b, etc., navigate to the appropriate simulation, appropriate sheet, and variable you want to change, change it, then save the config file. **DO NOT CHANGE THE FILE NAME OR MOVE IT**. Example: