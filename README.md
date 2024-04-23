These are Python scripts and data which will allow one to perform various simulations to explore cooperativity between viruses

**Set up**
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
           **You must navigate to HCMV_Dose_Infection and activate the environment every time you want to run the code**

  ****Running a simulation**** \n
  A config file, called "config.yml" is provided which supplies key parameters to the Python scripts that control the simulation. **DO NOT MOVE THIS FILE**
  This config file can be opened and modified with most text editors.

  To choose which simulation you want to run