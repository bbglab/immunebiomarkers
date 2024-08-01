import sys
import os

### Get current directory, must be run from launch pad folder 
wd = os.getcwd() + '/launch_pad/'
prep_cmd = 'python prepare_the_engines.py'

### Step 0 - run exhaustive
print("Run CPI Exhaustive Analyses..")
os.chdir(wd + '/0_exhaustive/')
os.system(prep_cmd) 
os.system('./0_go_exhaustive.sh') 

### Step 1 - xgboost study
print("Run XGboost study..")
os.chdir(wd + '/1_prediction_study/')
cmd1 = './1_go_study.sh'
os.system(cmd1)

### Step 2 - validation 
print("Run validation analyses..")
os.chdir(wd + '/2_validation/')
os.system(prep_cmd)
os.system('./2_go_validate.sh')

### Step 3 - figures
print("Create figures..") 
os.chdir(wd + '/3_figures/')
os.system(prep_cmd)
os.system('./3_make_figures.sh')
