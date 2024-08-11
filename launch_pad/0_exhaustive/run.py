import sys
import os

### Step 0 - run exhaustive
print("Run CPI Exhaustive Analyses..")
os.system('python prepare_the_engines.py') 
os.system('./0_go_exhaustive.sh') 
