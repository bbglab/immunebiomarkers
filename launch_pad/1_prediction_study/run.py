import sys
import os

### Step 0 - run exhaustive
print("Run prediction study..")
os.system('python prepare_the_engines.py') 
os.system('./1_go_study.sh') 
