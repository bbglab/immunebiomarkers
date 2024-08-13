import sys
import os

### Step 0 - run exhaustive
print("Run validation..")
os.system('python prepare_the_engines.py') 
os.system('./3_make_figures.sh') 
