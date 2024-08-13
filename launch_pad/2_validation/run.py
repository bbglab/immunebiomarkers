import sys
import os

### Step 0 - run exhaustive
print("Run validation..")
os.system('python prepare_the_engines.py') 
os.system('./2_go_validate.sh') 
