import sys
import os

I_DIR = f'{os.path.dirname(os.path.dirname(os.getcwd()))}/0_analyses/0_exhaustive/'
O_DIR = os.getcwd() + '/jobs/'

for i in os.listdir(I_DIR):
    if "ipynb_checkpoints" not in i:
    	cmd = 'jupyter nbconvert --to script ' + I_DIR + i + ' --output ' + O_DIR + i.split(".ipynb")[0]
    	os.system(cmd) 
