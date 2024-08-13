import sys
import os

I_DIR = f'{os.path.dirname(os.path.dirname(os.getcwd()))}'
O_DIR = os.getcwd() + '/jobs/'

### Run supplement data prep
tmp_dirs = [I_DIR + '/0_analyses/3_supplement/', I_DIR + '/1_figures/supplement/', I_DIR + '/1_figures/0_make_plot_collections/'] + [ I_DIR + '/1_figures/figure' + str(j+1) + '/' for j in range(1,5)]

for j in tmp_dirs:
    for i in os.listdir(j):
        if 'checkpoints' not in i and 'README' not in i:
            print(i)
            cmd = 'jupyter nbconvert --to script ' + j + i + ' --output ' + O_DIR + i.split(".ipynb")[0]
            os.system(cmd)