#%%
# ms-python.python added
import os

try:
    os.chdir(os.path.join(os.getcwd(), "merrycrispr/notebooks"))
    print(os.getcwd())
except:
    pass


#%%
import pandas as pd
import numpy as np
import pyfaidx
import regex

spacers = (
    "/Users/milessmith/workspace/mc_human_files/pml_spacers_to_off_target_score.fa"
)
otrf = "/Users/milessmith/workspace/mc_human_files/bowtie_results.csv"

#%%
