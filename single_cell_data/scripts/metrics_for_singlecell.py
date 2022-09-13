import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import glob
mypath = "metrics"
all_files = glob.glob("metrics/*.csv")

sample_metrics = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    sample_metrics.append(df)

frame = pd.concat(sample_metrics, axis=0, ignore_index=True)
frame.to_csv("results/single_cell_cellranger_diagnostics.csv")





