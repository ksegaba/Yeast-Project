"""
Code to make pairs of environment-specific fitness traits to run multi-trait 
models on.
"""
__author__ = "Kenia Segura Abá"


import os
import pandas as pd
from itertools import permutations


os.chdir("/mnt/home/seguraab/Shiu_Lab/Project/")


if __name__ == "__main__":
    # read in fitness data
    df = pd.read_csv("Data/Peter_2018/pheno.csv")

    # make permutations of trait pairs
    pairs = list(permutations(list(df.columns[1:]), 2))
    pairs = list(map(tuple, map(sorted, pairs))) # sort the tuples
    pairs = set(pairs) # keep only the unique tuples 595 out of 1190

    # write to file
    with open("../tmp/data/env_fitness_pairs.txt", "w") as f:
        for pair in pairs:
            f.write(f"{pair[0]}\t{pair[1]}\n")