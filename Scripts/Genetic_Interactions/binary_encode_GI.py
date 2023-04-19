"""
One-hot encoding of composite GI features into floats 
for downstream modeling and analysis.

Written by: Kenia Segura Abá
"""
import datatable as dt
import pandas as pd

def encode_GI(X, ID, save_name):
    X = X.to_pandas()
    
    # Binary encoding
    X = pd.get_dummies(X)
    X = X.to_dict(orient="records")

    # Write to file
    X = dt.Frame(X)
    X["ID"] = ID
    X.to_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/%s"%save_name)
    return X

def main():
    # Read in data 
    X1 = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All_1.csv", header=True)
    ID = X1[:,1] # subset isolate IDs
    X1 = X1[:,2:] # Subset only composite features (no instance IDs)

    encode_GI(X1, ID, "geno_GI_All_1_binary.csv")
    del X1

    X2 = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All_2.csv", header=True)
    encode_GI(X2, ID, "geno_GI_All_2_binary.csv")
    del X2

    X3 = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All_3.csv", header=True)
    encode_GI(X3, ID, "geno_GI_All_3_binary.csv")
    del X3

if __name__ == "__main__":
    main()
