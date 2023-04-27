#!/usr/bin/python3
"""
Pairplot of factors and RF test performance after feature selection
with regression line.
"""

import scipy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Read in data
    X = pd.read_csv('snps_RF_FS_4factors_factors_features.csv', index_col=0) # same for all data types
    snp = pd.read_csv('/mnt/home/seguraab/Shiu_Lab/Project/Results/RESULTS_RF_SNPs_FS.txt', sep='\t')
    orf = pd.read_csv('/mnt/home/seguraab/Shiu_Lab/Project/Results/RESULTS_RF_ORFs_FS.txt', sep='\t')
    cnv = pd.read_csv('/mnt/home/seguraab/Shiu_Lab/Project/Results/RESULTS_RF_CNVs_FS.txt', sep='\t')
    pc = pd.read_csv('/mnt/home/seguraab/Shiu_Lab/Project/Results/RESULTS_RF_PCs_sorted.txt', sep='\t')

    # Merge datasets
    X1 = snp[["cond", "r2_test"]].merge(X, left_on="cond", right_index=True)
    X1["Type"] = "SNPs"
    X2 = orf[["cond", "r2_test"]].merge(X, left_on="cond", right_index=True)
    X2["Type"] = "ORFs"
    X3 = cnv[["cond", "r2_test"]].merge(X, left_on="cond", right_index=True)
    X3["Type"] = "CNVs"
    X4 = pc[["cond", "r2_test"]].merge(X, left_on="cond", right_index=True)
    X4["Type"] = "PCs"

    # LM Plots
    # X_plot = pd.concat([X1, X2, X3, X4])
    for df in [X1, X2, X3, X4]:
        p = sns.pairplot(df, x_vars=["h2", "var", "median", "severity"],
                y_vars="r2_test", aspect=1, kind="reg", height=3) # hue="cond" # hue="Type"
        p.set(ylim=(0,1))
        # handles = p._legend_data.values()
        # labels = p._legend_data.keys()
        # p.fig.legend(handles=handles, labels=labels, loc="lower center", ncol=7,
        #              bbox_to_anchor=(0.5, -0.05))
        # plt.tight_layout()
        
        # Annotate linear equation for r2 vs h2 plot
        # dir(p.fig) # get attributes
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x=df.h2, y=df.r2_test)
            # x=p.fig.get_children()[1].get_lines()[0].get_xdata(),
            # y=p.fig.get_children()[1].get_lines()[0].get_ydata()) # these data do not match r2 and h2
        p.fig.get_children()[1].text(.4, .75, str(f"y = {slope:.2f} x + {intercept:.1f}\nPCC = {r_value:.2f}, R2={r_value**2:.2f}\nP-value = {p_value:.4f}"))
        print(f"Line for r2_test vs h2: y = {slope:.2f} x + {intercept:.1f}\nPCC = {r_value:.2f}\nP-value = {p_value:.4f}")

        # Annotate linear equation for r2 vs var plot
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x=df['var'], y=df.r2_test)
        p.fig.get_children()[2].annotate(f"y = {slope:.2f} x + {intercept:.1f}\nPCC = {r_value:.2f}, R2={r_value**2:.2f}\nP-value = {p_value:.4f}",
                                         xy=(.04, .75))
        print(f"Line for r2_test vs var: y = {slope:.2f} x + {intercept:.1f}\nPCC = {r_value:.2f}, R2={r_value**2:.2f}\nP-value = {p_value:.4f}")

        # Annotate linear equation for r2 vs med plot
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x=df['median'], y=df.r2_test)
        p.fig.get_children()[3].annotate(f"y = {slope:.2f} x + {intercept:.1f}\nPCC = {r_value:.2f}, R2={r_value**2:.2f}\nP-value = {p_value:.4f}",
                                         xy=(.25, .75))
        print(f"Line for r2_test vs med: y = {slope:.2f} x + {intercept:.1f}\nPCC = {r_value:.2f}, R2={r_value**2:.2f}\nP-value = {p_value:.4f}")

        # Annotate linear equation for r2 vs severity plot
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x=df.severity, y=df.r2_test)
        p.fig.get_children()[4].annotate(f"y = {slope:.2f} x + {intercept:.1f}\nPCC = {r_value:.2f}, R2={r_value**2:.2f}\nP-value = {p_value:.4f}",
                                         xy=(.15, .75))
        print(f"Line for r2_test vs med: y = {slope:.2f} x + {intercept:.1f}\nPCC = {r_value:.2f}, R2={r_value**2:.2f}\nP-value = {p_value:.4f}")

        # Save figure
        plt.savefig(f"RF_FS_{df.Type.values[0]}_4factors_lmplot.pdf")
