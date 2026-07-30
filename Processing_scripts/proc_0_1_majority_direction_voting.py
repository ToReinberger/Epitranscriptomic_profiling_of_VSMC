import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns

plot_histo = True

for file in glob.glob(r"direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4\diffmod.table"):
    path_out = "/".join(file.split("\\")[:-1])
    print(file)
    print(path_out)

    data = pd.read_csv(file)
    data["temp"] = data["kmer"] + "_" + data["mod_assignment"]

    print(data.head(5))
    print(out := data[["kmer", "mod_assignment"]].groupby(["kmer", "mod_assignment"]).value_counts().reset_index())

    if plot_histo:
        # Pivot the table to get "higher" and "lower" counts in separate columns
        pivot_df = out.pivot(index='kmer', columns='mod_assignment', values='count')
        # Calculate the ratio
        pivot_df['ratio'] = np.log2(pivot_df['higher'] / pivot_df['lower'])
        # Reset index if needed
        result = pivot_df.reset_index()
        result.kmer = result.kmer.apply(lambda x: x.replace("T", "U"))
        print(result.sort_values(by='ratio', ascending=False))
        # for kmer in result.kmer:
        #     print(result[result.kmer == kmer])
        #     print(result[result.kmer == kmer[::-1]])
        result["middle_base"] = result.kmer.apply(lambda x: x[2])
        result["middle_bases"] = result.kmer.apply(lambda x: x[1:4])
        plt.style.use('ggplot')
        plt.figure(figsize=(9, 8))
        plt.suptitle("C2 vs. C3")
        # dim = int(len(result["middle_bases"].unique())**0.5)
        dim = 2
        row_idx, col_idx = 0, 0
        for b, tbl in result.groupby("middle_base"):
            plt.subplot2grid((dim, dim), (row_idx, col_idx))
            row_idx += 1
            if row_idx == dim:
                col_idx += 1
                row_idx = 0
            sns.histplot(tbl.ratio.values, kde=True, bins=20)
            plt.xlabel("log2(higher/lower)")
            # plt.hist(tbl.ratio.values, bins=100)
            plt.title(b)
        plt.subplots_adjust(bottom=0.2, hspace=0.4, wspace=0.3)
        plt.savefig("direct_RNA_seq/0_PLOTS/Xpore_mode/C3_vs_C4_NEW.png", dpi=300)
        plt.show()

    out2 = out.sort_values(by=["kmer", "count"], ascending=False).drop_duplicates(subset="kmer", keep="first")
    out2["temp"] = out2["kmer"] + "_" + out2["mod_assignment"]
    data = data[data["temp"].isin(out2["temp"])].drop("temp", axis=1)
    # data.to_csv(path_out + "/majority_direction_kmer_diffmod_TR_temp.table", sep="\t", index=False)
    print(data)
