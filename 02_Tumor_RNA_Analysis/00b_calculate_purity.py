from puree import *

# combined_data_path = "/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts/data/processed/rna/asc_merged_with_gtex.txt"

# p = PUREE()
# purities_and_logs = p.get_output(combined_data_path, "HGNC")
# print(purities_and_logs)

df = pd.read_csv("/Users/hoyin/Downloads/gene_tpm_v10_bladder.gct", sep="\t", skiprows=2)
#print(df)
df = df.drop_duplicates(subset=["Description"]).set_index("Description").drop(columns=["Name"])
#df = df.T
df = df.rename_axis(None, axis=1)
print(df)

df.to_csv("/Users/hoyin/Downloads/gene_tpm_v10_bladder_processed.tsv", sep="\t")