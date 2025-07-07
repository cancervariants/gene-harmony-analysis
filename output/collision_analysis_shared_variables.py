from pathlib import Path
import pandas as pd

ROOT_PATH = Path(__file__).resolve().parent
# ROOT_PATH = '/Users/rsaxs014/Desktop/gene-harmony-analysis'

merged_alias_ap_collision_hgnc_df = pd.read_csv(ROOT_PATH / "merged_alias_ap_collision_hgnc_df.csv")

ap_record_set_hgnc = set(merged_alias_ap_collision_hgnc_df["HGNC_ID"])
ap_record_count_hgnc = len(ap_record_set_hgnc)

ap_ambiguous_symbol_set_hgnc = set(merged_alias_ap_collision_hgnc_df["collision"])
ap_ambiguous_symbol_count_hgnc = len(ap_ambiguous_symbol_set_hgnc)

merged_alias_ap_collision_ensg_df = pd.read_csv(ROOT_PATH / "merged_alias_ap_collision_ensg_df.csv")

ap_record_set_ensg = set(merged_alias_ap_collision_ensg_df["ENSG_ID"])
ap_record_count_ensg = len(ap_record_set_ensg)

ap_ambiguous_symbol_set_ensg = set(merged_alias_ap_collision_ensg_df["collision"])
ap_ambiguous_symbol_count_ensg = len(ap_ambiguous_symbol_set_ensg)

merged_alias_ap_collision_ncbi_df = pd.read_csv(ROOT_PATH / "merged_alias_ap_collision_ncbi_df.csv")

ap_record_set_ncbi = set(merged_alias_ap_collision_ncbi_df["NCBI_ID"])
ap_record_count_ncbi = len(ap_record_set_ncbi)

ap_ambiguous_symbol_set_ncbi = set(merged_alias_ap_collision_ncbi_df["collision"])
ap_ambiguous_symbol_count_ncbi = len(ap_ambiguous_symbol_set_ncbi)



merged_alias_aa_collision_hgnc_df = pd.read_csv(ROOT_PATH / "merged_alias_aa_collision_hgnc_df.csv")

aa_record_set_hgnc = set(merged_alias_aa_collision_hgnc_df["HGNC_ID"])
aa_record_count_hgnc = len(aa_record_set_hgnc)

aa_ambiguous_symbol_set_hgnc = set(merged_alias_aa_collision_hgnc_df["collision"])
aa_ambiguous_symbol_count_hgnc = len(aa_ambiguous_symbol_set_hgnc)

merged_alias_aa_collision_ensg_df = pd.read_csv(ROOT_PATH / "merged_alias_aa_collision_ensg_df.csv")

aa_record_set_ensg = set(merged_alias_aa_collision_ensg_df["ENSG_ID"])
aa_record_count_ensg = len(aa_record_set_ensg)

aa_ambiguous_symbol_set_ensg = set(merged_alias_aa_collision_ensg_df["collision"])
aa_ambiguous_symbol_count_ensg = len(aa_ambiguous_symbol_set_ensg)

merged_alias_aa_collision_ncbi_df = pd.read_csv(ROOT_PATH / "merged_alias_aa_collision_ncbi_df.csv")

aa_record_set_ncbi = set(merged_alias_aa_collision_ncbi_df["NCBI_ID"])
aa_record_count_ncbi = len(aa_record_set_ncbi)

aa_ambiguous_symbol_set_ncbi = set(merged_alias_aa_collision_ncbi_df["collision"])
aa_ambiguous_symbol_count_ncbi = len(aa_ambiguous_symbol_set_ncbi)



mini_hgnc_df = pd.read_csv(ROOT_PATH / "mini_hgnc_df.csv")

gene_record_set_hgnc = set(mini_hgnc_df["HGNC_ID"])
gene_record_count_hgnc = len(gene_record_set_hgnc)

mini_ensg_df = pd.read_csv(ROOT_PATH / "mini_ensg_df.csv")

gene_record_set_ensg = set(mini_ensg_df["ENSG_ID"])
gene_record_count_ensg = len(gene_record_set_ensg)

mini_ncbi_df = pd.read_csv(ROOT_PATH / "mini_ncbi_df.csv")

gene_record_set_ncbi = set(mini_ncbi_df["NCBI_ID"])
gene_record_count_ncbi = len(gene_record_set_ncbi)



primary_symbol_set_hgnc = set(mini_hgnc_df['gene_symbol'])
primary_symbol_count_hgnc = len(primary_symbol_set_hgnc)
primary_symbol_count_hgnc

alias_symbol_set_hgnc = set(mini_hgnc_df['alias_symbol'])
alias_symbol_count_hgnc = len(alias_symbol_set_hgnc)

total_symbol_count_hgnc = len(set(mini_hgnc_df['gene_symbol']) | set(mini_hgnc_df['alias_symbol']))

primary_symbol_set_ensg = set(mini_ensg_df['gene_symbol'])
primary_symbol_count_ensg = len(primary_symbol_set_ensg)

alias_symbol_set_ensg = set(mini_ensg_df['alias_symbol'])
alias_symbol_count_ensg = len(alias_symbol_set_ensg)

total_symbol_count_ensg = pd.concat([mini_ensg_df["gene_symbol"], mini_ensg_df["alias_symbol"]]).dropna().nunique()

primary_symbol_set_ncbi = set(mini_ncbi_df['gene_symbol'])
primary_symbol_count_ncbi = len(primary_symbol_set_ncbi)

alias_symbol_set_ncbi = set(mini_ncbi_df['alias_symbol'])
alias_symbol_count_ncbi = len(alias_symbol_set_ncbi)

total_symbol_count_ncbi = pd.concat([mini_ncbi_df["gene_symbol"], mini_ncbi_df["alias_symbol"]]).dropna().nunique()