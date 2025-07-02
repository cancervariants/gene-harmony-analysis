from pathlib import Path
import pandas as pd

ROOT_PATH = Path(__file__).resolve().parent
# ROOT_PATH = '/Users/rsaxs014/Desktop/gene-harmony-analysis'

merged_alias_ap_collision_hgnc_df = pd.read_csv(ROOT_PATH / "merged_alias_ap_collision_hgnc_df.csv")

ap_record_set_hgnc = set(merged_alias_ap_collision_hgnc_df["HGNC_ID"])
ap_record_count_hgnc = len(ap_record_set_hgnc)

merged_alias_ap_collision_ensg_df = pd.read_csv(ROOT_PATH / "merged_alias_ap_collision_ensg_df.csv")

ap_record_set_ensg = set(merged_alias_ap_collision_ensg_df["ENSG_ID"])
ap_record_count_ensg = len(ap_record_set_ensg)

merged_alias_ap_collision_ncbi_df = pd.read_csv(ROOT_PATH / "merged_alias_ap_collision_ncbi_df.csv")

ap_record_set_ncbi = set(merged_alias_ap_collision_ncbi_df["NCBI_ID"])
ap_record_count_ncbi = len(ap_record_set_ncbi)

merged_alias_aa_collision_hgnc_df = pd.read_csv(ROOT_PATH / "merged_alias_aa_collision_hgnc_df.csv")

aa_record_set_hgnc = set(merged_alias_aa_collision_hgnc_df["HGNC_ID"])
aa_record_count_hgnc = len(aa_record_set_hgnc)

merged_alias_aa_collision_ensg_df = pd.read_csv(ROOT_PATH / "merged_alias_aa_collision_ensg_df.csv")

aa_record_set_ensg = set(merged_alias_aa_collision_ensg_df["ENSG_ID"])
aa_record_count_ensg = len(aa_record_set_ensg)

merged_alias_aa_collision_ncbi_df = pd.read_csv(ROOT_PATH / "merged_alias_aa_collision_ncbi_df.csv")

aa_record_set_ncbi = set(merged_alias_aa_collision_ncbi_df["NCBI_ID"])
aa_record_count_ncbi = len(aa_record_set_ncbi)

mini_hgnc_df = pd.read_csv(ROOT_PATH / "mini_hgnc_df.csv")

gene_record_set_hgnc = set(mini_hgnc_df["HGNC_ID"])
gene_record_count_hgnc = len(gene_record_set_hgnc)

mini_ensg_df = pd.read_csv(ROOT_PATH / "mini_ensg_df.csv")

gene_record_set_ensg = set(mini_ensg_df["ENSG_ID"])
gene_record_count_ensg = len(gene_record_set_ensg)

mini_ncbi_df = pd.read_csv(ROOT_PATH / "mini_ncbi_df.csv")

gene_record_set_ncbi = set(mini_ncbi_df["NCBI_ID"])
gene_record_count_ncbi = len(gene_record_set_ncbi)