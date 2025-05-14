# gene-harmony-analysis

## Background

Human gene symbols are regulated and follow guidelines established by HGNC. All genes are designated an authoritative symbol (also known as a primary gene symbol), a descriptive name, and an HGNC identification number. Although primary gene symbols are monitored to be unique, alias symbols are not. 

Aliases are additional gene symbols and short descriptions that are used synonymously for the gene and/or any associated gene products. Gene symbols are curated from use in databases, experimental results, and literature. Primary gene symbols and aliases play a crucial role in referencing genes across publications, medical records, and data collections. 

Our preliminary research uncovered collisions, or instances where a single gene symbol was used for multiple different genes. We have categorized them into two kinds: 

a) Alias-primary collisions, which are gene symbols that are used as a primary gene symbol and an alias. The primary gene symbol KRAS is an alias in addition to a primary gene symbol. 

b) Alias-alias collisions are gene symbols that represent an alias of multiple genes. The gene symbol VH is an alias for 35 genes in the NCBI database. 

Out of the 45,646 genes in the HGNC database, 1.40% (637) had alias-primary collisions and 6.82% (3,113) had alias-alias collisions. 
The Ensembl database, which has 41,068 genes, was found to have alias-primary collisions in 1.65% (678) of genes and alias-alias collisions in 6.16% (2,530) of genes. 

The NCBI database, with 45,390 genes, had 3.72% (1,689) and 13.25% (6,013) of genes with alias-primary and alias-alias collisions respectively, illustrating the prevalence of ambiguity that challenges the aggregation of genomic knowledge. 

The [1_alias_primary_collision_analysis](./analysis/1_alias_primary_collision_analysis.ipynb) and [2_alias_alias_collision_analysis](./analysis/alias-alias_collision_analysis.ipynb) Jupyter notebooks show the analyses to get these values

![collision_graphic][def]

## Purpose

The difficulties associated with resolving ambiguity and ensuring accurate understanding of gene symbols restrict the rate of clinical decision-making and contribute to confusion in gene knowledge aggregation. The gene nomenclature system would be most effective if it is unambiguous with a tool to take existing knowledgebase entries as inputs to resolve. 

This curated collection of gene-symbol relationship data will be a foundation for disambiguating gene symbols.
A gene concept with all of the gene symbols that represent it:
![gene_symbol_relationship_graphicASPM](https://github.com/user-attachments/assets/1bd389e3-dd60-4690-a1f3-71bbe6e5868f)

To resolve an amgiuous symbol, these relationships would provide the neccessary context:
![gene_symbol_relationship_graphicASP](https://github.com/user-attachments/assets/c7af96c6-f12b-4b6a-aa16-774111f8c0b7)


# Notebook Dependencies
|   | Name of Notebook                            | Prerequisite Notebook(s) | Input files                            | Notes  |
|---|---------------------------------------------|--------------------------|----------------------------------------|---|
| 1 | alias_primary_collision_analysis            | none                     | ensg_biomart_gene20240626.txt          |   |
|   |                                             |                          | hgnc_biomart_gene20240626.txt          |   |
|   |                                             |                          | Homo_sapiens.gene_info20240627         |   |
| 2 | alias_alias_collision_analysis              | 1                        | none                                   |   |
| 3 | alias_alias_collision_distribution_analysis | 2, 1                     | none                                   |   |
| 4 | symbol_capture_generation                   | 1                        | ensg_mart_export_dros_murin_ortho.txt  | takes longer than an hour to run  |
|   |                                             |                          | ortholog_set_1_df.txt                  |   |
|   |                                             |                          | …                                      |   |
|   |                                             |                          | ortholog_set_10_df.txt                 |   |
| 5 | symbol_capture_analysis                     | 4, 1                     | none                                   | one cell needs to run overnight  |
| 6 | sqlite_symbol_capture_transformation        | 4, 1                     | ensg_biomart_gene20240626.txt          |   |
|   |                                             |                          | hgnc_biomart_gene20240626.txt          |   |
|   |                                             |                          | Homo_sapiens.gene_info20240627         |   |
|   |                                             |                          | ortholog_set_1_df.txt                  |   |
|   |                                             |                          | …                                      |   |
|   |                                             |                          | ortholog_set_10_df.txt                 |   |
| 7 | ambiguous_symbol_distribution_analysis      | 2, 1                     | none                                   |   |
| 8 | concordance_via_networkx_analysis           | 6, 4, 1                  | none                                   |   |
| 9 | concordance_via_upsetplot_analysis          | 6, 4, 1                  | none                                   |   |
| 10 | dgidb_gene_content_analysis      | 2, 1                     | dgidb_genes_JUNE.tsv                                   |   |
| 11 | dgidb_query_analysis           | 10, 2, 1                  | log_data.xlsx                                    |   |    

# Notebook Contents
**1_alias_primary_collision_analysis**
- How many ambiguous symbols resulting from alias-primary collisions are in each database (xxxx_alias_primary_collision_set)
- How many genes are involved in alias-primary collisions in each database
- How many and which genes are involved in alias-primary collisions in all 3 databases (common_ap_collisions)
- How many genes are involved in alias-primary collisions across all 3 databases

**2_alias_alias_collision_analysis**
- How many unique gene records are there in each database (xxxx_gene_id_set)
- How many ambiguous symbols resulting from alias-alias collisions are in each database 

**3_alias_alias_collision_distribution_analysis**
- How many unique primary gene symbols are in each database (xxxx_gene_symbol_count)
- How many primary symbols appear in all 3 databases (all_sources_unique_primary_symbol_count)
- How many unique primary symbols are found between all 3 databases (bw_all_sources_unique_primary_symbol_count)
- How many unique alias symbols there are per database and across all 3
- How many alias symbols appear in all 3 databases
- How many genes are involved in alias-alias collisions in each database
- How many and which genes are involved in alias-primary collisions in all 3 databases (all_sources_aa_collision_genes)
- How many genes are involved in alias-alias collisions across all 3 databases

- Per each database: How many genes are the ambiguous gene symbols, resulting from alias-alias collisions, being shared between
- How many gene concept to symbol relationships there are
  
# How can you help?

Contributing information on collisions that you come across will help collect data on the collisions that would be most impactful to resolve as well as increasing the data available for developing resolution strategies for downstream tool development.

## Contact Information

For any feedback, questions, or conversation, please make an issue.


[def]: https://github.com/cancervariants/gene-harmony-analysis/assets/109570522/91425d67-0884-4fbc-83ab-e7cfd8bd57bd
