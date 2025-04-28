# gene-harmony-analysis

## Background

Human gene symbols are regulated and follow guidelines established by HGNC. All genes are designated an authoritative symbol (also known as a primary gene symbol), a descriptive name, and an HGNC identification number. Although primary gene symbols are monitored to be unique, alias symbols are not. 

Aliases are additional gene symbols and short descriptions that are used synonymously for the gene and/or any associated gene products. Gene symbols are curated from use in databases, experimental results, and literature. Primary gene symbols and aliases play a crucial role in referencing genes across publications, medical records, and data collections. 

Our preliminary research uncovered collisions, or instances where a single gene symbol was used for multiple different genes. We have categorized them into two kinds: 

a) Alias-primary collisions, which are gene symbols that are used as a primary gene symbol and an alias. The primary gene symbol KRAS is an alias in addition to a primary gene symbol. 

b) Alias-alias collisions are gene symbols that represent an alias of multiple genes. The gene symbol VH is an alias for 35 genes in the NCBI database. 

Out of the 43,164 genes in the HGNC database, 483 (1.12%) had alias-primary collisions and 2,084 (4.83%) had alias-alias collisions. 
The Ensembl database, which has 40,353 genes, was found to have alias-primary collisions in 218 (0.54%) of genes and alias-alias collisions in 3,680 (9.12%) of genes. 

The NCBI database, which had the largest number of genes- 75,346, had 1,712 (2.27%) and 5,670 (7.53%) of genes with alias-primary and alias-alias collisions respectively, illustrating the prevalence of ambiguity that challenges the aggregation of genomic knowledge. 

The [alias-primary_collision_analysis](./alias-primary_collisions/alias-primary_collision_analysis.ipynb) and [alias-alias_collision_analysis](./alias-alias_collisions/alias-alias_collision_analysis.ipynb) Jupyter notebooks show the analyses to get these values

![collision_graphic](https://github.com/cancervariants/gene-harmony-analysis/assets/109570522/91425d67-0884-4fbc-83ab-e7cfd8bd57bd)

## Purpose

The difficulties associated with resolving ambiguity and ensuring accurate understanding of gene symbols restrict the rate of clinical decision-making and contribute to confusion in gene knowledge aggregation. The gene nomenclature system would be most effective if it is unambiguous with a tool to take existing knowledgebase entries as inputs to resolve. 

This curated collection of alias data will be a foundation for disambiguating gene symbols.

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
| 5 | symbol_capture_analysis                     | 4, 1                     | none                                   |   |
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

# How can you help?

Contributing information on collisions that you come across will help collect data on the collisions that would be most impactful to resolve as well as increasing the data available for developing resolution strategies for downstream tool development.

## Contact Information

For any feedback, questions, or conversation, please make an issue.
