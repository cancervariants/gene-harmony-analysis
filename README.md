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

The [total_alias_overlap](./total_alias_overlap.ipynb) Jupyter notebook shows the analysis to get these values (the notebook will be condensed and made more efficient)

![collision_graphic](https://github.com/cancervariants/gene-harmony-analysis/assets/109570522/91425d67-0884-4fbc-83ab-e7cfd8bd57bd)

## Purpose

The difficulties associated with resolving ambiguity and ensuring accurate understanding of gene symbols restrict the rate of clinical decision-making and contribute to confusion in gene knowledge aggregation. The gene nomenclature system would be most effective if it is unambiguous with a tool to take existing knowledgebase entries as inputs to resolve. 

This curated collection of collision data will be a foundation for disambiguating gene symbols.

# How can you help?

Contributing information on collisions that you come across will help collect data on the collisions that would be most impactful to resolve as well as increasing the data available for developing resolution strategies for downstream tool development.

1. In the collision records folders, there are collision records that are completed (but can always be updated) and a [blank sample record](./sample_collision_record.yaml) to use as a template.
2. The [contributing documentation](./CONTRIBUTING.md) explains the different features that are included in a collision record. 
3. To propose an update or make a new collision record create a personal fork to the repo. New collision records should be YAML files named after the collision in either the [alias-alias](./alias-alias_collision_records) or [alias-primary](./alias-primary_collision_records) collision records folders. Once created, the review process can be started with the creation of a pull request. 

## Contact Information

For any feedback, questions, or conversation, please make an issue.
