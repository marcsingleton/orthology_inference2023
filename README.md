# IDREvoDevo

This is the repo for the IDREvoDevo project. Its aim is to identify evolutionarily conserved properties of intrinsically disordered regions (IDRs) in specific proteins. Ideally IDRs across different proteins would share similar patterns of conserved properties, allowing them to be grouped into classes similar to the classification systems used for structured domains, such as Pfam. These analyses focus specifically on proteins found in genomes within the genus *Drosophila*.

## Project Organization

At the highest level, this project is organized into the following components:

```
IDREvoDevo/
	├── analysis/
	├── bin/
	├── data/
	├── src/
	└── README.md
```

Only `analysis/` and `src/`, which together contain all code written for this project, are explicitly tracked by Git. `bin/` contains third-party programs or code used in this project. Though this directory is not tracked by Git, scripts may reference it by path, so it is included here for completeness. Similarly, `data/`, which contains all the raw data used in this project, is not tracked by Git.

`analysis/` contains only directories, which serve to group related analyses. Some directories are "orphaned" and no longer contribute to any recent or ongoing analyses, but are included here for completeness. Currently it contains the following entries:
- `brownian/`: Application of Brownian motion model to EggNOGv5 orthologs segmented by IUPRED
- `EggNOGv5_validation/`: Analyses to validate *Drosophila* orthologs downloaded from the [EggNOG v5.0](http://eggnog5.embl.de/) database
- `ortho_cluster1`: Initial attempt at constructing orthologous groups using genome annotations largely from [Yang *et al.*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6305970/)
- `ortho_cluster2`: Construction of orthologous groups using all *Drosophila* genome annotations obtained from [NCBI](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/all/)
- `ortho_cluster3`: Final set of orthologous groups removing genomes which clustered poorly in `ortho_cluster2/`
- `ortho_MSA/`: Creation of multiple sequence alignments from orthologous groups generated in `ortho_cluster3/`
- `ortho_tree/`: Set of orthologous groups using all genomes in `ortho_cluster3/` and an outgroup *S. lebanonensis*; this directory exists solely to calculate a phylogenetic tree for the species in `ortho_cluster3/`
- `ortho_search/`: Scripts to run and parse BLAST searches for all genomes used in `ortho_cluster2/`
- `seg_condiv/`: Clustering of features calculated from diverged segments in EggNOGv5 orthologs
- `seg_idr`: Clustering of features calculated from disordered segments identified by IUPRED in EggNOGv5 orthologs
- `TF_CF_lists`: Analyses to parse and deduplicate lists of transcription factors and cofactors found in [Stampfel *et al.*](https://pubmed.ncbi.nlm.nih.gov/26550828/) and [Hens *et al.*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3929264/)
