# Orthology Inference

This is the repo for the Orthology Inference project. Its aim is to identify orthologous groups of proteins across 33 *Drosophila* genomes. Like many other ortholoy inference pipelines, we use a graph-based method of clustering reciprocal best hits from BLAST searches between pairs of genomes. The chief innovations, however, are 1) the extension of commonly used triangle clustering to its generalization of *k*-clique percolation and 2) the implementation of two phylo-HMMs to identity poorly supported regions and missing data in the alignments.

## Project Organization

At the highest level, this project is organized into the following components:

```
orthology_inference/
	├── analysis/
	├── bin/
	├── data/
	├── src/
	└── README.md
```

Only `analysis/` and `src/`, which together contain all code written for this project, are explicitly tracked by Git. `bin/` contains third-party programs or code used in this project. Though this directory is not tracked by Git, scripts may reference it by path, so it is included here for completeness. Similarly, `data/`, which contains all the raw data used in this project, is not tracked by Git.

`analysis/` contains only directories, which serve to group related analyses. Some directories are "orphaned" and no longer contribute to any recent or ongoing analyses, but are included here for completeness. Currently, it contains the following entries:
- `brownian/`: Application of Brownian motion model to orthologs from NCBI annotations, among other phylogenetic analyses
- `evosim/`: Simulations of evolution of alignments from reconstructed ancestors
- `GOpred/`: Prediction of GO terms associated with proteins using rates of feature evolution
- `ortho_cluster1/`: Construction of orthologous groups using all *Drosophila* genome annotations obtained from [NCBI](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/all/)
- `ortho_cluster2/`: Final set of orthologous groups removing genomes which clustered poorly in `ortho_cluster1/`
- `ortho_MSA/`: Creation of multiple sequence alignments from orthologous groups generated in `ortho_cluster2/`
- `ortho_tree/`: Set of orthologous groups using all genomes in `ortho_cluster2/` and an outgroup *S. lebanonensis*; this directory exists solely to calculate a phylogenetic tree for the species in `ortho_cluster2/`
- `ortho_search/`: Scripts to run and parse BLAST searches for all genomes used in `ortho_cluster1/`
- `TF_CF_ids/`: Analyses to parse and deduplicate lists of transcription factors and cofactors found in [Stampfel *et al.*](https://pubmed.ncbi.nlm.nih.gov/26550828/) and [Hens *et al.*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3929264/)

## Dependencies
Nearly all code is written in Python and run with version 3.9.12. The remaining code is written in Bash shell scripts. The following Python libraries were used.

|Library|Version|
|---|---|
|homomorph|0.3.0|
|matplotlib|3.5.1|
|NetworkX|2.7.1|
|NumPy|1.22.3|
|pandas|1.4.1|
|SciPy|1.9.1|
|scikit-bio|0.5.7|
|TensorFlow|2.7.0|

Scikit-Bio attempts to import some deprecated warnings from scipy.stats during import, so these lines were commented out to ensure compatibility.