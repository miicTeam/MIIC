# MIIC
MIIC (Multivariate Information Inductive Causation) learns a large class of causal or non-causal graphical models from purely observational data, while including the effects of unobserved latent variables, commonly found in many datasets. Starting from a complete graph, the method iteratively removes dispensable edges, by uncovering significant information contributions from indirect paths, and assesses edge-specific confidences from randomization of available data. The remaining edges are then oriented based on the signature of causality in observational data. This approach can be applied on a wide range of datasets and provide new biological insights on regulatory networks from single cell expression data, genomic alterations during tumor development and co-evolving residues in protein structures.

## References
Verny L., Sella N., Affeldt S., Singh PP., Isambert H.; Learning causal networks with latent variables from multivariate information in genomic data;  PLoS Comput. Biol., 2017.

## Useful links
This code is included in the MIIC web server, available at [https://miic.curie.fr](https://miic.curie.fr)


## Getting Started
Download or clone the project 

## Prerequisites
In order to compile MIIC executables you need to have c++ compiler and pthreads installed in your machine. 
In order to run MIIC you need the R environment (see https://www.r-project.org/) for help and download. You also need to install some R packages: MASS, getopt, plotrix, methods, igraph, ppcor, bnlearn 

## Running the test
```
cd ../miic/common
Rscript miic.R -i ../data/alarm1000samples.txt -l ../data/alarmLayout.txt -t ../data/alarmTrueEdges.txt -s ../data/alarmStateOrder.tsv -o ../data/resultExec -n
```
## Documentation
You can find the complete documentation inside the "doc" folder.

## Authors
- Verny Louis
- Sella Nadir
- Séverine Affeldt
- Hervé Isambert

## License
These materials are provided under a Creative Commons License license.
