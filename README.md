# Metabolite correlation networks reveal complex phenotypes of adaptation-driving mutations

This repository contains all the data and source code used to produce the results presented in the paper "Metabolite correlation networks reveal complex phenotypes of adaptation-driving mutations". 

## About

This paper is a result of a collaboration between the Bershtein Lab and the Pilosof Lab. We set out to understand how mutations reshape the organization of metabolism and, in turn, affect cellular fitness. Such insight is essential for linking molecular changes to systems-level metabolic adaptation. To address these questions, we engineered a metabolically suboptimal _E._ _coli_ strain by replacing the native _metK_ gene with an ortholog from _U._ _urealyticum_ and subjected it to adaptive laboratory evolution. Using theory and methodology inspired by Network Ecology, we constructed correlation-based metabolite interaction networks from untargeted LC–MS data (2,118 metabolites) to explore how adaptive mutations reorganize metabolic connectivity.

## Abstract

Living organisms are organized into hierarchical levels of increasing complexity. As mutational effects propagate through these levels, multiple phenotypes are produced. However, identifying which phenotypes influence fitness and drive selection remains a key challenge in understanding genotype-phenotype-fitness relationships. Here, we demonstrate that mutation-induced structural changes in metabolite correlation networks—an organizational level with emergent properties shaped by the cumulative effects of multiple biochemical reactions and regulatory mechanisms—constitute complex phenotypes linking mutations to bacterial fitness. We engineered a metabolically suboptimal _E._ _coli_ strain by replacing the _metK_ gene encoding methionine adenosyltransferase (MAT) with an ortholog from _U._ _urealyticum_ and subjected it to laboratory evolution. Analysis of correlation networks constructed from 2,118 untargeted metabolites revealed that adaptive mutations enhanced the strain’s fitness by reducing network density, removing node hubs, and extensively rewiring connectivity. These changes yielded smaller, more cohesive, and better-interconnected network clusters. Moreover, evolution shifted the node representing S-adenosylmethionine (SAM), the product of the MAT-catalyzed reaction, from a peripheral role to a key connector between network clusters with a significantly increased betweenness. None of the accumulated mutations potentially driving this transition directly influenced MAT activity or SAM metabolism, indicating that the shift in SAM’s role is an adaptive phenotype emerging from the metabolic system’s complexity. Targeted metabolomics of key nodes displaying similar network transitions unveiled additional metabolites involved in SAM-related pathways. We propose the construction and analysis of metabolite correlation networks as an experimental and analytical framework for mapping genotype-phenotype-fitness relationships and exploring the mechanisms of metabolic adaptation.

## Folders and files description

#### 1. source_code.R

This is the main R script used to perform all analyses.
Each step of the workflow is documented with in-line comments for clarity and reproducibility.

#### 2. Data/

This folder contains all input data and intermediate results used in the analysis:

  - LC-MS_results.xlsx — LC–MS results used in this study. These data were generated using Compound Discoverer 3.3 SP2 (Thermo Fisher Scientific).

  - multi_networks_results/ — Contains precomputed Infomap modularity analysis results.
    You can reproduce these results using the provided script, but the precomputed files are included to save processing time.

  - PCLRC_filtered_edges/ — Contains CSV files of PCLRC-filtered network edges, generated with the R script from Di Cesare et al (2022).

## System requirements:

R Programming language: 4.2.1

The script requires the following R packages, with the used versions specified:

* tidyverse - 2.0.0
* readxl - 1.4.3
* Hmisc - 4.7-1
* ggplot2 - 3.5.1
* igraph - 1.5.1
* data.table - 1.14.2
* reshape2 - 1.4.4
* dplyr - 1.1.4
* influential - 2.2.9
* emln - 0.0.1
* infomapecology - 2.0.1
* magrittr - 2.0.3
* combinat - 0.0-8
* khroma - 1.12.0
* Infomap (MacOS version): Version 1.7.1

Ensure these packages are installed and loaded before running the script.

## Run Time and Computational Considerations
The script can be executed on the original dataset used in the analysis. However, two steps are computationally intensive and are recommended to be executed on an HPC cluster:

  - Step 6: Generating Shuffled Networks and Comparing Jaccard Indexes. In the original analysis, 1,000 shuffles were used. However, in the provided script, the number of shuffles has been reduced to 10 to improve time efficiency. This adjustment enables the step to run on a standard desktop computer but at the cost of reduced result resolution.

  - Step 8: Modularity Detection using the infomapecology Package. To save time, precomputed results are available and can be loaded directly. Users not interested in recomputing this step can skip its execution.

These parts of the script were tested on the HPC OS version- Oracle Linux Server 8.7.






