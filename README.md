# vh_coexpression_networks

A network is a mathematical model to illustrate relationships between entities. It consist of a collection of nodes connected by links.

One relevant type of network in the field of biology are co-expression networks, where nodes are macromoecules (ARNs, proteins...), and interactions between them represent similar changes in their abundances throughout the time (expression profiles). These networks contain  information about global changes in the transcriptome or proteome in response to a perturbation, and the "similarity" between expression profiles (the existence or absence of connection and its weight) can be quantified by means of several methods whose study will be part of this work.

The main objective of this project is to systematize the construction and  study of Virus-Host co-expression networks based on proteomics data.

In order to achieve that, we have focused on three separate goals:

* Comparative of the performance of different network construction methods based on similarity.  
* Characterization of the structure, properties and dynamics of Virus-Host co-expression networks, with emphasis on the interconnection schema between Viral and Host networks and how it impacts the centrality of their respective nodes.
* Modular and functional analysis of the Host networks, to highlight some of the biological processes modulated by the viruses upon cell infection.

## Data Folder

The raw input data with protein relative abundances for the two viruses of study is located in the following files:
* Herpes Simplex Virus 1 ([data_HSV1_edgeR.txt](data/data_HSV1_edgeR.txt)) [1]
* Vaccinia ([data_VACV_edgeR_2reps.txt](data/data_VACV_edgeR_2reps.txt)) [2]

The data folder contains also the results of the pipeline to create and analyze the networks:
* [allNets.RData](data/allNets.RData): R data file containing all the networks built for the two viruses, by means of Pearson Correlation Coefficient, Graphical Gaussian Models and Concordance Coefficient [3]. 
* [allProteins.RData](data/allProteins.RData): Processed input data. 
* [HSV1_results_proportionality.RData](HSV1_results_proportionality.RData): Created network object(with the Concordance Coefficient), nodes properties, interactions and other relevant data frames for the HSV1 virus. 
* [VACV_results_proportionality.RData](VACV_results_proportionality.RData): Created network object(with the Concordance Coefficient), nodes properties, interactions and other relevant data frames for the VACV virus.

These files are created by the coexpr_net_main.R script and they are imported in coexpr_net_pictures.R.

Node properties and interactions have been saved as csv  files that can be imported in other network analysis tools such as Gephi, Cytoscape or Neo4j, for further analysis. Community-based functional enrichment results are stored in csv files as well.

## Code Folder

The code folder contains the following R scripts:
* [coexpr_net_functions.R](data/coexpr_net_functions.R): Helper functions. 
* [coexpr_net_main.R](data/coexpr_net_main.R): This script processes the input data, creates and analyzes the Virus-Host co-expression network, for the virus and method of choice. 
    * Some of the variables that can be adjusted before running the script are:
        * *virus*: one of "VACV", "HSV1"
        * *corr_method*: one of "pearson", "proportionality"(for the Concordance Coefficient), "gaussian"
        * *FDR_threshold*: to extract the significant differerentially expressed proteins 
        * *LogFC_threshold*: to extract the significant differerentially expressed proteins
        * *enrichment_threshold*: to filter functional annotation results
        * *lambda*: Penalty parameter in Gaussian Graphical Model, larger lambda results 
                        in faster calculation and a sparser graph

    * The outputs produced are the following:
        * Inline basic network properties and pictures
        * Functional annotation of communities in a csv file
        * Node properties and relations in csv format
        * Network object exported in .R file

* [coexpr_net_pictures.R](data/coexpr_net_pictures.R): This script produces several visualizations of network properties.

Information about the R packages versions required to reproduce this work can be found in [sessionInfo.txt](sessionInfo.txt).

## How to leverage this code to analyze a new dataset from a different virus. 

1. Modify the loadProteins function inside [coexpr_net_functions.R](data/coexpr_net_functions.R) to create a data frame from the new input data. The data frame must have the following columns:
UniProtID, GeneID, Origin(Virus|Host), time_sample_1,..., time_sample_n, FDR_1,..., FDR_n, logFC_1,...,logFC_n

2. In same places of the code, variables take some value or another depending on the virus analyzed. This is done by means of R switch functions and they need to be adjusted shall a new virus be introduced. 

## References

[1]  Soh, T. K., Davies, C., Muenzner, J., Hunter, L. M., Barrow, H. G., Connor, V., Bouton, C. R., Smith, C., Emmott, E., Antrobus, R., Graham, S. C., Weekes, M. P., & Crump, C. M. (2020). Temporal Proteomic Analysis of Herpes Simplex Virus 1 Infection Reveals Cell-Surface Remodeling via pUL56-Mediated GOPC Degradation. Cell reports, 33(1), 108235. https://doi.org/10.1016/j.celrep.2020.108235

[2] Soday, L., Lu, Y., Albarnaz, J. D., Davies, C., Antrobus, R., Smith, G. L., & Weekes, M. P. (2019). Quantitative Temporal Proteomic Analysis of Vaccinia Virus Infection Reveals Regulation of Histone Deacetylases by an Interferon Antagonist. Cell reports, 27(6), 1920–1933.e7. https://doi.org/10.1016/j.celrep.2019.04.042

[3] Erb, I., & Notredame, C. (2016). How should we measure proportionality on relative gene expression data?. Theory in biosciences, 135(1-2), 21–36. https://doi.org/10.1007/s12064-015-0220-8



