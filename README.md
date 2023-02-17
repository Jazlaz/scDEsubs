# scDEsubs
Package: scDEsubs
Version: 1.0
Date: 2021-09-22
Title: scDEsubs: an R package for flexible identification of
        differentially expressed subpathways using single cell RNA-seq 
        expression experiments
Author: Giannis Aslanis
Maintainer: Giannis Aslanis <giannisas84@gmail.com>
Description: scDEsubs is a network-based systems biology package that
        extracts disease-perturbed subpathways within a pathway network
        as recorded by RNA-seq experiments. It contains an extensive
        and customizable framework covering a broad range of operation
        modes at all stages of the subpathway analysis, enabling a
        case-specific approach. The operation modes refer to the
        pathway network construction and processing, the subpathway
        extraction, visualization and enrichment analysis with regard
        to various biological and pharmacological features. Its
        capabilities render it a tool-guide for both the modeler and
        experimentalist for the identification of more robust
        systems-level biomarkers for complex diseases.
Depends: R (>= 3.3), locfit
SystemRequirements:
License: GPL-3
Repository: Bioconductor
NeedsCompilation: no
LazyLoad: yes
Imports: graph, igraph, RBGL, circlize, limma, edgeR, MAST,
        stats, grDevices, graphics, pheatmap, utils, ggplot2,
        Matrix, jsonlite, tools, methods
Suggests: RUnit, BiocGenerics, knitr
VignetteBuilder: knitr
biocViews: SystemsBiology, GraphAndNetwork, Pathways, KEGG,
        GeneExpression, NetworkEnrichment, Network, RNASeq,
        DifferentialExpression, Normalization, ImmunoOncology
