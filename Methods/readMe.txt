This GitHub repository contains both the functional code and the analytical workflows used to implement sample-specific network inference methods on TCGA-BLCA data. A variety of computational approaches have been developed to infer sample-specific protein–protein interaction networks (ssPPIs) from transcriptomic data, which can be broadly categorized into two classes: reference-dependent methods (leveraging an external reference, typically normal tissue) and reference-independent methods (deriving networks without an explicit reference).

Reference-dependent methods

(1) SSN constructs individual-specific networks by quantifying differential Pearson correlation coefficients (ΔPCC) between a reference group and a perturbed network formed by adding a single sample. Because the original publication did not provide code, we developed our own Python implementation of SSN.

(2) iENA extends the SSN framework by incorporating higher-order edge-level statistics (e.g., variance and covariance of correlations) to identify significant perturbations. We used the R code released with the original iENA publication.

(3) SSP identifies disease-specific perturbations by comparing each sample’s gene-expression rank profile to a stable reference network derived from normal tissue samples. We used the R implementation provided by the SSP authors.

Reference-independent methods

(1) PPIXpress (PPIX) filters a global PPI network by retaining only interactions supported by gene expression above a data-driven threshold, typically estimated across a cohort. We used the Java implementation from the PPIXpress authors.

(2) LIONESS (LION) infers personalized biological networks through a leave-one-out strategy: it compares a full population network (using all samples) with networks reconstructed after excluding each sample in turn. The difference quantifies each sample’s contribution, enabling estimation of its specific network. We used the lionessR package to run this algorithm.

(3) CSN builds cell-specific networks from scRNA-seq data by assessing gene–gene dependency patterns within each cell while using the global distribution of all cells as background. We implemented our own R version of CSN based on the original MATLAB code.

(4) SWEET constructs single-sample networks by weighting samples according to expression similarity (via PCC) and quantifying sample-specific effects through differential analysis between population and perturbed networks. We used the Python implementation provided by the SWEET authors.
