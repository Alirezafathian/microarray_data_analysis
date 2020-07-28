Microarray Data Analysis
============================
R scripts for analysis of microarray data.

Microarrays data used to measure gene expression to compare expression of a set of genes from a cell maintained in a particular condition to the same set of genes from a reference cell maintained under normal conditions [[1]](#1).
We used microarray data to compare the expression of two different sets of genes from different cells maintained in particular conditions.

## Example

we used GEO dataset <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52509">GSE13601</a> as an example. In this dataset there are expression data from cigarette smoke-treated mice at 4 and 6 months of age [[1]](#1).


## Requirements
limma,
Biobase,
GEOquery,
pheatmap,
ggplot2,
plyr,
pheatmap	

## References
<a id="2">[1]</a>
John-Schuster G, Hager K, Conlon TM, Irmler M et al. Cigarette smoke-induced iBALT mediates macrophage activation in a B cell-dependent manner in COPD. Am J Physiol Lung Cell Mol Physiol 2014 Nov 1;307(9):L692-706.
