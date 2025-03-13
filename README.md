# Auto-GO: Reproducible, Robust and High Quality Ontology Enrichment Visualizations

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/autoGO)](https://cran.r-project.org/package=autoGO) 
[![](https://cranlogs.r-pkg.org/badges/autoGO)](https://CRAN.R-project.org/package=autoGO)
[![DOI](https://img.shields.io/badge/DOI-10.1109%2FBIBM58861.2023.10385722-red)](https://doi.org/10.1109/BIBM58861.2023.10385722)

---

Auto-GO is a framework that enables automated, high quality Gene Ontology enrichment analysis visualizations. It also features a handy wrapper for Differential Expression analysis around the `DESeq2`. The whole framework is structured in different, independent functions, in order to let the user decide which steps of the analysis to perform and which plot to produce.


## Main Features

#### Differential Expression Analysis:
- `deseq_analysis()`: A wrapper around `DESeq2` for performing differential expression analysis.
- `volcanoplot()`: Generates volcano plots to visualize DEGs.

#### Enrichment Analysis:
- `autoGO()`: Automates Gene Ontology enrichment analysis.

#### Visualization:
- `lolliGO()`: Creates lollipop plots for GO enrichment.
- `barplotGO()`: Creates barplots for GO enrichment.
- `heatmapGO()`: Creates heatmaps for GO enrichment.

#### Single Sample GSEA:
- `ssgseawrapper()`: A wrapper for single-sample GSEA.




## Installation

Install the released version of autoGO from CRAN:

```R
# Install from CRAN
install.packages("autoGO")
```

Or install the development version from GitHub with:

```R
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install AutoGO from GitHub
devtools::install_github("mpallocc/auto-go")
```


## Usage

For a step-by-step tutorial, check out the vignettes:
  
```R
browseVignettes("autoGO")
```

## Output Structure

The figure below illustrates the output structure of the autoGO package.

![](../tree/develop/vignettes/imgs/tree-structure.png)

By default, the output folder is named `results` (unless specified otherwise using the `where_results` and `outfolder` parameters).

You can find a detailed description of the output folders in this ![file](../blob/develop/output-structure.md).




---



## Citation

If you use AutoGO in your research, please cite the following publication:
I. Grassucci et al., "Auto-GO: Reproducible, Robust and High Quality Ontology Enrichment Visualizations," 2023 IEEE International Conference on Bioinformatics and Biomedicine (BIBM), Istanbul, Turkiye, 2023, pp. 3638-3641, doi: 10.1109/BIBM58861.2023.10385722. 


  
  
  
  
  
  