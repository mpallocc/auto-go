# Output Structure

By default, the output folder for the autoGO results is named `results` (unless specified otherwise using the `where_results` and `outfolder` parameters).

`where_results` is the working directory where the output folder will be created.

`outfolder` is the name of the folder in which you want to save results. 

---

Here’s a breakdown of the output structure:


```
results/
├── deseq_norm_data.txt            # Normalized counts from Deseq2
├── deseq_vst_data.txt             # VST-transformed normalized counts
├── CONTROL_vs_TREATMENT/          # A folder for each comparison performed
│   ├── all_res.tsv                # Differential expression analysis table
│   ├── filtered_DEGs/             # Folder with filtered results and plots
│   │   ├── volcano.png            # Volcano plot of DEGs
│   │   ├── up_genes/              # Up-regulated DEGs 
│   │   │   ├── enrichment_tables/ # GO enrichment data
│   │   │   ├── plots/             # GO enrichment plots 
│   │   ├── down_genes/            # Down-regulated DEGs 
│   │   │   ├── enrichment_tables/ # GO enrichment data
│   │   │   ├── plots/             # GO enrichment plots 
│   │   ├── up_down_genes/         # All DEGs
│   │   │   ├── enrichment_tables/ # GO enrichment data
│   │   │   ├── plots/             # GO enrichment plots 
│   ├── ssgsea/                    # Results from ssGSEA analysis
│   │   ├── tables/                # ssGSEA analysis tables
│   │   ├── plots/                 # ssGSEA related plots
├── ComparisonsHeatmap/            # Heatmaps for GO enrichment of multiple comparisons
  │   ├── up_genes/                 
  │   ├── down_genes/ 
  │   ├── up_down_genes/ 
```


**Explanation of the Outputs:**
  
- Normalized Counts: 
   - `deseq_norm_data.txt`: Contains normalized counts from the DESeq2 analysis.
   - `deseq_vst_data.txt`: Contains variance-stabilized transformed (VST) counts from DESeq2.


- Comparison Folders: 
   For each comparison (e.g., `CONTROL_vs_TREATMENT`), a dedicated folder is created, which includes:
     - `all_res.tsv`: A table with the complete results of the differential expression analysis.
     - `volcano.png`: A volcano plot.
     - `filtered_DEGs/`: A folder containing DEGs filtered according to log2FC and p-value defined thresholds. It also includes visualizations such as volcano plots.
       - `up_genes/`: Results and visualizations for up-regulated DEGs.
       - `down_genes/`: Results and visualizations for down-regulated DEGs.
       - `up_down_genes/`: Results and visualizations for both up and down-regulated genes.


- Enrichment Results:
  Each DEG subset (up_genes, down_genes, up_down_genes) has its own:
  - enrichment_tables/: Folder containing GO Enrichment results.
  - plots/: Folder containing visualizations (barplots and lolliplots)


- Comparisons Heatmap:
  - The `ComparisonsHeatmap/` folder contains heatmaps summarizing GO enrichment results across multiple comparisons. Subfolders are organized by DEG type (`up_genes`, `down_genes`, `up_down_genes`).


- ssGSEA: 
  - If ssGSEA analysis is performed, results are saved in the `ssgsea/` folder, including results files and specific plots.

---

![](../develop/vignettes/imgs/tree-structure.png)
*Visual Representation of an Example of the Output Structure.*