# 16s and ITS analysis with nanopore sequencing

## This repository contains scripts and workflows for analyzing 16S rRNA and ITS amplicon sequencing data generated using Oxford Nanopore Technologies (ONT) platforms. The analysis pipeline includes steps for quality control, taxonomic classification, and visualization of microbial communities.

## The differential abundance analysis was updated from traditional non-parametric tests (Mann-Whitney with relative abundance/rarefaction) to ANCOM-BC2 (Analysis of Compositions of Microbiomes with Bias Correction).

Why ANCOM-BC2?

Compositional Bias: Mathematically handles differences in sequencing depth across libraries.

Raw Data: Does not require and does not use rarefied data, preserving 100% of the biological information (crucial for samples with high read count variation, e.g., ~800k to 7.5M).

Rigorous FDR: Applies strict mathematical corrections (q-value) to avoid False Positives.

## Evironment Setup (Conda)
To avoid dependency conflicts between complex Python and R packages, this pipeline uses two distinct virtual environments.
### R Enviroment (Statistical Analysis)
```bash
# 1. Create the environment with the base packages
conda create -n env_ancombc -c conda-forge -c bioconda bioconductor-ancombc bioconductor-phyloseq r-tidyverse -y

# 2. Activate the environment
conda activate env_ancombc

# 3. Install the auxiliary dependency required by ANCOM-BC2
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cran.us.r-project.org"); BiocManager::install("microbiome", update=FALSE, ask=FALSE)'
```
### Python Enviroment (Data Processing)
Python Environment (Data Cleaning and Plots)
The standard project environment used in Jupyter Notebook (16s-nanopore or similar) containing:
- pandas
- matplotlib
- seaborn
- numpy

## 2. Pipeline execution
### Step 1: Quality Control
The first step in the analysis pipeline is to perform quality control on the raw sequencing data. This includes filtering out low-quality reads, trimming adapters, and removing chimeric sequences. `nanocomp` is a useful tool for assessing the quality of nanopore sequencing data, providing metrics such as read length distribution, quality scores, and yield. In next, filter reads based on quality scores and length thresholds to ensure that only high-quality reads are retained for downstream analysis, using `chopper` tools.
### Step 2: Taxonomic Classification
After quality control, the next step is to classify the reads into taxonomic groups. This can be done using tools like `kma`, which align the reads to reference databases and assign taxonomic labels based on sequence similarity. The output of this step is a table of taxonomic assignments for each read, which can be used for downstream analysis.
### Step 3: Differential Abundance Analysis
Once the taxonomic assignments are obtained, the next step is to perform differential abundance analysis to identify taxa that are significantly different between experimental groups. This can be done using ANCOM-BC2, which accounts for compositional bias and provides rigorous statistical testing. The output of this step is a list of taxa that are differentially abundant between groups, along with associated p-values and q-values.
### Step 4: Visualization
Finally, the results of the analysis can be visualized using various plotting libraries in Python or R. Common visualizations include bar plots, heatmaps, and ordination plots, which can help to interpret the results and identify patterns in the microbial communities.

## Important Notes
Rarefaction vs. Differential Abundance
Differential Abundance (ANCOM-BC2): Should NEVER use rarefied data. The R script strictly uses the clean matrix (otu_abundance_matrix_CLEAN.csv) with absolute counts (Raw Counts).

Alpha and Beta Diversity: MUST use rarefied data. For these ecological analyses, the matrix must be leveled by the sample with the lowest viable number of reads (rarefaction_depth = 816500, based on the CTL-2 sample) to ensure a fair comparison without excluding any samples.

### Expected Directory structure
```
raw_data/
    ├── 1_fastq/
    │   ├── sample.fastq
                ├── 2. renamed_fastq/
                            └── sample_renamed.fastq
                                        │                 
                                        downstream/
                                        ├── 1_qc_raw/
                                        │   └── nano_comp_report/
                                        ├── 2_trimmed/
                                        │   └── sample_trimmed.fastq
                                                ├── 3_filtered/
                                                │   └── chopper/
                                                ├── 4_kma_taxa/
                                                │   └── ncbi/
                                                │       └── otu_abundance_matrix_CLEAN.csv
                                                ├── 5_statistics/
                                                │   ├── kda/
                                                │   │   └── metadata_filtered.csv
                                                │   └── ancombc/
                                                │       ├── ancombc_vs_CTL.csv
                                                │       └── ancombc_pairwise_tratamentos.csv
                                                └── 6_plots/
                                                    ├── Gaphs/
                                                    └── Figures/
```
## References
1. Mandal, S., Van Treuren, W., White, R. A., Eggesbø, M., Knight, R., & Peddada, S. D. (2015). Analysis of composition of microbiomes: a novel method for studying microbial composition. Microbiome, 3(1), 1-11. https://doi.org/10.1186/s40168-015-0172-3
2. Lin, H., Peddada, S. D., & Chen, J. (2014). Analysis of compositions of microbiomes with bias correction. Nature Communications, 5(1), 1-10. https://doi.org/10.1038/ncomms4647
3. Lin, H., Peddada, S. D., & Chen, J. (2022). ANCOM-BC2: A method for differential abundance analysis of microbiome data with bias correction. Bioinformatics, 38(12), 3200-3208. https://doi.org/10.1093/bioinformatics/btac204
4. MCKINNEY, Wes. Data Structures for Statistical Computing in Python. In: PYTHON IN SCIENCE CONFERENCE. Anais... Austin, Texas: 2010. Disponível em: <https://doi.curvenote.com/10.25080/Majora-92bf1922-00a>. Acesso em: 24 jun. 2026 
