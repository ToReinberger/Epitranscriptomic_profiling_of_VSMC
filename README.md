# Epitranscriptomic profiling of VSMC phenotypes

## Project overview

This repository contains the Python scripts used for data processing, statistical analysis, and generation of the tables and figures presented in the study “Epitranscriptomic profiling of VSMC phenotypes reveals uridine modifications linked to post-transcriptional regulation.”

The study investigates epitranscriptomic changes associated with phenotypic transitions of human vascular smooth muscle cells (VSMCs), with a particular focus on uridine-centered RNA modifications and their potential involvement in post-transcriptional regulation. The provided scripts process the underlying experimental data and generate the corresponding processed datasets, summary tables, statistical results, and graphical outputs used in the manuscript and supplementary material.

This repository is intended to support transparency and reproducibility of the computational analyses performed in the study.

---

## Analysis scripts

The script names are prefixed with the number of the corresponding main figure in the publication. Scripts beginning with `proc_` perform general preprocessing used by several downstream analyses.

---

## Figure 1: Characterization of VSMC phenotypes

#### `Fig1_plot_pca_clustermap.py`

1. **Input:** Gene-level Illumina RNA-seq expression matrices and differential-expression result tables for the five VSMC conditions.
2. **Process:** Filters genes according to expression and variability, performs principal-component analysis and hierarchical clustering, and identifies gene clusters with condition-specific expression patterns.
3. **Output:** PCA plots, clustered expression heatmaps, gene-cluster plots, and gene lists used for functional enrichment analysis.

#### `Fig1_plot_enrichment_results.py`

1. **Input:** Functional-enrichment result tables generated for differentially expressed genes or condition-associated gene clusters.
2. **Process:** Filters and ranks enriched terms from MSigDB, KEGG, Reactome, and Gene Ontology, removes redundant terms, and formats pathway names for visualization.
3. **Output:** Processed enrichment tables and bar plots summarizing the biological pathways associated with the VSMC phenotypes.

---

## Figure 2: Distribution and sequence context of RNA modifications

#### `2_plot_DMR_distribution_CLEAN.py`

1. **Input:** Filtered xPore results annotated with the modified central base and relative transcript position.
2. **Process:** Assigns modification sites to the 5′-UTR, coding sequence, or 3′-UTR and calculates their distribution separately for A-, C-, G-, and U-centered sites.
3. **Output:** Stacked histograms and density plots showing the transcript-region distribution of differentially modified positions.

#### `2_plot_DMR_motifs_C3_C4(1).py`

1. **Input:** Filtered xPore DMR tables containing 5-mer sequences, modification assignments, and C3-versus-C4 effect directions.
2. **Process:** Groups sites by central nucleotide, modification assignment, and direction of change and calculates position-specific nucleotide frequencies.
3. **Output:** Sequence-motif visualizations for modification-associated 5-mers detected between C3 and C4.

#### `2_check_dmr_kmer_enrichment.py`

1. **Input:** Filtered xPore DMRs and transcript sequences used to determine the background frequency of all possible RNA 5-mers.
2. **Process:** Counts DMR-associated and background occurrences of each 5-mer, tests enrichment using contingency-based statistics, and applies false-discovery-rate correction.
3. **Output:** Tab-separated k-mer enrichment tables, intermediate JSON count files, and enrichment plots showing significantly overrepresented motifs.

#### `2_plot_GYUYY_motif_density.py`

1. **Input:** Significant C3-versus-C4 DMRs annotated with their 5-mer sequence and replicate-specific modification rates.
2. **Process:** Selects sites matching the GYUYY motif and compares their modification-rate distributions between C3 and C4.
3. **Output:** Density plots summarizing condition-specific modification rates at GYUYY-associated sites.

#### `2_plot_DMR_transcripts_per_gene.py`

1. **Input:** Gene-annotated xPore DMRs, transcript-level RNA-seq results, and transcript-to-gene annotation tables.
2. **Process:** Integrates DMRs with transcript expression and genomic annotation to examine how modification sites are distributed across alternative transcripts of individual genes.
3. **Output:** Transcript-level plots and summary visualizations showing DMR occurrence across transcripts and genes.

#### `2_plot_modRates_50_protein_level_counts.py`

1. **Input:** A gene-level master table integrating DMR counts, transcript properties, RNA abundance, and protein abundance.
2. **Process:** Calculates modification-site densities and examines relationships between DMR burden, transcript features, RNA expression, and protein abundance.
3. **Output:** Distribution, correlation, and exploratory multivariable plots describing gene-level characteristics associated with RNA modifications.

---

## Figure 3: Association of RNA modifications with gene and protein abundance

#### `3_plot_DMR_genes_overview_matrix.py`

1. **Input:** Selected genes and gene-annotated xPore DMRs containing transcript positions, central bases, effect directions, and modification assignments.
2. **Process:** Summarizes the location and direction of A-, C-, G-, and U-centered DMRs across selected genes.
3. **Output:** Gene-by-modification overview matrices showing DMR positions within the 5′-UTR, coding sequence, and 3′-UTR.

#### `3_plot_Gene_Level_Protein_level.py`

1. **Input:** A gene-level master table containing DMR counts, RNA-seq results, proteomics results, transcript lengths, and abundance measurements.
2. **Process:** Compares RNA- and protein-level changes with DMR occurrence and evaluates correlations for genes with different modification burdens.
3. **Output:** Scatter plots, correlation plots, and gene-level summary figures linking RNA modification patterns to RNA and protein abundance.

#### `3_run_and_plot_enrichment.py`

1. **Input:** Genes selected according to their DMR burden and a background set of expressed genes.
2. **Process:** Performs Enrichr-based functional enrichment using the expressed genes as the analysis background and selects representative terms from relevant gene-set libraries.
3. **Output:** Excel tables containing enrichment results and graphical summaries of pathways associated with DMR-affected genes.

#### `3_multivariat_analyzes_CLEAN.py`

1. **Input:** The integrated gene-level DMR, RNA-expression, and protein-abundance master table.
2. **Process:** Uses regression and multivariable models to test whether DMR occurrence is associated with differential RNA or protein abundance while accounting for transcript number, transcript length, and baseline abundance.
3. **Output:** Model coefficients, odds ratios, confidence intervals, significance values, and variable-importance summaries.

---

## Figure 4: Expression of genes affected by RNA modifications

#### `4_plot_pca_clustermap_ViolinPlots_RNAmod_genes.py`

1. **Input:** Illumina RNA-seq expression matrices, differential-expression results, and gene lists defined by RNA-modification status.
2. **Process:** Filters and normalizes expression data, performs clustering, and compares the expression of DMR-associated genes across the VSMC conditions.
3. **Output:** PCA plots, clustered heatmaps, and violin plots showing condition-dependent expression patterns of RNA-modification-associated genes.

---

## Figure 5: Poly(A)-tail analyses

#### `5_get_polyA_peaks.py`

1. **Input:** Read-level poly(A)-tail length estimates from Nanopore direct RNA sequencing for each sample and transcript.
2. **Process:** Filters reads by quality, estimates the dominant poly(A)-tail-length peak for each transcript using kernel-density estimation, combines sample-level results, and calculates transcript- and gene-level summaries.
3. **Output:** Per-sample poly(A)-peak tables, a combined transcript-level matrix, and gene-level poly(A)-tail summary tables.

#### `5_transcriptID2gene_names_biomart.py`

1. **Input:** Tables containing Ensembl transcript identifiers.
2. **Process:** Queries Ensembl BioMart in batches to map transcript identifiers to Ensembl gene identifiers and external gene names.
3. **Output:** Transcript-to-gene annotation tables and input datasets supplemented with gene names.

#### `5_plot_logFC_PolyA_diff.py`

1. **Input:** Gene-level Poly(A)-tail differences combined with the DMR, RNA-expression, and proteomics master table.
2. **Process:** Selects genes according to DMR burden and differential abundance and tests correlations between Poly(A)-tail changes and RNA- or protein-level fold changes.
3. **Output:** Regression and correlation plots comparing Poly(A)-tail-length differences with RNA and protein abundance changes.

---

## Figure 6: RNA secondary-structure analyses

#### `6_get_RNAfold_params_and_plot.py`

1. **Input:** Filtered xPore modification sites, Ensembl transcript sequences, and sequence windows surrounding the selected DMR positions.
2. **Process:** Extracts local transcript sequences, predicts RNA secondary structures using ViennaRNA and related folding tools, calculates structural energies, assigns nucleotides to structural elements using `forgi`, and compares alternative structures generated from shifted sequence windows.
3. **Output:** Per-site RNA-folding parameter tables, predicted dot-bracket structures, structural-element annotations, similarity matrices, minimum-free-energy profiles, and RNA secondary-structure plots.

#### `6_run_RNA_element_enrichment.py`

1. **Input:** RNA-folding result tables containing DMR positions, predicted `forgi` structural-element annotations, transcript-region assignments, central nucleotides, and xPore statistics.
2. **Process:** Determines the structural-element pattern surrounding each DMR, counts its occurrence within the corresponding transcript-sequence background, and tests whether specific RNA structural contexts are enriched at modification sites using Fisher’s exact and hypergeometric tests with multiple-testing correction.
3. **Output:** Summary tables containing structural-element frequencies, enrichment ratios, odds ratios, raw and adjusted significance values, and motif-specific enrichment results for coding sequences and 3′-UTRs.

#### `6_plot_RNAfold_enrichment.py`

1. **Input:** Summary tables generated by `6_run_RNA_element_enrichment.py`, including structural-element frequencies, enrichment odds ratios, and adjusted significance values.
2. **Process:** Selects significantly enriched RNA structural motifs and compares their frequencies and enrichment strengths across transcript regions and central nucleotide classes.
3. **Output:** Bar plots showing the proportion and over-representation of DMRs within specific RNA secondary-structure elements.

#### `6_plot_rnafold_from_RNAplfold_WSL_adjusted.py`

1. **Input:** Transcript FASTA sequences, RNAplfold base-pairing and unpaired-probability outputs, selected DMR coordinates, and predicted RNA structural elements.
2. **Process:** Parses local base-pair probabilities, selects compatible RNA base pairs, constructs secondary-structure graphs, assigns paired or unpaired probabilities to individual nucleotides, and highlights the sequence context surrounding selected modification sites.
3. **Output:** RNA secondary-structure and arc-diagram plots displaying local base pairing, structural accessibility, nucleotide positions, and selected DMR sites.

---

## Figure 7: miRNA binding and RNA-structure analyses

#### `7_map_miRNA_binding_sites.py`

1. **Input:** Filtered xPore DMRs, transcript sequences, TargetScan miRNA-family information, and small-RNA-seq expression results.
2. **Process:** Extracts sequence windows surrounding DMRs, searches for complementary miRNA seed sequences, maps candidate miRNAs to modification sites, and integrates miRNA-expression information.
3. **Output:** Tables of predicted miRNA-binding sites overlapping or surrounding DMRs, including candidate miRNAs and their expression characteristics.

#### `7_LitSearch_miRNA_genes.py`

1. **Input:** Candidate miRNAs, genes, authors, or predefined literature-search terms.
2. **Process:** Queries PubMed and retrieves publication metadata, abstracts, identifiers, and related information for literature-based assessment of candidate regulatory interactions.
3. **Output:** Structured Excel, XML, or JSON files containing publications associated with the selected miRNAs or genes.

#### `7_plot_rnafold_from_RNAplfold_WSL_adjusted.py`

1. **Input:** Transcript sequences, selected DMR positions, candidate miRNA sequences, and RNAplfold or RNAfold base-pairing results.
2. **Process:** Parses predicted RNA structures, maps base-pair probabilities, aligns miRNA and target sequences, and highlights DMRs and predicted miRNA-interaction sites.
3. **Output:** RNA secondary-structure network plots showing local structural accessibility, DMR positions, and predicted miRNA–target interactions.

---

## Data-processing workflow

The `proc_0_*` scripts prepare and integrate the xPore, transcript-annotation, RNA-expression, Poly(A)-tail, and proteomics datasets used by the downstream figure scripts. They are intended to be run approximately in numerical order.

### `proc_0_0_get_UTR_cDNA_GeneName_from_Biomart.py`

1. **Input:** xPore result tables containing Ensembl transcript identifiers.
2. **Process:** Queries Ensembl BioMart to retrieve transcript-level cDNA sequences, transcript lengths, genomic coordinates, 5′-UTR and 3′-UTR boundaries, Ensembl gene identifiers, and external gene names. The retrieved annotation can then be merged with the xPore results.
3. **Output:** Transcript annotation tables containing cDNA sequences, UTR coordinates, transcript lengths, gene identifiers, and gene names, as well as annotated xPore result tables.

### `proc_0_1_majority_direction_voting.py`

1. **Input:** Raw xPore differential-modification results containing a `higher` or `lower` modification assignment for each tested 5-mer position.
2. **Process:** Counts the two possible modification assignments for each 5-mer and determines the predominant assignment. Positions whose assignment does not agree with the majority assignment of the corresponding 5-mer are removed.
3. **Output:** A filtered xPore table in which each 5-mer is represented by a consistent predominant modification direction. The script can additionally generate diagnostic distributions of the `higher`-to-`lower` assignment ratio.

### `proc_0_2_FILTER_Xpore_data.py`

1. **Input:** xPore results after majority-direction filtering, including replicate-specific modification rates, coverage values, k-mer identities, effect sizes, and statistical test results.
2. **Process:** Applies the selected quality criteria, including missing-value, replicate-support, coverage, modification-rate difference, and adjusted-p-value filters. It also standardizes RNA k-mers by converting thymidine to uridine.
3. **Output:** High-confidence xPore DMR tables used as the basis for downstream annotation, integration, and figure generation.

### `proc_0_2_add_gene_name_to_transcript_id.py`

1. **Input:** Filtered xPore tables containing Ensembl transcript identifiers and a local transcript annotation reference table.
2. **Process:** Removes transcript-version suffixes and maps transcript identifiers to transcript names and Ensembl gene identifiers. Missing annotations are queried individually through the Ensembl REST API.
3. **Output:** xPore result tables supplemented with transcript names and gene identifiers.

### `proc_transcriptID2gene_names_biomart.py`

1. **Input:** An xPore table containing Ensembl transcript identifiers.
2. **Process:** Queries Ensembl BioMart in batches to map transcript identifiers to Ensembl gene identifiers, external gene names, and gene-name sources.
3. **Output:** A reusable transcript-to-gene mapping table and an xPore table supplemented with gene-level annotation.

### `proc_0_3_map_genome_cdna2xpore.py`

1. **Input:** A filtered xPore table and a transcript-reference table containing genomic coordinates, cDNA positions, UTR boundaries, transcript sequences, gene identifiers, and gene names.
2. **Process:** Standardizes transcript identifiers and joins the transcript-reference annotation to each xPore position.
3. **Output:** An annotated xPore table containing genomic, transcript, cDNA, UTR, sequence, and gene-level information for each tested modification site.

### `proc_0_4_CREATE_MASTER_TABLE_GENE_level.py`

1. **Input:** Gene-annotated xPore DMRs, ONT transcript- and gene-level expression tables, ONT differential-expression results, proteomics abundance and differential-abundance tables, and HGNC cross-reference information.
2. **Process:** Aggregates DMRs from individual transcript positions to the gene level, calculates total and nucleotide-specific DMR counts, counts affected transcripts and expressed isoforms, and integrates these measures with ONT gene expression, differential-expression statistics, protein abundance, differential-protein statistics, transcript lengths, and identifier mappings.
3. **Output:** A gene-level master table containing DMR burden, transcript characteristics, RNA abundance, RNA fold changes, protein abundance, protein fold changes, significance values, and gene/protein identifiers. This table is used by several Figure 2 and Figure 3 analyses.

### `proc_0_5_CREATE_MASTER_TABLE_Transcript_level.py`

1. **Input:** Transcript-annotated xPore DMRs, Illumina transcript-expression data, ONT transcript-expression data, transcript-level Poly(A)-tail summaries, and UTR annotations.
2. **Process:** Aggregates DMR counts and modification rates per transcript, separately summarizes `higher` and `lower` modification assignments and 3′-UTR-associated sites, and integrates these measures with Illumina and ONT RNA abundance, transcript length, UTR length, gene annotation, and Poly(A)-tail length.
3. **Output:** TSV and Excel transcript-level master tables containing DMR counts, nucleotide-specific DMR counts, median modification rates, condition differences, RNA abundance, transcript features, and Poly(A)-tail measurements.

### Recommended processing order

```text
Raw xPore output
    |
    v
proc_0_1_majority_direction_voting.py
    |
    v
proc_0_2_FILTER_Xpore_data.py
    |
    +--> proc_0_0_get_UTR_cDNA_GeneName_from_Biomart.py
    |
    +--> proc_0_2_add_gene_name_to_transcript_id.py
    |        or
    +--> proc_transcriptID2gene_names_biomart.py
    |
    v
proc_0_3_map_genome_cdna2xpore.py
    |
    +--> proc_0_4_CREATE_MASTER_TABLE_GENE_level.py
    |
    +--> proc_0_5_CREATE_MASTER_TABLE_Transcript_level.py
    |
    v
Figure-specific analysis and plotting scripts
```

---

## Input data

The analyses use processed outputs from:

- Illumina bulk RNA sequencing
- Oxford Nanopore direct RNA sequencing
- xPore differential-modification analysis
- Nanocount and DESeq2
- Poly(A)-tail estimation
- Quantitative proteomics
- Ensembl, HGNC and TargetScan reference resources

Raw sequencing data are not included in this repository. Data availability
and accession information are provided in the associated publication.

Before running the scripts, update the input and output paths defined near the
bottom of each script. The original analysis was performed using a
project-specific directory structure.

The scripts reproduce the processed data, statistical analyses and figure
panels used in the publication. Some scripts contain parameters selected for
the final manuscript analyses, including DMR thresholds, minimum coverage,
adjusted-p-value cutoffs, Poly(A)-tail filtering criteria and RNA-folding
window sizes.

Results may differ slightly across package or reference-database versions.
The original analysis used the Ensembl release and software versions specified
in the accompanying environment file.

---

## Requirements

Python ≥ 3.10

Main dependencies:
* pandas
* numpy
* scipy
* statsmodels
* scikit-learn
* matplotlib
* seaborn
* pybiomart
* ViennaRNA
* forgi
* networkx

``` pip install -r requirements.txt ```

---

When using this code, please cite:

Reinberger T et al. Epitranscriptomic profiling of VSMC
phenotypes reveals uridine modifications linked to post-transcriptional
regulation. [Journal, year, DOI].


## License

This repository is distributed under the [MIT/GPL-3.0] License.

## Contact

For questions concerning the analysis, please contact:

Tobias Reinberger  
Institute for Cardiogenetics,\
Universität zu Lübeck \
tobias.reinberger@uni-luebeck.de

Inken Wohlers  \
Research Center Borstel  \
iwohlers@fz-borstel.de


