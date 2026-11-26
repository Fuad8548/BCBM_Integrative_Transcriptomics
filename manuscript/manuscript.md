# Molecular Adaptation and Tumour–Microenvironment Interactions in Breast Cancer Brain Metastasis: An Integrative Transcriptomic Investigation

## Abstract

### Background 
Breast cancer brain metastasis (BCBM) represents a complex metastatic process involving tumour adaptation and interactions with the brain microenvironment. The transcriptional programmes associated with this process and their relationship to the cellular composition of metastatic tissue remain incompletely characterized.   

### Methods
We performed an integrative transcriptomic investigation combining paired bulk gene-expression analyses of primary breast tumours and brain metastases with independent validation, pathway enrichment, and single-nucleus transcriptomic characterization of brain metastatic specimens. Differential expression analyses incorporated paired modelling and a purity-adjusted model. Cross-dataset integration was used to identify directionally consistent transcriptional changes. Single-nucleus data were used to characterize cellular composition, extracellular matrix (ECM)-associated gene expression, and inter-sample variation in ECM-associated transcriptional scores.   

### Results
Bulk transcriptomic analyses identified substantial expression differences between primary breast tumours and brain metastases. Metastatic samples exhibited increased expression of oxidative phosphorylation and respiratory metabolism-associated programmes alongside reduced expression of extracellular matrix organization, collagen-associated processes, and selected immune and cytokine-associated pathways. These patterns were supported by pathway enrichment and cross-dataset directional comparisons. Single-nucleus analysis revealed a heterogeneous, tumour-dominated cellular landscape with contributions from myeloid, lymphoid, stromal, vascular, pericyte, astrocytic, and other populations. ECM-associated transcripts were distributed across multiple cellular populations, with stromal cells showing prominent expression of several core ECM-associated genes. ECM-associated scores varied continuously across metastatic samples and were positively associated with stromal abundance. Samples with higher ECM-associated scores exhibited differences in cellular composition and expression of selected ECM-associated genes.   

### Conclusion
Our integrated analyses identify coordinated metabolic, ECM-associated, and microenvironmental transcriptional patterns in breast cancer brain metastases. The bulk and single-nucleus results suggest that these patterns reflect complementary aspects of metastatic tissue biology, including tumour-associated transcriptional adaptation and multicellular contributions to ECM-associated expression. These findings provide a transcriptomic framework for further investigation of the biological processes associated with breast cancer brain metastasis.

## 1. Introduction

### 1.1 Breast cancer brain metastasis and the cerebral microenvironment
Brain metastasis represents a clinically challenging manifestation of breast cancer, contributing substantially to morbidity and mortality among patients with advanced disease. The brain presents a distinctive metastatic environment characterized by specialized vascular structures, neural cell populations, unique metabolic demands, and complex interactions between tumour cells and resident or recruited cells. Successful metastatic colonization therefore requires more than dissemination from the primary tumour; it involves adaptation to the biological constraints and opportunities of the cerebral microenvironment.   

The metastatic process encompasses multiple stages, including tumour-cell dissemination, vascular transit, interactions with the blood–brain barrier, extravasation, survival within brain tissue, and subsequent metastatic growth. These processes involve reciprocal interactions between cancer cells and components of the brain microenvironment, including astrocytes, microglia, endothelial cells, pericytes, and extracellular matrix-associated structures. Understanding the transcriptional characteristics of breast cancer brain metastases may therefore provide insight into tumour adaptation and the cellular context in which metastatic lesions develop.

### 1.2 Metabolic adaptation in metastatic tumour cells
Metastatic progression is accompanied by changes in cellular metabolism that can support survival, proliferation, stress adaptation, and biosynthetic requirements. Although cancer metabolism is frequently discussed in the context of increased glycolysis, tumour metabolic phenotypes are not restricted to a single pathway. Mitochondrial metabolism, oxidative phosphorylation, electron transport, and alternative metabolic processes may contribute differently depending on tumour context, nutrient availability, oxygenation, and cellular state.

The brain is a metabolically active organ, and metastatic tumour cells growing within cerebral tissue may encounter distinct metabolic conditions compared with those in their primary tissue environment. Transcriptional alterations involving mitochondrial respiration and oxidative metabolism may consequently represent relevant features of metastatic adaptation. However, transcriptomic evidence of metabolic pathway enrichment should be distinguished from direct measurements of metabolic flux, mitochondrial activity, or energetic output.   

Characterizing metabolic-associated transcriptional programmes in brain metastases may therefore help identify pathways associated with metastatic tissue adaptation and generate hypotheses for subsequent functional investigation. 

### 1.3 Extracellular matrix remodelling and the metastatic microenvironment
The extracellular matrix (ECM) is a complex network of structural and regulatory molecules that contributes to tissue organization, cellular adhesion, migration, signalling, and mechanical properties. ECM-associated components, including collagens, glycoproteins, proteoglycans, and matrix-regulatory molecules, participate in interactions between tumour cells and their surrounding environment. Changes in ECM-associated transcription may reflect alterations in stromal composition, cellular activation states, tissue organization, and interactions between different cellular populations.   

In metastatic tissues, ECM-associated programmes are not necessarily restricted to tumour cells. Stromal cells, vascular-associated populations, and other components of the microenvironment may contribute to the expression of matrix-associated genes. Consequently, bulk tissue measurements may reflect both transcriptional changes within individual cell populations and variation in the relative abundance of those populations.   

Integrating bulk transcriptomic analyses with single-nucleus data may provide a more informative framework for examining ECM-associated changes, distinguishing tissue-level transcriptional patterns from their potential cellular distribution and evaluating heterogeneity between metastatic specimens.

### 1.4 Cellular heterogeneity of breast cancer brain metastases
Brain metastatic specimens contain heterogeneous cellular populations that may include malignant cells, immune populations, stromal cells, vascular-associated cells, and brain-associated populations. Differences in cellular composition can influence bulk transcriptional profiles and complicate the interpretation of gene-expression changes. In addition, metastatic specimens may vary considerably in the relative abundance of individual populations and in their transcriptional characteristics.   

Single-nucleus RNA sequencing offers an approach to characterizing this cellular heterogeneity while enabling the examination of gene expression across distinct populations within metastatic tissue. The integration of single-nucleus transcriptomics with bulk expression analyses may help contextualize tissue-level findings, identify cellular populations associated with selected transcriptional programmes, and characterize inter-sample variation.   

Such integration is particularly relevant for ECM-associated transcription, where differences in cellular composition may contribute to observed variation in matrix-associated expression across specimens. 

### 1.5 Rationale and objectives of the present study
Despite the importance of tumour adaptation and microenvironmental interactions in brain metastasis, the relationship between tissue-level transcriptional changes, metabolic-associated programmes, ECM-associated expression, and cellular heterogeneity remains an important area of investigation.   

In this study, we conducted an integrative transcriptomic investigation of breast cancer brain metastasis by combining bulk gene-expression analyses, independent validation, pathway enrichment, and single-nucleus transcriptomic characterization.   

Our objectives were to:

1. Characterize transcriptional differences between primary breast tumours and brain metastases using paired differential-expression analysis.

2. Evaluate the robustness of identified transcriptional changes using purity-adjusted modelling and an independent validation cohort.

3. Identify biological pathways associated with brain metastasis, with particular attention to metabolic, ECM-associated, and immune-related programmes.

4. Characterize the cellular composition of brain metastatic specimens using single-nucleus transcriptomic data.

5. Investigate the cellular distribution of ECM-associated transcripts across major cell populations.

6. Assess inter-sample heterogeneity in ECM-associated transcription and its relationship to cellular composition.

7. Integrate bulk and single-nucleus findings into a coherent biological framework for understanding transcriptional adaptation and microenvironmental variation in breast cancer brain metastasis


## 2. Materials and Methods

### 2.1 Study design and analytical workflow
This study employed a retrospective, computational transcriptomic design combining bulk gene-expression and single-nucleus RNA sequencing (snRNA-seq) analyses to investigate molecular alterations associated with breast cancer brain metastasis (BCBM). The framework comprised two complementary branches:

- Bulk transcriptomic branch: Compared primary breast tumours and matched brain metastases, tested sensitivity to tumour purity, evaluated reproducibility in an independent paired validation cohort, and performed pathway enrichment and leading-edge analyses.   

- Single-nucleus transcriptomic branch: Characterized cellular composition across brain metastatic lesions, examined ECM-associated gene expression and detection prevalence across annotated cell populations, aggregated sample-by-cell-type pseudobulk profiles, scored ECM-associated transcriptional variation across samples, and evaluated associations between ECM scores and cellular composition. 

### 2.2 Data sources and cohort characteristics
Three publicly available Gene Expression Omnibus (GEO) datasets were analyzed:
- Discovery bulk cohort (`GSE184869`): 45 paired primary breast tumour and brain metastasis samples (90 samples total) spanning 16,714 protein-coding genes indexed by Ensembl gene identifiers.   
- Independent validation bulk cohort (`GSE125989`): An independent paired primary tumour and brain metastasis dataset (16 paired samples) used to evaluate reproducibility outside the discovery cohort.   
- Single-nucleus transcriptomic cohort (`GSE324453`): An annotated Seurat object comprising 56,268 nuclei, 35,259 features, and 28 brain metastatic samples across 12 annotated major cell populations.

### 2.3 Bulk transcriptomic data preprocessing
For GSE184869, the supplied expression matrix was obtained from GEO supplementary files as batch-corrected, $\log_2$-transformed, Trimmed Mean of M-values (TMM)-normalized counts per million (CPM), and these upstream preprocessing steps were retained as supplied. Sample and patient metadata were standardized, matched to the expression columns, and reordered accordingly. Integrity checks verified 90 matched samples representing 45 complete primary–metastasis pairs with no missing expression values. 

### 2.4 Differential gene expression analysis
Paired differential expression in GSE184869 was performed in R using limma (lmFit followed by empirical Bayes moderation via eBayes). In the discovery model (Model A), linear models were specified as ~ patientID + tumourID, with tumour state encoded as primary tumour (P) versus brain metastasis (M) and the contrast evaluated on the tumourIDM coefficient (positive $\log_2$ fold change indicating higher expression in brain metastases).   

To evaluate robustness to variation in tumour purity, a sensitivity analysis was conducted on GSE184869 samples with complete tumour-purity metadata and paired primary–metastasis representation. Purity values were parsed to numeric format (removing percentage symbols where present), and Model A (~ patientID + tumourID) was compared against purity-adjusted Model B (~ patientID + tumourID + purity) using limma and empirical Bayes moderation. Because both models used the same discovery dataset, Model B was treated as a sensitivity analysis rather than an independent validation. Multiple testing was controlled via Benjamini–Hochberg (BH) false discovery rate (FDR < 0.05), with $\vert{}\log_2\text{FC}\vert{} > 1$ used as a descriptive magnitude threshold. 

### 2.5 Independent validation cohort analysis
GSE125989 was analyzed separately to assess reproducibility of paired primary-versus-brain-metastasis expression shifts using FDR-adjusted significance and fold-change criteria.   

### 2.6 Consensus differential expression analysis
Differential expression outputs from the discovery model, purity-adjusted sensitivity model, and GSE125989 validation cohort were integrated across cohorts/models. Direction-consistent consensus genes were defined by concordant direction of expression change across all three analyses, and stringent consensus genes were identified as those reaching statistical significance across all three comparisons.

### 2.7 Functional pathway enrichment analysis
Genes from each bulk comparison were ranked by their limma moderated $t$-statistic (positive values toward brain metastasis, negative values toward primary tumour). For discovery GSEA, Ensembl IDs had version suffixes stripped and were mapped to Entrez Gene IDs using clusterProfiler and org.Hs.eg.db, retaining the entry with the largest absolute moderated $t$-statistic where multiple Ensembl IDs mapped to a single Entrez ID.   

Gene Ontology Biological Process (GO BP) and Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analyses were conducted using gseGO() and gseKEGG() across gene-set sizes of 10–500 with BH multiple-testing correction. Pathways were matched by database identifiers across the discovery, purity-adjusted, and validation GSEA runs; pathways achieving adjusted $p < 0.05$ and consistent normalized enrichment score (NES) direction across all three analyses were classified as robust consensus pathways (166 GO BP and 19 KEGG pathways).   

### 2.8 Leading-edge gene analysis
Leading-edge genes were extracted from the core_enrichment field for seven prespecified pathways spanning oxidative phosphorylation, aerobic respiration, ECM organization, collagen fibril organization, ECM–receptor interaction, integrin signalling, and cytokine–cytokine receptor interaction. Identifiers were harmonized via AnnotationDbi and org.Hs.eg.db across recognized symbols, aliases, and Entrez IDs (retaining unmapped raw IDs diagnostically while excluding them from mapped consensus comparisons). Leading-edge genes supported in at least two analyses were classified as consensus leading-edge genes, and those supported across all three were noted as shared across all analyses.

### 2.9 Single-nucleus RNA-seq data processing
The annotated Seurat object from GSE324453 (56,268 nuclei, 35,259 features, 28 brain metastatic samples) was analyzed along with supplied sample metadata and cellular annotations. Quality control distributions—including library-level and patient/sample-level metrics, UMI counts, detected genes, and mitochondrial transcript fractions — were inspected across samples and clusters to evaluate sequencing characteristics.   

### 2.10 Cell-type annotation and cellular composition analysis
Analyses used the 12 annotated major populations in GSE324453: tumour cells, proliferating tumour cells, myeloid cells, T/NK cells, stromal cells, pericytes, endothelial cells, neurons, oligodendrocyte-lineage cells, astrocytes, plasma cells, and B cells. For each of the 28 samples, cell-type proportions were calculated relative to the analyzed nuclei in that sample, and cohort-level mean proportions and standard deviations were computed across samples with equal sample weighting (rather than pooling all nuclei across samples). 

### 2.11 ECM-associated gene expression across cell populations
A 13-gene ECM-associated panel (COL6A3, COL14A1, COL15A1, AEBP1, SULF1, THBS2, THBS4, LUM, POSTN, MMP2, CTSK, CXCL12, and DPT) was examined across the 12 annotated populations. Detection prevalence (percentage of nuclei within a given cell population with detectable expression) and relative expression patterns were evaluated together as complementary metrics. 

Sample-by-cell-type pseudobulk aggregation was also performed by grouping nuclei by sample of origin and annotated cell population to reduce nuclear pseudoreplication and summarize expression at the biological-sample level.   

### 2.12 ECM signature scoring and sample stratification
A seven-gene core ECM signature (COL6A3, AEBP1, SULF1, THBS2, LUM, POSTN, and CTSK) was constructed to summarize sample-level ECM-associated transcriptional variation across the 28 brain metastatic samples, alongside complementary evaluation of the expanded 13-gene panel. Samples were ordered along the continuous ECM score distribution and stratified at the median score into high- and low-ECM analytical groups for descriptive and comparative analyses (without treating these strata as independently established biological subtypes).

### 2.13 Statistical analysis and visualization
All statistical analyses and figure generation were performed in R. Multiple testing in differential expression and pathway enrichment was controlled using BH FDR (< 0.05). Correlation and high- versus low-ECM group comparisons across the 28 metastatic samples were evaluated with multiple-testing correction where appropriate and interpreted cautiously given sample size and observational design. Main figures (Figures 1–6) and supplementary visualizations were generated directly from the corresponding analysis outputs. 

## 3. Results

### 3.1 Transcriptomic divergence between primary breast tumours and brain metastases
In the paired discovery cohort (GSE184869; 45 paired primary breast tumour and brain metastasis samples), differential expression analysis identified 5,247 significantly altered genes at $\text{FDR} < 0.05$, consisting of 2,816 upregulated and 2,431 downregulated genes in brain metastases. Applying the stricter filter of $\text{FDR} < 0.05$ and $\vert{}\log_2\text{FC}\vert{} > 1$ yielded 665 genes (146 upregulated and 519 downregulated).   Notable shifts in the discovery cohort included strong downregulation of HTRA3, multiple skeletal muscle-associated myosin genes (MYH1, MYH4, MYH13, MYH2, MYH8), CTSK, PODN, and LRRC15, alongside pronounced upregulation of GFAP in brain metastases. In the purity-adjusted sensitivity analysis (Model B), 242 genes remained significant at $\text{FDR} < 0.05$ (101 upregulated, 141 downregulated; 185 genes at $\vert{}\log_2\text{FC}\vert{} > 1$), with key alterations such as GFAP upregulation and HTRA3, CTSK, PODN, LRRC15, and MYH downregulation persisting.   

Independent validation in the paired GSE125989 cohort identified 87 unique genes (118 probes) meeting the stringent differential expression criteria (81 downregulated and 6 upregulated), with prominent reductions in ECM-related genes including AEBP1, COL1A1, COL1A2, GREM1, POSTN, SPON1, SPARC, and DPT. Across the discovery, purity-adjusted, and validation analyses, 165 genes showed consistent expression direction (99 downregulated and 66 upregulated), and 15 genes reached statistical significance across all three comparisons — all downregulated in brain metastases: AEBP1, DPT, POSTN, CTSK, COL15A1, THBS2, COL6A3, LUM, CXCL12, SULF1, COL14A1, MMP2, IGF1, PDGFRA, and THBS4.   (Refer to Figure 1 [Figure_1_Integrated_Transcriptomic_Landscape.png], illustrating the study workflow, discovery and purity-adjusted differential expression landscapes, and three-way consensus integration.) 

### 3.2 Functional pathway remodelling in BCBM
GSEA across the ranked differential expression profiles revealed coordinated pathway-level shifts between primary breast tumours and brain metastases, integrating into 166 robust GO BP pathways and 19 robust KEGG pathways significant with consistent NES direction across all three analyses.

#### 3.2.1 Increased mitochondrial energy-metabolism-associated transcription
Brain metastases showed consistent positive enrichment for mitochondrial energy-metabolism-associated transcriptional programmes across the discovery, purity-adjusted, and validation analyses:   
- Oxidative phosphorylation: $\text{NES} \approx +3.14$ (discovery), $+2.89$ (purity-adjusted), and $+2.34$ (validation).   
- Aerobic respiration: $\text{NES} \approx +2.98$ (discovery), $+2.67$ (purity-adjusted), and $+2.05$ (validation).   
- Respiratory electron transport / electron transport chain: Consistent positive enrichment across all three comparisons ($\text{NES} \approx +3.15, +2.62,$ and $+2.16$ for respiratory electron transport).   

Leading-edge analysis identified core pathway-contributing genes driving these metabolic enrichment peaks. In accordance with our analytical scope, these findings reflect increased mitochondrial energy-metabolism-associated transcription rather than direct proof of increased mitochondrial respiration or metabolic flux.   (Refer to Figure 2 [././results/figures/Figure_2/Figure_2_Functional_Pathway_Landscape.png], metabolic enrichment panel.)   

3.2.2 Reduced ECM and collagen-associated transcriptionConversely, ECM-related pathways demonstrated consistent negative enrichment in brain metastases across all three analyses:   ECM organization: $\text{NES} \approx -3.09$ (discovery), $-2.78$ (purity-adjusted), and $-2.84$ (validation).   Collagen fibril organization: $\text{NES} \approx -2.75$, $-2.39$, and $-2.84$.   Collagen metabolism: $\text{NES} \approx -2.81$, $-2.54$, and $-2.36$.   ECM assembly: $\text{NES} \approx -2.39$, $-1.95$, and $-2.51$.   Together with negative enrichment of ECM–receptor interaction and the downregulation of multiple ECM genes in the 15-gene consensus set, these results indicate reduced ECM-associated transcription at the bulk comparison level without directly measuring ECM protein deposition or matrix restructuring.   (Refer to Figure 2 [Figure_2_Functional_Pathway_Landscape.png], ECM enrichment panel.)  







