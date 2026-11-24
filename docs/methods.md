## 5. Materials and Methods

### 5.1. Study design and analytical framework

This study employed an integrative transcriptomic approach combining bulk transcriptomic and single-nucleus RNA sequencing analyses to investigate molecular alterations associated with breast cancer brain metastasis. The analytical framework consisted of two complementary branches. First, bulk transcriptomic data were used to identify differential gene expression between primary breast tumours and brain metastases, assess the robustness of these alterations through purity-adjusted analyses and independent validation, and characterize associated biological pathways. Second, single-nucleus transcriptomic data from brain metastatic specimens were analyzed to characterize cellular composition, investigate the cellular distribution of extracellular matrix (ECM)-associated gene expression, and evaluate sample-level heterogeneity in ECM-associated transcription.

Bulk transcriptomic analyses incorporated differential expression, validation, consensus integration, and functional enrichment approaches. Single-nucleus analyses included quality control, cellular composition profiling, cell-type-associated ECM expression analysis, pseudobulk characterization, ECM signature scoring, and assessment of associations between ECM-associated transcriptional scores and cellular abundance. Findings from the two analytical branches were integrated at the biological interpretation level to characterize metabolic transcriptional alterations and microenvironmental heterogeneity associated with breast cancer brain metastasis. As the datasets were not patient-matched, integration was performed across complementary findings rather than through direct patient-level comparisons.


### 5.2. Transcriptomic datasets

The bulk discovery analysis used the publicly available GSE184869 dataset, comprising paired primary breast tumour and brain metastasis samples from 45 patients. The expression matrix contained 16,714 protein-coding genes across 90 samples. Gene identifiers were provided as Ensembl gene IDs.

The expression data were obtained from the GEO supplementary files and associated study materials. The supplied expression matrix was already batch-corrected and processed using log2-transformed, TMM-normalized counts per million (CPM). Corresponding sample and patient metadata were obtained from the supplementary study metadata file.

Sample identifiers were standardized and matched between the expression matrix and metadata. Metadata records were reordered to correspond to the expression matrix columns. Data integrity checks confirmed 90 matched samples, 45 complete primary tumour–brain metastasis pairs, and no missing expression values.

The independent bulk validation dataset was GSE125989. The single-nucleus transcriptomic dataset was GSE324453. Their specific platforms, preprocessing procedures, and sample characteristics are described in the corresponding analytical subsections.


### 5.3 Bulk differential expression analysis

Differential expression between primary breast tumours and brain metastases in the GSE184869 discovery dataset was assessed using the limma package in R. The supplied log2-transformed, batch-corrected, TMM-normalized CPM expression matrix was used as the input.

A paired linear model was specified using patient identity and tumour type as explanatory variables. The design matrix was constructed using the formula `~ patientID + tumourID`, with tumour type encoded as a factor containing primary tumour (P) and brain metastasis (M) levels. The model was fitted using `lmFit`, followed by empirical Bayes moderation using `eBayes`.

The brain metastasis versus primary tumour contrast was evaluated using the `tumourIDM` coefficient. Consequently, positive log2 fold-change values represented higher expression in brain metastases, whereas negative values represented lower expression relative to primary tumours.

Gene-level significance was assessed using Benjamini–Hochberg-adjusted p-values, with an FDR threshold of 0.05. A more stringent effect-size criterion of absolute log2 fold change greater than 1 was applied in subsequent result summaries where specified.


### 5.4. Purity-adjusted differential expression analysis

To assess whether the observed transcriptional differences were sensitive to variation in tumour purity, a second paired differential expression model was fitted using the GSE184869 discovery cohort. Tumour-purity information was obtained from the supplementary study metadata. The purity variable was converted to numeric format, with percentage symbols removed where present.

Samples with missing tumour-purity information were excluded. Following this filtering, only patients with both primary tumour and brain metastasis samples were retained. Two models were then fitted to the resulting subset using the limma package.

Model A included patient identity and tumour type as explanatory variables (`~ patientID + tumourID`), whereas Model B additionally included tumour purity as a continuous covariate (`~ patientID + tumourID + purity`). Both models were fitted using `lmFit` and empirical Bayes moderation with `eBayes`. The brain metastasis versus primary tumour effect was extracted using the `tumourIDM` coefficient.

Benjamini–Hochberg-adjusted p-values were used to assess gene-level significance, with an FDR threshold of 0.05. Model A and Model B effect estimates were compared to examine changes in log2 fold change, retention of effect direction, and statistical significance following purity adjustment. Because both models were fitted using the same discovery cohort, the purity-adjusted analysis was treated as a sensitivity analysis rather than independent validation.


### 5.5. Independent validation and consensus integration

An independent bulk transcriptomic dataset, GSE125989, was analyzed to evaluate the reproducibility of differential expression findings identified in the discovery cohort. Differential expression between paired primary breast tumours and brain metastases was assessed within the validation dataset, and genes exhibiting significant expression differences were identified using false discovery rate-adjusted significance criteria and effect-size thresholds.

Consensus integration was performed by comparing differential expression findings across the discovery analysis, purity-adjusted analysis, and independent validation dataset. Genes exhibiting consistent directions of differential expression across the analytical datasets were prioritized as directionally concordant consensus genes. This integration approach was used to identify reproducible transcriptional patterns associated with brain metastases while distinguishing concordance across analytical models from independent cohort replication. As the discovery and purity-adjusted analyses used the same cohort, they were not considered independent sources of biological replication.


### 5.6 Functional pathway enrichment and leading-edge analysis

Gene set enrichment analysis (GSEA) was performed using ranked differential expression results to identify biological pathways associated with breast cancer brain metastasis. Genes were ranked by their moderated *t*-statistics, with positive values indicating higher expression in brain metastases and negative values indicating lower expression relative to primary breast tumours.

For the discovery analysis, Ensembl gene identifiers were stripped of version suffixes and mapped to Entrez Gene identifiers using `clusterProfiler` and `org.Hs.eg.db`. Where multiple Ensembl identifiers mapped to the same Entrez identifier, the entry with the largest absolute moderated *t*-statistic was retained. The resulting ranked gene list was used for Gene Ontology Biological Process (GO BP) and Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analyses using `gseGO()` and `gseKEGG()`, respectively. Gene sets containing 10–500 genes were considered, and Benjamini–Hochberg adjustment was applied to account for multiple testing.

Pathway-level reproducibility was assessed by integrating enrichment results from the discovery differential expression analysis, the purity-adjusted sensitivity analysis, and the independent validation cohort. Pathways were matched by their database identifiers. A pathway was classified as a robust consensus pathway when it showed an adjusted *p* value below 0.05 in all three analyses and maintained a consistent normalized enrichment score (NES) direction across them. This approach identified 166 robust GO biological processes and 19 robust KEGG pathways.

To investigate the genes contributing to selected biological signals, leading-edge genes were extracted from the `core_enrichment` field of the GSEA results for seven prespecified pathways covering oxidative phosphorylation, aerobic respiration, extracellular matrix organization, collagen fibril organization, ECM–receptor interaction, integrin signalling, and cytokine–cytokine receptor interaction. Identifier harmonization was performed using `AnnotationDbi` and `org.Hs.eg.db`, including the resolution of recognized gene symbols, aliases, and Entrez identifiers. Unmapped identifiers were retained as explicitly labelled raw identifiers in the diagnostic output but were excluded from mapped-gene consensus comparisons.

For each selected pathway, gene-level support was assessed across the discovery, purity-adjusted, and validation analyses. Genes detected in at least two analyses were classified as consensus leading-edge genes, while genes supported by all three analyses were classified as shared across all analyses. These analyses were used to assess the reproducibility of pathway-associated gene signals rather than to establish independent functional validation of individual genes.


### 5.7. Single-nucleus transcriptomic dataset and preprocessing

Single-nucleus transcriptomic data were obtained from GEO dataset GSE324453. An annotated Seurat object was used for downstream analysis. The object contained 56,268 nuclei, 35,259 features, and 28 samples. Cellular annotations represented 12 major populations: tumour cells, proliferating tumour cells, myeloid cells, T/NK cells, stromal cells, pericytes, endothelial cells, neurons, oligodendrocyte-lineage cells, astrocytes, plasma cells, and B cells.

The annotated cellular identities and sample-level metadata were used to characterize cellular composition and investigate cell-type-associated ECM transcription. Quality control assessment included examination of library-level and sample-level metrics, UMI counts, detected genes, and mitochondrial transcript proportions. Quality control distributions were visualized across samples and annotated cellular populations to assess variation in sequencing characteristics within the dataset.


### 5.9 ECM cellular-origin and expression analysis

To investigate the cellular distribution of extracellular matrix (ECM)-associated transcription, a predefined panel of ECM-related genes was examined across the major annotated cell populations in the GSE324453 single-nucleus dataset. The analysis focused on genes associated with ECM structure, organization, matrix-associated signaling, and remodeling, including COL6A3, COL14A1, COL15A1, AEBP1, SULF1, THBS2, THBS4, LUM, POSTN, MMP2, CTSK, CXCL12, and DPT.

Gene expression was evaluated separately within each annotated cell population to identify the cellular compartments in which individual ECM-associated transcripts were detected. Detection prevalence was calculated as the percentage of nuclei within a given cell population expressing the corresponding gene. This measure was used to characterize how widely each transcript was detected within a cell type.

Relative expression patterns were additionally examined across cell populations to identify cell-type-associated transcriptional signals. Expression patterns were visualized using heatmaps and related comparative plots. These analyses were used to distinguish broadly distributed ECM-associated transcription from signals concentrated in particular cellular compartments.

Detection prevalence and relative expression were interpreted as complementary measures. A transcript detected in a high proportion of nuclei may have a different expression profile from one detected in a smaller subset of nuclei. Accordingly, detection prevalence was not treated as a direct measure of transcript abundance.

The cellular distribution of ECM-associated transcripts was interpreted as evidence of potential cellular contributions to ECM-related transcriptional programs. However, transcript detection alone does not establish extracellular protein abundance, matrix deposition, or functional tissue remodeling.


### 5.9 ECM cellular-origin and expression analysis

To investigate the cellular distribution of extracellular matrix (ECM)-associated transcription, a predefined panel of ECM-related genes was examined across the major annotated cell populations in the GSE324453 single-nucleus dataset. The analysis focused on genes associated with ECM structure, organization, matrix-associated signaling, and remodeling, including COL6A3, COL14A1, COL15A1, AEBP1, SULF1, THBS2, THBS4, LUM, POSTN, MMP2, CTSK, CXCL12, and DPT.

Gene expression was evaluated separately within each annotated cell population to identify the cellular compartments in which individual ECM-associated transcripts were detected. Detection prevalence was calculated as the percentage of nuclei within a given cell population expressing the corresponding gene. This measure was used to characterize how widely each transcript was detected within a cell type.

Relative expression patterns were additionally examined across cell populations to identify cell-type-associated transcriptional signals. Expression patterns were visualized using heatmaps and related comparative plots. These analyses were used to distinguish broadly distributed ECM-associated transcription from signals concentrated in particular cellular compartments.

Detection prevalence and relative expression were interpreted as complementary measures. A transcript detected in a high proportion of nuclei may have a different expression profile from one detected in a smaller subset of nuclei. Accordingly, detection prevalence was not treated as a direct measure of transcript abundance.

The cellular distribution of ECM-associated transcripts was interpreted as evidence of potential cellular contributions to ECM-related transcriptional programs. However, transcript detection alone does not establish extracellular protein abundance, matrix deposition, or functional tissue remodeling.


### 5.10 Pseudobulk analysis

To account for the nested structure of the single-nucleus dataset, pseudobulk analysis was used to summarize gene expression at the biological-sample level rather than treating individual nuclei as independent replicates. Nuclei were grouped by their sample of origin and annotated cell population, and expression measurements were aggregated within each group to generate sample-level expression profiles.

This aggregation provides a framework for examining transcriptional patterns while reducing the risk of pseudoreplication arising from multiple nuclei originating from the same biological sample. Sample-level profiles were used as the units of observation for downstream comparisons and exploratory analyses, where applicable.

Pseudobulk results were interpreted in the context of cell population identity and between-sample variability. Because the number of nuclei recovered for each cell population may differ across samples, the availability and representation of each sample–cell population combination were considered when interpreting the results.

Pseudobulk expression was treated as a sample-level summary of transcription within annotated cell populations, rather than as a direct measurement of tissue-level gene expression or extracellular protein abundance.


### 5.11 ECM signature construction and scoring

An ECM-associated transcriptional signature was constructed to summarize sample-level variation in the expression of genes associated with extracellular matrix structure, organization, and remodeling. The signature was based on seven core ECM-associated genes: **COL6A3, AEBP1, SULF1, THBS2, LUM, POSTN, and CTSK**.

An expanded 13-gene panel was also examined to capture a broader range of ECM-associated and vascular-associated transcriptional signals. In addition to the seven core genes, this panel included **COL14A1, COL15A1, MMP2, CXCL12, THBS4, and DPT**.

Expression of the signature genes was evaluated across individual samples. A composite ECM signature score was then used to summarize the overall expression pattern of the core gene panel in each sample. The resulting scores enabled comparisons of ECM-associated transcriptional variation across the 28 brain-metastasis samples.

The core seven-gene panel was used for the principal sample-level ECM scoring analysis, whereas the expanded panel was used for complementary examination of additional ECM-associated transcripts. This distinction was maintained to avoid treating all genes in the broader panel as equally representative of the core stromal ECM signature.

The ECM score was interpreted as a relative measure of ECM-associated transcriptional activity within the analyzed dataset. It was not considered a direct measurement of extracellular matrix abundance, protein deposition, or tissue-level matrix remodeling.


### 5.12 Sample-level heterogeneity and statistical comparisons

Sample-level heterogeneity in ECM-associated transcription was investigated using the core seven-gene ECM signature scores calculated for the 28 brain-metastasis samples. The distribution of scores was examined to characterize the extent and pattern of inter-sample variation. Samples were ordered according to their signature scores to assess whether ECM-associated transcription followed a continuous distribution or suggested distinct groups.

For exploratory comparisons, samples were divided into high- and low-ECM groups based on the median signature score. This grouping was used to compare cellular composition and gene-level expression patterns between samples with relatively higher and lower ECM-associated transcription. The high- and low-score groups were treated as analytical subsets of a continuous distribution, not as independently established biological subtypes.

The relationship between the ECM signature score and stromal-cell proportion was examined to assess whether variation in ECM-associated transcription was associated with the relative representation of stromal nuclei across samples. Associations between the signature score and cellular composition were evaluated at the sample level.

Differences in cell-type proportions and individual ECM-associated gene expression between the high- and low-score groups were examined using statistical comparisons. Where multiple gene-level or cell-population comparisons were performed, adjusted significance values were considered to account for multiple testing.

All comparisons were interpreted as exploratory analyses of inter-sample variation. Associations between ECM scores and cellular proportions were not interpreted as evidence that changes in one directly caused changes in the other.


### 5.13 Statistical analysis and reproducibility

Statistical analyses were conducted using R. Differential expression and functional enrichment analyses were evaluated using false discovery rate (FDR)-adjusted significance values to account for multiple hypothesis testing. Unless otherwise specified, an FDR-adjusted \(p\)-value below 0.05 was considered statistically significant.

For differential expression analyses, genes were additionally examined using an absolute log2 fold-change threshold of greater than 1 when applying the more stringent effect-size criterion. Pathway enrichment results were interpreted using normalized enrichment scores (NES), with positive values indicating enrichment toward the brain-metastasis condition and negative values indicating enrichment toward the primary-tumour condition.

For single-nucleus analyses, cellular composition, ECM-associated expression, and sample-level signature scores were summarized using descriptive statistics and exploratory statistical comparisons. Multiple-testing correction was applied where appropriate. Correlation and group-comparison results were interpreted according to their statistical significance and biological context, with particular consideration of the limited number of biological samples.

Analytical workflows were organized into separate project directories for bulk transcriptomics, single-nucleus transcriptomics, and integrated results. Intermediate data, statistical outputs, figures, and analysis scripts were retained to support traceability of the computational workflow. Figures and summary tables were generated from the corresponding analysis outputs.

The study relied on publicly available transcriptomic datasets and computational analyses. No new biological specimens were collected for this investigation. All findings were interpreted within the limitations of the underlying datasets and the observational nature of the analyses.

