# Molecular Adaptation and Tumour–Microenvironment Interactions in Breast Cancer Brain Metastasis: An Integrative Transcriptomic Investigation

## Abstract

### Background

Breast cancer brain metastasis (BCBM) represents a clinically challenging manifestation of advanced breast cancer. Successful colonization of the brain requires tumour cells to adapt to a distinct tissue environment characterized by specialized metabolic demands, neural-associated cell populations, vascular structures, and a unique extracellular matrix (ECM). However, the transcriptional adaptations accompanying brain metastasis and the cellular contributions to the metastatic microenvironment remain incompletely characterized.

### Objectives

This study aimed to investigate transcriptional differences between primary breast tumours and brain metastases, identify reproducible metabolic and microenvironment-associated pathways, characterize the cellular composition of brain metastatic lesions, and examine the cellular origins and inter-sample heterogeneity of ECM-associated transcription.

### Methods

An integrative transcriptomic approach combined bulk RNA sequencing datasets with single-nucleus RNA sequencing data. Differential expression analysis was performed using paired primary breast tumour and brain metastasis samples from GSE184869. A tumour-purity-adjusted sensitivity analysis was conducted to assess the robustness of differential expression findings. An independent paired validation cohort, GSE125989, was used to evaluate reproducibility. Gene set enrichment analysis (GSEA), cross-analysis pathway integration, and leading-edge analysis were performed to identify reproducible biological programmes and contributing genes.

Single-nucleus transcriptomic data from GSE324453 were used to characterize cellular composition and investigate ECM-associated gene expression across annotated cell populations. A seven-gene core ECM signature was used to assess sample-level variation, its relationship with stromal-cell representation, and differences in cellular composition between analytically defined high- and low-scoring samples.

### Results

Bulk transcriptomic analyses identified increased expression of oxidative phosphorylation, aerobic respiration, and electron transport chain-associated pathways in brain metastases. Conversely, extracellular matrix organization, collagen-associated processes, ECM–receptor interactions, and selected immune-related pathways showed reduced expression relative to primary breast tumours. These patterns were examined using a paired discovery model, a purity-adjusted sensitivity model, and an independent validation cohort.

Cross-analysis integration identified 166 robust Gene Ontology Biological Process pathways and 19 robust KEGG pathways with statistical significance and consistent enrichment direction across all three analyses. Fifteen genes were significantly altered in the same direction across all three differential-expression analyses.

Single-nucleus analysis of 56,268 nuclei from 28 metastatic samples demonstrated a tumour-dominated cellular landscape with substantial representation of myeloid, lymphoid, stromal, vascular-associated, and neural-associated populations. ECM-associated transcription was distributed across several non-malignant cell populations, particularly stromal cells, with additional gene-specific contributions from pericytes, endothelial cells, and astrocytes.

The core ECM signature showed substantial variation across metastatic samples and was positively associated with stromal-cell proportion (Spearman rho = 0.623, false discovery rate = 0.00236). Comparisons between high- and low-scoring samples revealed differences in stromal representation and selected ECM-associated transcripts.

### Conclusion

Breast cancer brain metastases exhibit a transcriptional landscape characterized by increased mitochondrial energy-metabolism-associated expression and reduced expression of several ECM-associated and selected immune-related programmes relative to primary breast tumours. Within brain metastases, ECM-associated transcription remains detectable across multiple cellular compartments and varies between samples, partly in association with stromal-cell representation. These findings support a model of metabolic transcriptional adaptation accompanied by a heterogeneous, multicellular tumour–microenvironment landscape.

**Keywords:** Breast cancer brain metastasis; transcriptomics; single-nucleus RNA sequencing; oxidative phosphorylation; extracellular matrix; tumour microenvironment; stromal cells; metabolic adaptation.

---

# 1. Introduction

## 1.1 Breast cancer brain metastasis

Breast cancer is a heterogeneous malignancy capable of spreading to distant organs, including the brain. Brain metastasis is associated with substantial morbidity and presents major therapeutic challenges because the brain differs markedly from the primary tumour environment in its cellular composition, vascular organization, metabolic requirements, and extracellular matrix characteristics.

The development of brain metastases involves multiple biological processes, including tumour-cell dissemination, vascular interactions, adaptation to the brain environment, survival within a distinct tissue compartment, and interactions with resident and recruited non-malignant cells.

Although molecular characteristics of primary breast cancer have been extensively studied, the transcriptional changes accompanying brain metastatic progression remain an important area of investigation.

![./report_illustrations/introduction_1.1/1.jpg](https://www.biorender.com/template/breast-cancer-to-brain-metastasis)
![./report_illustrations/introduction_1.1/2.jpg](https://www.mdpi.com/2073-4409/11/8/1263)
![./report_illustrations/introduction_1.1/3.jpg](https://www.nature.com/articles/nrc.2016.25)
![./report_illustrations/introduction_1.1/4.jpg](https://www.mdpi.com/2073-4409/11/8/1263)
![./report_illustrations/introduction_1.1/5.jpg](https://www.mdpi.com/1422-0067/24/15/12034)
![./report_illustrations/introduction_1.1/6.jpg](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2024.1448089/full)
![./report_illustrations/introduction_1.1/7.jpg](https://link.springer.com/article/10.1007/s12609-023-00488-0)
![./report_illustrations/introduction_1.1/8.jpg](https://www.nature.com/articles/nrclinonc.2011.58)
![./report_illustrations/introduction_1.1/9.jpg](https://academic.oup.com/neuro-oncology/article/28/2/319/8275781)
![./report_illustrations/introduction_1.1/10.jpg](https://www.biorender.com/template/structure-of-the-btb)



## 1.2 Metabolic adaptation in the brain metastatic environment

Metastatic tumour cells encounter environmental conditions that differ from those of their primary tissue. Adaptation to these conditions may involve changes in nutrient utilization, mitochondrial metabolism, oxidative phosphorylation, and other energy-associated pathways.

The brain has distinctive metabolic requirements and a highly specialized cellular environment. Therefore, changes in mitochondrial energy-metabolism-associated transcription may represent an important component of the molecular landscape of breast cancer brain metastases.

However, increased expression of metabolic genes does not necessarily establish increased metabolic activity. Transcriptomic analyses can identify pathway-associated expression patterns, but functional measurements are required to determine actual respiratory activity, ATP production, and metabolic flux.

![./report_illustrations/introduction_1.2/1.jpg](https://bio.libretexts.org/Bookshelves/Introductory_and_General_Biology/General_Biology_2e_%28OpenStax%29/02%3A_Unit_II-_The_Cell/2.04%3A_Cellular_Respiration/2.4.05%3A_Oxidative_Phosphorylation)
![./report_illustrations/introduction_1.2/2.jpg]()
![./report_illustrations/introduction_1.2/3.jpg](https://www.cambridge.org/core/books/basic-physiology-for-anaesthetists/energy-balance/3B4A3DAB1C21402A39939DF91C7838BE)
![./report_illustrations/introduction_1.2/4.jpg](https://livephysics.com/infographics/biology-the-electron-transport-chain/)


## 1.3 Extracellular matrix and tumour–microenvironment interactions

The extracellular matrix is a complex network of structural and regulatory molecules that contributes to tissue organization, cellular adhesion, signalling, migration, and tissue homeostasis.

ECM-associated genes, including collagens, matricellular proteins, proteases, and matrix-associated regulatory molecules, are expressed by different cellular populations depending on tissue context.

In metastatic cancer, changes in ECM-associated transcription may reflect alterations in cellular composition, gene expression within individual populations, or both.

Importantly, bulk transcriptomic measurements cannot independently determine which cell populations contribute to an observed expression pattern. Single-cell and single-nucleus approaches can provide additional cellular resolution and help identify potential cellular sources of selected transcriptional programmes.


![./report_illustrations/introduction_1.3/1.jpg](https://uen.pressbooks.pub/anatomyphysiology/chapter/features-of-connective-tissues/)
![./report_illustrations/introduction_1.3/2.jpg](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1274474/full)
![./report_illustrations/introduction_1.3/3.jpg](https://cordis.europa.eu/project/id/101138496/reporting/fr)
![./report_illustrations/introduction_1.3/4.jpg](https://link.springer.com/article/10.1186/s13046-024-03236-z)
![./report_illustrations/introduction_1.3/5.jpg](https://link.springer.com/article/10.1186/s13046-024-03236-z)
![./report_illustrations/introduction_1.3/6.jpg](https://www.eurekalert.org/multimedia/1054662)
![./report_illustrations/introduction_1.3/7.jpg](https://link.springer.com/chapter/10.1007/978-3-030-98950-7_5)
![./report_illustrations/introduction_1.3/8.jpg](https://pmc.ncbi.nlm.nih.gov/articles/PMC11452922/)
![./report_illustrations/introduction_1.3/9.jpg](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2024.1363004/full)
![./report_illustrations/introduction_1.3/10.jpg](https://www.nature.com/articles/s41389-017-0011-9)


## 1.4 The need for integrative transcriptomic analysis

Bulk RNA sequencing provides a broad view of transcriptional differences across tissue samples but combines signals from malignant and non-malignant cellular populations.

Single-nucleus RNA sequencing provides an opportunity to investigate cellular heterogeneity and identify cell populations associated with specific transcriptional programmes.

Combining these approaches provides complementary information:

- Bulk transcriptomics identifies tissue-level transcriptional differences between primary tumours and brain metastases.
- Single-nucleus transcriptomics characterizes the cellular landscape within brain metastases.
- Cross-analysis integration evaluates whether biological pathways and genes show reproducible patterns across analytical models.
- Sample-level signature analysis investigates variation between individual metastatic samples.

This study applies an integrative strategy to investigate metabolic adaptation, ECM-associated transcription, and cellular heterogeneity in breast cancer brain metastasis.

---

# 2. Research Objectives

## 2.1 General objective

To investigate the transcriptional adaptations and tumour–microenvironment interactions associated with breast cancer brain metastasis using an integrative bulk and single-nucleus transcriptomic approach.

## 2.2 Specific objectives

1. Identify differentially expressed genes between primary breast tumours and brain metastases using paired bulk transcriptomic data.
2. Assess the robustness of differential expression findings after adjustment for tumour purity.
3. Validate transcriptional changes using an independent paired breast cancer brain metastasis dataset.
4. Identify reproducible biological pathways associated with metabolic adaptation, ECM organization, and immune-related processes.
5. Characterize the cellular composition of brain metastatic samples using single-nucleus transcriptomic data.
6. Investigate the cellular distribution of selected ECM-associated genes.
7. Construct and evaluate an ECM-associated transcriptional signature across metastatic samples.
8. Examine the relationship between ECM signature variation and cellular composition.
9. Develop an integrated biological model connecting tissue-level transcriptional changes with cellular heterogeneity within brain metastases.

---

# 3. Materials and Methods

## 3.1 Study design

This investigation used a retrospective, computational transcriptomic design integrating publicly available bulk and single-nucleus gene-expression datasets.

The bulk transcriptomic branch compared primary breast tumours with matched brain metastases, incorporated a tumour-purity-adjusted sensitivity analysis, and evaluated findings in an independent validation cohort.

The single-nucleus branch investigated cellular composition, ECM-associated gene expression, and sample-level heterogeneity within brain metastatic lesions.

The two branches were integrated at the biological interpretation level. They were not patient-matched and were not used for direct patient-level correlation.

## 3.2 Bulk transcriptomic datasets

### Discovery cohort: GSE184869

The discovery dataset, GSE184869, consisted of 45 paired primary breast tumour and brain metastasis samples, representing 90 tissue samples.

The supplied expression matrix contained 16,714 protein-coding genes identified by Ensembl gene identifiers.

The expression data were supplied as batch-corrected, log2-transformed, TMM-normalized counts per million (CPM). These preprocessing steps were therefore inherited from the supplied dataset rather than independently repeated as part of the differential-expression analysis.

### Independent validation cohort: GSE125989

GSE125989 was used as an independent validation dataset comprising 16 paired primary tumour and brain metastasis samples.

The validation analysis evaluated the reproducibility of transcriptional changes observed in the discovery cohort.

## 3.3 Paired differential-expression analysis

Differential expression analysis was performed using the limma framework.

The paired discovery model incorporated patient identity and tumour type:

`~ patientID + tumourID`

The contrast of interest represented brain metastasis relative to primary breast tumour.

Empirical Bayes moderation was applied to the fitted model. Multiple testing was addressed using the Benjamini–Hochberg procedure.

Genes with an adjusted p value below 0.05 were considered statistically significant. A more stringent descriptive threshold of absolute log2 fold change greater than 1 was also used to identify genes with larger expression differences.

Positive log2 fold changes indicated higher expression in brain metastases, whereas negative values indicated lower expression.

## 3.4 Tumour-purity-adjusted sensitivity analysis

A second differential-expression model incorporated tumour purity as a continuous covariate:

`~ patientID + tumourID + purity`

This model was used to investigate whether observed expression differences were sensitive to variation in tumour purity.

The purity-adjusted model was fitted using limma with empirical Bayes moderation.

Results were compared with the original paired model in terms of effect direction, effect magnitude, and statistical significance.

Because both models used the same discovery cohort, the purity-adjusted analysis was considered a sensitivity analysis rather than an independent validation experiment.

## 3.5 Independent validation and consensus integration

The GSE125989 dataset was analyzed separately to evaluate reproducibility of differential-expression findings.

Direction-consistent genes were identified by comparing the direction of expression changes across the discovery, purity-adjusted, and validation analyses.

Genes significant in all three analyses were identified as a more stringent subset of reproducible differential-expression findings.

## 3.6 Functional pathway enrichment

Gene set enrichment analysis was performed using ranked differential-expression statistics.

The discovery analysis used moderated t-statistics to rank genes. Ensembl identifiers were stripped of version suffixes and mapped to Entrez identifiers using clusterProfiler and org.Hs.eg.db.

When multiple Ensembl identifiers mapped to the same Entrez identifier, the entry with the largest absolute moderated t-statistic was retained.

Gene Ontology Biological Process and KEGG pathway enrichment were conducted using gseGO and gseKEGG, respectively. Gene-set sizes between 10 and 500 were considered, with Benjamini–Hochberg adjustment.

Pathway reproducibility was assessed by integrating results from discovery, purity-adjusted, and validation analyses. Pathways were matched by database identifiers.

Robust consensus pathways were defined by statistical significance in all three analyses and consistent normalized enrichment score direction.

Leading-edge analysis was performed on seven selected pathways representing metabolic, ECM-associated, and immune-related biological processes.

## 3.7 Single-nucleus transcriptomic dataset

The single-nucleus branch used the annotated Seurat object from GSE324453.

The analyzed dataset contained:

- 56,268 nuclei.
- 35,259 features.
- 28 metastatic samples.
- 12 annotated major cellular populations.

These populations included tumour cells, proliferating tumour cells, myeloid cells, T/NK cells, stromal cells, pericytes, endothelial cells, neurons, oligodendrocyte-lineage cells, astrocytes, plasma cells, and B cells.

The supplied annotated object was used as the basis for downstream cellular composition and ECM-associated analyses.

## 3.8 Cellular composition analysis

Cellular composition was evaluated using sample-level cell population proportions.

For each sample, the proportion of each annotated cellular population was calculated relative to the analyzed nuclei.

Mean proportions were subsequently calculated across the 28 samples, with each sample contributing equally to the average.

Variability was represented using the standard deviation of sample-level proportions.

These values describe the relative composition of recovered and analyzed nuclei rather than absolute cell abundance within the original tissue.

## 3.9 ECM-associated gene expression and cellular origin

Thirteen ECM-associated genes were investigated:

COL6A3, COL14A1, COL15A1, AEBP1, SULF1, THBS2, THBS4, LUM, POSTN, MMP2, CTSK, CXCL12, and DPT.

Expression patterns were examined across annotated cellular populations to identify potential cellular sources of ECM-associated transcription.

Detection prevalence was used to characterize the proportion of nuclei expressing individual genes.

The analysis focused on cellular distribution and gene-specific expression patterns rather than direct evidence of ECM protein deposition or extracellular matrix architecture.

## 3.10 Pseudobulk analysis

Sample-by-cell-type pseudobulk analysis was included to investigate aggregated expression patterns across cellular compartments.

Pseudobulk analysis provides a sample-level framework for examining gene-expression patterns within annotated cell populations.

The precise aggregation and normalization parameters should be reported according to the finalized computational workflow.

## 3.11 ECM signature construction and scoring

A seven-gene core ECM-associated signature was constructed using:

- COL6A3
- AEBP1
- SULF1
- THBS2
- LUM
- POSTN
- CTSK

An expanded 13-gene panel was also considered, incorporating additional vascular-associated and ECM-associated genes.

The core signature was used for the principal sample-level heterogeneity analysis.

An updated scoring workflow was used for the final ECM score analysis.

## 3.12 Sample-level heterogeneity analysis

ECM scores were examined across metastatic samples to characterize inter-sample variation.

The relationship between ECM scores and stromal-cell proportion was evaluated using Spearman correlation.

Samples were divided into high- and low-scoring groups using the median ECM score for descriptive and comparative analyses.

Differences in cellular composition and individual ECM-associated gene expression were evaluated between the two groups, with multiple-testing correction applied where appropriate.

The high- and low-scoring groups were analytical categories derived from a continuous distribution, not independently established biological subtypes.

## 3.13 Statistical analysis and reproducibility

Differential-expression and pathway analyses incorporated multiple-testing correction using the Benjamini–Hochberg procedure.

Statistical findings were interpreted according to their effect sizes, adjusted p values, consistency across analyses, and biological context.

Analyses were performed in R using established Bioconductor and single-cell transcriptomic packages.

Exact software versions, package versions, and the statistical procedures used in individual single-nucleus comparisons should be documented from the finalized computational environment.

---

# 4. Results

## 4.1 Integrated transcriptomic landscape of breast cancer brain metastasis

**Figure 1. Integrated transcriptomic landscape of breast cancer brain metastasis.**

The paired discovery analysis identified extensive transcriptional differences between primary breast tumours and brain metastases.

In the discovery cohort, 5,247 genes were statistically significant at an adjusted p value below 0.05, comprising 2,816 upregulated and 2,431 downregulated genes.

Applying the additional absolute log2 fold-change threshold greater than 1 identified 665 genes, including 146 upregulated and 519 downregulated genes.

The purity-adjusted analysis identified 242 statistically significant genes, comprising 101 upregulated and 141 downregulated genes. Under the more stringent fold-change threshold, 185 genes remained significant.

Several genes showed consistent direction across the original and purity-adjusted analyses.

Among the notable downregulated genes were HTRA3, MYH1, MYH4, MYH13, MYH2, MYH8, CTSK, PODN, and LRRC15. GFAP showed increased expression.

The independent validation cohort identified 118 significant probes under the specified statistical and fold-change thresholds, corresponding to 87 unique genes.

Notable downregulated validation genes included AEBP1, COL1A1, GREM1, COL1A2, GAS1, SPON1, POSTN, SPARC, IGKC, and DPT.

Cross-analysis integration identified 165 genes with consistent expression direction across the discovery, purity-adjusted, and validation analyses. Of these, 99 were consistently downregulated and 66 consistently upregulated.

Fifteen genes were statistically significant in all three analyses, all showing reduced expression in brain metastases:

AEBP1, DPT, POSTN, CTSK, COL15A1, THBS2, COL6A3, LUM, CXCL12, SULF1, COL14A1, MMP2, IGF1, PDGFRA, and THBS4.

These findings demonstrate reproducible transcriptional differences between primary breast tumours and brain metastases, with a prominent contribution from ECM-associated genes.

## 4.2 Functional pathway landscape

**Figure 2. Functional pathway landscape of breast cancer brain metastasis.**

Pathway enrichment analysis revealed two prominent biological patterns: increased mitochondrial energy-metabolism-associated transcription and reduced expression of ECM-associated and selected immune-related pathways.

### 4.2.1 Metabolic pathway enrichment

Oxidative phosphorylation showed positive enrichment across the discovery, purity-adjusted, and validation analyses, with normalized enrichment scores of approximately +3.14, +2.89, and +2.34, respectively.

Aerobic respiration also showed positive enrichment, with corresponding NES values of approximately +2.98, +2.67, and +2.05.

Electron transport chain-associated processes demonstrated similarly consistent positive enrichment.

These results indicate that mitochondrial energy-metabolism-associated transcription was reproducibly elevated in brain metastases relative to primary breast tumours.

### 4.2.2 ECM-associated pathway suppression

ECM organization showed negative enrichment across the three analyses, with NES values of approximately −3.09, −2.78, and −2.84.

Collagen fibril organization and collagen metabolism also demonstrated consistent negative enrichment.

ECM assembly, ECM–receptor interaction, and integrin-associated pathways showed similar patterns.

These findings support a broad reduction in the expression of ECM-associated transcriptional programmes in brain metastases.

### 4.2.3 Immune and cytokine-associated pathways

Selected immune-related pathways, including cytokine–cytokine receptor interaction, integrin-associated processes, and negative regulation of leukocyte migration, showed reduced enrichment.

The results suggest that the transcriptional landscape of brain metastases differs from that of primary tumours in selected immune-associated and cell-interaction programmes.

However, these results should not be interpreted as direct evidence of reduced immune-cell abundance or impaired immune function.

## 4.3 Cellular composition of brain metastases

**Figure 3. Single-nucleus cellular composition of breast cancer brain metastases.**

Single-nucleus analysis included 56,268 nuclei from 28 metastatic samples.

Tumour cells represented the dominant cellular population, accounting for 32,887 nuclei, followed by proliferating tumour cells with 7,157 nuclei.

Other major populations included myeloid cells, T/NK cells, stromal cells, pericytes, endothelial cells, neurons, oligodendrocyte-lineage cells, astrocytes, plasma cells, and B cells.

When mean cellular proportions were calculated across samples, tumour cells represented approximately 52.38% of analyzed nuclei, followed by:

- Proliferating tumour cells: 12.57%.
- Myeloid cells: 11.28%.
- Stromal cells: 5.16%.
- T/NK cells: 5.03%.
- Pericytes: 3.28%.
- Endothelial cells: 3.16%.

The remaining populations each represented smaller average proportions.

These findings demonstrate that brain metastatic samples are predominantly composed of tumour cells but retain a heterogeneous microenvironment containing immune, stromal, vascular-associated, and neural-associated populations.

The reported mean proportions represent averages across individual samples rather than proportions calculated from the pooled nuclei of all samples.

## 4.4 Cellular distribution of ECM-associated transcription

**Figure 4. Cellular distribution and detection prevalence of ECM-associated genes.**

The 13-gene ECM-associated panel demonstrated distinct expression patterns across the annotated cellular populations.

Stromal cells showed the broadest contribution to the core ECM-associated transcriptional landscape.

Among stromal cells, COL6A3 showed the highest detection prevalence at approximately 96.0%, followed by THBS2 at 81.7%, SULF1 at 73.7%, and AEBP1 at 73.1%.

LUM, POSTN, MMP2, CTSK, COL14A1, CXCL12, and other genes showed varying levels of detection.

Pericytes demonstrated additional ECM-associated expression, including COL14A1, AEBP1, and COL6A3.

Endothelial cells showed relatively prominent expression of COL15A1 and CXCL12, supporting a vascular-associated component of the ECM transcriptional landscape.

Astrocytes exhibited selected ECM-associated transcripts, particularly THBS4, SULF1, and AEBP1.

The detection prevalence of individual genes differed substantially between cell populations.

One notable observation was the low detection prevalence of DPT in stromal cells despite its relatively elevated standardized expression in selected comparisons. This highlights the importance of interpreting expression intensity alongside detection frequency.

Overall, the results indicate that ECM-associated transcription is distributed across multiple non-malignant cell populations, with stromal cells representing a major contributor.

## 4.5 ECM-associated transcriptional heterogeneity across metastatic samples

**Figure 5. Sample-level ECM-associated transcriptional heterogeneity.**

The seven-gene core ECM signature was used to characterize variation across 28 metastatic samples.

ECM scores showed substantial inter-sample variation, with a standard deviation of approximately 0.7697 and a median score of approximately 0.1423.

The distribution was continuous rather than clearly separated into distinct groups.

### 4.5.1 Relationship between ECM score and stromal-cell representation

A positive association was observed between the ECM signature score and stromal-cell proportion.

Spearman correlation demonstrated a rho of 0.623, with an adjusted p value of 0.00236.

This indicates that samples with higher ECM-associated transcriptional scores tended to contain a greater proportion of stromal cells.

However, the association does not establish that changes in stromal-cell abundance are the sole determinant of ECM-associated transcription.

### 4.5.2 Cellular composition differences between high- and low-scoring samples

Samples were divided into high- and low-scoring groups using the median ECM score.

The high-scoring group demonstrated a higher mean stromal-cell proportion than the low-scoring group.

Stromal representation increased from approximately 1.20% in the low-scoring group to 6.79% in the high-scoring group, corresponding to a difference of 5.59 percentage points.

This difference was statistically significant after multiple-testing correction.

Other cellular populations showed variable differences, although not all remained statistically significant after adjustment.

These findings support a relationship between ECM-associated transcriptional scores and stromal-cell representation.

### 4.5.3 Gene-level differences between high- and low-scoring samples

Several individual ECM-associated genes showed higher expression in the high-scoring group.

POSTN, CXCL12, SULF1, CTSK, LUM, THBS2, AEBP1, and COL6A3 showed statistically significant differences after multiple-testing correction.

Other genes showed smaller differences or did not demonstrate statistically significant separation.

DPT displayed a relatively large difference in standardized expression, but its low detection prevalence warrants cautious interpretation.

The observed differences indicate that the ECM-associated signature reflects variation across several genes rather than a uniform change in every member of the expanded ECM panel.

## 4.6 Integrated biological model

**Figure 6. Integrative biological model of breast cancer brain metastasis.**

The findings from bulk and single-nucleus transcriptomic analyses were integrated into a conceptual model connecting tissue-level transcriptional differences with cellular heterogeneity within brain metastases.

At the bulk tissue level, brain metastases demonstrated increased mitochondrial energy-metabolism-associated transcription and reduced expression of several ECM-associated and selected immune-related programmes relative to primary breast tumours.

At the single-nucleus level, metastatic samples showed a tumour-dominated cellular composition with contributions from stromal, vascular-associated, immune, and neural-associated populations.

ECM-associated transcription was distributed across multiple non-malignant cell populations, particularly stromal cells, with additional gene-specific contributions from pericytes, endothelial cells, and astrocytes.

Furthermore, sample-level ECM-associated transcription varied continuously and was partly associated with stromal-cell representation.

Together, these findings provide an integrated description of metabolic transcriptional adaptation and a heterogeneous multicellular ECM-associated landscape in breast cancer brain metastases.

---

# 5. Discussion

## 5.1 Overview of the principal findings

This study provides an integrative transcriptomic characterization of breast cancer brain metastasis by combining bulk and single-nucleus gene-expression analyses.

Four principal findings emerged.

First, brain metastases exhibited increased mitochondrial energy-metabolism-associated transcription relative to primary breast tumours.

Second, ECM-associated and selected immune-related pathways showed reduced bulk transcriptional enrichment.

Third, single-nucleus analysis revealed that ECM-associated transcription was distributed across multiple non-malignant cellular populations, particularly stromal and vascular-associated compartments.

Fourth, ECM-associated transcription varied substantially between metastatic samples and was positively associated with stromal-cell representation.

These observations provide complementary insights into the transcriptional landscape of brain metastasis.

## 5.2 Metabolic transcriptional adaptation in brain metastases

The consistent enrichment of oxidative phosphorylation, aerobic respiration, and electron transport chain-associated pathways represents one of the major findings of this investigation.

The consistency of these signals across the discovery model, purity-adjusted analysis, and independent validation cohort strengthens confidence that the observed pattern is not restricted to a single analytical model.

Mitochondrial energy metabolism is particularly relevant in the context of metastatic adaptation because tumour cells must survive within tissue environments that differ in nutrient availability, cellular interactions, and metabolic conditions.

The observed enrichment suggests that mitochondrial energy-metabolism-associated transcription is a prominent feature of the brain metastatic state.

Nevertheless, the biological significance of this finding requires careful interpretation.

Transcriptomic enrichment does not demonstrate increased mitochondrial respiration or establish that metastatic cells rely exclusively on oxidative phosphorylation.

Additional functional studies, including oxygen-consumption measurements, mitochondrial activity assays, and metabolomic profiling, would be required to establish whether the transcriptional pattern corresponds to altered metabolic activity.

## 5.3 Reduced ECM-associated transcription at the tissue level

The reduction in ECM organization, collagen-associated pathways, ECM–receptor interactions, and integrin-related processes indicates that the transcriptional landscape of brain metastases differs substantially from that of primary breast tumours.

The consistent downregulation of several ECM-associated genes across the discovery, purity-adjusted, and validation analyses supports the reproducibility of this observation.

However, bulk expression measurements combine transcriptional signals from multiple cellular populations.

Consequently, reduced bulk ECM-associated expression may reflect changes in the abundance of ECM-producing cells, altered expression within specific populations, or a combination of these processes.

The present analysis does not directly establish reduced ECM protein deposition, diminished collagen abundance, or changes in matrix stiffness.

Protein-level and spatial studies will be necessary to investigate these possibilities.

## 5.4 Multicellular contributions to the ECM-associated transcriptional landscape

Single-nucleus analysis provided additional cellular resolution by identifying several populations expressing ECM-associated genes.

Stromal cells demonstrated the broadest representation of the core ECM-associated transcriptional programme.

Pericytes and endothelial cells contributed additional gene-specific expression patterns, suggesting that vascular-associated populations may also participate in the transcriptional organization of the metastatic microenvironment.

Astrocytes showed selected ECM-associated transcripts, indicating that neural-associated populations may contribute to the broader transcriptional landscape.

These findings support a multicellular interpretation of ECM-associated transcription within brain metastases.

Importantly, the data identify potential cellular sources of ECM-associated transcripts rather than directly demonstrating the deposition, organization, or functional activity of the corresponding extracellular proteins.

## 5.5 Relationship between ECM-associated transcription and stromal representation

The positive association between ECM signature scores and stromal-cell proportions suggests that cellular composition contributes to variation in ECM-associated transcription across metastatic samples.

The high- and low-scoring comparisons further demonstrated a substantial difference in stromal representation.

However, the relationship is unlikely to be explained solely by cellular abundance.

Variation in gene expression within individual cell populations may also contribute to the observed differences.

The present analysis does not fully separate cell-composition effects from cell-intrinsic transcriptional variation.

Future analyses incorporating cell-type-specific pseudobulk modelling and appropriate sample-level statistical designs could help distinguish these contributions.

## 5.6 Inter-sample heterogeneity rather than discrete ECM subtypes

The continuous distribution of ECM signature scores suggests that metastatic samples vary along a quantitative transcriptional gradient.

Although high- and low-scoring groups were useful for comparative analysis, the data do not establish the existence of two distinct ECM-defined biological subtypes.

This distinction is important because dichotomizing a continuous score can obscure intermediate states and exaggerate apparent differences between groups.

The findings are therefore more appropriately interpreted as evidence of inter-sample heterogeneity in ECM-associated transcription.

## 5.7 Reconciling bulk and single-nucleus findings

An important aspect of this study is the apparent difference between reduced ECM-associated expression in bulk tissue comparisons and detectable ECM-associated transcription across multiple cell populations within brain metastases.

These observations are not contradictory because the two analyses address different biological questions.

Bulk transcriptomics compares primary breast tumours with brain metastases and measures aggregate tissue-level expression.

Single-nucleus transcriptomics examines the cellular composition and gene-expression patterns within brain metastatic samples.

Thus, reduced ECM-associated transcription at the tissue-comparison level can coexist with localized expression of selected ECM-associated genes in particular cellular populations within metastatic lesions.

This distinction highlights the value of integrating bulk and single-nucleus approaches when investigating complex tumour microenvironments.

## 5.8 Biological implications and future directions

The findings generate several hypotheses for future investigation.

First, the increased mitochondrial energy-metabolism-associated transcription suggests that mitochondrial function may be an important component of the brain metastatic phenotype. Functional metabolic experiments could determine whether this transcriptional pattern corresponds to altered respiratory activity.

Second, the multicellular distribution of ECM-associated transcription suggests that interactions between tumour cells, stromal cells, and vascular-associated populations may contribute to the metastatic microenvironment.

Third, the association between ECM scores and stromal representation suggests that tissue composition and cell-intrinsic expression should be investigated together.

Finally, spatial transcriptomic and protein-level analyses could determine whether the cellular distribution of ECM-associated transcripts corresponds to localized matrix deposition and organization.

Together, these directions could extend the present computational findings toward mechanistic understanding.

---

# 6. Limitations

Several limitations should be considered when interpreting the findings.

## 6.1 Publicly available retrospective datasets

The investigation relied on publicly available datasets generated by different research studies.

Differences in sample collection, sequencing platforms, preprocessing procedures, and annotation strategies may introduce technical variation.

## 6.2 Limited functional validation

The study used transcriptomic measurements and did not directly assess mitochondrial respiration, ATP production, metabolic flux, ECM protein abundance, collagen deposition, or matrix mechanical properties.

Therefore, biological interpretations should remain at the transcriptional level unless supported by independent experimental evidence.

## 6.3 Bulk tissue cellular composition

Bulk transcriptomic measurements combine signals from malignant and non-malignant populations.

Although tumour-purity adjustment was performed, this does not fully resolve cellular-composition effects or identify the cellular origins of individual expression changes.

## 6.4 Single-nucleus RNA sequencing limitations

Single-nucleus sequencing captures nuclear transcripts rather than the complete cellular transcriptome.

Detection of individual genes may therefore be influenced by transcript localization, sequencing depth, nuclear RNA abundance, and technical dropout.

The cellular proportions obtained from the analyzed nuclei should not be interpreted as absolute tissue-cell abundance.

## 6.5 Incomplete cellular resolution of ECM biology

The expression of an ECM-associated gene in a particular cell population does not establish that the corresponding protein is secreted, deposited, or functionally active in the extracellular environment.

Protein-level and spatial validation are required to establish these relationships.

## 6.6 ECM signature limitations

The seven-gene core signature provides a focused representation of selected ECM-associated transcriptional features.

It does not capture the entire complexity of ECM biology.

The observed score may also be influenced by cellular composition, gene-specific expression, and the characteristics of the analyzed single-nucleus dataset.

## 6.7 Interpretation of high- and low-scoring groups

The high- and low-scoring groups were defined using the median ECM score for comparative purposes.

They should not be interpreted as independently validated molecular subtypes.

## 6.8 Cohort and patient-level integration

The bulk and single-nucleus datasets were not patient-matched.

Consequently, the study cannot establish direct patient-level relationships between bulk differential expression and single-nucleus cellular composition.

The integration was performed at the level of biological interpretation rather than matched-sample correlation.

---

# 7. Future Directions

The findings suggest several opportunities for further investigation.

### 7.1 Functional validation of metabolic adaptation

- Measure mitochondrial oxygen-consumption rates.
- Evaluate respiratory chain activity.
- Investigate ATP production and metabolic flux.
- Examine mitochondrial-associated transcriptional programmes using independent experimental systems.

### 7.2 Spatial characterization of ECM-associated populations

- Apply spatial transcriptomic approaches to brain metastatic tissue.
- Investigate the localization of stromal and vascular-associated populations.
- Examine whether ECM-associated transcripts correspond to localized protein deposition.

### 7.3 Protein-level ECM characterization

- Evaluate selected collagen and matricellular proteins using immunohistochemistry or immunofluorescence.
- Examine ECM abundance and distribution.
- Investigate matrix organization and mechanical properties.

### 7.4 Cell-type-specific transcriptional modelling

- Perform more detailed pseudobulk analyses.
- Investigate cell-intrinsic expression changes within stromal, pericyte, endothelial, and tumour populations.
- Separate cellular abundance effects from within-cell transcriptional differences.

### 7.5 Clinical and biological validation

Future studies involving larger independent cohorts could evaluate whether metabolic and ECM-associated transcriptional patterns are associated with clinical outcomes, metastatic characteristics, or treatment response.

Such analyses would require appropriately matched clinical information and independent validation.

---

# 8. Conclusion

This integrative transcriptomic investigation identifies coordinated metabolic and tumour–microenvironment-associated transcriptional patterns in breast cancer brain metastasis.

Bulk transcriptomic analyses demonstrate increased mitochondrial energy-metabolism-associated expression alongside reduced expression of several ECM-associated and selected immune-related pathways compared with primary breast tumours.

Single-nucleus analysis reveals that ECM-associated transcription is distributed across multiple cellular populations within brain metastases, particularly stromal cells, with additional contributions from pericytes, endothelial cells, and astrocytes.

Furthermore, ECM-associated transcription varies substantially between metastatic samples and is partly associated with stromal-cell representation.

The integration of tissue-level and cellular-resolution transcriptomic findings supports a model in which metabolic transcriptional adaptation coexists with a heterogeneous, multicellular ECM-associated landscape in breast cancer brain metastases.

Although functional and spatial validation remain necessary, the findings provide a framework for understanding the relationship between metabolic adaptation, cellular composition, and ECM-associated transcription in the brain metastatic microenvironment.

---

# 9. Data Availability

The publicly available datasets analyzed in this study are:

- Discovery bulk transcriptomic dataset: GEO accession GSE184869.
- Independent bulk validation dataset: GEO accession GSE125989.
- Single-nucleus transcriptomic dataset: GEO accession GSE324453.

The computational analyses were conducted using R and relevant Bioconductor and single-cell analysis packages.

Project-specific scripts, processed results, and figures are maintained within the project repository, subject to the applicable data-use and redistribution conditions.

# 10. Figure Legends

### Figure 1. Integrated transcriptomic landscape of breast cancer brain metastasis.

(A) Overview of the integrative study design, including bulk differential expression, purity-adjusted analysis, independent validation, pathway integration, and single-nucleus analysis. (B) Differential-expression landscape from the paired discovery model. (C) Differential-expression landscape after adjustment for tumour purity. (D) Integration of direction-consistent differential-expression findings across discovery, purity-adjusted, and validation analyses.

### Figure 2. Functional pathway landscape of breast cancer brain metastasis.

(A) Global pathway enrichment landscape across the discovery, purity-adjusted, and validation analyses. (B) Enrichment of mitochondrial energy-metabolism-associated pathways, including oxidative phosphorylation and aerobic respiration. (C) Reduced enrichment of ECM organization, collagen-associated, and ECM interaction pathways. (D) Selected immune- and cytokine-associated pathway patterns.

### Figure 3. Cellular composition of breast cancer brain metastases.

Single-nucleus cellular composition across 28 metastatic samples and 56,268 nuclei. Sample-level cellular proportions are summarized across the annotated populations, with mean proportions calculated using equal sample weighting and variability represented by the standard deviation across samples.

### Figure 4. Cellular distribution of ECM-associated transcription in brain metastases.

(A) Cellular expression landscape of selected ECM-associated genes across annotated populations. (B) Detection prevalence of individual ECM-associated genes within selected cellular compartments. The analysis identifies stromal, pericyte, endothelial, and astrocyte-associated expression patterns.

### Figure 5. Sample-level ECM-associated transcriptional heterogeneity.

(A) Distribution of the core ECM signature across metastatic samples. (B) Relationship between ECM score and stromal-cell proportion. (C) Cellular composition comparison between high- and low-scoring samples. (D) Gene-level differences in ECM-associated expression between the two groups.

### Figure 6. Integrative biological model of breast cancer brain metastasis.

Conceptual summary of the principal findings, illustrating increased mitochondrial energy-metabolism-associated transcription, reduced bulk ECM-associated transcription, the multicellular distribution of ECM-associated gene expression, and inter-sample variation in ECM-associated transcription. The model summarizes associations and does not imply experimentally established causal mechanisms.

# 11. References and Bibliographic Completion

The final reference list should include:

1. The original publication associated with GSE184869.
2. The original publication associated with GSE125989.
3. The original publication associated with GSE324453.
4. The limma methodology reference.
5. The clusterProfiler methodology reference.
6. The Seurat methodology reference, where applicable.
7. Relevant literature on breast cancer brain metastasis, oxidative phosphorylation, ECM biology, and tumour–microenvironment interactions.

Full author names, publication titles, journal information, years, and DOI identifiers should be retrieved and checked against the original publications before final submission.

---

**End of project report draft**