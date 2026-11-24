# Results

## 1. Transcriptomic divergence between primary breast tumours and brain metastases

To characterize transcriptional changes associated with breast cancer brain metastasis (BCBM), paired primary breast tumour and brain metastatic samples were analyzed using a differential expression framework accounting for patient-specific variation. In the discovery cohort (GSE184869; 45 paired samples), 5,247 genes were significantly differentially expressed between primary tumours and brain metastases at a false discovery rate (FDR) of 0.05, comprising 2,816 upregulated and 2,431 downregulated genes in brain metastases. Applying a more stringent threshold of FDR < 0.05 and absolute log₂ fold change > 1 identified 665 genes, including 146 upregulated and 519 downregulated genes.

Several of the most prominent transcriptional changes involved genes associated with extracellular matrix biology and tissue-specific expression. HTRA3, multiple skeletal muscle-associated myosin genes, CTSK, PODN, and LRRC15 showed marked downregulation in brain metastases, whereas GFAP exhibited strong upregulation. These changes indicate substantial transcriptional divergence between primary breast tumours and their brain metastatic counterparts.

To assess whether these differences were robust to tumour purity, a sensitivity analysis incorporating purity estimates into the differential expression model was performed. Although the number of significant genes decreased substantially following purity adjustment, important transcriptional changes persisted. HTRA3, CTSK, multiple MYH genes, LRRC15, and PODN remained downregulated, while GFAP remained strongly upregulated. This persistence suggests that the observed transcriptional differences are not entirely explained by variation in tumour purity.

Independent validation using the paired cohort GSE125989 identified 87 unique genes meeting the stringent differential expression criteria, including 81 downregulated and 6 upregulated genes. The validation cohort showed prominent reductions in ECM-associated genes, including AEBP1, COL1A1, COL1A2, GREM1, POSTN, SPON1, SPARC, and DPT.

Integration of the discovery, purity-adjusted, and validation analyses identified 165 direction-consistent genes, including 99 downregulated and 66 upregulated genes. Fifteen genes reached statistical significance across all three analytical comparisons, all showing downregulation in brain metastases. These included several ECM-associated genes, such as AEBP1, COL14A1, COL15A1, COL6A3, CTSK, CXCL12, DPT, LUM, MMP2, POSTN, SULF1, THBS2, THBS4, IGF1, and PDGFRA.

Together, these analyses demonstrate extensive transcriptomic divergence between primary breast tumours and brain metastases, with a reproducible pattern of downregulation involving multiple ECM-associated genes. The persistence of these changes in purity-adjusted and independent validation analyses supports the robustness of the identified transcriptional alterations.


## 2. Functional pathway remodeling in breast cancer brain metastases

To characterize the biological processes associated with the transcriptional differences identified between primary breast tumours and brain metastases, pathway enrichment analysis was performed using the consensus differential expression framework. The resulting enrichment landscape revealed coordinated changes involving mitochondrial energy metabolism, extracellular matrix organization, and immune-associated signaling.

### 2.1. Enhanced mitochondrial energy-metabolism programs in brain metastases

Brain metastases exhibited positive enrichment of mitochondrial energy-metabolism pathways, particularly oxidative phosphorylation, aerobic respiration, and respiratory electron transport. These pathways showed consistent positive enrichment across the discovery, purity-adjusted, and validation analyses, supporting a reproducible metabolic-associated transcriptional shift.

Oxidative phosphorylation showed normalized enrichment scores (NES) of 3.14, 2.89, and 2.34 across the three analytical comparisons. Aerobic respiration was similarly enriched, with NES values of 2.98, 2.67, and 2.05, while respiratory electron transport showed NES values of 3.15, 2.62, and 2.16. The consistency of these enrichment patterns suggests that mitochondrial energy-metabolism-associated transcription is a prominent feature of the brain metastatic phenotype.

These findings indicate a shift toward enhanced expression of genes associated with mitochondrial respiration and oxidative energy metabolism in brain metastases. However, pathway enrichment reflects coordinated transcriptomic changes and does not directly establish increased respiratory activity, ATP production, or mitochondrial function.

### 2.2. Reduced extracellular matrix and collagen-associated transcription

In contrast to the metabolic pathways, extracellular matrix-related pathways showed consistent negative enrichment in brain metastases. ECM organization, collagen fibril organization, collagen metabolism, and ECM assembly were among the prominent negatively enriched processes across the discovery, purity-adjusted, and validation analyses.

ECM organization exhibited NES values of −3.09, −2.78, and −2.84 across the respective comparisons. Collagen fibril organization showed NES values of −2.75, −2.39, and −2.84, while collagen metabolism showed NES values of −2.81, −2.54, and −2.36. ECM assembly was also consistently negatively enriched, with NES values of −2.39, −1.95, and −2.51.

The concordant suppression of multiple related ECM and collagen-associated pathways, together with the downregulation of individual ECM-associated genes in the differential expression analyses, indicates a reproducible reduction in bulk ECM-associated transcription in brain metastases relative to primary breast tumours.

### 2.3. Reduced immune-associated and cytokine signaling programs

Pathway enrichment analysis also revealed negative enrichment of several immune-associated processes, including cytokine signaling, integrin-related pathways, and leukocyte migration-associated programs. Cytokine–receptor interaction showed NES values of −2.86, −1.44, and −1.98 across the three analytical comparisons. Integrin-related signaling showed NES values of −2.59, −1.74, and −1.87, while leukocyte migration-associated pathways were also negatively enriched, with NES values of −2.45, −1.71, and −2.09.

These findings indicate that brain metastases are associated with broad changes in immune-related and cell-interaction transcriptional programs. However, because bulk transcriptomic measurements reflect both expression changes and cellular composition, these pathway-level differences cannot be attributed exclusively to altered signaling activity within tumour cells.

### 2.4. Integrated pathway-level interpretation

Collectively, the pathway enrichment results reveal three major features of the brain metastatic transcriptomic landscape: increased mitochondrial energy-metabolism-associated expression, reduced ECM and collagen-associated expression, and reduced immune-associated signaling. These coordinated changes provide a functional context for the differential expression results and motivate subsequent single-nucleus analyses to investigate the cellular composition and distribution of ECM-associated transcription within brain metastases.


## 3. Cellular composition of breast cancer brain metastases
This section transitions from bulk transcriptomics to the single-nucleus dataset. Its purpose is to characterize the cellular landscape in which the ECM-associated transcriptional signals occur.

### 3.1. Single-nucleus dataset overview
The single-nucleus transcriptomic dataset GSE324453 comprised 56,268 nuclei from 28 breast cancer brain metastasis samples, annotated into 12 cellular populations.

The dataset was dominated by tumour cells, with additional proliferating tumour cells and multiple microenvironmental populations, including myeloid, stromal, T/NK, vascular, perivascular, and brain-associated cell types.

### 3.2. Tumour-dominated cellular composition
Tumour cells represented the largest population across the dataset, accounting for a mean of 52.38% of nuclei per sample. Proliferating tumour cells constituted an additional 12.57%, while myeloid cells represented 11.28%. Stromal and T/NK populations contributed mean proportions of 5.16% and 5.03%, respectively.

Other populations were present at comparatively lower mean proportions, including pericytes (3.28%), endothelial cells (3.16%), oligodendrocyte-lineage cells (2.31%), plasma cells (1.81%), neurons (1.70%), astrocytes (1.01%), and B cells (0.32%).

**Table 1. Mean cellular composition across samples**

|   Cell population    | Mean proportion |
| :------------------: | :-------------: |
|     Tumour cells     |     52.38%      |
| Tumour cells (prol.) |     12.57%      |
|    Myeloid cells     |     11.28%      |
|    Stromal cells     |      5.16%      |
|      T/NK cells      |      5.03%      |
|      Pericytes       |      3.28%      |
|  Endothelial cells   |      3.16%      |
|         ODG          |      2.31%      |
|     Plasma cells     |      1.81%      |
|       Neurons        |      1.70%      |
|      Astrocytes      |      1.01%      |
|       B cells        |      0.32%      |

The cellular landscape was therefore characterized by a predominant tumour compartment alongside a diverse collection of non-tumour populations, providing cellular context for the transcriptomic programs identified in the bulk analysis.

### 3.3. Inter-sample variability in cellular composition
Although tumour cells constituted the predominant cellular population overall, their relative abundance varied across individual brain metastasis samples. The sample-level stacked composition plot and cellular composition heatmap revealed differences in the relative contributions of tumour, proliferating tumour, and microenvironmental populations among samples.

This variability indicates that the single-nucleus cohort does not represent a uniform cellular landscape. Individual brain metastases differ in their relative cellular composition, including the abundance of stromal, myeloid, vascular, lymphoid, and brain-associated populations.

### 3.4. Integrated interpretation
Together, the single-nucleus analysis demonstrates that breast cancer brain metastases are predominantly composed of tumour cells but retain a heterogeneous microenvironment containing multiple non-tumour cellular populations. Although stromal cells represent a relatively small fraction of nuclei on average, their presence alongside vascular, perivascular, immune, and brain-associated populations provides a cellular framework for investigating ECM-associated transcription within the metastatic tissue.

The observed variability in cellular composition also highlights an important consideration for interpreting bulk transcriptomic data, in which gene-expression differences may reflect both changes in transcriptional activity and differences in the relative abundance of cellular populations.


### ECM-associated transcription is predominantly localized to non-tumour cellular compartments

To determine which cellular populations contribute to the ECM-associated transcriptional signature observed in breast cancer brain metastases, we examined the expression patterns and detection prevalence of 13 ECM-associated genes across the annotated cell populations. The analysis revealed a predominantly non-tumour-associated distribution, with stromal cells exhibiting the most prominent and broadly represented ECM-associated transcriptional signature.

Stromal cells showed relatively high expression and broad detection of several core ECM-associated genes, including **COL6A3, THBS2, SULF1, AEBP1, LUM, POSTN, and CTSK**. Among these, COL6A3 was detected in 96.0% of stromal cells, followed by THBS2 (81.7%), SULF1 (73.7%), and AEBP1 (73.1%). LUM, POSTN, and MMP2 were also represented across substantial fractions of the stromal population. These findings identify stromal cells as the principal cellular compartment associated with the ECM transcriptional signature in this dataset.

In contrast, tumour and proliferating tumour populations generally showed lower relative expression and limited detection of most ECM-associated genes, with detection frequencies predominantly below approximately 15.5%. This distribution suggests that the ECM-associated transcriptional signal observed across brain metastasis samples is largely associated with non-tumour cellular populations rather than being uniformly expressed by the tumour compartment.

Several additional populations exhibited gene-specific ECM-associated expression patterns, indicating that the signature was not exclusively restricted to stromal cells. Pericytes showed relatively prominent expression of COL14A1, AEBP1, and COL6A3, with detection frequencies of 36.8%, 43.9%, and 40.0%, respectively. Endothelial cells exhibited relatively prominent COL15A1, CXCL12, and MMP2 expression, detected in 56.3%, 41.3%, and 30.2% of cells, respectively. Astrocytes showed a distinct pattern characterized by THBS4, SULF1, AEBP1, and THBS2 expression, with THBS4 detected in 46.9% of astrocytes.

The complementary expression heatmap and detection-prevalence analysis further highlighted the distinction between relative expression intensity and the proportion of cells expressing each gene. Although some genes exhibited comparatively high expression within a particular population, their detection prevalence remained limited. For example, DPT showed relatively prominent stromal expression in the expression heatmap but was detected in only approximately 2.7% of stromal cells, suggesting that its signal may reflect a relatively restricted cellular subset rather than a broadly distributed stromal expression pattern.

Together, these findings indicate that ECM-associated transcription in breast cancer brain metastases is predominantly associated with stromal cells, accompanied by distinct gene-specific expression patterns in pericytes, endothelial cells, and astrocytes. This cellular distribution provides a context for interpreting the ECM-associated signal observed in bulk transcriptomic analyses, although transcriptional expression alone does not establish extracellular protein abundance, deposition, or functional remodelling of the tissue matrix.


### Sample-level heterogeneity of ECM-associated transcription

Having established the cellular distribution of ECM-associated transcription, we next examined its variation across individual brain metastasis samples. An ECM-associated transcriptional score calculated using a core set of stromal ECM-associated genes demonstrated substantial inter-sample variability, indicating that the ECM-associated transcriptional landscape was not uniform across the cohort. The score showed a continuous distribution rather than a clear separation into distinct molecular subgroups.

The sample-level ECM score was positively associated with stromal cell abundance (Spearman’s ρ = 0.623, FDR = 0.00236), indicating that samples with higher ECM-associated transcriptional scores generally contained a greater proportion of stromal cells. Nevertheless, the association did not establish that stromal abundance alone accounted for the observed variation, suggesting that differences in cellular composition and transcriptional characteristics may both contribute to the sample-level ECM signal.

To further characterize this heterogeneity, samples were stratified into high- and low-ECM-score groups. The high-ECM group exhibited a significantly greater stromal proportion than the low-ECM group (6.79% versus 1.20%; difference = 5.59 percentage points, FDR = 0.00583). Endothelial proportions also showed a modest increase in the high-ECM group, although this difference did not remain statistically significant after multiple-testing correction. Differences in myeloid and pericyte proportions were similarly suggestive but did not reach adjusted statistical significance. Conversely, the proportion of tumour cells was lower in the high-ECM group, although the difference did not meet the adjusted significance threshold. These observations indicate that elevated ECM-associated transcription was accompanied most consistently by increased stromal representation, while differences in other cellular populations were less conclusive.

Gene-level comparisons further characterized the transcriptional differences between the high- and low-ECM groups. POSTN, CXCL12, SULF1, CTSK, LUM, THBS2, AEBP1, and COL6A3 exhibited significant differences between groups after multiple-testing correction. POSTN, CXCL12, SULF1, and CTSK showed particularly prominent differences, while several additional ECM-associated genes demonstrated comparatively modest changes. These findings indicate that the sample-level ECM-associated signature reflected variation across multiple constituent genes rather than a uniform change in every ECM-associated transcript.

Although DPT exhibited a comparatively large difference between the high- and low-ECM groups, its interpretation requires caution because its detection was limited in the stromal population. Moreover, the core ECM score was constructed using a predefined set of ECM-associated genes; therefore, comparisons of these same genes between score-defined groups are descriptive and are not independent validation of the score.

Collectively, these findings demonstrate substantial sample-level heterogeneity in ECM-associated transcription across breast cancer brain metastases. This heterogeneity was associated with stromal abundance and differences in selected ECM-associated genes, supporting a continuum of ECM-associated transcriptional variation rather than the presence of clearly defined ECM molecular subtypes.

### Integrative interpretation of bulk and single-nucleus transcriptomic findings

Integration of the bulk transcriptomic and single-nucleus analyses revealed complementary aspects of breast cancer brain metastasis-associated transcriptional remodeling. Bulk transcriptomic comparisons between primary breast tumours and brain metastases identified prominent differences in metabolic, extracellular matrix-associated, and immune-related pathways, while single-nucleus analyses provided cellular context for interpreting these changes within brain metastatic samples.

Across bulk transcriptomic comparisons, brain metastases exhibited increased expression of mitochondrial energy-metabolism-associated programs, particularly oxidative phosphorylation, aerobic respiration, and electron transport chain pathways. In contrast, ECM organization, collagen-associated processes, and several immune- and cytokine-related pathways showed reduced enrichment. These findings characterize a transcriptional shift associated with the metastatic tissue context, although pathway enrichment does not directly establish metabolic activity or functional pathway output.

Single-nucleus analyses demonstrated that brain metastatic samples were predominantly composed of tumour cells but also contained stromal, myeloid, vascular, lymphoid, and neural-associated populations. Within this cellular landscape, ECM-associated transcription was predominantly localized to stromal cells, with additional gene-specific expression patterns observed in pericytes, endothelial cells, and astrocytes. This cellular distribution provides a context for interpreting the reduced ECM-associated signal observed in bulk comparisons, indicating that bulk transcriptional differences may reflect both changes in expression and differences in the relative representation of ECM-associated cellular populations.

The single-nucleus analysis further revealed substantial variation in ECM-associated transcription between individual brain metastasis samples. Higher ECM-associated scores were associated with greater stromal abundance and differences in the expression of selected ECM-associated genes. Rather than defining distinct ECM molecular subtypes, these findings support a continuous spectrum of ECM-associated transcriptional variation across the analysed samples.

Taken together, the integrated analyses highlight three interconnected features of breast cancer brain metastases: altered metabolic-associated transcription, differences in ECM- and immune-related transcription relative to primary breast tumours, and spatially contextualized ECM-associated transcription within a heterogeneous metastatic cellular environment. These observations provide an integrative transcriptomic framework for understanding the molecular and cellular heterogeneity of breast cancer brain metastases, while remaining associative and not establishing causal mechanisms of metastatic adaptation.








