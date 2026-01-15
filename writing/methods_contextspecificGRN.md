# Methods: Computational Pipeline for Constructing the Gene Regulatory Network

## Overview

We previously have developed a computational pipeline to infer putative regulatory relationships between transcription factors (TFs) and their target genes by leveraging single-cell multimodal data. The pipeline constructs a gene × TF matrix that quantifies the regulatory potential of each TF on each gene by linking three orthogonal data sources: (i) correlations between gene expression and cis-regulatory peak accessibility at the metacell level, (ii) TF motif binding predictions from sequence-based analysis, and (iii) in-silico ChIP-seq scores that capture the functional relationship between TF expression and chromatin accessibility dynamics. The pipeline is implemented as a Snakemake workflow, ensuring reproducibility and modularity across all computational steps.

## Pipeline Architecture

The pipeline consists of four sequential modules executed in the following order:

1. **Gene-Peak Correlation Module** (`gp_corr`): Computes correlations between gene expression and peak accessibility
2. **Peak-TF Binding Module** (`peak_tf`): Predicts TF binding sites within accessible chromatin regions
3. **In-Silico ChIP Module** (`compute_ins_chip`): Integrates TF expression with peak accessibility to compute functional binding scores
4. **Gene × TF Assembly Module** (`gene_x_tf`): Constructs the final regulatory matrix through matrix multiplication

Each module is described in detail below.

---

## Module 1: Gene-Peak Correlation Analysis

Gene-peak correlations were computed using the SEACells framework<sup>3</sup>. Briefly, SEACells aggregates single-cell data into metacells—groups of highly similar cells that reduce sparsity while preserving biological heterogeneity. For each gene, candidate cis-regulatory peaks within ±100 kb of the transcript boundaries were identified using genomic annotations (hg38). Pearson correlation coefficients were then computed between normalized, log-transformed gene expression and peak accessibility across metacells, with statistical significance assessed via two-tailed *P*-values.

For downstream GRN construction, we applied additional filtering criteria: gene-peak associations were retained only if they exhibited positive correlation (*r* > 0) and met the significance threshold (*P* < 0.1). Genes with fewer than two significantly associated peaks were excluded to ensure robust regulatory inference.

---

## Module 2: TF Motif Binding Prediction

TF binding site predictions were generated using FIMO (Find Individual Motif Occurrences) from the MEME Suite<sup>1</sup>. Peak sequences (150 bp centered on summit) were scanned against TF position weight matrices from the CIS-BP database<sup>2</sup>. For each peak *p* and TF *t*, we converted the FIMO-reported *P*-value to a motif binding score: *F*<sub>*p*,*t*</sub> = −log<sub>10</sub>(*P*<sub>FIMO</sub>). Results were assembled into a peak × TF matrix **F**.

---

## Module 3: In-Silico ChIP-seq Score Computation

### Rationale

While sequence-based motif predictions identify potential TF binding sites, they do not capture the functional activity of TF binding in specific cellular contexts. To address this limitation, we compute in-silico ChIP-seq scores<sup>3</sup> that integrate three signals: (i) the presence of a TF binding motif, (ii) the chromatin accessibility of the peak, and (iii) the correlation between TF expression and peak accessibility across cell states.

### Computational Procedure

For each TF *t*, candidate peaks were restricted to those with predicted binding sites (*F*<sub>*p*,*t*</sub> > 0). Using the SEACells metacell framework<sup>3</sup>, we computed Spearman rank correlations (*ρ*<sub>*p*,*t*</sub>) between TF expression and peak accessibility across metacells. Only positive correlations were retained to focus on activating regulatory relationships (*ρ*<sup>+</sup><sub>*p*,*t*</sub> = max(*ρ*<sub>*p*,*t*</sub>, 0)).

The in-silico ChIP-seq score for each peak-TF pair integrates three components:

*I*<sub>*p*,*t*</sub> = *ρ*<sup>+</sup><sub>*p*,*t*</sub> × *A*<sup>max</sup><sub>*p*</sub> × *F*<sup>norm</sup><sub>*p*,*t*</sub>

where *A*<sup>max</sup><sub>*p*</sub> is the maximum accessibility of peak *p* across all metacells, and *F*<sup>norm</sup><sub>*p*,*t*</sub> is the FIMO score normalized per TF. This multiplicative formulation ensures high scores are assigned only to TF-peak pairs satisfying all three criteria: strong positive TF-accessibility correlation, high chromatin accessibility, and high-confidence motif match.

---

## Module 4: Gene × TF Matrix Assembly

The final gene × TF regulatory matrix integrates gene-peak associations from Module 1 with in-silico ChIP scores from Module 3.

### Computational Procedure

Gene-peak associations are filtered to retain only those with positive correlation (*r* > 0) and statistical significance (*P* < 0.1). Genes with fewer than two associated peaks are excluded.

The gene × TF matrix **G** is computed by aggregating in-silico ChIP scores across all cis-regulatory peaks associated with each gene:

**G**[*g*, *t*] = Σ<sub>*p* ∈ *P*<sub>*g*</sub></sub> *I*<sub>*p*,*t*</sub>

where *P*<sub>*g*</sub> denotes the set of significantly associated peaks for gene *g*, and *I*<sub>*p*,*t*</sub> is the in-silico ChIP score for peak *p* and TF *t*. The resulting matrix contains continuous values representing the cumulative regulatory evidence from all linked peaks.

### Interpretation

The resulting gene × TF matrix encodes putative regulatory relationships. A non-zero entry **G**[*g*, *t*] indicates that TF *t* has predicted regulatory potential on gene *g*, supported by: (i) at least one accessible chromatin peak within 100 kb of the gene, (ii) a predicted TF binding site within that peak, (iii) positive correlation between TF expression and peak accessibility, and (iv) significant correlation between peak accessibility and gene expression. Higher values indicate stronger evidence for regulation, reflecting either stronger TF-peak associations, more accessible chromatin, or multiple regulatory elements contributing to the relationship.

---

## Pipeline Parameters Summary

| Module | Parameter | Default Value | Description |
|--------|-----------|---------------|-------------|
| Gene-Peak Correlation | `transcript_span` | 100,000 bp | Distance threshold for peak-to-gene assignment |
| Gene-Peak Correlation | `min_corr` | 0 (exclusive) | Minimum Pearson correlation coefficient |
| Gene-Peak Correlation | `max_pval` | 0.1 (exclusive) | Maximum *P*-value for significance |
| Gene-Peak Correlation | `min_peaks` | 2 (inclusive) | Minimum peaks per gene for inclusion |
| Peak-TF Binding | `width` | 150 bp | Peak sequence width for motif scanning |
| Peak-TF Binding | FIMO threshold | *P* < 0.0001 | Default FIMO significance threshold |
| In-Silico ChIP | Correlation type | Spearman | Rank-based correlation for robustness |
| In-Silico ChIP | Correlation filter | > 0 | Retain only positive correlations |

---

## Software Implementation

The pipeline is implemented as a Snakemake workflow ensuring reproducibility and parallelization. Key dependencies include SEACells<sup>3</sup> for metacell aggregation and correlation computation, FIMO from the MEME Suite<sup>1</sup> for motif scanning, and the CIS-BP database<sup>2</sup> for TF position weight matrices. The in-silico ChIP score formulation and gene × TF matrix assembly extend the SEACells framework with novel algorithmic contributions for GRN inference.

---

## Data Flow Summary

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         INPUT DATA                                          │
├─────────────────────────────────────────────────────────────────────────────┤
│  SEACell RNA-seq     SEAcell ATAC-seq    GTF Annotation      TF Motif DB    │
│  (expression)        (accessibility)      (gene coords)     (CIS-BP PWMs)   │
└────────┬─────────────────────┬─────────────────┬────────────────┬──────────-┘
         │                     │                 │                │
         ▼                     ▼                 ▼                ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 1: Gene-Peak Correlation                                            │
│  ───────────────────────────────                                            │
│  • Assign peaks to genes (±100 kb)                                          │
│  • Compute Pearson correlations                                             │
│  • Filter by correlation and P-value                                        │
│  • Require minimum peaks per gene                                           │
│                                                                             │
│  Output: gp_corrs (gene, peak, cor, pval)                                   │
└────────────────────────────────┬────────────────────────────────────────────┘
                                 │
         ┌───────────────────────┼───────────────────────┐
         │                       │                       │
         ▼                       │                       ▼
┌────────────────────────┐       │       ┌────────────────────────────────────┐
│  MODULE 2: Peak-TF     │       │       │  MODULE 3: In-Silico ChIP          │
│  Binding Prediction    │       │       │  ────────────────────────          │
│  ─────────────────     │       │       │  • TF expression ↔ peak            │
│  • Extract sequences   │       │       │    accessibility correlation       │
│  • FIMO motif scan     │       │       │  • Filter positive correlations    │
│  • Build peak×TF       │       │       │  • Integrate: corr × acc × FIMO    │
│    matrix              │       │       │                                    │
│                        │       │       │  Output: InSilicoChip matrix       │
│  Output: FIMO matrix   │───────┼──────▶│                                    │
└────────────────────────┘       │       └─────────────────┬──────────────────┘
                                 │                         │
                                 ▼                         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 4: Gene × TF Matrix Assembly                                        │
│  ───────────────────────────────────                                        │
│  • Filter gene-peak associations (r > 0, P < 0.1)                           │
│  • Aggregate in-silico ChIP scores per gene                                 │
│  • Matrix multiplication: (gene × peak) × (peak × TF) = (gene × TF)        │
│                                                                             │
│  Output: geneXTF matrix (continuous values)                                 │
└────────────────────────────────┬────────────────────────────────────────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         OUTPUT                                              │
├─────────────────────────────────────────────────────────────────────────────┤
│  Gene × TF regulatory matrix                                                │
│  • Rows: Genes with ≥2 significantly associated peaks                       │
│  • Columns: Transcription factors with motif predictions                    │
│  • Values: Aggregated in-silico ChIP scores                                 │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Limitations

### Modeling of Transcriptional Repression

Following the in-silico ChIP methodology<sup>4</sup>, the pipeline retains only TF-peak pairs with positive correlations between TF expression and peak accessibility. This design choice reflects the assumption that TF binding at functional regulatory elements increases local chromatin accessibility—a relationship that holds for transcriptional activators but is inverted for repressors, where binding often leads to chromatin compaction and reduced accessibility.

Consequently, the gene × TF matrix captures putative *activating* regulatory relationships. Canonical repressors (e.g., REST, polycomb complex members) and dual-function TFs acting in repressive contexts are not modeled by this framework. For comprehensive GRN inference including repression, we recommend complementing this analysis with perturbation-based methods or repressive histone mark data (H3K27me3, H3K9me3).

---

## References

1. MEME Suite: Bailey, T.L. et al. The MEME Suite. *Nucleic Acids Res.* 43, W39–W49 (2015).
2. CIS-BP Database: Weirauch, M.T. et al. Determination and inference of eukaryotic transcription factor sequence specificity. *Cell* 158, 1431–1443 (2014).
3. SEACells: Persad, S. et al. SEACells infers transcriptional and epigenomic cellular states from single-cell genomics data. *Nat. Biotechnol.* 41, 1746–1757 (2023).
4. In-silico ChIP-seq: Argelaguet, R. et al. Decoding gene regulation in the mouse embryo using single-cell multi-omics. *bioRxiv* 2022.06.15.496239 (2022).
