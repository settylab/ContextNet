## Methods

### ContextNet: Integration of Context-Specific Gene Regulatory Networks with Cell-Cell Communication Inference

#### Overview

ContextNet is a computational framework that extends the NicheNet^1^ cell-cell communication inference method by integrating context-specific gene regulatory networks (GRNs) derived from single-cell multiomics data. Standard cell-cell communication inference methods rely on curated prior knowledge networks that aggregate regulatory relationships across diverse biological contexts. These generic networks often suffer from bias in the input data, and may fail to capture tissue-specific, disease-specific, or developmental stage-specific regulatory relationships. ContextNet addresses this limitation by enabling seamless incorporation of custom GRNs constructed from matched scRNA-seq and scATAC-seq data, thereby grounding ligand-target predictions in experimentally observed regulatory relationships.

#### GRN Import and Preprocessing

ContextNet accepts GRNs as a gene × transcription factor (TF) weight matrix where rows represent target genes, columns represent TFs, and values represent regulatory weights quantifying the strength of each TF-target relationship. This matrix representation captures the output of multiomics GRN inference methods that jointly analyze chromatin accessibility and gene expression to identify high-confidence regulatory edges.

Upon import, regulatory weights undergo per-TF normalization to ensure comparability across TFs with different weight distributions. For each TF column, weights are divided by the maximum value in that column, yielding normalized weights in the range [0, 1] where the highest-confidence target for each TF receives a weight of 1.0. This normalization preserves the relative ranking of targets within each TF's regulon while enabling meaningful comparisons across TFs. Edges with zero weights are removed, and an optional minimum weight threshold can be applied to retain only high-confidence regulatory relationships.

#### Network Integration with NicheNet Prior Knowledge

ContextNet integrates the custom GRN with NicheNet's curated prior knowledge networks, which comprise: (1) a ligand-receptor network defining known ligand-receptor interactions, and (2) weighted signaling and regulatory networks capturing intracellular signal transduction from receptors to downstream genes. NicheNet's prior networks are organized into a ligand-receptor signaling component (lr_sig) encoding receptor-to-downstream signaling, and a gene regulatory component (gr) encoding TF-to-target relationships.

During integration, the custom GRN completely replaces NicheNet's gr component, yielding a ligand-target matrix derived entirely from context-specific regulatory relationships. NicheNet's ligand-receptor and signaling networks (lr_sig) are preserved unchanged, ensuring that known ligand-receptor interactions and proximal signaling events remain informed by curated knowledge while downstream regulatory cascades reflect the specific biological context.

#### Ligand-Target Matrix Construction

Following network integration, ContextNet constructs a ligand-target regulatory potential matrix using the Personalized PageRank (PPR) algorithm^2,3^, consistent with NicheNet's established methodology. The resulting ligand-target matrix $M$ has dimensions $|targets| \times |ligands|$, where entry $M_{i,j}$ represents the regulatory potential of ligand $j$ affecting target gene $i$. Higher values indicate stronger predicted regulatory influence. A ligand-target filtering (LTF) cutoff is applied to remove low-confidence predictions below a threshold (default derived from NicheNet hyperparameter optimization).

#### Output and Downstream Analysis

ContextNet produces a ContextNetwork object containing the ligand-target regulatory potential matrix along with the integrated signaling and regulatory networks. This output is directly compatible with downstream NicheNet and MultiNicheNet analyses, including ligand activity prediction, ligand-target link inference, and sender-receiver prioritization. The `as_nichenet_input()` function exports the integrated networks in standard NicheNet format, enabling seamless integration with existing analysis workflows.

#### Implementation

ContextNet is implemented as an R package compatible with R version 4.1.0 and above. The package leverages nichenetr for PPR-based ligand-target matrix construction and depends on standard data manipulation packages (dplyr, reshape2, Matrix). All S4 class definitions include validity checking to ensure data integrity. The package is available at [repository URL].

#### References

1. Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. *Nat. Methods* **17**, 159–162 (2020).

2. Page, L., Brin, S., Motwani, R. & Winograd, T. The PageRank citation ranking: Bringing order to the web. Technical Report, Stanford InfoLab (1999).

3. Haveliwala, T. H. Topic-sensitive PageRank: A context-sensitive ranking algorithm for web search. *IEEE Trans. Knowl. Data Eng.* **15**, 784–796 (2003).
