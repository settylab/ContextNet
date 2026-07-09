# Generate the small public smoke-test fixture used by Tutorial 01.
#
# The fixture is intentionally synthetic. It is not a biological reference
# dataset; it exists so new users can verify that ContextNet, nichenetr, Seurat,
# and the NicheNet-style prior files are installed and wired together.

suppressPackageStartupMessages({
    library(Seurat)
})

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
                        error = function(e) NA_character_)
if (!is.na(script_path)) {
    repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                               mustWork = FALSE)
} else {
    repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
    repo_root <- normalizePath(getwd(), mustWork = TRUE)
}

out_dir <- file.path(repo_root, "inst", "extdata", "contextnet_smoke")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)

genes <- c(
    "TGFB1", "FGF2", "TGFBR1", "FGFR1", "SMAD3", "FOS",
    "SERPINE1", "CDKN1A", "MMP1", "JUNB", "GAPDH", "ACTB"
)

# Context-specific gene regulatory network (GRN): target genes x TFs.
gene_x_tf <- matrix(
    0,
    nrow = 4,
    ncol = 2,
    dimnames = list(c("SERPINE1", "CDKN1A", "MMP1", "JUNB"), c("SMAD3", "FOS"))
)
gene_x_tf["SERPINE1", "SMAD3"] <- 1.0
gene_x_tf["CDKN1A", "SMAD3"] <- 0.9
gene_x_tf["MMP1", "FOS"] <- 1.0
gene_x_tf["JUNB", "FOS"] <- 0.8
utils::write.csv(gene_x_tf, file.path(out_dir, "gene_x_tf.csv"))

# Minimal NicheNet-compatible priors.
lr_network <- data.frame(
    from = c("TGFB1", "FGF2"),
    to = c("TGFBR1", "FGFR1"),
    database = "contextnet-smoke",
    stringsAsFactors = FALSE
)
saveRDS(lr_network, file.path(out_dir, "lr_network_human_21122021.rds"))

lr_sig <- data.frame(
    from = c("TGFB1", "TGFBR1", "TGFBR1", "FGF2", "FGFR1", "FGFR1"),
    to = c("TGFBR1", "SMAD3", "FOS", "FGFR1", "FOS", "SMAD3"),
    weight = c(1.0, 1.0, 0.2, 1.0, 1.0, 0.1),
    stringsAsFactors = FALSE
)
gr <- data.frame(
    from = c("SMAD3", "SMAD3", "FOS", "FOS"),
    to = c("SERPINE1", "CDKN1A", "MMP1", "JUNB"),
    weight = c(0.1, 0.1, 0.1, 0.1),
    stringsAsFactors = FALSE
)
saveRDS(list(lr_sig = lr_sig, gr = gr),
        file.path(out_dir, "weighted_networks_nsga2r_final.rds"))

# Tiny Seurat object with two sender populations, one receiver population, and
# one reference population. LIG1 is designed to recover the receiver target
# program more strongly than FGF2.
n_cells <- 80
cell_type <- rep(c("sender_a", "sender_b", "receiver", "reference"), each = 20)
counts <- matrix(
    rpois(length(genes) * n_cells, lambda = 1),
    nrow = length(genes),
    dimnames = list(genes, paste0("cell", seq_len(n_cells)))
)
counts["TGFB1", cell_type == "sender_a"] <- rpois(20, lambda = 8)
counts["FGF2", cell_type == "sender_b"] <- rpois(20, lambda = 8)
counts["TGFBR1", cell_type == "receiver"] <- rpois(20, lambda = 6)
counts["FGFR1", cell_type == "reference"] <- rpois(20, lambda = 6)
counts[c("SERPINE1", "CDKN1A"), cell_type == "receiver"] <-
    counts[c("SERPINE1", "CDKN1A"), cell_type == "receiver"] + 6

meta <- data.frame(
    cell_type = cell_type,
    condition = ifelse(cell_type == "receiver", "program", "baseline"),
    stringsAsFactors = FALSE,
    row.names = colnames(counts)
)

seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    project = "ContextNetSmoke"
)
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
saveRDS(seurat_obj, file.path(out_dir, "seurat_obj.rds"))

message("Wrote ContextNet smoke-test fixture to: ", out_dir)
