## =====================================================================
## 01_prepare_tables.R
##
## Le o objeto phyloseq limpo e gera, para cada nivel taxonomico
## (phylum, class, order, family, genus), tres tabelas com a coluna
## 'Lineage' (taxonomia separada por ";") e uma coluna por amostra/grupo:
##
##   *_samples_rel.tsv : abundancia relativa (%) por amostra     (dados brutos)
##   *_groups_rel.tsv  : media da abundancia relativa (%) por grupo (dados brutos)
##   *_groups_rare.tsv : contagens rarefeitas somadas por grupo   (presenca/ausencia)
##
## Regras:
##  - Abundancia usa os dados BRUTOS (proporcoes), sem descartar reads.
##  - Presenca/ausencia (Venn, UpSet, riqueza) usa dados RAREFEITOS na mesma
##    profundidade, para que a profundidade de cada amostra nao pese.
##  - ASV sem classificacao em um nivel recebe o nome "unclassified_<ultimo
##    nome conhecido>" (ex.: unclassified_Gaiellales). Assim a linhagem fica
##    unica e nenhum read se perde.
## =====================================================================

source("00_config.R")
suppressMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
})

## ---- 1. Objeto phyloseq e metadata -----------------------------------

ps <- readRDS(ps_limpo)

meta <- data.frame(sample_data(ps), check.names = FALSE) %>%
  rownames_to_column("sample")

if (!"cave" %in% names(meta)) stop("A metadata precisa da coluna 'cave'.")

meta$grupo   <- sub("_16S", "", meta$cave)             # ex.: jane_16S_chuva -> jane_chuva
meta$caverna <- sub("_.*$", "", meta$grupo)            # jane
meta$estacao <- ifelse(grepl("chuva", meta$grupo), "chuva",
                ifelse(grepl("seca",  meta$grupo), "seca", NA))

extra <- setdiff(unique(meta$grupo), ordem_grupos)
if (length(extra) > 0) {
  stop("Grupos que nao estao em 'ordem_grupos' (00_config.R): ",
       paste(extra, collapse = ", "))
}

grupos <- intersect(ordem_grupos, unique(meta$grupo))
meta   <- meta[order(match(meta$grupo, grupos), meta$sample), ]
print(table(meta$grupo)[grupos])

write.table(meta, file.path(dir_tab, "metadata_grupos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 2. Filtro de seguranca -------------------------------------------

ps <- subset_taxa(ps, Kingdom %in% c("Bacteria", "Archaea") &
                      (is.na(Order)  | Order  != "Chloroplast") &
                      (is.na(Family) | Family != "Mitochondria"))

stopifnot(setequal(sample_names(ps), meta$sample))

## ---- 3. Rarefacao (somente para presenca/ausencia) --------------------

prof <- sort(sample_sums(ps))
print(prof)
if (min(prof) < prof_rare) {
  stop("Ha amostras com menos de ", prof_rare, " reads; ajuste 'prof_rare' em 00_config.R")
}

ps_rare <- rarefy_even_depth(ps, sample.size = prof_rare, rngseed = semente,
                             replace = FALSE, trimOTUs = TRUE, verbose = FALSE)

## ---- 4. Funcoes auxiliares ---------------------------------------------

niveis      <- c(phylum = 2, class = 3, order = 4, family = 5, genus = 6)
colunas_tax <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Linhagem de cada ASV ate o nivel pedido, preenchendo NA com "unclassified_"
monta_lineage <- function(ps, rank) {
  tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  m   <- as.matrix(tax[, colunas_tax[1:niveis[[rank]]], drop = FALSE])
  for (j in seq_len(ncol(m))) {
    falta <- is.na(m[, j])
    if (any(falta)) {
      anterior <- if (j == 1) rep("Root", nrow(m)) else m[, j - 1]
      anterior <- anterior[falta]
      m[falta, j] <- ifelse(startsWith(anterior, "unclassified_"),
                            anterior, paste0("unclassified_", anterior))
    }
  }
  setNames(apply(m, 1, paste, collapse = ";"), rownames(tax))
}

# Soma as ASVs que compartilham a mesma linhagem: linhas = linhagens, colunas = amostras
agrega <- function(ps, lin) {
  mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) mat <- t(mat)
  rowsum(mat, group = lin[rownames(mat)])
}

rel <- function(m) sweep(m, 2, colSums(m), "/") * 100

media_grupo <- function(m, meta) {
  sapply(grupos, function(g) rowMeans(m[, meta$sample[meta$grupo == g], drop = FALSE]))
}

soma_grupo <- function(m, meta) {
  sapply(grupos, function(g) rowSums(m[, meta$sample[meta$grupo == g], drop = FALSE]))
}

salva <- function(m, rank, tipo) {
  m   <- m[rowSums(m) > 0, , drop = FALSE]
  out <- data.frame(Lineage = rownames(m), m, check.names = FALSE, row.names = NULL)
  arq <- file.path(dir_tab, paste0(marcador, "_", rank, "_", tipo, ".tsv"))
  write.table(out, arq, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  salvo:", basename(arq), "-", nrow(out), "linhagens\n")
}

## ---- 5. Gera as tabelas para cada nivel ---------------------------------

for (r in names(niveis)) {
  cat("\n== Nivel:", r, "==\n")

  # Dados brutos -> abundancia relativa
  m_raw <- agrega(ps, monta_lineage(ps, r))
  m_rel <- rel(m_raw)
  nao_class <- sum(m_raw[grepl("unclassified_", rownames(m_raw)), ]) / sum(m_raw)
  cat("  fracao dos reads sem classificacao neste nivel:", round(nao_class, 3), "\n")

  salva(m_rel[, meta$sample, drop = FALSE], r, "samples_rel")
  salva(media_grupo(m_rel, meta),           r, "groups_rel")

  # Dados rarefeitos -> presenca/ausencia por grupo, controlando diferenças
  # na profundidade de sequenciamento entre as amostras.
  m_rare <- agrega(ps_rare, monta_lineage(ps_rare, r))
  salva(soma_grupo(m_rare, meta), r, "groups_rare")
}

cat("\nPronto. Tabelas em:", dir_tab, "\n")
