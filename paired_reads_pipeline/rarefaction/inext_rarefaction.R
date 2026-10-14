# Rarefação por diversidade e consolidação por família/gênero
#
# Este script combina a análise de diversidade com a agregação taxonômica em
# níveis superiores (famílias e gêneros). Ele gera gráficos de rarefação e salva
# objetos com a matriz resumida para uso posterior.

library(phyloseq)
library(iNEXT)
library(ggplot2)
library(patchwork)

ps <- readRDS("ps-filT270.rds")
prof_rare <- 30000
path_pipeline <- "aggregated"
dir.create("plots", showWarnings = FALSE)
dir.create(path_pipeline, showWarnings = FALSE)

# Função para gerar painéis de rarefação por diversidade de ordem q.
rarefaction_plot <- function(ps, output, amostras = NULL, type = 1, nboot = 50) {

  # 1. Matriz de abundância (amostras nas linhas)
  mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) mat <- t(mat)
  if (!is.null(amostras)) mat <- mat[amostras, , drop = FALSE]

  # 2. Conferência básica
  print(rowSums(mat))

  # 3. Converter para lista (removendo zeros)
  lista <- lapply(rownames(mat), function(s) {
    x <- mat[s, ]
    x[x > 0]
  })
  names(lista) <- rownames(mat)

  # 4. Rodar iNEXT para diferentes ordens de diversidade.
  richness <- iNEXT(lista, q = 0, datatype = "abundance", nboot = nboot)
  shannon  <- iNEXT(lista, q = 1, datatype = "abundance", nboot = nboot)
  simpson  <- iNEXT(lista, q = 2, datatype = "abundance", nboot = nboot)

  # 5. Gerar plots
  n <- length(lista)

  faz_plot <- function(x, titulo) {
    ggiNEXT(x, type = type) +
      geom_vline(
        xintercept = prof_rare,
        linetype = "dashed",
        colour = "black",
        linewidth = 0.8
      ) +
      scale_shape_manual(values = rep(c(16, 17, 15, 18), length.out = n)) +
      ggtitle(titulo) +
      theme_bw() +
      guides(colour = guide_legend(nrow = 4), fill = "none", shape = "none")
  }

  richness_plot <- faz_plot(richness, "Richness (q = 0)")
  shannon_plot  <- faz_plot(shannon,  "Shannon (q = 1)")
  simpson_plot  <- faz_plot(simpson,  "Simpson (q = 2)")

  # 6. Cria painel com uma única legenda.
  painel <- (richness_plot | shannon_plot | simpson_plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  # 7. Salva painel em PDF.
  ggsave(file.path("plots", paste0(output, "_painel.pdf")),
         painel, width = 14, height = 5)

  # 8. Retorna os objetos gerados.
  return(list(painel = painel, richness = richness_plot,
              shannon = shannon_plot, simpson = simpson_plot))
}

# Todas as amostras.
res_todas <- rarefaction_plot(ps, "16S_todas_amostras")

# Uma figura por caverna.
meta <- data.frame(sample_data(ps))
for (cv in unique(meta$cave)) {
  rarefaction_plot(ps, paste0("16S_", cv),
                   amostras = rownames(meta)[meta$cave == cv])
}

# ----------------------------
# Agregação por família
# ----------------------------
rank_cols <- c("Phylum", "Class", "Order", "Family")
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

# Rótulo = linhagem completa; inclui os níveis acima para evitar nomes repetidos.
lin <- apply(tax[, rank_cols], 1, function(x) {
  ok <- which(!is.na(x))
  if (length(ok) == 0) return("Unclassified")
  last <- max(ok)
  if (last == 4) paste(x, collapse = "|")
  else paste0(paste(x[1:last], collapse = "|"), "|unclassified")
})

mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) mat <- t(mat)
lin <- lin[colnames(mat)]

fam_mat <- rowsum(t(mat), group = lin)

ps_fam <- phyloseq(otu_table(fam_mat, taxa_are_rows = TRUE),
                   sample_data(ps))

stopifnot(all(sample_sums(ps_fam) == sample_sums(ps)))

saveRDS(ps_fam, file.path(path_pipeline, "ps_familia.rds"))
saveRDS(fam_mat, file.path(path_pipeline, "fam_mat.rds"))

# ----------------------------
# Agregação por gênero
# ----------------------------
rank_cols_g <- c("Phylum", "Class", "Order", "Family", "Genus")
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

lin <- apply(tax[, rank_cols_g], 1, function(x) {
  ok <- which(!is.na(x))
  if (length(ok) == 0) return("Unclassified")
  last <- max(ok)
  if (last == 5) paste(x, collapse = "|")
  else paste0(paste(x[1:last], collapse = "|"), "|unclassified")
})

mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) mat <- t(mat)
lin <- lin[colnames(mat)]

gen_mat <- rowsum(t(mat), group = lin)

ps_gen <- phyloseq(otu_table(gen_mat, taxa_are_rows = TRUE),
                   sample_data(ps))

stopifnot(all(sample_sums(ps_gen) == sample_sums(ps)))

saveRDS(ps_gen, file.path(path_pipeline, "ps_genus.rds"))
saveRDS(gen_mat, file.path(path_pipeline, "gen_mat.rds"))

res_todas <- rarefaction_plot(ps_gen, "16S_todas_amostras_genus")
