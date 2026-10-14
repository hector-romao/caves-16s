## =====================================================================
## 06_heat_tree.R - heat trees (metacoder) por grupo
##
## A arvore taxonomica e a mesma para todos os grupos (layout fixo com
## set.seed); so muda a cor dos nos, que mostra a abundancia relativa
## media do grupo. Nos claros = pouco ou nenhum read naquele grupo.
##
## Pacote: install.packages("metacoder")
## =====================================================================

source("00_config.R")
library(metacoder)



## ---- Parametros que voce pode editar ------------------------------------

rank_heat <- "order"    # "family" e legivel; "genus" gera arvores enormes
log_trans <- FALSE       # TRUE aplica log1p nas abundancias antes de colorir

plot_heat_tree <- function(grupo, rank = rank_heat, log_trans = FALSE, salvar = TRUE) {

  dados   <- le_tabela(rank, "groups_rel")
  colunas <- setdiff(names(dados), "Lineage")
  if (!grupo %in% colunas) stop("Grupo nao encontrado: ", grupo)

  if (log_trans) dados[colunas] <- log1p(dados[colunas])

  taxa <- parse_tax_data(dados, class_cols = "Lineage", class_sep = ";")

  # Remove linhagens sem nenhuma abundancia
  sem_reads <- rowSums(taxa$data$tax_data[, colunas]) == 0
  taxa_filt <- filter_obs(taxa, data = "tax_data", !sem_reads, drop_taxa = TRUE)

  # Soma das abundancias dos descendentes em cada taxon
  taxa_filt$data$tax_abund <- calc_taxon_abund(taxa_filt, "tax_data", cols = colunas)

  # Gradiente do cinza claro ate a cor do grupo
  pal <- colorRampPalette(c("#ededed", cores_grupos[[grupo]]))

  set.seed(1)   # layout igual em todas as figuras
  p <- heat_tree(taxa_filt,
                 node_label = taxon_names,
                 node_color = taxa_filt$data$tax_abund[[grupo]],
                 node_color_axis_label = "Mean relative abundance (%)",
                 layout = "davidson-harel",       # era "davidson-harel"
                 initial_layout = "reingold-tilford",
                 node_color_range = pal(3),
                 title = grupo)

  if (salvar) {
    arq <- file.path(dir_fig, paste0("heat_tree_", marcador, "_", rank, "_", grupo, ".pdf"))
    pdf(arq, width = 8, height = 8)
    print(p)
    dev.off()
  }
  p
}

## ---- Um heat tree por grupo -------------------------------------------------

for (g in ordem_grupos) {
  plot_heat_tree(g, rank = rank_heat, log_trans = log_trans)
}
