## =====================================================================
## 02_barplot.R - barras de abundancia relativa (top N taxons + Other)
##
## por = "amostra": uma barra por amostra (21), agrupadas por caverna/estacao
## por = "grupo"  : uma barra por grupo (media das amostras do grupo)
## =====================================================================

source("00_config.R")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})

meta <- read.delim(file.path(dir_tab, "metadata_grupos.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)

bar_plot <- function(rank, por = "amostra", top_n = 20, salvar = TRUE) {

  tipo  <- if (por == "amostra") "samples_rel" else "groups_rel"
  dados <- le_tabela(rank, tipo)

  # Nome do taxon = ultimo elemento da linhagem (como no script original)
  dados <- dados %>%
    mutate(Taxon = str_extract(Lineage, "[^;]+$")) %>%
    select(-Lineage)

  dados_long <- dados %>%
    pivot_longer(-Taxon, names_to = "Sample", values_to = "Abundance") %>%
    group_by(Taxon, Sample) %>%                       # junta nomes repetidos
    summarise(Abundance = sum(Abundance), .groups = "drop")

  top <- dados_long %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Abundance), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    slice_head(n = top_n) %>%
    pull(Taxon)

  dados_plot <- dados_long %>%
    mutate(Taxon = ifelse(Taxon %in% top, Taxon, "Other")) %>%
    group_by(Sample, Taxon) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(RelAbundance = Abundance / sum(Abundance)) %>%
    ungroup()

  dados_plot$Taxon <- factor(dados_plot$Taxon, levels = c(top, "Other"))

  # Ordem e rotulos das barras (primeiro grupo no topo)
  if (por == "amostra") {
    ord <- meta$sample
    rot <- paste0(meta$grupo, " | ", sub("_L001$", "", sub("^16S-", "", meta$sample)))
  } else {
    ord <- ordem_grupos[ordem_grupos %in% unique(dados_plot$Sample)]
    rot <- ord
  }
  dados_plot$Sample <- factor(dados_plot$Sample, levels = rev(ord), labels = rev(rot))

  cores <- c(setNames(viridis(length(top), option = "turbo"), top), Other = "grey80")

  p <- ggplot(dados_plot, aes(x = Sample, y = RelAbundance, fill = Taxon)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = cores) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme_minimal(base_family = "serif") +
    labs(title = paste0(toupper(marcador), " - ", rank),
         y = "Relative abundance", x = NULL, fill = "Taxa")

  if (salvar) {
    arq <- file.path(dir_fig, paste0("barplot_", marcador, "_", rank, "_", por, ".pdf"))
    ggsave(arq, p, width = 9, height = if (por == "amostra") 8 else 5, device = cairo_pdf)
  }
  p
}

## ---- Um grafico por nivel e por tipo (amostra / grupo) -----------------

for (por in c("amostra", "grupo")) {
  plots <- lapply(ranks, bar_plot, por = por)
  names(plots) <- ranks

  painel <- wrap_plots(plots, ncol = 2)
  ggsave(file.path(dir_fig, paste0("barplot_", marcador, "_painel_", por, ".pdf")),
         painel, width = 16, height = 22, device = cairo_pdf)
}




###############################################################################
#####################Realizar o plot por amostra##############################
###############################################################################

barplot_amostra <- function(amostra, rank = "family", top_n = 15, salvar = TRUE) {

  dados <- le_tabela(rank, "samples_rel")
  if (!amostra %in% names(dados)) stop("Amostra nao encontrada: ", amostra)

  grupo <- meta$grupo[meta$sample == amostra]

  # Abundancia relativa (%) de cada taxon nesta amostra, ordenada
  d <- dados %>%
    transmute(Taxon = str_extract(Lineage, "[^;]+$"), Abundance = .data[[amostra]]) %>%
    group_by(Taxon) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(Abundance))

  # Top N desta amostra + o restante agrupado em "Other"
  top <- slice_head(d, n = top_n)

  outros <- tibble(
    Taxon = "Other",
    Abundance = sum(d$Abundance[-seq_len(nrow(top))], na.rm = TRUE)
  )

  plot_df <- bind_rows(top, outros)

  # Definir os níveis dos táxons
  plot_df$Taxon <- factor(
    plot_df$Taxon,
    levels = plot_df$Taxon
  )

  # Criar uma cor para cada táxon
  cores <- c(
    setNames(
      viridis::viridis(nrow(top), option = "turbo"),
      as.character(top$Taxon)
    ),
    Other = "grey80"
  )

  # Gráfico de barras empilhadas
  p <- ggplot(plot_df, aes(x = amostra, y = Abundance, fill = Taxon)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = cores) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_minimal(base_family = "serif") +
    labs(
      title = paste0(amostra, " (", grupo, ")"),
      subtitle = paste0(toupper(marcador), " - ", rank),
      x = NULL,
      y = "Relative abundance (%)",
      fill = "Taxon"
    )

  print(p)


  if (salvar) {
    arq <- file.path(dir_fig, "barplot", "por_amostra",
                     paste0("barplot_", marcador, "_", rank, "_", amostra, ".pdf"))
    dir.create(dirname(arq), recursive = TRUE, showWarnings = FALSE)
    ggsave(arq, p, width = 8, height = 6, device = cairo_pdf)
  }
  p
}

# Uma amostra:
barplot_amostra("16S-J1_S34_L001", rank = "family")

# Todas as 21 amostras, nos 5 niveis taxonomicos (105 PDFs):
for (a in meta$sample) {
  for (r in ranks) barplot_amostra(a, rank = r)
}
