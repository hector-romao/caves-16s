## =====================================================================
## 05_alluvial.R - graficos aluviais Phylum -> Class -> Genus
##
## alluvial_plot()    : fluxo ponderado pela abundancia relativa media do grupo
## alluvial_riqueza() : fluxo onde cada genero presente vale 1 (riqueza),
##                      usando as contagens rarefeitas e ignorando "unclassified_*"
##
## Pacote: install.packages("ggalluvial")
## =====================================================================

source("00_config.R")
suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggalluvial)
  library(ggplot2)
  library(patchwork)
})

## ---- Parametros que voce pode editar ------------------------------------

# Com ~50 mil ASVs ha centenas de generos; manter so os N mais abundantes
# deixa o grafico legivel. NULL = todos.
top_n_generos <- 30

# Mesma cor para o mesmo filo em todos os graficos e paineis
todos_filos <- sort(unique(sapply(strsplit(le_tabela("genus", "groups_rel")$Lineage, ";"),
                                  `[`, 2)))
cores_filo  <- setNames(scales::hue_pal()(length(todos_filos)), todos_filos)

## ---- Preparo dos dados de UM grupo ---------------------------------------

prepara <- function(tipo, grupo, top_n = NULL, excluir_nao_class = FALSE) {
  dados <- le_tabela("genus", tipo)
  if (!grupo %in% names(dados)) stop("Grupo nao encontrado: ", grupo)

  dados <- dados[dados[[grupo]] > 0, c("Lineage", grupo)]
  names(dados)[2] <- "Abundance"
  if (excluir_nao_class) dados <- dados[!grepl("unclassified_", dados$Lineage), ]

  dados <- dados %>%
    separate(Lineage, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
             sep = ";") %>%
    group_by(Phylum, Class, Genus) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(Abundance))

  if (!is.null(top_n)) dados <- slice_head(dados, n = top_n)

  dados$Phylum <- factor(dados$Phylum, levels = unique(dados$Phylum))
  dados$Class  <- factor(dados$Class,  levels = unique(dados$Class))
  dados$Genus  <- factor(dados$Genus,  levels = unique(dados$Genus))
  dados
}

## ---- Desenho do grafico --------------------------------------------------

desenha <- function(dados_agg, titulo, ylab, log_trans = FALSE, tam_texto = 3.5) {
  p <- ggplot(dados_agg,
              aes(axis1 = Phylum, axis2 = Class, axis3 = Genus, y = Abundance)) +
    geom_alluvium(aes(fill = Phylum), width = 1/1200) +
    geom_stratum(width = 1/100, fill = "#1c0e11", color = "white", linewidth = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              check_overlap = TRUE, hjust = 0, size = tam_texto, family = "serif") +
    scale_x_discrete(limits = c("Phylum", "Class", "Genera"), expand = c(.05, .05)) +
    scale_fill_manual(values = cores_filo) +
    labs(title = titulo, y = ylab) +
    theme_minimal()

  if (log_trans) p <- p + scale_y_log10()
  p
}

## ---- Versao por abundancia ------------------------------------------------

alluvial_plot <- function(grupo, log_trans = FALSE, top_n = top_n_generos,
                          salvar = TRUE, tam_texto = 3.5) {

  dados_agg <- prepara("groups_rel", grupo, top_n)

  titulo <- paste0(grupo, if (!is.null(top_n)) paste0(" (top ", top_n, " genera)") else "")
  p <- desenha(dados_agg, titulo, "Mean relative abundance (%)", log_trans, tam_texto)

  if (salvar) {
    arq <- file.path(dir_fig, paste0("alluvial_", marcador, "_", grupo,
                                     if (log_trans) "_log10" else "", ".pdf"))
    ggsave(arq, p, width = 16, height = 14, device = cairo_pdf)
  }
  p
}

## ---- Versao por riqueza (cada genero conta 1) ------------------------------

alluvial_riqueza <- function(grupo, salvar = TRUE, tam_texto = 2.5) {

  dados_agg <- prepara("groups_rare", grupo, top_n = NULL, excluir_nao_class = TRUE) %>%
    mutate(Abundance = 1)

  p <- desenha(dados_agg, paste0(grupo, " - genus richness"), "Genera (count)",
               FALSE, tam_texto)

  if (salvar) {
    arq <- file.path(dir_fig, paste0("alluvial_riqueza_", marcador, "_", grupo, ".pdf"))
    ggsave(arq, p, width = 14, height = 18, device = cairo_pdf)
  }
  p
}

## ---- Gera tudo ---------------------------------------------------------------

# Um grafico por grupo (abundancia e riqueza)
for (g in ordem_grupos) {
  alluvial_plot(g)
  alluvial_riqueza(g)
}

# Painel com todos os grupos (abundancia)
plots <- lapply(ordem_grupos, alluvial_plot, salvar = FALSE, tam_texto = 2)
painel <- wrap_plots(plots, ncol = 5)
ggsave(file.path(dir_fig, paste0("alluvial_", marcador, "_painel.pdf")),
       painel, width = 40, height = 20, device = cairo_pdf)
