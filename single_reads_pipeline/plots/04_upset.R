## =====================================================================
## 04_upset.R - UpSet plot: linhagens compartilhadas entre os grupos
##
## Funciona como o Venn, mas aguenta muitos conjuntos (aqui, os 10 grupos).
## Presenca/ausencia vem das contagens rarefeitas somadas por grupo.
##
## Pacote: install.packages("UpSetR")
## =====================================================================

source("00_config.R")
library(UpSetR)

## ---- Parametros que voce pode editar ------------------------------------

ranks_upset <- c("family", "genus")
remover_nao_classificados <- TRUE     # ignora linhagens "unclassified_*"

plot_upsets <- function(rank, sets = ordem_grupos, nintersects = 30, sufixo = "") {

  dados <- le_tabela(rank, "groups_rare")
  if (remover_nao_classificados) dados <- dados[!grepl("unclassified_", dados$Lineage), ]

  rownames(dados) <- dados$Lineage
  dados <- dados[, sets, drop = FALSE]

  # Converter para presenca (1) ou ausencia (0)
  dados_binarios <- dados
  dados_binarios[] <- lapply(dados, function(x) as.integer(x > 0))
  dados_binarios <- dados_binarios[rowSums(dados_binarios) > 0, , drop = FALSE]

  # Com 10 conjuntos, mostrar intersecoes vazias so poluiria o grafico,
  # por isso empty.intersections fica desligado (no original estava "on").
  p <- UpSetR::upset(
    dados_binarios,
    nsets        = length(sets),
    nintersects  = nintersects,
    sets         = sets,
    order.by     = "freq",
    keep.order   = TRUE,
    show.numbers = "yes",
    point.size   = 3,
    text.scale   = c(2, 2, 2, 2, 1.5, 2),
    sets.bar.color = unname(cores_grupos[sets])
  )

  arq <- file.path(dir_fig, paste0("upset_", marcador, "_", rank, sufixo, ".pdf"))
  cairo_pdf(arq, width = 12, height = 8)
  print(p)
  dev.off()
  cat(">>> Salvo em:", arq, "\n")
}

## ---- Todos os grupos ----------------------------------------------------------

for (r in ranks_upset) plot_upsets(r)

## ---- Exemplo: so os grupos com estacao conhecida -------------------------------
## (descomente se quiser)
# com_estacao <- grep("chuva|seca", ordem_grupos, value = TRUE)
# for (r in ranks_upset) plot_upsets(r, sets = com_estacao, sufixo = "_com_estacao")
