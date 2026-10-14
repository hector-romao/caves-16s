## =====================================================================
## 03_venn.R - diagramas de Venn (2 a 4 conjuntos) por nivel taxonomico
##
## Cada "conjunto" e um grupo (caverna x estacao); os elementos sao as
## linhagens com contagem rarefeita > 0 naquele grupo. Os numeros mostram
## quantas linhagens sao exclusivas de cada grupo e quantas sao compartilhadas.
##
## Pacote: install.packages("VennDiagram")
## =====================================================================

source("00_config.R")
suppressMessages(library(VennDiagram))
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

## ---- Parametros que voce pode editar ------------------------------------

ranks_venn <- c("family", "genus")

# TRUE = ignora linhagens "unclassified_*" (nao sao taxons de verdade)
remover_nao_classificados <- TRUE

# Cada item e um grupo de 2 a 4 colunas a comparar. O nome do item vira o
# nome do arquivo .png. dores e samambaia nao entram nas comparacoes de
# estacao porque a pasta nao informa a estacao; elas aparecem no UpSet.
grupos_venn <- list(
  Jane_chuva_vs_seca   = c("jane_chuva",   "jane_seca"),
  Tarzan_chuva_vs_seca = c("tarzan_chuva", "tarzan_seca"),
  Timor_chuva_vs_seca  = c("timor_chuva",  "timor_seca"),
  Val_chuva_vs_seca    = c("val_chuva",    "val_seca"),
  Cavernas_chuva       = c("jane_chuva", "tarzan_chuva", "timor_chuva", "val_chuva"),
  Cavernas_seca        = c("jane_seca",  "tarzan_seca",  "timor_seca",  "val_seca")
)

## ---- Funcao que gera UM diagrama -----------------------------------------

gerar_venn <- function(tabela, rank, nome_grupo, colunas) {

  if (!length(colunas) %in% 2:4) stop("Cada grupo precisa ter de 2 a 4 colunas.")
  faltando <- setdiff(colunas, names(tabela))
  if (length(faltando) > 0) stop("Colunas nao encontradas: ", paste(faltando, collapse = ", "))

  conjuntos <- lapply(colunas, function(col) tabela$Lineage[tabela[[col]] > 0])
  names(conjuntos) <- colunas

  arquivo <- file.path(dir_fig, paste0("venn_", marcador, "_", rank, "_", nome_grupo, ".png"))

  venn.diagram(
    x               = conjuntos,
    filename        = arquivo,
    imagetype       = "png",
    resolution      = 300,
    fill            = unname(cores_grupos[colunas]),
    alpha           = 0.5,
    cex             = 1.3,
    cat.cex         = 1.2,
    cat.fontface    = "bold",
    margin          = 0.08,
    main            = paste0(nome_grupo, " (", rank, ")"),
    main.cex        = 1.4,
    disable.logging = TRUE
  )
  cat(">>> Diagrama salvo em:", arquivo, "\n")
}

## ---- Gera todos os diagramas ------------------------------------------------

for (rank in ranks_venn) {
  tab <- le_tabela(rank, "groups_rare")
  if (remover_nao_classificados) tab <- tab[!grepl("unclassified_", tab$Lineage), ]
  for (nome in names(grupos_venn)) {
    gerar_venn(tab, rank, nome, grupos_venn[[nome]])
  }
}
