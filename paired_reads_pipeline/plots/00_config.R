## =====================================================================
## 00_config.R - configuracao comum a todos os scripts de graficos
##
## Todos os scripts fazem source("00_config.R"), entao rode-os com o
## diretorio de trabalho na pasta onde estes arquivos estao
## (ex.: setwd("paired_reads_pipeline/plots") a partir da raiz do projeto).
##
## Ordem de execucao:
##   01_prepare_tables.R  (uma vez, gera as tabelas)
##   02_barplot.R, 03_venn.R, 04_upset.R, 05_alluvial.R, 06_heat_tree.R
## =====================================================================

# O objeto e os diretórios de saída podem ser informados por variáveis de
# ambiente para permitir execução em outra máquina sem editar os scripts.
ps_limpo <- Sys.getenv(
  "CAVES_16S_PHYLOSEQ",
  unset = "/media/lgbio-nas1/hectorromao/fungi_caves/3.dada2/rds_files/ps-filT270.rds"
)

dir_tab <- file.path(getwd(), "tables_paired")
dir_fig <- file.path(getwd(), "figures_paired")

marcador  <- "16s"     # prefixo dos arquivos
prof_rare <- 30000     # profundidade da rarefacao (so para presenca/ausencia)
semente   <- 123

ranks <- c("phylum", "class", "order", "family", "genus")

# Grupos = valores da coluna 'cave' da metadata, sem o sufixo "_16S".
# dores e samambaia ficam sem estacao enquanto voce nao confirmar.
ordem_grupos <- c("dores",
                  "jane_chuva", "jane_seca",
                  "samambaia",
                  "tarzan_chuva", "tarzan_seca",
                  "timor_chuva", "timor_seca",
                  "val_chuva", "val_seca")

# Uma cor por caverna; chuva = tom escuro, seca = tom claro
cores_grupos <- c(dores        = "#8c564b",
                  jane_chuva   = "#1f77b4", jane_seca   = "#9ecae1",
                  samambaia    = "#2ca02c",
                  tarzan_chuva = "#d62728", tarzan_seca = "#ff9896",
                  timor_chuva  = "#ff7f0e", timor_seca  = "#ffbb78",
                  val_chuva    = "#9467bd", val_seca    = "#c5b0d5")

dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

options(bitmapType = "cairo")   # PNG sem X11 no servidor

# Le uma das tabelas geradas por 01_prepare_tables.R
#   tipo: "samples_rel" = abundancia relativa (%) por amostra
#         "groups_rel"  = media da abundancia relativa (%) por grupo
#         "groups_rare" = contagens rarefeitas somadas por grupo (presenca/ausencia)
le_tabela <- function(rank, tipo) {
  arq <- file.path(dir_tab, paste0(marcador, "_", rank, "_", tipo, ".tsv"))
  if (!file.exists(arq)) stop("Tabela nao encontrada: ", arq,
                              "  (rode 01_prepare_tables.R antes)")
  read.delim(arq, check.names = FALSE, stringsAsFactors = FALSE)
}
