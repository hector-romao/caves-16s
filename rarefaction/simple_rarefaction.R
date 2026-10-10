# Verificação e preparação dos dados para rarefação
#
# Este script valida a consistência das tabelas finais de ASVs, constrói a
# matriz de abundância e gera uma curva de rarefação para avaliar se o número
# de reads por amostra é adequado para comparações ecológicas.

library("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(phyloseq)
library(dplyr)

# Arquivos de saída do DADA2.
path_rds <- "/media/lgbio-nas1/hectorromao/fungi_caves/3.dada2/rds_files"

# ----------------------------
# Validação dos dados
# ----------------------------
# Confere se todas as amostras foram processadas corretamente antes do cálculo da rarefação.
dir <- "~/fungi_caves/4.rarefaction/dataset_270-TruncLen"
files <- list.files(dir, pattern = "^clean-.*_asv_taxonomia\\.csv$", full.names = TRUE)
length(files)

dados <- map_dfr(files, read_csv, show_col_types = FALSE)
n_distinct(dados$sample)

# Contagens inteiras e não negativas.
all(dados$abundance == round(dados$abundance)); min(dados$abundance)

# Cada ASV_ID corresponde a uma única sequência e vice-versa.
dados %>% distinct(ASV_ID, sequence) %>% count(ASV_ID) %>% filter(n > 1)
dados %>% distinct(ASV_ID, sequence) %>% count(sequence) %>% filter(n > 1)

# Cada sequência possui uma única taxonomia.
dados %>% distinct(sequence, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  count(sequence) %>% filter(n > 1)

# Verifica se sobraram mitocôndrias ou cloroplastos.
dados %>% filter(Family == "Mitochondria" | Order == "Chloroplast") %>% nrow()

# Profundidade de leitura de cada amostra após limpeza.
dados %>% group_by(sample) %>% summarise(reads = sum(abundance)) %>% arrange(reads)

# ----------------------------
# Preparação do objeto phyloseq
# ----------------------------
# Matriz de abundância: amostras nas linhas, ASVs nas colunas, zeros preenchidos.
mat <- dados %>%
  select(sample, ASV_ID, abundance) %>%
  pivot_wider(names_from = ASV_ID, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Tabela de taxonomia: uma linha por ASV.
tax <- dados %>%
  distinct(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  column_to_rownames("ASV_ID") %>%
  as.matrix()

# Metadados da amostra.
meta <- read.table("/media/lgbio-nas1/hectorromao/fungi_caves/4.rarefaction/metadata.tsv",
                  row.names = 1, header = TRUE)
all(rownames(mat) %in% rownames(meta))

ps <- phyloseq(otu_table(mat, taxa_are_rows = FALSE),
               tax_table(tax),
               sample_data(meta))
ps
sort(sample_sums(ps))
saveRDS(ps, file.path(ps, "ps-filT270.rds"))

# ----------------------------
# Curva de rarefação
# ----------------------------
# Produz a curva de riqueza observada em função da profundidade de leitura.
rc <- rarecurve(mat, step = 2000, tidy = TRUE)
meta_df <- data.frame(sample_data(ps)) %>% tibble::rownames_to_column("sample")
rc <- left_join(rc, meta_df, by = c("Site" = "sample"))

ggplot(rc, aes(Sample, Species, group = Site, colour = cave)) +
  geom_line() +
  geom_vline(xintercept = 30000, linetype = 2) +
  labs(x = "Reads por amostra", y = "ASVs observadas", colour = "Estação")
