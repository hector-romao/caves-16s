library("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(vegan)

# Saídas gráficas desta etapa.
dir.create("plots", showWarnings = FALSE)

#Checagem dos dados
dir <- file.path("..", "data", "asvs")
files <- list.files(dir, pattern = "^clean-.*_asv_taxonomia\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("Nenhum CSV de ASV encontrado em: ", dir)
length(files)                          # deve ser 21

dados <- map_dfr(files, read_csv, show_col_types = FALSE)
n_distinct(dados$sample)               # deve ser 21

# Contagens inteiras e não negativas
all(dados$abundance == round(dados$abundance)); min(dados$abundance)

# Cada ASV_ID corresponde a uma única sequência (e vice-versa)?
dados %>% distinct(ASV_ID, sequence) %>% count(ASV_ID)   %>% filter(n > 1)  # vazio
dados %>% distinct(ASV_ID, sequence) %>% count(sequence) %>% filter(n > 1)  # vazio

# Cada sequência tem uma única taxonomia?
dados %>% distinct(sequence, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  count(sequence) %>% filter(n > 1)                                         # vazio

# Sobrou mitocôndria ou cloroplasto?
dados %>% filter(Family == "Mitochondria" | Order == "Chloroplast") %>% nrow()  # 0

# Profundidade de cada amostra depois da limpeza
dados %>% group_by(sample) %>% summarise(reads = sum(abundance)) %>% arrange(reads)




#########################Rarefação#######################################

# Matriz de abundância: amostras nas linhas, ASVs nas colunas, zeros preenchidos
mat <- dados %>%
  select(sample, ASV_ID, abundance) %>%
  pivot_wider(names_from = ASV_ID, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Tabela de taxonomia: uma linha por ASV
tax <- dados %>%
  distinct(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  column_to_rownames("ASV_ID") %>%
  as.matrix()

# Metadados: crie um CSV com uma linha por amostra
meta <- read.delim("metadata.tsv", row.names = 1, check.names = FALSE)
if (!setequal(rownames(mat), rownames(meta))) {
  stop("As amostras das tabelas de ASVs e dos metadados nao correspondem.")
}

ps <- phyloseq(otu_table(mat, taxa_are_rows = FALSE),
               tax_table(tax),
               sample_data(meta))
ps
sort(sample_sums(ps))

saveRDS(ps, "ps-filT270.rds")

#Rarefaction

rc <- rarecurve(mat, step = 2000, tidy = TRUE)    # colunas: Site, Sample, Species
meta_df <- data.frame(sample_data(ps)) %>% tibble::rownames_to_column("sample")
rc <- left_join(rc, meta_df, by = c("Site" = "sample"))

rarefaction_plot <- ggplot(rc, aes(Sample, Species, group = Site, colour = cave)) +
  geom_line() +
  geom_vline(xintercept = 80000, linetype = 2) +
  labs(x = "Reads por amostra", y = "ASVs observadas", colour = "Caverna")
ggsave(file.path("plots", "simple_rarefaction_single.pdf"),
       rarefaction_plot, width = 10, height = 7)
