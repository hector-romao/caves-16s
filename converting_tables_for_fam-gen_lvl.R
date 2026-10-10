# Conversão de matrizes por família e gênero
#
# Este script reorganiza a tabela de abundância em um formato mais útil para
# exportação por amostra, separando os níveis taxonômicos e salvando um CSV por
# amostra na pasta correspondente.

library(tidyverse)

# ----------------------------
# Agrupamento por família
# ----------------------------
# A matriz fam_mat costuma conter as famílias nas linhas e as amostras nas colunas.
# O código transforma esse formato em uma tabela longa e exporta um arquivo CSV por amostra.
longa <- as.data.frame(fam_mat) %>%
  rownames_to_column("Lineage") %>%
  pivot_longer(-Lineage, names_to = "sample", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  separate(Lineage, into = c("Phylum", "Class", "Order", "Family"),
           sep = "\\|", fill = "right", remove = FALSE) %>%
  mutate(across(Phylum:Family, ~ na_if(.x, "unclassified"))) %>%
  select(sample, Lineage, abundance, Phylum, Class, Order, Family) %>%
  arrange(sample, desc(abundance))

dir.create("familia/por_amostra", recursive = TRUE, showWarnings = FALSE)

longa %>%
  group_by(sample) %>%
  group_walk(~ write.csv(.x %>% mutate(sample = .y$sample),
                         file.path("familia/por_amostra",
                                   paste0("clean-", .y$sample, "_familia.csv")),
                         row.names = FALSE))


# ----------------------------
# Agrupamento por gênero
# ----------------------------
# O mesmo processo é repetido para gen_mat, resultando em tabelas por gênero.
library(tidyverse)

longa <- as.data.frame(gen_mat) %>%
  rownames_to_column("Lineage") %>%
  pivot_longer(-Lineage, names_to = "sample", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  separate(Lineage, into = c("Phylum", "Class", "Order", "Family", "Genus"),
           sep = "\\|", fill = "right", remove = FALSE) %>%
  mutate(across(Phylum:Family, ~ na_if(.x, "unclassified"))) %>%
  select(sample, Lineage, abundance, Phylum, Class, Order, Family) %>%
  arrange(sample, desc(abundance))

dir.create("genero/por_amostra", recursive = TRUE, showWarnings = FALSE)

longa %>%
  group_by(sample) %>%
  group_walk(~ write.csv(.x %>% mutate(sample = .y$sample),
                         file.path("genero/por_amostra",
                                   paste0("clean-", .y$sample, "_genero.csv")),
                         row.names = FALSE))
