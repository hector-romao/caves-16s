# Pipeline DADA2 para amplicons de 16S com truncamento em 270 pb
#
# Este script executa a etapa de qualidade, filtragem, denoising, merge de reads,
# remoção de quimeras e atribuição taxonômica. Ele é a base do processamento
# das amostras e produz os arquivos finais em formato CSV e RDS.

library(dada2)

# Diretórios usados na análise.
# Importante: estes caminhos estão fixados para um armazenamento externo.
# Ajuste para o ambiente local antes de executar em outra máquina.
path_reads <- "/media/lgbio-nas1/hectorromao/fungi_caves/2.trimming/trimmed_paired"
path_pipeline <- "/media/lgbio-nas1/hectorromao/fungi_caves/3.dada2/trimmed_paired"

fnFs <- sort(list.files(path_reads, pattern = "_R1.fastq.gz", full.names = TRUE, recursive = TRUE))
fnRs <- sort(list.files(path_reads, pattern = "_R2.fastq.gz", full.names = TRUE, recursive = TRUE))
sample.names <- sub("_R1.fastq.gz", "", basename(fnFs))

# Perfil de qualidade: define onde cortar
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

# Filtragem de leituras com parâmetros conservadores.
# A opção truncLen = c(270, 270) foi escolhida para manter qualidade e reduzir perdas.
filtFs <- file.path(path_pipeline, "filtered", paste0(sample.names, "_F_filt2.fastq.gz"))
filtRs <- file.path(path_pipeline, "filtered", paste0(sample.names, "_R_filt2.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                     truncLen = c(270, 270), compress = TRUE, multithread = TRUE)
head(out)
tail(out)

# Estimação de erros para cada extremo.
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Denoising (inferência de ASVs)
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Verificação do objeto gerado pela denoising.
dadaFs[[1]]
dadaRs[[1]]

# Merge de reads pareadas.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                      minOverlap = 12, maxMismatch = 2, verbose = TRUE)

# Constrói a tabela de abundância de ASVs.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspeciona o tamanho das sequências.
table(nchar(getSequences(seqtab)))

# Remove quimeras e checa a proporção de reads mantidas.
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim) / sum(seqtab)

# Verificando amostras e tamanho de sequência.
seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% 399:434]
sum(seqtab2) / sum(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab2, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
sum(seqtab.nochim) / sum(seqtab2)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                    "merged", "lengthfilt", "nonchim")
rownames(track) <- sample.names
track

# Atribuição de taxonomia usando bases SILVA.
path_db_species <- "/media/hd16-databases/silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
path_db_genus <- "/media/hd16-databases/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz"
path_db_assignSpecies <- "/media/hd16-databases/silva/silva_v138.2_assignSpecies.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, path_db_genus, multithread = TRUE)
taxa <- addSpecies(taxa, path_db_assignSpecies)

saveRDS(taxa, file.path(path_rds, "taxa.rds"))

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# ----------------------------
# Exportação para tabelas por amostra
# ----------------------------
# Cria um arquivo CSV por amostra contendo abundância e taxonomia das ASVs.
library(dplyr)
library(tidyr)
library(tibble)

seqs <- colnames(seqtab.nochim)
asv_ids <- sprintf("ASV%04d", seq_along(seqs))

tax_df <- as.data.frame(taxa) %>%
  rownames_to_column("sequence") %>%
  mutate(ASV_ID = asv_ids[match(sequence, seqs)]) %>%
  select(ASV_ID, everything())

abund_long <- as.data.frame(seqtab.nochim) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sequence", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  left_join(tax_df, by = "sequence") %>%
  select(sample, ASV_ID, abundance, Kingdom, Phylum, Class,
         Order, Family, Genus, Species, sequence) %>%
  arrange(sample, desc(abundance))

out_dir <- "/media/lgbio-nas1/hectorromao/fungi_caves/3.dada2/trimmed_paired/filtro_270"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lista_amostras <- split(abund_long, abund_long$sample)

for (nome in names(lista_amostras)) {
  write.csv(lista_amostras[[nome]],
            file = file.path(out_dir, paste0(nome, "_asv_taxonomia.csv")),
            row.names = FALSE)
}

salvar_tabelas <- function(seqtab, taxa_obj, out_dir) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  seqs <- colnames(seqtab)
  asv_ids <- sprintf("ASV%04d", seq_along(seqs))

  tax_df <- as.data.frame(taxa_obj) %>%
    rownames_to_column("sequence") %>%
    mutate(ASV_ID = asv_ids[match(sequence, seqs)]) %>%
    select(ASV_ID, everything())

  niveis <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (n in setdiff(niveis, colnames(tax_df))) tax_df[[n]] <- NA_character_

  abund_long <- as.data.frame(seqtab) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "sequence", values_to = "abundance") %>%
    filter(abundance > 0) %>%
    left_join(tax_df, by = "sequence") %>%
    select(sample, ASV_ID, abundance, all_of(niveis), sequence) %>%
    arrange(sample, desc(abundance))

  lista_amostras <- split(abund_long, abund_long$sample)
  for (nome in names(lista_amostras)) {
    write.csv(lista_amostras[[nome]],
              file = file.path(out_dir, paste0(nome, "_asv_taxonomia.csv")),
              row.names = FALSE)
  }

  otu_tab <- seqtab
  colnames(otu_tab) <- asv_ids
  write.csv(otu_tab, file.path(out_dir, "tabela_abundancia_ASV.csv"))
  write.csv(tax_df, file.path(out_dir, "tabela_taxonomia_ASV.csv"), row.names = FALSE)

  message("Salvas ", length(lista_amostras), " tabelas em: ", out_dir)
}

# Salvamento de objetos intermediários para reuso posterior.
path_rds <- "/media/lgbio-nas1/hectorromao/fungi_caves/3.dada2/rds_files"
saveRDS(errF, file.path(path_rds, "errF-filT270.rds"))
saveRDS(errR, file.path(path_rds, "errR-filT270.rds"))
saveRDS(dadaFs, file.path(path_rds, "dadaFs-filT270.rds"))
saveRDS(dadaRs, file.path(path_rds, "dadaRs-filT270.rds"))
saveRDS(mergers, file.path(path_rds, "mergers-filT270.rds"))
saveRDS(seqtab.nochim, file.path(path_rds, "seqtab.nochim_filT270.rds"))
