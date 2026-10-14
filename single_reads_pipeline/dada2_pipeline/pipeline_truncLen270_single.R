## =====================================================================
## DADA2 para amostras SINGLE-END (sem par)
##
##   R1 : 16S-K                (so tem R1, comeca no primer forward 341F)
##   R2 : 16S-L, 16S-M, 16S-N  (so tem R2, comeca no primer reverse 805R)
##
## Por que R1 e R2 sao processados SEPARADOS:
##  - R1 e R2 cobrem extremidades diferentes do amplicon, e R2 esta na fita
##    oposta. A mesma bacteria gera sequencias diferentes em cada um, entao
##    nao dá para colocar os dois numa mesma tabela de ASVs.
##  - Cada orientacao tem o seu proprio modelo de erro (learnErrors).
##  - Nao ha mergePairs, e o filtro de tamanho 399:434 (usado nos pares
##    fundidos) nao se aplica: todas as sequencias ficam com o tamanho do
##    truncLen.
##
## Saidas (em out_base/R1 e out_base/R2):
##   <amostra>_asv_taxonomia.csv, tabela_abundancia_ASV.csv, tabela_taxonomia_ASV.csv
## =====================================================================

library(dada2)
library(dplyr)
library(tidyr)
library(tibble)

## ---- 1. Caminhos e parametros ------------------------------------------

path_reads    <- "/media/lgbio-nas1/hectorromao/fungi_caves/2.trimming/trimmed_single"
path_pipeline <- "/media/lgbio-nas1/hectorromao/fungi_caves/3.dada2/trimmed_single"
path_rds      <- "/media/lgbio-nas1/hectorromao/fungi_caves/3.dada2/rds_files"
out_base      <- file.path(path_pipeline, "filtro_270")

path_db_genus         <- "/media/hd16-databases/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz"
path_db_assignSpecies <- "/media/hd16-databases/silva/silva_v138.2_assignSpecies.fa.gz"

# Depois do cutadapt, R1 tem no maximo ~283 pb e R2 ~279 pb. Reads mais curtos
# que o truncLen sao DESCARTADOS. Em single-end nao ha restricao de overlap,
# entao escolha o corte so pela qualidade (veja os PDFs quality_R1/R2.pdf).
tl <- c(R1 = 270, R2 = 270)

mt <- 16    # nucleos (servidor compartilhado: evite TRUE, que usa todos)

dir.create(path_pipeline, showWarnings = FALSE, recursive = TRUE)
dir.create(path_rds,      showWarnings = FALSE, recursive = TRUE)
dir.create(out_base,      showWarnings = FALSE, recursive = TRUE)

## ---- 2. Funcao que roda o pipeline para UMA orientacao ------------------

processa_single <- function(orient, tl) {

  stopifnot(orient %in% c("R1", "R2"))

  # --- 2.1 Arquivos desta orientacao ---
  padrao <- paste0("_", orient, "\\.fastq\\.gz$")
  fns <- sort(list.files(path_reads, pattern = padrao, full.names = TRUE, recursive = TRUE))
  if (length(fns) == 0) stop("Nenhum arquivo ", orient, " encontrado em ", path_reads)

  sample.names <- sub(padrao, "", basename(fns))
  names(fns) <- sample.names
  cat("\n=====", orient, "-", length(fns), "amostra(s):",
      paste(sample.names, collapse = ", "), "=====\n")

  # --- 2.2 Perfil de qualidade (salvo em PDF; ajuda a escolher o truncLen) ---
  pdf(file.path(path_pipeline, paste0("quality_", orient, ".pdf")), width = 10, height = 8)
  print(plotQualityProfile(fns))
  dev.off()

  # --- 2.3 Filtragem ---
  dir_filt <- file.path(path_pipeline, paste0("filtered_", orient))
  dir.create(dir_filt, showWarnings = FALSE, recursive = TRUE)
  filt <- file.path(dir_filt, paste0(sample.names, "_", orient, "_filt.fastq.gz"))
  names(filt) <- sample.names

  # Single-end: truncLen e maxEE sao numeros unicos (nao vetores c(F, R))
  out <- filterAndTrim(fns, filt, truncLen = tl[[orient]],
                       maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
                       compress = TRUE, multithread = mt)
  rownames(out) <- sample.names
  print(out)

  # Amostras sem nenhum read apos o filtro nao geram arquivo e quebrariam o dada()
  keep <- out[, "reads.out"] > 0
  if (!all(keep)) {
    warning("Sem reads apos a filtragem (removidas): ",
            paste(sample.names[!keep], collapse = ", "))
  }
  if (!any(keep)) stop("Nenhuma amostra ", orient, " passou na filtragem. Reduza o truncLen.")
  filt <- filt[keep]

  # --- 2.4 Modelo de erro ---
  err <- learnErrors(filt, multithread = mt)
  pdf(file.path(path_pipeline, paste0("errors_", orient, ".pdf")), width = 8, height = 8)
  print(plotErrors(err, nominalQ = TRUE))
  dev.off()

  # --- 2.5 Inferencia de ASVs ---
  dd <- dada(filt, err = err, multithread = mt)
  # Com UMA amostra o dada() devolve um objeto e nao uma lista (caso do K)
  if (inherits(dd, "dada")) dd <- list(dd)
  names(dd) <- names(filt)

  # --- 2.6 Tabela de sequencias e quimeras ---
  seqtab <- makeSequenceTable(dd)
  cat("Dimensao da tabela (amostras x sequencias):", dim(seqtab), "\n")
  print(table(nchar(getSequences(seqtab))))      # deve concentrar-se no truncLen

  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = mt, verbose = TRUE)
  cat("Fracao dos reads mantida apos quimeras:",
      round(sum(seqtab.nochim) / sum(seqtab), 3), "\n")

  # --- 2.7 Rastreio dos reads ---
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out[names(filt), , drop = FALSE],
                 sapply(dd, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoised", "nonchim")
  print(track)

  # --- 2.8 Taxonomia ---
  # R2 esta na fita oposta ao gene: tryRC = TRUE testa tambem o reverso-complementar
  tryrc <- (orient == "R2")
  taxa <- assignTaxonomy(seqtab.nochim, path_db_genus, multithread = mt, tryRC = tryrc)
  taxa <- addSpecies(taxa, path_db_assignSpecies, tryRC = tryrc)

  taxa.print <- taxa
  rownames(taxa.print) <- NULL
  print(head(taxa.print))
  cat("Fracao das ASVs sem Familia:", round(mean(is.na(taxa[, "Family"])), 3), "\n")

  # --- 2.9 Objetos do R ---
  sufixo <- paste0("single-", orient, "-filT", tl[[orient]], ".rds")
  saveRDS(err,           file.path(path_rds, paste0("err-",           sufixo)))
  saveRDS(dd,            file.path(path_rds, paste0("dada-",          sufixo)))
  saveRDS(seqtab.nochim, file.path(path_rds, paste0("seqtab.nochim-", sufixo)))
  saveRDS(taxa,          file.path(path_rds, paste0("taxa-",          sufixo)))
  saveRDS(track,         file.path(path_rds, paste0("track-",         sufixo)))

  # --- 2.10 Tabelas em CSV ---
  salvar_tabelas(seqtab.nochim, taxa, file.path(out_base, orient), orient)

  invisible(list(seqtab.nochim = seqtab.nochim, taxa = taxa, track = track))
}

## ---- 3. Funcao que salva as tabelas (uma por amostra + gerais) -----------

salvar_tabelas <- function(seqtab, taxa_obj, out_dir, orient) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # IDs curtos COM o prefixo da orientacao (ASV_R1_0001, ASV_R2_0001) para nao
  # confundir com os IDs da analise pareada (ASV0001), que sao outras sequencias
  seqs    <- colnames(seqtab)
  asv_ids <- sprintf("ASV_%s_%04d", orient, seq_along(seqs))

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
  write.csv(tax_df,  file.path(out_dir, "tabela_taxonomia_ASV.csv"), row.names = FALSE)

  message("Salvas ", length(lista_amostras), " tabelas em: ", out_dir)
}

## ---- 4. Rodar -------------------------------------------------------------

res_R1 <- processa_single("R1", tl)   # K
res_R2 <- processa_single("R2", tl)   # L, M, N

res_R1$track
res_R2$track