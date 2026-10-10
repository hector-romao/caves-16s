# Script de controle de qualidade e etapa de trimming
#
# Este script roda FastQC, MultiQC e seqkit para avaliar a qualidade bruta das
# amostras, remover primers/adaptadores e registrar o status de reads não pareadas.

work_dir="/media/lgbio-nas1/hectorromao/fungi_caves"

# Análise de qualidade das amostras com FastQC e MultiQC.
while read -r line; do
    fastqc "$work_dir/0.raw_data/$line"/*fastq* --noextract --nogroup -t 16 -o "$work_dir/1.QC/$line"
done < "$work_dir/caves.list"

multiqc "$work_dir/1.QC/**" -o "$work_dir/multiqc" -n multiqc_report

# Verificando o número de reads por amostra com seqkit.
seqkit stats "$work_dir/0.raw_data/"*/*fastq* -T -a -j 16 > "$work_dir/1.QC/seqkit_stats.tsv"

# ----------------------------
# Remoção de primers e adaptadores
# ----------------------------
work_dir="/media/lgbio-nas1/hectorromao/fungi_caves"
out=$work_dir/2.trimming/trimmed
shopt -s nullglob

mkdir -p "$out"
: > "$out/unpaired.txt"

for dir in "$work_dir"/0.raw_data/*/; do
  cave=$(basename "$dir")
  mkdir -p "$out/$cave"

  # Pareados: para cada R1, localiza o R2 correspondente e remove primers.
  for R1 in "$dir"*_R1_001.fastq*; do
    R2=${R1/_R1_001.fastq/_R2_001.fastq}

    if [[ ! -f "$R2" ]]; then
      echo "SEM R2: $R1" >> "$out/unpaired.txt"
      continue
    fi

    base=$(basename "$R1")
    base=${base%%_R1_001.fastq*}

    cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC \
      --discard-untrimmed -j 16 \
      -o "$out/$cave/${base}_R1.fastq.gz" \
      -p "$out/$cave/${base}_R2.fastq.gz" \
      "$R1" "$R2" \
      > "$out/$cave/${base}.cutadapt.log" < /dev/null
  done

  # Reads R2 que não têm par em R1 são registradas para revisão manual.
  for R2 in "$dir"*_R2_001.fastq*; do
    R1=${R2/_R2_001.fastq/_R1_001.fastq}
    [[ -f "$R1" || -f "$R1.gz" || -f "${R1%.gz}" ]] || echo "SEM R1: $R2" >> "$out/unpaired.txt"
  done
done

cat "$out/unpaired.txt"

# Exceções observadas em amostras com leitura única.
# SEM R1: ...
# SEM R2: ...

raw=/media/lgbio-nas1/hectorromao/fungi_caves/0.raw_data
out=/media/lgbio-nas1/hectorromao/fungi_caves/2.trimming/trimmed_single
mkdir -p "$out"

# K: só tem R1 -> primer forward.
cutadapt -g CCTACGGGNGGCWGCAG \
  --discard-untrimmed -j 16 \
  -o "$out/16S-K_S224_L001_R1.fastq.gz" \
  "$raw/samambaia_16S/16S-K_S224_L001_R1_001.fastq" \
  > "$out/16S-K_S224_L001_R1.cutadapt.log"

# L, M, N: só têm R2 -> primer reverse.
for f in "$raw/samambaia_16S/16S-L_S225_L001_R2_001.fastq" \
         "$raw/samambaia_16S/16S-M_S226_L001_R2_001.fastq" \
         "$raw/funil_16S/16S-N_S227_L001_R2_001.fastq"; do
  base=$(basename "$f" _R2_001.fastq)
  cutadapt -g GACTACHVGGGTATCTAATCC \
    --discard-untrimmed -j 16 \
    -o "$out/${base}_R2.fastq.gz" \
    "$f" > "$out/${base}_R2.cutadapt.log"
done

# Verificação da perda de reads após o corte de primers e adaptadores.
seqkit stats "$work_dir/2.trimming/trimmed_paired/"*/*fastq*gz -T -a -j 16 > "$work_dir/2.trimming/seqkit_stats_trimmed_P.tsv"
seqkit stats "$work_dir/2.trimming/trimmed_single/"*fastq*gz -T -a -j 16 > "$work_dir/2.trimming/seqkit_stats_trimmed_S.tsv"

