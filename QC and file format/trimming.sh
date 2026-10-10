# Script de trimming para reads pareadas
#
# Esta rotina percorre cada diretório de amostras, identifica pares R1/R2,
# remove primers e adaptadores com cutadapt e salva os arquivos processados em
# uma pasta específica para o conjunto pareado.

work_dir="/media/lgbio-nas1/hectorromao/fungi_caves"
out=$work_dir/2.trimming/trimmed
shopt -s nullglob

mkdir -p "$out"
: > "$out/unpaired.txt"

for dir in "$work_dir"/0.raw_data/*/; do
  cave=$(basename "$dir")
  mkdir -p "$out/$cave"

  # 1) Pares: para cada R1, procura o R2 correspondente.
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

  # 2) R2 sem par em R1 são registrados para revisão manual.
  for R2 in "$dir"*_R2_001.fastq*; do
    R1=${R2/_R2_001.fastq/_R1_001.fastq}
    [[ -f "$R1" || -f "$R1.gz" || -f "${R1%.gz}" ]] || echo "SEM R1: $R2" >> "$out/unpaired.txt"
  done
done

cat "$out/unpaired.txt"
