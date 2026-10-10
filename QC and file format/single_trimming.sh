# Script de trimming para reads sem par
#
# Algumas amostras apresentaram apenas um dos pares (R1 ou R2). Neste caso, a
# remoção de primers é feita de forma individual para evitar perda extra de reads.

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
