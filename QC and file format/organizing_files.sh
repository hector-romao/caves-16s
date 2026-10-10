# Renomeação de reads em lote para padronizar nomes de amostras
#
# Este script reformata os arquivos FASTQ removendo nomes genéricos e incluindo o
# nome da caverna, do marcador e do ponto amostrado. Ele pode operar em modo
# simulação (DRY=1) ou de execução real (DRY=0).

work_dir="/media/lgbio-nas1/hectorromao/fungi_caves"

#!/usr/bin/env bash
# Simulação:  ./renomear_fastq.sh
# Executar:   DRY=0 ./renomear_fastq.sh
DRY=${DRY:-1}
cd ~/fungi_caves/0.raw_data || exit 1

# Tabela de conversão: código da amostra -> ponto amostrado.
declare -A PONTO=(
  [K]="P2_P4_P10" [L]="P1"      [M]="P3"      [N]="F1_F3"
  [O]="F2"        [P]="T1-T5"   [J1]="SD1_SD5" [J2]="SD2-SD4"
  [J3]="SV1-SV4"  [J4]="SV5"    [J5]="SJ1-SJ4" [J6]="SJ5"
  [J7]="SM1"      [J8]="SM2_SM3" [J9]="ST1-ST3" [J10]="ST4"
)

for dir in */; do
  dir=${dir%/}
  IFS=_ read -r cave marker season <<< "$dir"
  cave=${cave^}; marker=${marker^^}

  for f in "$dir"/16S-*.fastq*; do
    [[ -e $f ]] || continue
    base=$(basename "$f")

    # Extrai ID da amostra, R1/R2 e extensão com regex.
    if [[ $base =~ ^16S-(.+)_S[0-9]+_L001_(R[12])_001\.(fastq(\.gz)?)$ ]]; then
      id=${BASH_REMATCH[1]}; sentido=${BASH_REMATCH[2]}; ext=${BASH_REMATCH[3]}
    else
      echo "PULADO (nome fora do padrão): $f" >&2; continue
    fi

    # Busca o ponto; se o ID não estiver na tabela, mantém o ID e avisa.
    ponto=${PONTO[$id]:-$id}
    [[ -z ${PONTO[$id]} ]] && echo "AVISO: '$id' não está na tabela ($f)" >&2

    partes=("$cave" "$marker")
    [[ -n $season ]] && partes+=("$season")
    partes+=("$ponto" "$sentido")
    novo="$(IFS=_; echo "${partes[*]}").$ext"

    if (( DRY )); then
      echo "$f  ->  $dir/$novo"
    else
      mv -n -v "$f" "$dir/$novo"
      printf '%s\t%s\n' "$f" "$dir/$novo" >> mapa_renomeacao.tsv
    fi
  done
done

