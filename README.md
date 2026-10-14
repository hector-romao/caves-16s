# Análise de comunidades microbianas de cavernas — 16S

Este repositório reúne o processamento de amplicons 16S e as análises de
composição e diversidade de comunidades microbianas. O fluxo está organizado em
duas estratégias independentes: **reads pareadas** e **reads single-end**.
Resultados de R1 e R2 single-end são mantidos separados, pois representam
extremidades diferentes do amplicon e não devem ser combinados como se fossem
pares.

## Análises incluídas

- Controle de qualidade, organização dos FASTQ e remoção de primers com
  `cutadapt`.
- Inferência de ASVs, remoção de quimeras e atribuição taxonômica com `DADA2`.
- Tabelas de abundância relativa por amostra e por grupo em cinco níveis:
  filo, classe, ordem, família e gênero.
- Curvas de rarefação de ASVs e curvas de diversidade iNEXT para q = 0
  (riqueza), q = 1 (Shannon) e q = 2 (Simpson).
- Comparação da composição entre cavernas/estações com barras, diagramas de
  Venn, UpSet, gráficos aluviais e heat trees.
- Exportação das abundâncias agregadas por família e gênero para arquivos CSV
  por amostra.

## Ferramentas

| Etapa | Ferramentas principais |
|---|---|
| QC e preparação de FASTQ | FastQC, MultiQC, SeqKit, Cutadapt |
| Processamento de amplicons | R, DADA2, SILVA |
| Objeto e manipulação de microbioma | phyloseq, tidyverse |
| Rarefação e diversidade | vegan, iNEXT, patchwork |
| Visualização | ggplot2, viridis, VennDiagram, UpSetR, ggalluvial, metacoder |

Instale as ferramentas de linha de comando e os pacotes R no ambiente em que a
análise será executada. Os scripts não instalam pacotes automaticamente. Para
as etapas de QC e trimming, use Bash em Linux/WSL; os scripts R podem ser
executados com R ou RStudio.

## Organização do repositório

```text
.
├── README.md
├── .gitignore
├── converting_tables_for_fam-gen_lvl.R
├── QC and file format/
│   ├── QC.sh
│   ├── organizing_files.sh
│   ├── trimming.sh
│   ├── single_trimming.sh
│   └── seqkit_stats_trimmed_*.tsv
├── paired_reads_pipeline/
│   ├── dada2_pipeline/
│   │   └── pipeline_truncLen270.R
│   ├── data/
│   │   ├── asvs/
│   │   ├── familia/por_amostra/
│   │   └── genero/por_amostra/
│   ├── plots/
│   │   ├── 00_config.R ... 06_heat_tree.R
│   │   ├── tables_paired/
│   │   └── figures_paired/
│   └── rarefaction/
│       ├── metadata.tsv
│       ├── simple_rarefaction.R
│       ├── inext_rarefaction.R
│       └── plots/
└── single_reads_pipeline/
    ├── dada2_pipeline/
    │   └── pipeline_truncLen270_single.R
    ├── data/
    │   ├── asvs/
    │   ├── familia/por_amostra/
    │   └── genero/por_amostra/
    ├── plots/
    │   ├── 00_config.R ... 06_heat_tree.R
    │   ├── tables/
    │   └── figures/
    └── rarefaction/
        ├── metadata.tsv
        ├── simple_rarefaction.R
        ├── inext_rarefaction.R
        └── plots/
```

## Fluxo de processamento

### 1. QC, organização e trimming

Os scripts em `QC and file format/` verificam as leituras, padronizam nomes,
identificam arquivos sem par e removem primers. `trimming.sh` processa pares
R1/R2; `single_trimming.sh` trata separadamente bibliotecas single-end. Os
scripts contêm caminhos do ambiente original: ajuste as variáveis de diretório
e confirme os nomes dos FASTQ antes de executá-los. Preserve os arquivos brutos.

### 2. Processamento DADA2

`paired_reads_pipeline/dada2_pipeline/pipeline_truncLen270.R` filtra e trunca
reads pareadas, aprende os modelos de erro, denoisa cada extremidade, combina
os pares, remove quimeras e atribui taxonomia. O comprimento de truncamento é
270 pb por extremidade; a tabela final mantém sequências com o comprimento
esperado após o merge.

`single_reads_pipeline/dada2_pipeline/pipeline_truncLen270_single.R` processa
R1 e R2 em execuções separadas, cada qual com seu modelo de erro e atribuição
taxonômica. Não há etapa de merge para essas leituras. Os caminhos dos FASTQ,
dos resultados e das bases SILVA estão declarados no início dos scripts e
devem ser ajustados ao ambiente local.

As tabelas CSV em `data/asvs/` são a entrada das etapas de rarefação e
visualização. As tabelas exportadas por amostra em `data/familia/` e
`data/genero/` contêm abundâncias agregadas por linhagem taxonômica.

### 3. Rarefação e diversidade

Execute cada análise a partir de sua própria pasta `rarefaction/`, depois de
confirmar a correspondência entre `metadata.tsv` e as amostras em `data/asvs/`:

```bash
cd paired_reads_pipeline/rarefaction
Rscript simple_rarefaction.R
Rscript inext_rarefaction.R
```

Repita em `single_reads_pipeline/rarefaction/` para as leituras single-end.
`simple_rarefaction.R` valida as contagens e os metadados, constrói o objeto
`phyloseq` e salva a curva de riqueza observada. `inext_rarefaction.R` produz
curvas iNEXT para as amostras completas e por caverna e também salva objetos
agregados por família e gênero. A linha de referência de profundidade é 30.000
reads para o conjunto pareado e 80.000 para o single-end.

Os arquivos de figura são salvos em `rarefaction/plots/`; objetos agregados
ficam em `rarefaction/aggregated/`. A profundidade indicada é um parâmetro de
referência visual e não substitui a avaliação da cobertura amostral.

### 4. Tabelas e visualizações taxonômicas

Os scripts em `plots/` devem ser executados a partir da pasta correspondente.
Defina `CAVES_16S_PHYLOSEQ` com o caminho do objeto `phyloseq` limpo produzido
pela análise DADA2 antes de gerar as tabelas:

```bash
cd paired_reads_pipeline/plots
export CAVES_16S_PHYLOSEQ=/caminho/para/ps-filT270.rds
Rscript 01_prepare_tables.R
Rscript 02_barplot.R
Rscript 03_venn.R
Rscript 04_upset.R
Rscript 05_alluvial.R
Rscript 06_heat_tree.R
```

Para single-end, use `single_reads_pipeline/plots/` e o objeto `phyloseq`
correspondente. O arquivo `00_config.R` centraliza o objeto de entrada, os
diretórios de saída, a profundidade usada para presença/ausência, a semente,
os níveis taxonômicos, os grupos e as cores.

`01_prepare_tables.R` grava três tabelas TSV por nível taxonômico:

| Sufixo | Conteúdo | Uso |
|---|---|---|
| `samples_rel` | Abundância relativa (%) por amostra, sem rarefação | Barras por amostra |
| `groups_rel` | Média da abundância relativa (%) por grupo | Barras, aluviais e heat trees |
| `groups_rare` | Contagens após rarefação somadas por grupo | Presença/ausência, Venn e UpSet |

As abundâncias permanecem não rarefeitas. Para as comparações de
presença/ausência, as amostras são rarefeitas sem reposição na mesma
profundidade; a semente e essa profundidade são configuradas em `00_config.R`.
As linhagens sem classificação são mantidas nas tabelas, com rótulos
`unclassified_*`, e podem ser excluídas das figuras de interseção.

### 5. Figuras produzidas

| Script | Análise | Arquivos de saída |
|---|---|---|
| `02_barplot.R` | Top táxons e categoria `Other`; barras por amostra e por grupo | `barplot_*.pdf` |
| `03_venn.R` | Táxons compartilhados/exclusivos em comparações de 2–4 grupos | `venn_*.png` |
| `04_upset.R` | Interseções entre vários grupos | `upset_*.pdf` |
| `05_alluvial.R` | Fluxo filo → classe → gênero, ponderado por abundância e por riqueza | `alluvial_*.pdf` |
| `06_heat_tree.R` | Abundância relativa projetada na árvore taxonômica | `heat_tree_*.pdf` |
| `simple_rarefaction.R` | Riqueza observada em função do número de reads | `simple_rarefaction_*.pdf` |
| `inext_rarefaction.R` | Diversidade q = 0, 1 e 2, por amostra e caverna | `*_painel.pdf` |

As figuras são gravadas nas pastas `plots/figures_paired/` (pareadas),
`plots/figures/` (single-end) e `rarefaction/plots/` (rarefação). Os scripts
permitem ajustar níveis taxonômicos, número de táxons exibidos, grupos,
transformação logarítmica e filtros de linhagens não classificadas.

## Requisitos de entrada e reprodutibilidade

- Os CSVs de ASV devem conter `sample`, `ASV_ID`, `abundance` e os níveis
  taxonômicos usados pelo script.
- `metadata.tsv` deve ter uma linha por amostra; a primeira coluna deve conter
  os mesmos identificadores da coluna `sample`. A coluna `cave` define os
  grupos usados nas análises.
- Para a preparação das tabelas de gráficos, o objeto `phyloseq` deve conter
  contagens brutas, taxonomia e metadados compatíveis. O caminho pode ser
  fornecido por `CAVES_16S_PHYLOSEQ`; o valor padrão aponta para o ambiente
  original.
- QC, trimming, bases de referência e arquivos FASTQ brutos não são
  autocontidos neste repositório. Ajuste os caminhos absolutos existentes nos
  scripts antes de reproduzir essas etapas.
- As análises de abundância e presença/ausência usam entradas distintas de
  propósito; não substitua as contagens brutas pelas rarefeitas nos gráficos
  de abundância.

## Arquivos da versão anterior

O repositório remoto já continha os diretórios `dada2_pipeline/`, `data/` e
`rarefaction/` na raiz. Eles são mantidos para preservar a versão anterior; as
análises pareadas e single-end organizadas e documentadas neste README ficam
nas pastas `paired_reads_pipeline/` e `single_reads_pipeline/`.

## Observação sobre os dados single-end

As bibliotecas single-end têm poucas amostras e podem representar tanto R1
quanto R2. A orientação é mantida na identificação e no processamento das ASVs;
por isso, os resultados dessa análise não devem ser somados aos resultados de
reads pareadas nem entre orientações sem uma justificativa metodológica.
