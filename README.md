# Projeto de análise de comunidades microbianas em cavernas

Este workspace reúne os scripts, tabelas e dados usados para a análise de amplicons 16S, desde a etapa de qualidade das leituras até a produção de matrizes por família e gênero e a rarefação da diversidade.

## Visão geral

A organização do projeto foi pensada em etapas sequenciais:

1. Aquisição e organização dos arquivos brutos
2. Controle de qualidade e trimming de primers/adaptadores
3. Processamento com DADA2
4. Agrupamento taxonômico por família e gênero
5. Análise de diversidade por rarefação

## Estrutura das pastas

```text
code_availability/
├── README.md
├── converting_tables_for_fam-gen_lvl.R
├── dada2_pipeline/
│   └── pipeline_truncLen270.R
├── data/
│   ├── asvs/
│   │   └── ... arquivos CSV de ASV + taxonomia
│   ├── familia/
│   │   └── por_amostra/
│   │       └── ... arquivos por amostra, nível família
│   └── genero/
│       └── por_amostra/
│           └── ... arquivos por amostra, nível gênero
├── QC and file format/
│   ├── QC.sh
│   ├── organizing_files.sh
│   ├── single_trimming.sh
│   ├── trimming.sh
│   ├── seqkit_stats_trimmed_P.tsv
│   ├── seqkit_stats_trimmed_S.tsv
│   └── ...
├── rarefaction/
│   ├── inext_rarefaction.R
│   ├── metadata.tsv
│   ├── simple_rarefaction.R
│   └── ...
└── ...
```

## Descrição dos principais diretórios

### `dada2_pipeline/`

Contém a pipeline principal de processamento de reads com `dada2`. O script `pipeline_truncLen270.R` realiza:

- leitura dos arquivos FASTQ pares
- perfil de qualidade
- filtragem e truncamento
- aprendizagem de erros
- denoising
- merge de reads pareadas
- construção da tabela de ASVs
- remoção de quimeras
- atribuição taxonômica
- exportação de CSVs por amostra

### `data/`

Armazena os dados finais em tabelas organizadas por nível taxonômico:

- `asvs/`: tabelas com ASVs e taxonomia associada
- `familia/por_amostra/`: abundância por família, em arquivos separados por amostra
- `genero/por_amostra/`: abundância por gênero, em arquivos separados por amostra

### `QC and file format/`

Agrupa os scripts usados para:

- validação de qualidade (`QC.sh`)
- renomeação e organização dos arquivos FASTQ (`organizing_files.sh`)
- trimming de primers e adaptadores para reads pareadas (`trimming.sh`)
- trimming de reads não pareadas (`single_trimming.sh`)
- estatísticas de sequenciamento em TSV (`seqkit_stats_*`)

### `rarefaction/`

Contém scripts e metadados para avaliação da riqueza e diversidade microbiana:

- `simple_rarefaction.R`: validação da matriz e curvas simples de rarefação
- `inext_rarefaction.R`: análise de diversidade com `iNEXT` e agregação por família/gênero
- `metadata.tsv`: informações amostrais para associação com os dados de microbioma

## Fluxo de análise recomendado

1. Organizar e renomear os arquivos FASTQ
2. Executar os scripts de QC e trimming
3. Rodar a pipeline de DADA2
4. Exportar as tabelas por amostra
5. Gerar matrizes por família/gênero com `converting_tables_for_fam-gen_lvl.R`
6. Realizar a rarefação e comparações de diversidade

## Observações importantes

- Muitos scripts foram escritos para trabalhar com caminhos absolutos em um ambiente específico de armazenamento/servidor. Se você for reproduzir o projeto em outro computador, ajuste os caminhos em cada script antes de executar.
- O projeto usa nomes específicos de arquivos e convenções de amostras que devem ser respeitadas para evitar erros de leitura.
- Os dados finais em CSV são mais fáceis de inspecionar e utilizar em análise estatística e visualização.

## Dicas de prática

- Mantenha os scripts dentro das pastas temáticas e evite misturar arquivos brutos com arquivos processados.
- Sempre preserve uma cópia dos dados originais antes de rodar modificações em nomes ou processamento.
- Use o README como referência central da organização do projeto para facilitar manutenção por outras pessoas.
