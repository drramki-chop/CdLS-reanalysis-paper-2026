# CdLS Reanalysis — Analysis Code

This directory contains the custom scripts and pipeline configuration files used
for the analyses described in the manuscript. It is intended to accompany the
deposited code repository referenced in the manuscript's **Code availability**
statement; the underlying sequencing/expression data referenced by these scripts
are deposited separately on Zenodo/FigShare (see **Data availability**).

## Overview

| Analysis step | Tool (version) | Script(s) | Purpose |
|---|---|---|---|
| Read alignment, small-variant, CNV & SV calling | DRAGEN v3.9 | `run_dragen.sh` | Aligns exome reads, calls SNVs/indels, CNVs, and structural variants (via DRAGEN's bundled Manta-based SV caller) from BAM input |
| Read-depth based CNV calling | CNVpytor v1.3.1 | `combined_scripts_for_paper.sh` | Single-sample CNV calling from CRAM read-depth signal |
| Phenotype-driven variant prioritization | Exomiser CLI v14.0.0 | `run_exomizer.sh`, `CDL-111-00Pexome-analysis.yml` | Prioritizes trio exome variants against patient HPO terms and inheritance models |
| RNA-seq expression outlier detection | OUTRIDER (Bioconductor) | `outrider.R` | Detects aberrantly expressed genes from raw RNA-seq count matrices |
| RNA-seq splicing outlier detection | LeafCutterMD | *(standard workflow, no custom script — see below)* | Detects per-sample aberrant splicing events (outlier intron-excision clusters) from RNA-seq junction reads |

---

## 1. Read alignment, variant calling & SV calling — DRAGEN (`run_dragen.sh`)

SLURM submission script that runs the Illumina DRAGEN v3.9 Bio-IT platform on a
sample BAM, performing alignment, small-variant calling, CNV calling, and
structural-variant (SV) calling in a single pass. SV calling is enabled via
`--enable-sv true`, which invokes DRAGEN's integrated Manta-based SV caller —
this is the source of the Manta-derived SV calls referenced in the manuscript;
no standalone Manta invocation was required.

```bash
#!/bin/bash
#
#SBATCH --job-name=testDragen
#SBATCH --export=NONE
#SBATCH --cpus-per-task=48
#SBATCH --partition=dgddragenq
#SBATCH --nodes=1

ulimit -n 65535

/opt/edico/bin/dragen -r /staging/human/reference/broad_hg38_reference/dragen/ -f \
	--bam-input CDL-111.bam \
	--enable-map-align true \
	--enable-map-align-output true \
	--output-directory ./dragen_output/ \
	--intermediate-results-dir /staging/intermediate.results/ \
	--output-file-prefix CDL-111 \
	--pair-by-name true \
	--enable-cnv true \
	--cnv-enable-self-normalization true \
	--enable-variant-caller true \
	--enable-sv true \
	--output-format cram --enable-vcf-compression true \
	--enable-sort true \
	--enable-duplicate-marking true
```

---

## 2. Read-depth based CNV calling — CNVpytor (`combined_scripts_for_paper.sh`)

Single-sample CNV calling pipeline using CNVpytor v1.3.1, run on aligned CRAM
files. Steps proceed from read-depth import through histogram generation,
partitioning, CNV calling, and filtered export to VCF using the thresholds
shown (`Q0_range`, `p_range`, `dG_range`).

```bash
conda install bioconda::cnvpytor==1.3.1

# CNVpytor v1.3.1 - Single Sample CNV Calling Pipeline

# Step 1 - Read depth import
cnvpytor -v d -j 32 \
  -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
          chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
          chr20 chr21 chr22 chrX chrY chrM \
  -T /path/to/hg38.fasta \
  -root SAMPLE.pytor \
  -rd SAMPLE.cram

# Step 2 - Read depth histogram
cnvpytor -root SAMPLE.pytor -his 500

# Step 3 - Partition
cnvpytor -root SAMPLE.pytor -partition 500

# Step 4 - CNV calling
cnvpytor -root SAMPLE.pytor -call 500

# Step 5 - Filter and export
cnvpytor -root SAMPLE.pytor -view 500 <<ENDL
set Q0_range 0 0.5
set p_range 0 0.01
set dG_range 10000 inf
set print_filename SAMPLE.vcf
print calls
ENDL
```

---

## 3. Phenotype-driven variant prioritization — Exomiser (`run_exomizer.sh`, `CDL-111-00Pexome-analysis.yml`)

`run_exomizer.sh` is a SLURM array job that runs Exomiser CLI v14.0.0 across a
batch of samples (one analysis YAML per array index, listed in
`sampletorun.txt`). `CDL-111-00Pexome-analysis.yml` is a representative
per-sample analysis configuration: it specifies the trio VCF/PED, the proband's
HPO terms, allowed inheritance modes with their maximum allele-frequency
cutoffs, the population frequency and pathogenicity sources used for
filtering, and the prioritization steps (frequency/pathogenicity/inheritance
filters, OMIM and hiPhive prioritizers) and output formats.

### `run_exomizer.sh`

```bash
#!/bin/bash

#SBATCH --export=ALL 
#SBATCH -c 10
#SBATCH --mem-per-cpu=20G
#SBATCH -t 5:000:00
#SBATCH --array=0-116
#SBATCH --output=Array_test.%A_%a.out

SAMPLE_LIST=($(<sampletorun.txt))

YMLFILE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

module load Java/17.0.6.lua

java -jar exomiser-cli-14.0.0.jar --analysis ${YMLFILE}
```

### `CDL-111-00Pexome-analysis.yml` (representative per-sample configuration)

```yaml
analysis:
    genomeAssembly: hg38
    vcf: CDL-111-00P.trio.vcf
    ped: CDL-111-00P.ped
    proband: CDL-111-00P
    hpoIds: ["HP:0001770", "HP:0001377", "HP:0200055", "HP:0004209", "HP:0000347", "HP:0000664", "HP:0002553", "HP:0000687"]
    inheritanceModes: {
      AUTOSOMAL_DOMINANT: 0.1,
      AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,
      AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
      X_DOMINANT: 0.1,
      X_RECESSIVE_HOM_ALT: 0.1,
      X_RECESSIVE_COMP_HET: 2.0,
    }
    analysisMode: PASS_ONLY
    frequencySources: [
        THOUSAND_GENOMES,
        TOPMED,
        UK10K,

        ESP_AFRICAN_AMERICAN, ESP_EUROPEAN_AMERICAN, ESP_ALL,

        EXAC_AFRICAN_INC_AFRICAN_AMERICAN, EXAC_AMERICAN,
        EXAC_SOUTH_ASIAN, EXAC_EAST_ASIAN,
        EXAC_FINNISH, EXAC_NON_FINNISH_EUROPEAN,
        EXAC_OTHER,

        GNOMAD_E_AFR,
        GNOMAD_E_AMR,
        GNOMAD_E_EAS,
        GNOMAD_E_FIN,
        GNOMAD_E_NFE,
        GNOMAD_E_OTH,
        GNOMAD_E_SAS,

        GNOMAD_G_AFR,
        GNOMAD_G_AMR,
        GNOMAD_G_EAS,
        GNOMAD_G_FIN,
        GNOMAD_G_NFE,
        GNOMAD_G_OTH,
        GNOMAD_G_SAS
    ]
    pathogenicitySources: [ REVEL, MVP ]
    steps: [
        failedVariantFilter: { },
        variantEffectFilter: {
          remove: [
              FIVE_PRIME_UTR_EXON_VARIANT,
              FIVE_PRIME_UTR_INTRON_VARIANT,
              THREE_PRIME_UTR_EXON_VARIANT,
              THREE_PRIME_UTR_INTRON_VARIANT,
              NON_CODING_TRANSCRIPT_EXON_VARIANT,
              NON_CODING_TRANSCRIPT_INTRON_VARIANT,
              CODING_TRANSCRIPT_INTRON_VARIANT,
                UPSTREAM_GENE_VARIANT,
                DOWNSTREAM_GENE_VARIANT,
                INTERGENIC_VARIANT,
                REGULATORY_REGION_VARIANT
            ]
        },
        frequencyFilter: {maxFrequency: 2.0},
        pathogenicityFilter: {keepNonPathogenic: true},
        inheritanceFilter: {},
        omimPrioritiser: {},
        hiPhivePrioritiser: {},
    ]
outputOptions:
    outputContributingVariantsOnly: false
    numGenes: 0
    outputFileName: CDL-111-00P-PASS-ONLY
    outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]
```

---

## 4. RNA-seq expression outlier detection — OUTRIDER (`outrider.R`)

R script (adapted from the OUTRIDER Bioconductor vignette) that loads the raw
gene-level RNA-seq count matrix (`cdls_rnaseq_counts.txt`), filters
low-expressed genes, fits the OUTRIDER autoencoder model, and exports
significant aberrant-expression outliers (Benjamini-Hochberg adjusted
p ≤ 0.05) to `outrider_significant_outliers.csv`.

```r
# ==============================================================================
# OUTRIDER Expression Outlier Detection Pipeline
# Adapted from the Bioconductor Vignette
# ==============================================================================

# Load required libraries
library(OUTRIDER)
library(GenomicFeatures)

# ------------------------------------------------------------------------------
# Step 1: Input Data Preparation
# ------------------------------------------------------------------------------
ctsFile <- as.data.frame(data.table::fread("cdls_rnaseq_counts.txt"))
rownames(ctsFile) <- paste0(ctsFile$V1,"_",ctsFile$V2)
ctsFile.matrix <- as.matrix(ctsFile[,3:140])
ods <- OutriderDataSet(countData = ctsFile.matrix)

# ------------------------------------------------------------------------------
# Step 2: Filtering Low-Expressed Genes (Crucial for Autoencoder Performance)
# ------------------------------------------------------------------------------
ods <- filterExpression(ods, 
                        minCounts = TRUE, 
                        filterGenes = TRUE, 
                        fpkmCutoff = 1)

# ------------------------------------------------------------------------------
# Step 3: Run the Core OUTRIDER Pipeline
# ------------------------------------------------------------------------------
ods <- OUTRIDER(ods)
res <- results(ods, padjCutoff = 0.05)

# Export results matrix for your manuscript supplements
write.csv(as.data.frame(res), file = "outrider_significant_outliers.csv", row.names = FALSE)
```

---

## 5. RNA-seq splicing outlier detection — LeafCutterMD (standard workflow)

Splicing outliers were called with **LeafCutterMD**, the per-sample outlier
extension of LeafCutter, run with its default/standard workflow as documented
on the LeafCutter website (https://davidaknowles.github.io/leafcutter/ and
https://github.com/davidaknowles/leafcutter) — no custom modifications were
made to the published pipeline. The three standard stages are:

1. **Convert aligned BAM files to intron-junction files** with the bundled
   `bam2junc.sh` helper (which wraps `regtools junctions extract`), and list
   the resulting `.junc` files for clustering.
2. **Cluster intron-excision events** across all samples with
   `leafcutter_cluster_regtools.py`, producing a per-cluster intron-usage
   count matrix.
3. **Run LeafCutterMD** on the count matrix with `leafcutterMD.R`, which fits
   the per-cluster outlier model and reports, for every intron cluster and
   sample, a p-value and an effect size — corresponding to this study's
   `leafcutter_outlier_clusterPvals.txt.gz` and
   `leafcutter_outlier_effSize.txt.gz` output files.

```bash
# --- Step 1: Convert BAM files to intron-junction files (regtools-based) ---
for bamfile in samples/*.bam
do
    echo "Converting $bamfile to $bamfile.junc"
    sh scripts/bam2junc.sh $bamfile $bamfile.junc
done

# List the junction files to be clustered together
ls samples/*.junc > juncfiles.txt

# --- Step 2: Intron clustering across samples ---
python leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -m 50 \
    -l 500000 \
    -o cdls_outlier

# --- Step 3: Run LeafCutterMD outlier detection ---
Rscript --vanilla leafcutterMD.R \
    --num_threads 4 \
    -o cdls_outlier \
    cdls_outlier_perind_numers.counts.gz
```

The two key outputs from Step 3 — per-cluster outlier p-values
(`cdls_outlier_pVals.txt`) and per-cluster effect sizes
(`cdls_outlier_effSize.txt`) — were renamed/compressed for deposition as
`leafcutter_outlier_clusterPvals.txt.gz` and `leafcutter_outlier_effSize.txt.gz`
respectively (see Data availability).

---

## Repository contents

```
scripts/
├── README.md                          (this file)
├── run_dragen.sh                       DRAGEN alignment / variant / CNV / SV calling
├── combined_scripts_for_paper.sh       CNVpytor read-depth CNV calling pipeline
├── run_exomizer.sh                     Exomiser CLI batch (SLURM array) submission script
├── CDL-111-00Pexome-analysis.yml       Representative per-sample Exomiser analysis configuration
└── outrider.R                          OUTRIDER RNA-seq expression outlier detection
```
