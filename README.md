# PRO-seq Pausing Index Calculator (`calculate_pausing_index.R`)

Repository for R code designed to compute the **Pausing Index (PI)** across 50 bp sliding/tiled windows across whole gene bodies for both **positive (+)** and **negative (-)** strand genes using PRO-seq (Precision Run-On Sequencing) bigWig files and BED annotations.

---

## Overview

The **Pausing Index (PI)** measures RNA Polymerase II density in a specific genomic window relative to the rest of the gene body. This script divides genes into **50 bp non-overlapping/sliding windows** and calculates:
1. **Window-level Read Depth** (RPM and estimated raw read counts).
2. **Pausing Index (PI)** for each 50 bp window relative to the remainder of the gene body.
3. **Maximum Read Depth & Coordinates** within each 50 bp window.
4. **Functional Annotation**: Tags windows as **TSS** (Transcription Start Site region: within 500 bp of TSS) or **GB** (Gene Body).

---

## Prerequisites & Dependencies

Ensure the following R packages are installed before running the script:

```R
install.packages(c("purrr", "parallel", "dplyr"))

# Install GenomicRanges from Bioconductor or source
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
```

### R Packages Required:
* **`GenomicRanges` / `rtracklayer`**: Imports and manages genomic ranges/BED files.
* **`parallel`**: Enables multi-core parallel processing (`mclapply`).
* **`purrr`**: Fast data structure manipulation (`map_df`).
* **`dplyr`**: Data filtering and transformation.

---

## Input Requirements

| Input Type | Description | File Format / Type | Example / Notes |
| :--- | :--- | :--- | :--- |
| **Gene Annotations** | Standard BED file containing full gene coordinates, chromosome, strand, and gene names/IDs. | `.bed` | Must contain standard BED fields (Chr, Start, End, Name, Strand). |
| **PRO-seq BigWig Files** | BigWig track files representing normalized genomic read signal/depth. | `.bw` or `.bigWig` | Separate or combined stranded tracks. |
| **Sample Sheet/Table** | A dataframe mapping sample identifiers to bigWig file paths. | `data.frame` in R | Hardcoded or read via script (`PROseq.bw.file`). |

### Configuration Parameters in Script:
* `options(mc.cores = 60)`: Set the number of CPU cores for parallel processing.
* `/path/to/gene_bed_file.bed`: Path to the input gene BED file.
* `/path/to/bigWig_plus`: Path to strand-specific bigWig file(s).

---

## Pausing Index Formula & Calculation Methodology

For each 50 bp window $i$ of a gene of length $L$:

Pausing Index (PI) = (Reads RPM in Window i)/50 / ((Total Gene Reads RPM - Reads RPM in Window i) / (L-50))

* **TSS Region Flag**: Assigned as **`TSS`** if the window falls within the first 500 bp from the start of transcription ($TSS is TSS+500 for `+` strand; TSS-500 is TSS for `-` strand). Otherwise tagged as **`GB`** (Gene Body).
* **Strand handling**: Separately evaluates `+` and `-` strand genes to maintain correct transcription orientation, coordinate indexing, and TSS region flagging.

---

## Output File Specification

The script generates separate CSV files for each sample and strand:
* Positive Strand Output: `<SampleID>_posstr_pauses_allwindows.csv`
* Negative Strand Output: `<SampleID>_negstr_pauses_allwindows.csv`

### Output CSV Schema:

| Column Name | Data Type | Description |
| :--- | :--- | :--- |
| **`Chr`** | `String` | Chromosome name (e.g., `chr1`). |
| **`Start`** | `Integer` | 1-based start coordinate of the 50 bp window. |
| **`End`** | `Integer` | 1-based end coordinate of the 50 bp window. |
| **`Strand`** | `String` | Gene strand (`+` or `-`). |
| **`Gene.Name`** | `String` | Unique gene identifier from the BED file. |
| **`Reads.RPM`** | `Numeric` | Total normalized PRO-seq signal (Reads Per Million) in the 50 bp window. |
| **`Read.Count`** | `Integer` | Estimated raw read count computed via single-read RPM scaling. |
| **`PI`** | `Numeric` | Calculated Pausing Index for the 50 bp window. |
| **`MaxRead_RPM`** | `Numeric` | Maximum single-base PRO-seq signal intensity within the window. |
| **`MaxReadCoord`** | `Integer` | Genomic coordinate corresponding to the peak signal within the window. |
| **`TSS.or.GB`** | `String` | Classification tag: `TSS` (within 500 bp of TSS) or `GB` (Gene Body). |

---

## How to Run

1. Edit the script paths in `calculate_pausing_index.R`:
   ```R
   # 1. Set your target thread count
   options(mc.cores = 16)

   # 2. Update path to gene BED file
   genes <- import("/your/actual/path/to/genes.bed")

   # 3. Update paths to sample BigWig files
   PROseq.bw.file = data.frame(
     Path = c("/path/to/Sample1.bw", "/path/to/Sample2.bw"),
     Sample = c("Sample1", "Sample2")
   )
   ```
2. Run from command line:
   ```bash
   Rscript calculate_pausing_index.R
   ```
