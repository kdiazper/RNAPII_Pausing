# PRO-seq Pausing Index Calculator (`calculate_pausing_index.R`)

A high-performance R pipeline designed to compute the **Pausing Index (PI)** across 50 bp windows across whole gene bodies for both **positive (+)** and **negative (-)** strand genes using PRO-seq (Precision Run-On Sequencing) bigWig files and BED annotations, followed by cohort-level pause site calling and filtering.

---

## 📌 Overview

The **Pausing Index (PI)** measures RNA Polymerase II density in a specific genomic window relative to the rest of the gene body. This script executes a two-stage workflow:

### Stage 1: Per-Sample 50 bp Window PI Calculation
1. **Gene Filtering**: Filters input annotations by strand (`+` and `-`) and retains genes $> 50\text{ bp}$ in length.
2. **Window-Level Read Depth**: Divides genes into 50 bp windows and calculates normalized signal (RPM) and estimated raw read counts.
3. **Pausing Index (PI) Calculation**: Computes the PI for each window relative to the remainder of the gene body.
4. **Peak Signal Identification**: Records the maximum single-base PRO-seq signal (`MaxRead_RPM`) and its exact genomic coordinate (`MaxReadCoord`) within each window.
5. **Region Classification**: Tags windows as **`TSS`** (promoter region: within $500\text{ nt}$ downstream of the TSS) or **`GB`** (gene body: downstream of $\text{TSS} + 500\text{ nt}$).

### Stage 2: Cohort Aggregation & Shared Pause Calling
1. **Window Aggregation**: Combines per-sample, per-strand window outputs and computes cross-sample summary statistics (`MeanPI`, `MedianPI`, `Median_PI_Window`), capping extreme PI values at $1000$ and replacing missing (`NA`) values with $0$.
2. **Pause Site Calling ($\text{PI} \ge 15$)**: Identifies candidate pause windows with $\text{PI} \ge 15$ present in $\ge 8$ individuals (`min_samples = 8`).
   * **Promoter-Proximal Pauses (`TSS`)**: Selects the single 50 bp window per gene with the highest median PI among shared candidate windows.
   * **Gene Body Pauses (`GB`)**: Retains all shared candidate windows downstream of the promoter region.
3. **Cohort Integration & Transformation**: Retains measurements across all individuals for the selected pause windows, counts per-gene promoter and gene body pauses per sample, and applies a $\log_{10}(\text{PI} + 1)$ transformation for downstream statistical analysis.

---

## 🛠️ Prerequisites & Dependencies

The script is self-contained for bigWig reading and Rle extraction). Ensure the following standard CRAN and Bioconductor R packages are installed:

```R
install.packages(c("purrr", "parallel", "dplyr", "tidyr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges"))
```

### R Packages Required:
* **`rtracklayer` & `GenomicRanges`**: Imports BED annotations and stranded bigWig signal tracks (`BigWigFile`, `import`, `RleList`).
* **`parallel`**: Enables multi-core parallel processing (`mclapply`).
* **`purrr`**: Fast list-to-dataframe conversion (`map_df`).
* **`dplyr` & `tidyr`**: Data grouping, filtering, joins, and missing-value handling.

---

## 📥 Input Requirements

| Input Type | Description | File Format / Type | Example / Notes |
| :--- | :--- | :--- | :--- |
| **Gene Annotations** | Standard BED file containing full gene coordinates, chromosome, strand, and unique gene names/IDs. | `.bed` | Must contain standard BED fields (`Chr`, `Start`, `End`, `Name`, `Strand`). Genes $\le 50\text{ bp}$ are automatically filtered out. |
| **PRO-seq BigWig Files** | Strand-specific bigWig track files representing normalized genomic read signal/depth (RPM). | `.plus.bw` & `.minus.bw` | Supply the `.plus.bw` path; the script automatically locates the matching `.minus.bw` file in the same directory. |
| **Sample Table** | A dataframe mapping sample identifiers to `.plus.bw` file paths. | `data.frame` in R | Defined in `PROseq.bw.file` within the script. |

### Configuration Parameters in Script:
* `options(mc.cores = 60)`: Number of CPU cores allocated for parallel processing.
* `/path/to/gene_bed_file.bed`: Path to the input gene BED file.
* `PROseq.bw.file`: Dataframe containing `Path` (ending in `.plus.bw`) and `Sample` IDs.
* `input_dir` & `output_dir`: Directories used in Stage 2 to read per-sample window CSVs and write combined/filtered outputs.
* `min_samples <- 8`: Minimum number of individuals required to call a shared pause window (default: `8`).

---

## 🧮 Pausing Index Formula & Methodology

For each 50 bp window $i$ of a gene of length $L$ ($L > 50$):

$$\text{Pausing Index (PI)} = \frac{\text{Reads RPM in Window } i / 50}{\left(\text{Total Gene Reads RPM} - \text{Reads RPM in Window } i\right) / (L - 50)}$$

* **TSS Region Flag**: Assigned as **`TSS`** if the window falls within the first $500\text{ nt}$ downstream of transcription start ($\text{TSS} \to \text{TSS}+500$ for `+` strand; $\text{TSS}-500 \leftarrow \text{TSS}$ for `-` strand). Otherwise tagged as **`GB`** (Gene Body).
* **Strand Handling**: Separately evaluates `+` and `-` strand genes, reversing base-level `Rle` signal vectors for `-` strand features so all windowing proceeds $5' \to 3'$.
* **Statistical Transformation**: Final called pause windows are transformed via $\text{PI}_{\log_{10}} = \log_{10}(\text{PI} + 1)$ to prevent $-\infty$ values when $\text{PI} = 0$.

---

## 📤 Output File Specifications

### 1. Per-Sample, Per-Strand Window Files (Stage 1)
* **Positive Strand Output**: `<SampleID>_posstr_pauses_allwindows.csv`
* **Negative Strand Output**: `<SampleID>_negstr_pauses_allwindows.csv`

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
| **`TSS.or.GB`** | `String` | Classification tag: `TSS` (within $500\text{ nt}$ of TSS) or `GB` (Gene Body). |

### 2. Combined & Filtered Cohort Files (Stage 2)
* **`All_windows_sites.csv`**: All windows combined across samples and strands, including `Sample`, `pause_id` (`Chr:MaxReadCoord:Strand`), `window_id` (`Chr:Start:End:Strand`), `PI_raw`, capped `PI` ($\le 1000$), `MeanPI`, `MedianPI`, and `Median_PI_Window`.
* **`Windows_PI15PerSample_MaxTSSOnly_Shared.csv`**: Final shared pause sites ($\text{PI} \ge 15$ in $\ge 8$ individuals; top median-PI promoter window per gene plus all shared gene body windows), including per-sample pause counts (`Count_Pause_TSS`, `Count_Pause_GB`) and `PI_log10`.

---

## 🚀 How to Run

1. Update the configuration paths and sample definitions in `calculate_pausing_index.R`:
   ```R
   # 1. Set your target CPU core count
   options(mc.cores = 16)

   # 2. Update path to gene BED file (in both + and - strand sections)
   genes <- import("/your/actual/path/to/genes.bed")

   # 3. Update paths to sample BigWig files (.minus.bw must exist in the same directory)
   PROseq.bw.file = data.frame(
     Path = c("/path/to/Sample1.plus.bw", "/path/to/Sample2.plus.bw"),
     Sample = c("Sample1", "Sample2")
   )

   # 4. Update Stage 2 input/output directories
   input_dir  <- "/path/to/allwindows_csv_dir"
   output_dir <- "/path/to/output_dir"
   ```
2. Run the script from the command line:
   ```bash
   Rscript calculate_pausing_index.R
   ```
