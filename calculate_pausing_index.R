# Calculate pausing index using 50 bp windows
library(rtracklayer)
library(GenomicRanges)
library(purrr)
library(parallel)
library(dplyr)
library(tidyr)
options(dplyr.width = Inf)

options(mc.cores=60)

# Function to extract raw read counts later
one.read.rpm = function(bigwig.plus.file) {
  y = bw.read.depth.stranded(bigwig.plus.file,
                             GRanges("chr1:1-100000000"))[[1]]
  min(y[y>0])
}

# Function to obtain RLeList from bigWig files
get.rle = function(rle, region){
  r = as.list(rle[region])
  i = (as.character(strand(region)) == "-")
  r[i] = lapply(r[i], rev)
  names(r) = names(region)
  r
}

# Function to read bigWig files as RLeList
bw.read.depth = function(bigWigName, region) {
  bigWig = BigWigFile(bigWigName)
  rle = suppressWarnings(import(bigWig, which = region, as = "RleList"))
  get.rle(rle, region)
}

# Function to read in stranded bigWig files
bw.read.depth.stranded = function(bigwig.plus.file, region, abs.value = T) {
  bigWigPlusFile = bigwig.plus.file
  bigWigMinusFile = sub("\\.plus\\.bw$", ".minus.bw", bigWigPlusFile)
  r = as.list(rep(NA, length(region)))
  plus = (as.character(strand(region)) != "-")
  r[plus] = bw.read.depth(bigWigPlusFile, region[plus])
  r[!plus] = bw.read.depth(bigWigMinusFile, region[!plus])
  if (abs.value) 
    r = lapply(r, abs)
  names(r) = names(region)
  r
}

######################## Loop Over + Strand Genes ######################## 
#Read bed file and FPKM GRanges object
genes <- import("/path/to/gene_bed_file.bed")
genes <- genes[which(strand(genes) == "+"),]
genes <- genes[width(genes) > 50]

# Load every sample information - ensure matching .minus.bw files exist within the same directories as .plus.bw files
PROseq.bw.file = data.frame(Path = c(paste0("/path/to/bigWig.plus.bw"),
                                    paste0("/path/to/bigWig.plus.bw")),
                            Sample = c("SampleID_1","SampleID_2"))

genelist <- list()

write.genelist <- function(gene_name) {
  # Length of gene
  len <- width(genes[which(genes$name == gene_name),]) 
  
  tss <- start(genes[which(genes$name == gene_name),]):(start(genes[which(genes$name == gene_name),])+500)
  
  bw <- bw.allgenes[gene_name]
  
  # Isolate reads for gene
  bw_v2 <- bw[[1]]
  
  # Computes sum of every window.size numbers in x.
  windowed.sum.1 = function(window.size) function(x) {
    # trim x to be a multiple of window.size in length
    n = length(x)
    x = x[ 1 : (n - (n%%window.size)) ]
    s = cumsum(c(Rle(0), x))
    num.windows = length(x) / window.size
    #print(num.windows)
    window.start = window.size * c(0:(num.windows-1)) + 1
    #head(window.start)
    s[window.start+window.size] - s[window.start]
  }
  windowed_sum <- windowed.sum.1(50)(bw_v2) 
  
  #--------- Start to Last 50 BP of Gene ###########
  numerator <- windowed_sum/50
  
  denominator <- (sum(bw_v2)-windowed_sum)/(len-50)
  
  # PI for windows in start to last 50 bp of gene
  pi_start <- numerator/denominator
  
  # Find start and ends
  coord_overall <- (start(genes[which(genes$name == gene_name),]):(end(genes[which(genes$name == gene_name),]-49)))
  coord <- coord_overall[seq(1,length(coord_overall),50)]
  
  # Create dataframe for windows from start to last 50 bp of gene
  gene_start <- data.frame("Chr"=chrom(genes[which(genes$name == gene_name),]), "Start"=coord, "End"= coord+49, "Strand"=strand(genes[which(genes$name == gene_name),]),
                           "Gene.Name"=gene_name, "Reads.RPM"=as.numeric(windowed_sum), "Read.Count"= round(as.numeric(windowed_sum)/one_read_rpm),
                           "PI"= as.numeric(pi_start))
  
  
  #--------- Last 50 BP of Gene ---------# 
  y <- Rle(rep(0, 49))
  bw_v3 <- c(bw_v2, y)
  wind50_sum_subsetgene_pluszero <- runsum(bw_v3,50)
  
  wind50_sum_subsetgene_pluszero_idx <- (length(wind50_sum_subsetgene_pluszero)-48):(length(wind50_sum_subsetgene_pluszero))
  wind50_sum_subsetgene_pluszero_subset <- wind50_sum_subsetgene_pluszero[wind50_sum_subsetgene_pluszero_idx]
  
  windowsize_end <- c(49:1)
  numerator_end <- wind50_sum_subsetgene_pluszero_subset/windowsize_end
  
  denominator_end <- (sum(bw_v2)-wind50_sum_subsetgene_pluszero_subset)/(len-windowsize_end)
  
  # PI for windows in last 50 bp of gene
  pi_end <- numerator_end/denominator_end
  
  # Find start and ends
  coord_end <- (end(genes[which(genes$name == gene_name),])-48):(end(genes[which(genes$name == gene_name),]))
  
  # Create dataframe for last 50 bp of gene
  gene_end <- data.frame("Chr"=chrom(genes[which(genes$name == gene_name),]), "Start"= coord_end, "End"= (end(genes[which(genes$name == gene_name),])), "Strand"=strand(genes[which(genes$name == gene_name),]),
                         "Gene.Name"=gene_name, "Reads.RPM"=as.numeric(wind50_sum_subsetgene_pluszero_subset), "Read.Count"= round(as.numeric(wind50_sum_subsetgene_pluszero_subset)/one_read_rpm),
                         "PI"= as.numeric(pi_end))
  end_of_start <- tail(gene_start,1)
  gene_end <- gene_end[which(gene_end$Start == (end_of_start$End)+1),]
  
  
  #--------- Whole Gene ---------# 
  gene_whole <- rbind(gene_start, gene_end)
  
  gene_subset <- gene_whole
  rm(gene_whole)
  
  if (nrow(gene_subset) >= 1) {
    maxcoord_vector <- rep(NA,nrow(gene_subset))
    max_vector <- rep(NA,nrow(gene_subset))
    gene_subset$start_tx <- gene_subset$Start-start(genes[mcols(genes)$name == gene_name])+1
    gene_subset$end_tx <- gene_subset$start_tx+49
    
    for (row_number in 1:nrow(gene_subset)){
      i = gene_subset[row_number, "start_tx"]
      idx <- which.max(bw_v3[i:(i+49)])
      maxcoord_vector[row_number] <- gene_subset$Start[row_number] + idx - 1
      max_vector[row_number] <- max(bw_v3[i:(i+49)])
    }
    
    gene_subset$MaxRead_RPM <- max_vector
    gene_subset$MaxReadCoord <- maxcoord_vector
    gene_subset$start_tx <- NULL
    gene_subset$end_tx <- NULL
    
    # Add TSS or GB Flag
    gene_subset$TSS.or.GB <- ifelse(gene_subset$Start %in% tss, "TSS", "GB")
    
    genelist[[gene_name]] <- gene_subset
    return(gene_subset)
    
  }
  else {
    print(paste0(gene_name," ","is skipped"))
  }
}

# Extract PI per window per gene per sample
for (sample_name in PROseq.bw.file$Sample) {
  one_read_rpm <- one.read.rpm(PROseq.bw.file[PROseq.bw.file$Sample == sample_name, "Path"])
  
  bw.allgenes <- bw.read.depth.stranded(PROseq.bw.file[PROseq.bw.file$Sample == sample_name, "Path"],
                                        genes)
  names(bw.allgenes) <- genes$name
  
  print(paste0(Sys.time()," ", sample_name))
  
  # For every gene in every sample
  genelist <- mclapply(genes$name, write.genelist)
  
  genelist_df <- map_df(genelist, ~as.data.frame(.x), .id="id")
  genelist_df$id <- NULL
  genelist_df <- subset(genelist_df, !is.na(genelist_df$Read.Count))
  genelist_df$.x <- NULL
  
  print(paste0(Sys.time()," ", "writing windows"))
  
  write.csv(genelist_df, file = paste0("/path/to/",sample_name,"_posstr_pauses_allwindows.csv"), row.names = F)
}


######################## Loop Over - Strand Genes ######################## 
#Read bed file and FPKM GRanges object
genes <- import("/path/to/gene_bed_file.bed")
genes <- genes[which(strand(genes) == "-"),]
genes <- genes[width(genes) > 50]

#Load every sample information - ensure matching .minus.bw files exist within the same directories as .plus.bw files
PROseq.bw.file = data.frame(Path = c(paste0("/path/to/bigWig.plus.bw"),
                                    paste0("/path/to/bigWig.plus.bw")),
                            Sample = c("SampleID_1","SampleID_2"))
genelist <- list()

write.genelist <- function(gene_name) {
  
  # Length of gene
  len <- width(genes[which(genes$name == gene_name),]) 
  
  tss <- end(genes[which(genes$name == gene_name),]):(end(genes[which(genes$name == gene_name),])-500)
  
  bw <- bw.allgenes[gene_name]
  
  # Isolate reads for gene
  bw_v2 <- bw[[1]]
  
  # Computes sum of every window.size numbers in x.
  windowed.sum.1 = function(window.size) function(x) {
    # trim x to be a multiple of window.size in length
    n = length(x)
    x = x[ 1 : (n - (n%%window.size)) ]
    s = cumsum(c(Rle(0), x))
    num.windows = length(x) / window.size
    #print(num.windows)
    window.start = window.size * c(0:(num.windows-1)) + 1
    #head(window.start)
    s[window.start+window.size] - s[window.start]
  }
  windowed_sum <- windowed.sum.1(50)(bw_v2) 
  
  
  #--------- Start to Last 50 BP of Gene ---------# 
  numerator <- windowed_sum/50
  
  denominator <- (sum(bw_v2)-windowed_sum)/(len-50)
  
  # PI for windows in start to last 50 bp of gene
  pi_start <- numerator/denominator
  
  # Find start and ends
  coord_overall <- end(genes[which(genes$name == gene_name),]):(start(genes[which(genes$name == gene_name),]-49))
  coord <- coord_overall[seq(1,length(coord_overall),50)]
  
  # Create dataframe for windows from start to last 50 bp of gene
  gene_start <- data.frame("Chr"=chrom(genes[which(genes$name == gene_name),]), "Start"=coord-49, "End"= coord, "Strand"=strand(genes[which(genes$name == gene_name),]),
                           "Gene.Name"=gene_name, "Reads.RPM"=as.numeric(windowed_sum), "Read.Count"= round(as.numeric(windowed_sum)/one_read_rpm),
                           "PI"= as.numeric(pi_start))
  
  
  #--------- Last 50 BP of Gene ---------# 
  y <- Rle(rep(0, 49))
  bw_v3 <- c(bw_v2, y)
  wind50_sum_subsetgene_pluszero <- runsum(bw_v3,50)
  
  wind50_sum_subsetgene_pluszero_idx <- (length(wind50_sum_subsetgene_pluszero)-48):(length(wind50_sum_subsetgene_pluszero))
  wind50_sum_subsetgene_pluszero_subset <- wind50_sum_subsetgene_pluszero[wind50_sum_subsetgene_pluszero_idx]
  
  windowsize_end <- c(49:1)
  numerator_end <- wind50_sum_subsetgene_pluszero_subset/windowsize_end
  
  denominator_end <- (sum(bw_v2)-wind50_sum_subsetgene_pluszero_subset)/(len-windowsize_end)
  
  # PI for windows in last 50 bp of gene
  pi_end <- numerator_end/denominator_end
  
  # Find start and ends
  coord_end <- (start(genes[which(genes$name == gene_name),])+48):(start(genes[which(genes$name == gene_name),]))
  
  # Create dataframe for last 50 bp of gene
  gene_end <- data.frame("Chr"=chrom(genes[which(genes$name == gene_name),]), "Start"= (start(genes[which(genes$name == gene_name),])), "End"= coord_end, "Strand"=strand(genes[which(genes$name == gene_name),]),
                         "Gene.Name"=gene_name, "Reads.RPM"=as.numeric(wind50_sum_subsetgene_pluszero_subset), "Read.Count"= round(as.numeric(wind50_sum_subsetgene_pluszero_subset)/one_read_rpm),
                         "PI"= as.numeric(pi_end))
  end_of_start <- tail(gene_start,1)
  gene_end <- gene_end[which(gene_end$End == (end_of_start$Start)-1),]
  
  
  #--------- Whole Gene ---------# 
  gene_whole <- rbind(gene_start, gene_end)
  
  gene_subset <- gene_whole
  rm(gene_whole)
  gene_subset$end_tx <- end(genes[mcols(genes)$name == gene_name])-gene_subset$End+1
  gene_subset$start_tx <- gene_subset$end_tx-49
  
  if (nrow(gene_subset) >= 1) {
    maxcoord_vector <- rep(NA,nrow(gene_subset))
    max_vector <- rep(NA,nrow(gene_subset))
    gene_subset$end_tx <- end(genes[mcols(genes)$name == gene_name])-gene_subset$End+1
    gene_subset$start_tx <- gene_subset$end_tx-49
    
    for (row_number in 1:nrow(gene_subset)){
      i = gene_subset[row_number, "end_tx"]
      idx <- which.max(bw_v3[i:(i+49)])
      maxcoord_vector[row_number] <- gene_subset$End[row_number] - idx + 1
      max_vector[row_number] <- max(bw_v3[i:(i+49)])
    }
    
    gene_subset$MaxRead_RPM <- max_vector
    gene_subset$MaxReadCoord <- maxcoord_vector
    gene_subset$start_tx <- NULL
    gene_subset$end_tx <- NULL
    
    # Add TSS or GB Flag
    gene_subset$TSS.or.GB <- ifelse(gene_subset$End %in% tss, "TSS", "GB")
    
    genelist[[gene_name]] <- gene_subset
    return(gene_subset)    
  }
  else {print(paste0(gene_name," ","is skipped"))}
}

# Extract PI per window per gene per sample
for (sample_name in PROseq.bw.file$Sample) {
  
  one_read_rpm <- one.read.rpm(PROseq.bw.file[PROseq.bw.file$Sample == sample_name, "Path"])
  
  bw.allgenes <- bw.read.depth.stranded(PROseq.bw.file[PROseq.bw.file$Sample == sample_name, "Path"],
                                        genes)
  names(bw.allgenes) <- genes$name
  
  print(paste0(Sys.time()," ", sample_name))
  
  # For every gene in every sample
  genelist <- mclapply(genes$name, write.genelist)
  
  genelist_df <- map_df(genelist, ~as.data.frame(.x), .id="id")
  genelist_df$id <- NULL
  genelist_df <- subset(genelist_df, !is.na(genelist_df$Read.Count))
  genelist_df$.x <- NULL
  
  print(paste0(Sys.time()," ", "writing windows"))
  
  write.csv(genelist_df, file = paste0("/path/to/",sample_name,"_negstr_pauses_allwindows.csv"), row.names = F)
}


######### Combine PI Windows #########
input_dir  <- "/path/to/allwindows_csv_dir"
output_dir <- "/path/to/output_dir"

files <- list.files(path = input_dir, pattern = "_allwindows\\.csv$", full.names = TRUE)
sites <- lapply(files, read.csv)
print(paste0(Sys.time(), " Finished Reading All Files"))

# Extract sample IDs from filenames
sample_ids <- gsub("(_posstr|_negstr)_pauses_allwindows\\.csv$", "", basename(files))
names(sites) <- sample_ids

all_sites <- bind_rows(sites, .id = "Sample")
print(paste0(Sys.time(), " Converted List to DF"))

all_sites <- all_sites %>%
  mutate(
    pause_id  = paste(Chr, MaxReadCoord, Strand, sep = ":"),
    window_id = paste(Chr, Start, End, Strand, sep = ":"),
    PI_raw    = PI,
    PI        = ifelse(is.na(PI), 0, pmin(PI, 1000)) # cap at 1000, set NA to 0
  ) %>%
  group_by(Gene.Name, window_id) %>%
  mutate(
    MeanPI           = mean(PI_raw, na.rm = TRUE),
    MedianPI         = median(PI_raw, na.rm = TRUE),
    Median_PI_Window = median(PI)
  ) %>%
  ungroup()

print(paste0(Sys.time(), " Writing All Windows - ", nrow(all_sites), " rows total"))
write.csv(all_sites, file = file.path(output_dir, "All_windows_sites.csv"), quote = FALSE, row.names = FALSE)


######### Call Shared Pause Sites (PI >= 15 in >= 8 of 11 Individuals) ########
min_samples <- 8 # Minimum number of individuals required

# Promoter-proximal pauses
tss_selected_windows <- all_sites %>%
  filter(TSS.or.GB == "TSS", PI >= 15) %>%
  group_by(Gene.Name, window_id) %>%
  filter(n_distinct(Sample) >= min_samples) %>%
  summarise(Median_PI_Window = first(Median_PI_Window), .groups = "drop") %>%
  group_by(Gene.Name) %>%
  slice_max(Median_PI_Window, n = 1, with_ties = FALSE) %>%
  ungroup()

# Retain data across all individuals for the selected promoter-proximal pause windows
tss_shared_all <- all_sites %>%
  semi_join(tss_selected_windows, by = c("Gene.Name", "window_id")) %>%
  group_by(Sample, Gene.Name) %>%
  mutate(Count_Pause_TSS = n()) %>%
  ungroup()

# Get GB pauses/windows
gb_selected_windows <- all_sites %>%
  filter(TSS.or.GB == "GB", PI >= 15) %>%
  group_by(Gene.Name, window_id) %>%
  filter(n_distinct(Sample) >= min_samples) %>%
  summarise(.groups = "drop")

# Retain data across all individuals for the selected gene body pause windows
gb_shared_all <- all_sites %>%
  semi_join(gb_selected_windows, by = c("Gene.Name", "window_id")) %>%
  group_by(Sample, Gene.Name) %>%
  mutate(Count_Pause_GB = n()) %>%
  ungroup()

# Combine TSS and GB pause calls
tss_counts <- distinct(tss_shared_all, Sample, Gene.Name, Count_Pause_TSS)
gb_counts  <- distinct(gb_shared_all, Sample, Gene.Name, Count_Pause_GB)

all_sites_pi15 <- bind_rows(tss_shared_all, gb_shared_all) %>%
  select(-any_of(c("Count_Pause_TSS", "Count_Pause_GB"))) %>%
  left_join(tss_counts, by = c("Sample", "Gene.Name")) %>%
  left_join(gb_counts,  by = c("Sample", "Gene.Name")) %>%
  mutate(
    Count_Pause_TSS = tidyr::replace_na(Count_Pause_TSS, 0L),
    Count_Pause_GB  = tidyr::replace_na(Count_Pause_GB, 0L),
    # Add small pseudocount (+1) so PI == 0 does not become -Inf on log10 scale
    PI_log10        = log10(PI + 1)
  )

write.csv(
  all_sites_pi15,
  file = file.path(output_dir, "Windows_PI15PerSample_MaxTSSOnly_Shared.csv"),
  quote = FALSE,
  row.names = FALSE
)

