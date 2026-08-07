# Calculate pausing index using 50 bp windows
######################## Loop Over + Strand Genes ######################## 
library(bigWigUtils)
library(purrr)
library(parallel)
library(dplyr)

options(mc.cores=60)

#Read bed file and FPKM GRanges object
genes <- import("/path/to/gene_bed_file.bed")
genes <- genes[which(strand(genes) == "+"),]

# Function to extract raw read counts later
one.read.rpm = function(bigwig.plus.file) {
  y = bw.read.depth.stranded(bigwig.plus.file,
                             GRanges("chr1:1-100000000"))[[1]]
  min(y[y>0])
}

# Load every sample information
PROseq.bw.file = data.frame(Path = c(paste0("/path/to/bigWig_plus"),
                                    paste0("/path/to/bigWig_plus")),
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
library(bigWigUtils)
library(purrr)
library(parallel)
library(dplyr)

options(mc.cores=60)

#Read bed file and FPKM GRanges object
genes <- import("/path/to/gene_bed_file.bed")
genes <- genes[which(strand(genes) == "-"),]

# Function to extract raw read counts later
one.read.rpm = function(bigwig.plus.file) {
  y = bw.read.depth.stranded(bigwig.plus.file,
                             GRanges("chr1:1-100000000"))[[1]]
  min(y[y>0])
}

# Load every sample information
PROseq.bw.file = data.frame(Path = c(paste0("/path/to/bigWig_plus"),
                                    paste0("/path/to/bigWig_plus")),
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

