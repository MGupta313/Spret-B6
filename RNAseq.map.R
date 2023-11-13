# RNA-seq data QC/Mapping/post-processing pipeline. First script in the bulk RNAseq pipeline.

# Read in bistools; sources many important functions for this script
source("/home/boss_lab/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic
source("/home/boss_lab/Apps/bitbucket/sequencing-sample-qc/mapping_barplot.R") # Magic, plotting scripts
library("stringr") # str_split_1()

projectDir = "/home/boss_lab/Projects/Scharer_Projects/F1_hybrid_SpretB6/MZB_FoB_2023.11/RNA/" # Replace, must have '/' at the end

# Set working directory
homeDir = paste0(projectDir, "pipeline/")
setwd(homeDir)

# Manifest file
filesDir = homeDir # Replace, if necessary
filesFile = "RNAseq.sample.manifest.txt" # Replace, if necessary
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = T, as.is = T)

# Set max resources; ESB has 128 threads and 512G, allowing for multiple of these at the default of 32T/32GB to run concurrently
# Replace, adjust according to priority
threads = 32
sortMem = "32G"

# Replace, set STAR genome
# Choose one from available genomes, including: hg19, hg38, mm9, mm10, rn6, MacaM
starGenome = select_ref_STAR("mm10") # Replace

# Adapter trimming options
# Replace, choose one from: nextera, illumina
# "nextera" signifies Nextera Tn5 adapters, "illumina" signifies TRUSEQ adapters
seq_platform = "nextera"

# Picard call
picardCmd = "picard MarkDuplicates "

# Set output directory
outDir = paste0(projectDir, "data/")

# STAR call
starCmd = paste("STAR --runThreadN", threads, "--outSAMtype BAM Unsorted --outSAMunmapped Within KeepPairs --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --alignIntronMax 1000000 --outSAMmultNmax 1 --genomeDir", starGenome, "--outFileNamePrefix")

# Helper functions
is_single_end = function(fastq) { fq = get0(fastq, envir = .GlobalEnv); return(!is.null(fq) && !is.na(fq) && fq != "") }
is_paired_end = function(R1, R2) { return(is_single_end(R1) && is_single_end(R2)) }

# Name of STAR's output
starBamFile = "Aligned.out.bam"

for (i in 1:nrow(files)) {
  
  if (!files$include[i]) { next }
  
  time_start = Sys.time()
  print(paste(files$sample[i], time_start))
  
  print("Copy files to data dir") # Record where the data was copied from first
  # sampleDir = paste0(outDir, files$sample[i], "/")
  # if (!file.exists(paste0(sampleDir, files$fqMate1[i])) || !file.exists(paste0(sampleDir, files$fqMate2[i]))) {
  #   files$fqMate1_input[i] = paste0(files$dir[i], stringr::str_split_1(files$fqMate1[i], "\\|"), collapse = "|")
  #   files$fqMate2_input[i] = paste0(files$dir[i], stringr::str_split_1(files$fqMate2[i], "\\|"), collapse = "|")
  # }
  # mvFiles(files$fqMate1[i], files$dir[i], sampleDir)
  # files$dir[i] = mvFiles(files$fqMate2[i], files$dir[i], sampleDir)
  # 
  print("Format and concatenate fastq files")
  # files$fqMate1[i] = formatFastq(files$fqMate1[i], files$dir[i], paste0(files$sample[i], "_1"), threads = threads)
  # files$fqMate2[i] = formatFastq(files$fqMate2[i], files$dir[i], paste0(files$sample[i], "_2"), threads = threads)
  # 
  print("Fastqc")
  # system(paste0("fastqc ", files$dir[i],files$fqMate1[i], " -o ", files$dir[i]))
  # system(paste0("fastqc ", files$dir[i],files$fqMate2[i], " -o ", files$dir[i]))
  # 
  print("Cut 3 prime adapter sequences") # These files are only temporary and are removed after mapping
  # pelist = fqPECutadapt(files$fqMate1[i], files$fqMate2[i], files$dir[i], seq_platform, threads = threads)
  # fqMate1_trimmed = pelist[[1]]
  # fqMate2_trimmed = pelist[[2]]
  
  ##############################
  print("Step 5")
  print("Starting mapping")
  
  # Set genomes for mapping
  genomeStr = "mm10" # Replace
  B6Genome = "/home/boss_lab/Apps/genomes/STAR/mm10.ERCC"
  SpretGenome = "/home/boss_lab/Apps/genomes/species/spret-b6/SPRET_EiJ.STAR.genome"
  
  # Map to B6 genome
  print("Maternal B6 STAR mapping")
  
  #paired-end STAR alignment based on manifest columns
  B6sample = paste0(files$dir[i], files$sample[i], ".B6")
  genome = B6Genome
  
  if (!file.exists(paste0(files$dir[i], files$fqMate1[i], ".gz"))) {
    cmd1 = paste("gzip -k", paste0(files$dir[i], files$fqMate1[i]))
    cmd2 = paste("gzip -k", paste0(files$dir[i], files$fqMate2[i]))
    system(cmd1)
    system(cmd2)
  }
  
  fq1 = paste0(files$dir[i], files$fqMate1[i], ".gz")
  fq2 = paste0(files$dir[i], files$fqMate2[i], ".gz")
  
  starB6Cmd = paste("STAR --runMode alignReads --runThreadN 1 --genomeDir", genome, "--genomeLoad NoSharedMemory --readFilesIn", fq1, fq2, "--readFilesCommand zcat --outFileNamePrefix", B6sample, "--outSAMtype BAM SortedByCoordinate --outBAMcompression 6 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.06 --outFilterMatchNmin 30 --outSJfilterOverhangMin 30 10 10 10 --seedSearchStartLmax 30 --alignIntronMin 20 --alignIntronMax 20000 --alignEndsType Local")
  
  system(starB6Cmd)
  files$bamB6File = paste0(B6sample, ".bam")
  
  ##############################
  # Map to Spret genome
  print("Paternal Spret STAR mapping")
  
  #paired-end STAR alignment based on manifest columns
  Spretsample = paste0(files$dir[i], files$sample[i], ".Spret")
  genome = SpretGenome
  fq1 = paste0(files$dir[i], files$fqMate1[i], ".gz")
  fq2 = paste0(files$dir[i], files$fqMate2[i], ".gz")
  
  starSpretCmd = paste("STAR --runMode alignReads --runThreadN 1 --genomeDir", genome, "--genomeLoad NoSharedMemory --readFilesIn", fq1, fq2, "--readFilesCommand zcat --outFileNamePrefix", Spretsample, "--outSAMtype BAM SortedByCoordinate --outBAMcompression 6 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.06 --outFilterMatchNmin 30 --outSJfilterOverhangMin 30 10 10 10 --seedSearchStartLmax 30 --alignIntronMin 20 --alignIntronMax 20000 --alignEndsType Local")
  
  system(starSpretCmd)
  files$bamSpretFile = paste0(Spretsample, ".bam")
  
  
  # print("Sort bam")
  # files$bamFile[i] = sortBam(paste0(files$dir[i], starBamFile),
  #                            bamSortFile = paste0(files$dir[i], files$sample[i], ".sort"),
  #                            delBam = T, threads = threads, mem = sortMem)
  # 
  # print("Mark duplicates")
  # files$bamFile[i] = markDups(paste0(files$dir[i], files$bamFile[i]), picardCmd, delBam = T)
  # 
  # print("Get mapping stats")
  # cts = getBamCts(paste0(files$dir[i], files$bamFile[i]))
  # 
  # files$unmapped.reads[i] = cts[1]
  # files$mapped.reads[i] = cts[2]
  # files$unique.reads[i] = cts[3]
  # files$paired.reads[i] = cts[4]
  # 
  # # Calculate time to run this sample
  # files$timeTaken[i] = difftime(Sys.time(), time_start, units = "hours")
  # print(paste("Time taken for sample", files$sample[i], "(hrs):", files$timeTaken[i]))
  # 
  # Update manifest file and note that this sample finished mapping by setting 'include' to FALSE
  files$include[i] = FALSE
  write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)
}

# Update manifest file
files$include = TRUE
write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)

# ####################
# # Plot mapping stats
# 
# plot_mapping_stats(files)

#########################
# Print software versions

# Output R session
capture.output(sessionInfo(), file = paste0("Rsession.Info.", gsub("\\D", "", Sys.time()), ".txt"))

picardTool = "picard"
picardArgs = "MarkDuplicates --version"

# Capture versions of other tools
tools = c("skewer",
          "fastqc",
          "STAR",
          "samtools",
          picardTool)

args = c("-version",  # skewer
         "--version", # fastqc
         "--version", # STAR
         "--version", # samtools
         picardArgs)

names = c("skewer",
          "fastqc",
          "STAR",
          "samtools",
          "MarkDuplicates")

version = lapply(1:length(tools), FUN = function (x) {system2(tools[x], args[x], stdout=TRUE, stderr=TRUE)})
names(version) = names

capture.output(version, file = paste0("Software.versions.", gsub("\\D", "", Sys.time()), ".txt"))
