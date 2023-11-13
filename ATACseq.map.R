# ATAC-seq QC/Mapping/post-processing script; First script in the bulk ATACseq pipeline.

# Read in bistools; sources many functions important for this script
source("/home/boss_lab/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic
#source("/home/boss_lab/Apps/bitbucket/ATACseq/ATACseq.calc.TSSenrichment.R") # Magic, calcTSSe()
source("/home/boss_lab/Apps/bitbucket/sequencing-sample-qc/mapping_barplot.R") # Magic, plotting scripts
library("stringr") # str_split_1()

projectDir = "/home/boss_lab/Projects/Scharer_Projects/F1_hybrid_SpretB6/MZB_FoB_2023.11/ATAC/" # Replace, must have a '/' at the end

# Set working directory
homeDir = paste0(projectDir, "pipeline/")
setwd(homeDir)

# Manifest file
filesDir = homeDir # Replace, if necessary
filesFile = "ATACseq.sample.manifest.txt" # Replace, if necessary
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = TRUE, as.is = T)

# Set max resources; ESB has 128 threads and 512G, allowing for multiple of these at the default of 32T/32GB to run concurrently
# Replace, adjust according to priority
threads = 32
sortMem = "32G"

# Set genome, choose one from hg19, hg38, mm9, mm10
genomeStr = "mm10" # Replace
bowtieGenome = select_ref_bw2(genomeStr)
macsGenome = select_ref_macs(genomeStr)

# Set the output directory for all the samples
outDir = paste0(projectDir, "data/")

# Adapter trimming options
seq_platform = "nextera"  # Nextera Tn5 adapters

# Mapping options
bowtieCmd = "bowtie2"
bowtieOptions = paste("-k 1 -X2000 -t --quiet --mm --threads", threads, "-x", bowtieGenome)

# Tools calls
samCmd = "samtools"
picardCmd = "picard MarkDuplicates "
macsCMD = "macs2"

bamSortExt = ".sort"

# Chromosomes to exclude in sample bigwig and FRiP score
exChr = c("chrM|chrY")

# Transcript annotations for relevant genome
TxDB = select_ref_txdb(genomeStr)
txs = transcripts(TxDB)

# Helper functions
is_single_end = function(fastq) { fq = get0(fastq, envir = .GlobalEnv); return(!is.null(fq) && !is.na(fq) && fq != "") }
is_paired_end = function(R1, R2) { return(is_single_end(R1) && is_single_end(R2)) }

# Initialize read counts by chromosome variable
chrCounts = NA

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
  
  # # Set genomes for mapping
  # genomeStr = "mm10" # Replace
  # B6Genome = "/home/boss_lab/Apps/genomes/STAR/mm10.ERCC"
  # SpretGenome = "/home/boss_lab/Apps/genomes/species/spret-b6/SPRET_EiJ.STAR.genome"
  # 
  # # Map to B6 genome
  # print("Maternal B6 STAR mapping")
  # 
  # #paired-end STAR alignment based on manifest columns
  # B6sample = paste0(files$dir[i], files$sample[i], ".B6")
  # genome = B6Genome
  # 
  # if (!file.exists(paste0(files$dir[i], files$fqMate1[i], ".gz"))) {
  #   cmd1 = paste("gzip -k", paste0(files$dir[i], files$fqMate1[i]))
  #   cmd2 = paste("gzip -k", paste0(files$dir[i], files$fqMate2[i]))
  #   system(cmd1)
  #   system(cmd2)
  # }
  # 
  # fq1 = paste0(files$dir[i], files$fqMate1[i], ".gz")
  # fq2 = paste0(files$dir[i], files$fqMate2[i], ".gz")
  # 
  # starB6Cmd = paste("STAR --runMode alignReads --runThreadN 1 --genomeDir", genome, "--genomeLoad NoSharedMemory --readFilesIn", fq1, fq2, "--readFilesCommand zcat --outFileNamePrefix", B6sample, "--outSAMtype BAM SortedByCoordinate --outBAMcompression 6 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.06 --outFilterMatchNmin 22 --outFilterMatchNminOverLread 0.55 --outFilterScoreMinOverLread 0.55 --seedSearchStartLmax 30 --alignIntronMax 1 --alignEndsType Local")
  # 
  # system(starB6Cmd)
  # files$bamB6File[i] = paste0(B6sample, ".bam")
  
  ##############################
  # Map to Spret genome
  print("Step 6")
  print("Paternal Spret STAR mapping")
  
  # #paired-end STAR alignment based on manifest columns
  # Spretsample = paste0(files$dir[i], files$sample[i], ".Spret")
  # genome = SpretGenome
  # fq1 = paste0(files$dir[i], files$fqMate1[i], ".gz")
  # fq2 = paste0(files$dir[i], files$fqMate2[i], ".gz")
  # 
  # starSpretCmd = paste("STAR --runMode alignReads --runThreadN 1 --genomeDir", genome, "--genomeLoad NoSharedMemory --readFilesIn", fq1, fq2, "--readFilesCommand zcat --outFileNamePrefix", Spretsample, "--outSAMtype BAM SortedByCoordinate --outBAMcompression 6 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.06 --outFilterMatchNmin 22 --outFilterMatchNminOverLread 0.55 --outFilterScoreMinOverLread 0.55 --seedSearchStartLmax 30 --alignIntronMax 1 --alignEndsType Local")
  # 
  # system(starSpretCmd)
  # files$bamSpretFile[i] = paste0(Spretsample, ".bam")
  
  print("Step 7")
  print("Running lapels")
  
  # Run this before pylapels
  # It will change \tnM: to \tNM:
  files$bamSpretNMFile[i] = gsub(".out.bam",".NM.out.bam",files$bamSpretFile[i])
  SpretCmd = paste("samtools view -h", paste0(files$dir[i],files$bamSpretFile[i]), "| sed 's/\tnM:/\tNM:/' | samtools view -bS - >", paste0(files$dir[i],files$bamSpretNMFile[i]))
  system(SpretCmd)
  
  files$bamB6NMFile[i] = gsub(".out.bam",".NM.out.bam",files$bamB6File[i])
  B6Cmd = paste("samtools view -h", paste0(files$dir[i],files$bamB6File[i]), "| sed 's/\tnM:/\tNM:/' | samtools view -bS - >", paste0(files$dir[i],files$bamB6NMFile[i]))
  system(B6Cmd)
  
  # Running pylaplels on Spret aligned bams
  files$lapelsAlignedBam[i] = gsub(".out.bam",".lapels.sortedByName.bam",files$bamSpretFile[i])
  pylapelsCmd = paste("pylapels -n -o", paste0(files$dir[i],files$lapelsAlignedBam[i]), "/home/boss_lab/Apps/genomes/species/spret-b6/SPRET_EiJ.genome/SPRETmm10.mod", paste0(files$dir[i],files$bamSpretNMFile[i]))
  system(pylapelsCmd)
  
  # Need to sort spret bam files by co-ordinates again to mark duplicates
  files$lapelsAlignedCoordSortedBam[i] = gsub(".lapels.sortedByName.bam",".lapels.sortedByCoord.bam",files$lapelsAlignedBam[i])
  sortCmd = paste("samtools sort", paste0(files$dir[i],files$lapelsAlignedBam[i]), "-o", paste0(files$dir[i],files$lapelsAlignedCoordSortedBam[i]))
  system(sortCmd)
  
  print("Step 8")
  print("Mark duplicates")
  
  #### For B6 ####
  files$markDupsB6Bam[i] = markDups(paste0(files$dir[i], files$bamB6NMFile[i]), picardCmd, delBam = F)
  bamFile = paste0(files$dir[i], files$markDupsB6Bam[i])
  
  print("Get mapping stats")
  cts = getBamCts(bamFile)
  
  files$B6.unmapped.reads[i] = cts[1]
  files$B6.mapped.reads[i] = cts[2]
  files$B6.unique.reads[i] = cts[3]
  files$B6.paired.reads[i] = cts[4]
  
  #### For Spret ####
  files$markDupsSpretBam[i] = markDups(paste0(files$dir[i], files$lapelsAlignedCoordSortedBam[i]), picardCmd, delBam = F)
  bamFile = paste0(files$dir[i], files$markDupsSpretBam[i])
  
  print("Get mapping stats")
  cts = getBamCts(bamFile)
  
  files$Spret.unmapped.reads[i] = cts[1]
  files$Spret.mapped.reads[i] = cts[2]
  files$Spret.unique.reads[i] = cts[3]
  files$Spret.paired.reads[i] = cts[4]
  
  # print("Step 9")
  # print("Running suspenders")
  # 
  # # sort both spret and b6 dupMarked bams by query name before suspenders
  # files$markDupsB6SortedByNameBam[i] = gsub()
  # sortB6cmd = paste("samtools sort -n", files$markDupsB6Bam[i], "-o", files$markDupsB6SortedByNameBam[i])
  
  # print("Make BigWig")
  # files$bwFile[i] = paste0(files$dir[i], files$sample[i], ".rpm.bw")
  # bamToBigWig(bamFile, bwFile = files$bwFile[i], flag = scanBamFlag(isDuplicate = F), removeChrs = exChr, sigDigits = 2)
  # 
  # print("Call peaks with MACS2")
  # macsOutDir = paste0(files$dir[i], files$sample[i], ".MACS2.Peaks/")
  # files$PeakFile[i] = paste0(macsOutDir, files$sample[i], "_peaks.narrowPeak")
  # callMACS = paste(macsCMD, "callpeak -t", bamFile, "-f BAM -g", macsGenome,
  #                  "-n", files$sample[i], "--outdir", macsOutDir)
  # system(callMACS)
  # 
  # print("Get peaks count")
  # if (file.exists(paste0(macsOutDir, files$sample[i], "_summits.bed"))) {
  #   peaks = read.table(paste0(macsOutDir, files$sample[i], "_summits.bed"))
  #   files$MACS.peaks[i] = dim(peaks)[1]
  # } else {
  #   files$MACS.peaks[i] = 0
  # }
  # 
  # print("Calculate sample FRiP score") # Fraction of Reads In Peaks (FRiP) = sum(peakReads) / uniqueReads
  # # Read in peak file
  # bed = import.bed(files$PeakFile[i],
  #                  genome=genomeStr,
  #                  extraCols = c(signalValue = "numeric", # extraCols reads in the narrowPeak format
  #                                pValue = "numeric",
  #                                qValue = "numeric",
  #                                peak = "integer"))
  # # Remove peaks from certain chromosomes
  # bed = bed[!grepl(exChr, seqnames(bed))]
  # # Retrieve the bam's desired reads given the parameters
  # reads = readGAlignmentPairs(bamFile, param = ScanBamParam(flag = scanBamFlag(isDuplicate = F)))
  # # With "Union", all reads that overlap a peak by any amount are counted towards that peak's read count
  # so = summarizeOverlaps(bed, reads, mode = "Union", singleEnd = F, ignore.strand = T)
  # # Calculate sample frip score (sFRiP) and round to 3 significant digits
  # files$sFrip[i] = round(sum(assays(so)$counts) / files$unique.reads[i], digits = 3)
  # 
  # # TODO: Debug--why does this take ~5 hours and give na's back?
  # #print("Calculate TSS Enrichment")
  # #files$tsse[i] = round(calcTSSe(bamFile, txs), 2)
  # 
  # print("Remove any remaining temp files")
  # system(paste0("rm ", files$dir[i], "*.sam"))
  # system(paste0("rm ", files$dir[i], "*trimmed*"))
  # 
  # print("Count reads for each chromosome")
  # counts = readsByChr(bamFile)
  # chrCounts = rbind(chrCounts, counts[2,])
  # row.names(chrCounts)[nrow(chrCounts)] = files$sample[i]
  # 
  # # Calculate the time taken to run one sample
  # files$timeTaken[i] = difftime(Sys.time(), time_start, units = "hours")
  # print(paste("Time taken for sample", files$sample[i],  "(hrs):", files$timeTaken[i]))
  # 
  # Update manifest file and note that this sample finished mapping by setting 'include' to FALSE
  files$include[i] = FALSE
  write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)
}

# # Write read counts by chr to table and chrM percentage to manifest
# chrCounts = chrCounts[2:dim(chrCounts)[1], ]
# write.table(chrCounts, "reads.by.chrs.txt", quote = F, sep = "\t")
# 
# chrm = chrCounts[, "chrM"]/rowSums(chrCounts)
# chrm[is.na(chrm)] = 0 # Occurs if a sample was skipped in mapping, however, NA is placed in first position...
# files$chrM_percentage = chrm
# 
# # Note that all samples have been sequenced by resetting all includes to TRUE
files$include = TRUE
write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)
# 
# #######################################
# # Plot Mapping and Peak QC
# 
# plot_mapping_stats(files)
# 
# # Encode recommends a score greater than 20%
# # We use a Rule-of-Thumb of under 5% to 'fail' a sample
# plot_FRiP(files)
# 
# plot_peak_count(files)
# 
# # TODO: Revive after debugging TSSe
# # https://www.encodeproject.org/data-standards/terms/#enrichment
# # https://www.encodeproject.org/atac-seq/
# # hg38 metrics: < 5 concerning; 5-7 Acceptable; > 7 Ideal
# # mm10 metrics: < 10 concerning; 10-15 Acceptable; > 15 Ideal
# #plot_TSSe(files, genomeStr)

############################
# Software versions

# Output R session
capture.output(sessionInfo(), file = paste0("Rsession.Info.", gsub("\\D", "", Sys.time()), ".txt"))

picardTool = "picard"
picardArgs = "MarkDuplicates --version"

tools = c("fastqc",
          "skewer",
          "bowtie2",
          "samtools",
          picardTool,
          "macs2")

args = c("--version", # fastqc
         "-version",  # skewer
         "--version", # bowtie2
         "--version", # samtools
         picardArgs,  # Picard
         "--version"  # MACS2
)

names = c("fastqc",
          "skewer",
          "bowtie2",
          "samtools",
          "MarkDuplicates",
          "MACS2")

version = lapply(1:length(tools), FUN = function (x) {system2(tools[x], args[x], stdout=TRUE, stderr=TRUE)})
names(version) = names

capture.output(version, file = paste0("Software.versions.", gsub("\\D", "", Sys.time()), ".txt"))

