# This script generates read counts for each present with an ATAC peak

library(foreach)
library(doParallel)

group = "f" # Replace

homeDir = "/ATAC/analysis/snp/readCounts/"
setwd(homeDir)
rppm_file <- "/ATAC/coverage/MZB_FoB_2023.11.rppm.txt"
rppm_data <- read.table(rppm_file, header = TRUE, sep = "\t", quote="", comment.char="", check.names=F)

getReads = function(letter){
  v = paste0("SPRET_EiJ.snp.100k.subset_",group,letter) # Replace
  vcf_file <- paste0("/ATAC/analysis/snp/vcf/", v, ".vcf")
  print(vcf_file)  
  vcf <- VariantAnnotation::readVcf(vcf_file)
  
  # get which peaks have snps
  start_time <- Sys.time()
  present_snp = data.frame()
  
  for (i in 1:length(vcf)){
    snp <- vcf[i,]
    snp_chr <- paste0("chr",as.character(seqnames(snp)))
    peak_df <- subset(rppm_data, chr == snp_chr)
    selected_peak = which(peak_df$start <= end(snp) & peak_df$end >= start(snp))
    peak = peak_df[selected_peak, ]
    if(nrow(peak) > 0){
      peak$include = TRUE
      peak$snp_pos = start(snp)
      peak$ref = as.character(ref(snp)[[1]]) #snp@fixed@listData[["REF"]][[1]]
      if (length(as.character(alt(snp)[[1]])) > 1){
        peak$alt = paste(as.character(alt(snp)[[1]]), collapse = ",")
      }
      else{
        peak$alt = as.character(alt(snp)[[1]]) #snp@fixed@listData[["ALT"]]
      }
      peak$rsID = snp@rowRanges@ranges@NAMES
      present_snp = rbind(present_snp, peak)
    }  
  }
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(paste("loop1:",elapsed_time))
  
  snpPos = present_snp$snp_pos
  for (j in 1:8){
    sample = paste0("SCA00",j)
    bamFilename = paste0("/ATAC/data/SCA00",j, "/SCA00",j, ".sortedByCoordMerged.bam")
    bam <- BamFile(bamFilename)
    
    for (n in snpPos){
      snpRanges <- GRanges(seqnames = "chr1", ranges = IRanges(start = n, end = n))
      for (i in 1:3){
        sbp <- ScanBamParam(
          flag = scanBamFlag(isDuplicate = F),
          which = snpRanges,
          tag = "po",
          tagFilter = list(po = c(i))
        )
        reads_at_position <- readGAlignmentPairs(bam, param = sbp)
        num_reads <- length(reads_at_position)
        stout = paste(sample, "- For snp at", n, "po:i;", i, "# of reads = ", num_reads)
        print(stout)
        present_snp[present_snp$snp_pos==n, paste0(sample,"_po", i, "_reads")] = num_reads
      }
    }
  }
  
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(paste("loop2",elapsed_time))
  
  write.table(present_snp, file = paste0(homeDir,group,"/",v,".readCounts.txt"), sep = "\t", quote = F, row.names = F)
  print(paste("done for", v))
}

cores <- 4  # Number of cores to use
cl <- makeCluster(cores)
registerDoParallel(cl)

results = foreach(letter = letters[1:12], .combine = c) %dopar% {
  library(VariantAnnotation)
  library(Rsamtools)
  library(GenomicAlignments)
  getReads(letter)
}

