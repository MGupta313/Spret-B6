# README #

### What is this repository for? ###

* Quick summary: This repo contains all the necessary script to map F1 hybrid samples to both parents (B6 and Spret), call peaks, find peaks that contain Spret SNPs, and count number of reads coming from each parent mapping to each SNP present within an ATAC peak.
* Version: v1.0
* [About lapels and suspenders](https://github.com/holtjma/suspenders/wiki/Lapels-and-Suspenders-Pipeline)

### How do I get set up? ###

#### Files required: ####
* mm10.fa  
* spret.mod  
* SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf  
* SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf  
* STAR indexed mm10_spret genome  

#### Dependencies ####
* star  
* samtools  
* picard marDuplicates  
* macs2  
* pylapels  
* pysuspenders  
* VariantAnnotation  
* Rsamtools  
* GenomicAlignments  
* foreach  
* doParallel  

#### How to run scripts  ####
1. ATACseq.map.R  
2. plotMappingstats.R  
3. ATACseq.coverage.R  
4. makeSpretSNPvcf.sh  
5. getSnpReads.R  

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact