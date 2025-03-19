
# from paired_dedup_fragment_3852_sample1.bed to end motif profile

#xue gang

rm(list=ls())
workDir = "/mnt/raid7/xuegang/mission/end_motif/"
setwd(workDir)

library(tidyverse)
library(data.table)
library(future.apply)

options(future.globals.maxSize = Inf)


end_motif<-function(File,ref,kmer,outPath,runSample){
  
  fileName <- str_split(fileName, pattern = "3852_", simplify = T)[2]
  sampleName=str_split(fileName,pattern = "\\.",simplify = T)[1]
  
  cmd.sample=paste("echo",File,">>",runSample)
  system(cmd.sample)
  
  rawBed<-fread(File,header = F,sep = "\t",nThread=1)
  
  #filter
  filterBed<- rawBed[rawBed$V1 %in% paste0("chr", 1:22),] %>%
    dplyr::rename(chr=V1,start=V2,end=V3) %>%
    dplyr::mutate(fragment_len=end-start) %>%
    dplyr::filter(fragment_len<=1000) %>%
    dplyr::filter(fragment_len>=20) %>%
    dplyr::select(chr,start,end)
  
  left_motif<-filterBed %>%
    dplyr::mutate(start_kmer=start) %>%
    dplyr::mutate(end_kmer=start+kmer) %>%
    unite("read_region",chr,start,end,sep = "_",remove = F) %>%
    dplyr::mutate(score=rep(1,nrow(.))) %>%
    dplyr::mutate(strand=rep("+",nrow(.))) %>%
    dplyr::select(chr,start_kmer,end_kmer,read_region,score,strand) %>%
    dplyr::mutate_if(is.numeric, as.integer)
  
  right_motif<-filterBed %>%
    dplyr::mutate(start_kmer=end-kmer) %>%
    dplyr::mutate(end_kmer=end) %>%
    unite("read_region",chr,start,end,sep = "_",remove = F) %>%
    dplyr::mutate(score=rep(1,nrow(.))) %>%
    dplyr::mutate(strand=rep("-",nrow(.))) %>%
    dplyr::select(chr,start_kmer,end_kmer,read_region,score,strand) %>%
    dplyr::mutate_if(is.numeric, as.integer)
  
  
  filterBed_kmer<-bind_rows(left_motif,right_motif)
  
  rm(filterBed)
  gc()
  
  cmd_0<-paste("if [ ! -d",paste0(outPath,sampleName),"];then mkdir -p",paste0(outPath,sampleName),";fi")
  system(cmd_0)
  
  data.table::fwrite(filterBed_kmer,file = paste0(outPath,sampleName,"/filterBed_",kmer,"mer_",fileName),
                     sep = "\t",row.names = FALSE, col.names = FALSE,nThread=1)
  
  cmd_1=paste("cat",paste0(outPath,sampleName,"/filterBed_",kmer,"mer_",fileName),"| bedtools getfasta -fi",ref,"-bed stdin -s -name -tab >",
              paste0(outPath,sampleName,"/filterBed_",kmer,"mer_motif_",fileName))
  system(cmd_1)
  
  
  #The frequency of occurrence of each motif was calculated and normalized by the total number of ends among the sequenced fragments. 
  filterBed_kmer_motif_per<-fread(paste0(outPath,sampleName,"/filterBed_",kmer,"mer_motif_",fileName),header=F,nThread=1) %>% 
    dplyr::rename(position=V1,motif=V2) %>%
    dplyr::group_by(motif) %>% dplyr::summarise(motif_count = n(),.groups = "drop") %>%
    dplyr::mutate(motif_per=motif_count/nrow(filterBed_kmer))
  
  data.table::fwrite(filterBed_kmer_motif_per,file = paste0(outPath,sampleName,"/filterBed_",kmer,"mer_motif_per_",fileName),
                     sep = "\t",row.names = FALSE, col.names = T,nThread=1)
  
  rm(filterBed_kmer)
  gc()
  
  return(File)
  
  
}




############################################################################################################################


system("if [ ! -d /mnt/raid7/xuegang/mission/end_motif/run_code_result/end_kmer/mission1/5mC ];then mkdir -p /mnt/raid7/xuegang/mission/end_motif/run_code_result/end_kmer/mission1/5mC;fi")
system("if [ ! -d /mnt/raid7/xuegang/mission/end_motif/run_code_result/end_kmer/mission1/5hmC ];then mkdir -p /mnt/raid7/xuegang/mission/end_motif/run_code_result/end_kmer/mission1/5hmC;fi")


############################################################################################################################

plan(multisession,workers = 50)

# You should prepare the sample_bedPath as follows:
#eg : /mnt/raid8/bwa_mapping/cfDNA20240814/medip/mapping_result/sample1/paired_dedup_fragment_3852_sample1.bed

#4-mer
end_motif_5mC_result_4mer<-future_lapply(sampleInfo_5mC$sample_bedPath,
                                         end_motif,
                                         ref="/mnt/mem/reference/reference_bwa/GRCh37.p13.genome_spike-in_20181228.fa",
                                         kmer=4,
                                         outPath="/mnt/raid7/xuegang/mission/end_motif/run_code_result/end_kmer/mission1/5mC/",
                                         runSample="/mnt/raid7/xuegang/mission/end_motif/run_code_result/mission1_5mC_sample_4mer.txt",
                                         future.seed = 123)

print("Done for identifying 4mer-end motif of mission1 5mC samples !")


end_motif_5hmC_result_4mer<-future_lapply(sampleInfo_5hmC$sample_bedPath,
                                          end_motif,
                                          ref="/mnt/mem/reference/reference_bwa/GRCh37.p13.genome_spike-in_20181228.fa",
                                          kmer=4,
                                          outPath="/mnt/raid7/xuegang/mission/end_motif/run_code_result/end_kmer/mission1/5hmC/",
                                          runSample="/mnt/raid7/xuegang/mission/end_motif/run_code_result/mission1_5hmC_sample_4mer.txt",
                                          future.seed = 123)

print("Done for identifying 4mer-end motif of mission1 5hmC samples!")


#6-mer
end_motif_5mC_result_6mer<-future_lapply(sampleInfo_5mC$sample_bedPath,
                                         end_motif,
                                         ref="/mnt/mem/reference/reference_bwa/GRCh37.p13.genome_spike-in_20181228.fa",
                                         kmer=6,
                                         outPath="/mnt/raid7/xuegang/mission/end_motif/run_code_result/end_kmer/mission1/5mC/",
                                         runSample="/mnt/raid7/xuegang/mission/end_motif/run_code_result/mission1_5mC_sample_6mer.txt",
                                         future.seed = 123)

print("Done for identifying 6mer-end motif of mission1 5mC samples !")


end_motif_5hmC_result_6mer<-future_lapply(sampleInfo_5hmC$sample_bedPath,
                                          end_motif,
                                          ref="/mnt/mem/reference/reference_bwa/GRCh37.p13.genome_spike-in_20181228.fa",
                                          kmer=6,
                                          outPath="/mnt/raid7/xuegang/mission/end_motif/run_code_result/end_kmer/mission1/5hmC/",
                                          runSample="/mnt/raid7/xuegang/mission/end_motif/run_code_result/mission1_5hmC_sample_6mer.txt",
                                          future.seed = 123)

print("Done for identifying 6mer-end motif of mission1 5hmC samples!")

