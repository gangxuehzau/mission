#!/bin/bash

# your sample list
samples=("sample1" "sample2" "sample3")  


base_path="/mnt/raid8/bwa_mapping/cfDNA20240814/medip"
ref_path="/mnt/mem/huwen/script/Mapping/Database/GRCh37.p13.genome_spike-in_20181228.fa"
script_path="/mnt/mem/huwen/script/Mapping/tools"


for sample in "${samples[@]}"; do
    
    mkdir -p "$base_path/mapping_result/$sample"
    cd "$base_path/mapping_result/$sample" 

    $script_path/bwa mem -t 8 "$ref_path" \
        "$base_path/filter_data/$sample/${sample}_1.fq.gz" \
        "$base_path/filter_data/$sample/${sample}_2.fq.gz" > "${sample}.sam"


    $script_path/samtools view -@ 7 -b -S "${sample}.sam" > "${sample}.bam"
    
    $script_path/samtools view -@ 7 -bh -f 2 -F 1548 -q 30 "${sample}.bam" > "filtered_${sample}.bam"

 
    $script_path/java -jar $script_path/picard.jar SortSam \
        INPUT="filtered_${sample}.bam" \
        OUTPUT="sorted_filtered_${sample}.bam" \
        SORT_ORDER=coordinate

    
    $script_path/java -jar $script_path/picard.jar MarkDuplicates \
        REMOVE_DUPLICATES=true \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 \
        INPUT="sorted_filtered_${sample}.bam" \
        OUTPUT="dedup_sorted_filtered_${sample}.bam" \
        METRICS_FILE="dedup_sorted_filtered_${sample}.metrics"

    
    $script_path/samtools flagstat "${sample}.bam" > "${sample}.bam.flagstat"

    
    $script_path/samtools sort -n "dedup_sorted_filtered_${sample}.bam" -o "dedup_sorted_name_filtered_${sample}.bam"

    $script_path/samtools view -f 3 -F 3852 "dedup_sorted_name_filtered_${sample}.bam" -o "dedup_sorted_name_filtered_3852_${sample}.bam"

    $script_path/bedtools bamtobed -i "dedup_sorted_name_filtered_3852_${sample}.bam" -bedpe -mate1 | grep -v "WARNING" > "dedup_sorted_name_filtered_3852_${sample}.bed"

    awk -F '\t' -v OFS='\t' '{if($1!=$4) next; if($9=="+") {s=$2;e=$6} else {s=$5;e=$3} if (e>s) print $1,s,e,$7,$8,$9}' "dedup_sorted_name_filtered_3852_${sample}.bed" > "paired_dedup_fragment_3852_${sample}.bed"
done
