rule star_map:
    input:
        ref_dir = "VaDiR/refs/ucsc.hg19.fasta",
        fastqR1 = expand("{fastq_dir}/{sample}_R1.fastq.gz", 
                         sample = config['sample_name'], 
                         fastq_dir = config['fastq_dir']),
        fastqR2 = expand("{fastq_dir}/{sample}_R2.fastq.gz", 
                         sample = config['sample_name'],
                         fastq_dir = config['fastq_dir'])
    output:
        expand("mapped_reads/{sample}.sam", sample = config['sample_name'])
    log:
        "logs/bwa_rule.log"
    shell:
        """
        "$star" \
    --genomeDir "$star_index_dir" \
    --readFilesIn "$read_1" "$read_2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$pass2_dir/$sample_id." \
    --outSAMattributes NH HI AS nM NM MD \
    --outSAMtype BAM Unsorted \
    --twopassMode Basic \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --runThreadN $threads \
    > $output_dir/log.txt 2>&1 > {output} 2> {log} || (e=$?; cat {log}; exit $e)
        """
        
        
        
output_dir="$3"
sample_dir="$1"
sample_id="$2"
threads="$5"

read_1="$sample_dir/$sample_id$r1_fastq_file_suffix"
read_2="$sample_dir/$sample_id$r2_fastq_file_suffix"
log "Input read 1:  $read_1"
log "Input read 2:  $read_2"

pass2_dir="$output_dir/star-pass2"
log "Output directory:  $pass2_dir"
mkdir -p "$pass2_dir"


