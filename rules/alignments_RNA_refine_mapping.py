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
        star mem {input.ref_dir} {input.fastqR1} {input.fastqR2} > {output} 2> {log} || (e=$?; cat {log}; exit $e)
        """