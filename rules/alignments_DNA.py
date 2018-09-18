# VaDiR pipeline
# includes BWA for mapping and additional mapping "refinements"
# final output is recalibrated and aligned bam
# also Marking adapters, mapping, sorting, adding read groups
# log "Picard MarkDuplicates - Marking duplicates & creating index"
# dedupe_metrics="$output/$sample_dna_id.dedupe_metrics.txt"
# deduped_bam="$output/$sample_dna_id.deduped.bam"
# deduped_bai="$output/$sample_dna_id.deduped.bai"
# log "GATK BaseRecalibrator - Recalibrating base scores & generating recalibration report"
# recalibration_report="$output/$sample_dna_id.recal_report.grp"
# log "GATK PrintReads - Writing recalibrated scores to BAM"
# recalibrated_bam="$output/$sample_dna_id.refined.bam"

rule fastq_to_ubam:
    #Picard FastqToSam - Converting to uBAM
    # recomended java requirements
    # FASTQ2 = {input.fastqR2}
    input:    
        fastqR1= expand("input/fastq/{sample}_R1.fastq.gz", 
                         sample = config['sample_name']),
        fastqR2 = expand("input/fastq/{sample}_R2.fastq.gz", 
                         sample = config['sample_name'])
    params:
        sample_name = config['sample_name']
    output:
        output_ubam = expand("output/01_ubam/{sample}.ubam", sample = config['sample_name'])
    log:
        "output/logs/picard_fastq_to_ubam_rule.log"
    shell:
        """
         picard FastqToSam \
         FASTQ={input.fastqR1} \
         O={output.output_ubam} \
         SAMPLE_NAME={params.sample_name} \
         SORT_ORDER=queryname \
         TMP_DIR=output/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """

rule mark_adapters:
    # Recomended java settings: -d64 -Xmx16G
    input:    
        ubam= expand("output/01_ubam/{sample}.ubam", 
                         sample = config['sample_name']) 
    params:
        sample_name = config['sample_name']
    output:
        sam = expand("output/02_mark_adapters/{sample}.unaligned.sam", sample = config['sample_name']),
        adapter_metrics = expand("output/02_mark_adapters/{sample}.adapter_metrics.txt", sample = config['sample_name'])
    log:
        "output/logs/picard_mark_illumina_adapters.log"
    shell:
        """
         picard MarkIlluminaAdapters \
         I={input.ubam} \
         O={output.sam} \
         M={output.adapter_metrics} \
         TMP_DIR=output/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """

rule sam_to_interleaved_fastq:
    #  Recomended java parameters -d64 -Xmx16G
    input:    
        unaligned_sam = expand("output/02_mark_adapters/{sample}.unaligned.sam",  
                         sample = config['sample_name']) 
    params:
        sample_name = config['sample_name']
    output:
        interleaved_fastq = expand("output/03_interleaved_fastq/{sample}.interleaved.fastq", sample = config['sample_name'])
    log:
        "output/logs/picard_unaligned_sam_to_fastq.log"
    shell:
        """
         picard SamToFastq \
         I={input.unaligned_sam} \
         FASTQ={output.interleaved_fastq} \
         CLIPPING_ATTRIBUTE=XT \
         CLIPPING_ACTION=2 \
         INTERLEAVE=true \
         NON_PF=true \
         TMP_DIR=output/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """

rule DNA_bwa_map:
    #  "$bwa" mem -M -t $threads -p "$ref_genome_fasta" \
    input:
        reference = expand("input/refs/{reference}", 
                           reference = config['reference_DNA']), 
        interleaved_fastq = expand("output/03_interleaved_fastq/{sample}.interleaved.fastq", 
                         sample = config['sample_name']) 
    output:
        mapped_sam = expand("output/04_dna_mapped_reads/{sample}.sam", sample = config['sample_name'])
    log:
        "output/logs/bwa_rule.log"
    shell:
        """
        bwa mem -M -t 8 {input.reference} \
        {input.interleaved_fastq} \
        > {output.mapped_sam} 2> {log} || (e=$?; cat {log}; exit $e)
        """

rule merged_mapped_and_unampped_dna_reads:
    #Merge alignment data from a SAM or BAM with data in an unmapped BAM file. This tool produces a new SAM or BAM
    #file that includes all aligned and unaligned reads and also carries forward additional read attributes from the 
    #unmapped BAM (attributes that are otherwise lost in the process of alignment). 
    #The purpose of this tool is to use information from the unmapped BAM to fix up aligner output. 
    # Recomended java paremeters -d64 -Xmx16G   
    input:    
        mapped_sam = expand("output/04_dna_mapped_reads/{sample}.sam",  
                         sample = config['sample_name']),
        ubam = expand("output/01_ubam/{sample}.ubam",  
                         sample = config['sample_name']),
        reference = expand("input/refs/{reference}", 
                           reference = config['reference_DNA'])
    output:
        merged_mapped_and_unampped_dna_bam = expand("output/05_merged_mapped_and_unampped_dna_reads/{sample}.merged.bam", sample = config['sample_name'])
    log:
        "output/logs/merged_mapped_and_unampped_dna_reads.log"
    shell:
        """
         picard MergeBamAlignment \
         ALIGNED_BAM={input.mapped_sam} \
         UNMAPPED_BAM={input.ubam} \
         OUTPUT={output.merged_mapped_and_unampped_dna_bam} \
         R={input.reference} \
         CREATE_INDEX=true \
         ADD_MATE_CIGAR=true \
         CLIP_ADAPTERS=false \
         CLIP_OVERLAPPING_READS=true \
         INCLUDE_SECONDARY_ALIGNMENTS=true \
         MAX_INSERTIONS_OR_DELETIONS=-1 \
         PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
         ATTRIBUTES_TO_RETAIN=XS
         TMP_DIR=output/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """   

rule add_read_groups_to_bam:
    #  "$JAVA" -d64 -Xmx16g -jar "$picard_jar" AddOrReplaceReadGroups
    input:    
        bam = expand("output/05_merged_mapped_and_unampped_dna_reads/{sample}.merged.bam", 
                                                    sample = config['sample_name'])
    params:
        sample_name = config['sample_name']
    output:
        bam_with_read_groups = expand("output/06_add_read_groups_to_bam/{sample}.bam", sample = config['sample_name'])
    log:
        "output/logs/add_read_groups_to_bam.log"
    shell:
        """
        picard AddOrReplaceReadGroups \
        INPUT={input.bam} \
        OUTPUT={output.bam_with_read_groups} \
        SORT_ORDER=coordinate \
        RGID="params.sample_name" \
        RGLB=stranded \
        RGPL=illumina \
        RGPU="barcode" \
        RGSM={params.sample_name} \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=output/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """   

rule mark_duplicates:
    #  "$JAVA" -d64 -Xmx16g -jar "$picard_jar" MarkDuplicates
    input:    
        bam = expand("output/06_add_read_groups_to_bam/{sample}.bam", sample = config['sample_name'])
    output:
        deduped_bam = expand("output/07_mark_duplicates/{sample}.deduped.bam", sample = config['sample_name']),
        mark_duplicate_metrics = expand("output/07_mark_duplicates/{sample}.dedup_metrics.txt", sample = config['sample_name'])
    log:
        "output/logs/mark_duplicates.log"
    shell:
        """
        picard MarkDuplicates \
        INPUT={input.bam} \
        OUTPUT={output.deduped_bam} \
        CREATE_INDEX=true \
        METRICS_FILE={output.mark_duplicate_metrics} \
        TMP_DIR=output/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """   

rule base_recalibration_report:
    # "$JAVA" -d64 -Xmx8g -jar 
    input:    
        bam = expand("output/07_mark_duplicates/{sample}.deduped.bam", sample = config['sample_name']),
        known_1000g_indel_vcf = expand("input/refs/{indel_vcf_1000g}", 
                           indel_vcf_1000g = config['indel_vcf_1000g']),
        known_dbsnp_snp_vcf = expand("input/refs/{snp_vcf_dbsnp}", 
                           snp_vcf_dbsnp = config['snp_vcf_dbsnp']),
        reference = expand("input/refs/{reference}", 
                           reference = config['reference_DNA'])
    output:
        recalibration_report = expand(
            "output/08_base_recalibration_report/{sample}_recalibration_report.grp", 
            sample = config['sample_name'])
    params:
        threads = 8
    log:
        "output/logs/base_recalibration_report.log"
    shell:
        """
        ../../tools/gatk-4.0.8.1/gatk BaseRecalibrator \
        --reference {input.reference} \
        --input {input.bam} \
        --known-sites {input.known_1000g_indel_vcf} \
        --known-sites {input.known_dbsnp_snp_vcf} \
        --output {output.recalibration_report} \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """   

rule analyze_base_quality_covariates:
    # log "GATK AnalyzeCovariates - Generating comparison of base score accuracy before & after recalibration"
    # recalibration_plots="$output/$sample_dna_id.recal_plots.pdf"
    # "$JAVA" -d64 -Xmx8g -jar "$gatk_jar"
    #     -before "$recalibration_report" \
    #     -after "$post_recalibration_report" \
    #     -plots "$recalibration_plots" \
    #     >> $output/log.txt 2>&1 \
    #     || true  # Don't abort if this non-critical step fails
    input:
        recalibration_report = expand(
            "output/08_base_recalibration_report/{sample}_recalibration_report.grp", 
            sample = config['sample_name']),
        reference = expand("input/refs/{reference}", 
                           reference = config['reference_DNA'])
    output:
        recalibration_plots_pdf = expand("output/09_post_recalibration_report/{sample}_recalibration_report.pdf", 
                                         sample = config['sample_name'])
        
    log:
        "output/logs/post_recalibration_report.log"
    shell:
        """
        ../../tools/gatk-4.0.8.1/gatk AnalyzeCovariates \
        --before {input.recalibration_report} \
        --plots {output.recalibration_plots_pdf} \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """   



rule print_recalibrated_bam:
    # "$JAVA" -d64 -Xmx8g -jar "$gatk_jar" -T PrintReads \
    #     -R "$ref_genome_fasta" \
    #     -I "$deduped_bam" \
    #     -BQSR "$recalibration_report" \
    #     -o "$recalibrated_bam" \
    #     -nct $threads \
    #     >> $output/log.txt 2>&1
    input:
        recalibration_report = expand(
            "output/08_base_recalibration_report/{sample}_recalibration_report.grp", 
            sample = config['sample_name']),
        reference = expand("input/refs/{reference}", 
                           reference = config['reference_DNA'])
    output:
        recalibration_plots_pdf = expand("output/09_post_recalibration_report/{sample}_recalibration_report.pdf", 
                                         sample = config['sample_name'])
        
    log:
        "output/logs/post_recalibration_report.log"
    shell:
        """
        ../../tools/gatk-4.0.8.1/gatk AnalyzeCovariates \
        --before {input.recalibration_report} \
        --plots {output.recalibration_plots_pdf} \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """   



