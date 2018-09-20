# DNA mapping and refining alignments based on VaDiR

rule R01_fastq_to_ubam:
    #Picard FastqToSam - Converting to uBAM
    # recomended java requirements
    # FASTQ2 = {input.fastqR2}
    input:
        fastqR1= expand("{input_dir}/fastq/{sample}_R1.fastq.gz",
                         sample = config['sample_name'],
                         input_dir = config['input_dir']),
        fastqR2 = expand("{input_dir}/fastq/{sample}_R2.fastq.gz",
                         sample = config['sample_name'],
                         input_dir = config['input_dir'])
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir']
    output:
        output_ubam = expand("{output_dir}/01_ubam/{sample}.ubam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R01_picard_fastq_to_ubam_rule.log",
        output_dir = config['output_dir'])
    shell:
        """
         /opt/miniconda3/bin/picard FastqToSam \
         FASTQ={input.fastqR1} \
         O={output.output_ubam} \
         SAMPLE_NAME={params.sample_name} \
         SORT_ORDER=queryname \
         TMP_DIR={params.output_dir}/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R02_mark_adapters:
    # Recomended java settings: -d64 -Xmx16G
    input:
        ubam = expand("{output_dir}/01_ubam/{sample}.ubam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir']
    output:
        sam = expand("{output_dir}/02_mark_adapters/{sample}.unaligned.sam",
        sample = config['sample_name'],
        output_dir = config['output_dir']),
        adapter_metrics = expand(
            "{output_dir}/02_mark_adapters/{sample}.adapter_metrics.txt",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R02_picard_mark_illumina_adapters.log",
        output_dir = config['output_dir'])
    shell:
        """
         /opt/miniconda3/bin/picard MarkIlluminaAdapters \
         I={input.ubam} \
         O={output.sam} \
         M={output.adapter_metrics} \
         TMP_DIR={params.output_dir}/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R03_sam_to_interleaved_fastq:
    #  Recomended java parameters -d64 -Xmx16G
    input:
        unaligned_sam = expand(
        "{output_dir}/02_mark_adapters/{sample}.unaligned.sam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir']
    output:
        interleaved_fastq = expand(
        "{output_dir}/03_interleaved_fastq/{sample}.interleaved.fastq",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R03_picard_unaligned_sam_to_fastq.log",
        output_dir = config['output_dir'])
    shell:
        """
         /opt/miniconda3/bin/picard SamToFastq \
         I={input.unaligned_sam} \
         FASTQ={output.interleaved_fastq} \
         CLIPPING_ATTRIBUTE=XT \
         CLIPPING_ACTION=2 \
         INTERLEAVE=true \
         NON_PF=true \
         TMP_DIR={params.output_dir}/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R04_DNA_bwa_map:
    #  "$bwa" mem -M -t $threads -p "$ref_genome_fasta" \
    input:
        reference = expand("{input_dir}/refs/{reference}",
                           reference = config['reference_DNA'],
                           input_dir = config['input_dir']),
        interleaved_fastq = expand(
        "{output_dir}/03_interleaved_fastq/{sample}.interleaved.fastq",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    output:
        mapped_sam = expand("{output_dir}/04_dna_mapped_reads/{sample}.sam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R04_bwa_rule.log",
        output_dir = config['output_dir'])
    shell:
        """
        /opt/miniconda3/bin/bwa mem -M -t 8 {input.reference} \
        {input.interleaved_fastq} \
        > {output.mapped_sam} 2> {log} || (e=$?; cat {log}; exit $e)
        """
rule R05_merged_mapped_and_unampped_dna_reads:
    #Merge alignment data from a SAM or BAM with data in an unmapped BAM file.
    #This tool produces a new SAM or BAM
    #file that includes all aligned and unaligned reads and also carries forward
    #additional read attributes from the
    #unmapped BAM (attributes that are otherwise lost in the process of
    #alignment).
    #The purpose of this tool is to use information from the unmapped BAM to
    #fix up aligner output.
    # Recomended java paremeters -d64 -Xmx16G
    input:
        mapped_sam = expand("{output_dir}/04_dna_mapped_reads/{sample}.sam",
        sample = config['sample_name'],
        output_dir = config['output_dir']),
        ubam = expand("{output_dir}/01_ubam/{sample}.ubam",
                         sample = config['sample_name'],
                         output_dir = config['output_dir']),
        reference = expand("{input_dir}/refs/{reference}",
                           reference = config['reference_DNA'],
                           input_dir = config['input_dir'])
    params:
        output_dir = config['output_dir']
    output:
        merged_mapped_and_unampped_dna_bam = expand(
        "{output_dir}/05_merged_mapped_and_unampped_dna_reads/{sample}.merged.bam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R05_merged_mapped_and_unampped_dna_reads.log",
        output_dir = config['output_dir'])
    shell:
        """
         /opt/miniconda3/bin/picard MergeBamAlignment \
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
         TMP_DIR={params.output_dir}/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R06_add_read_groups_to_bam:
    #  "$JAVA" -d64 -Xmx16g -jar "$picard_jar" AddOrReplaceReadGroups
    input:
        bam = expand(
        "{output_dir}/05_merged_mapped_and_unampped_dna_reads/{sample}.merged.bam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir']
    output:
        bam_with_read_groups = expand(
        "{output_dir}/06_add_read_groups_to_bam/{sample}.bam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R06_add_read_groups_to_bam.log",
        output_dir = config['output_dir'])
    shell:
        """
        /opt/miniconda3/bin/picard AddOrReplaceReadGroups \
        INPUT={input.bam} \
        OUTPUT={output.bam_with_read_groups} \
        SORT_ORDER=coordinate \
        RGID="params.sample_name" \
        RGLB=stranded \
        RGPL=illumina \
        RGPU="barcode" \
        RGSM={params.sample_name} \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR={params.output_dir}/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R07_mark_duplicates:
    #  "$JAVA" -d64 -Xmx16g -jar "$picard_jar" MarkDuplicates
    input:
        bam = expand("{output_dir}/06_add_read_groups_to_bam/{sample}.bam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    params:
        output_dir = config['output_dir']
    output:
        deduped_bam = expand(
        "{output_dir}/07_mark_duplicates/{sample}.deduped.bam",
        sample = config['sample_name'],
        output_dir = config['output_dir']),
        mark_duplicate_metrics = expand(
        "{output_dir}/07_mark_duplicates/{sample}.dedup_metrics.txt",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R07_mark_duplicates.log",
        output_dir = config['output_dir'])
    shell:
        """
        /opt/miniconda3/bin/picard MarkDuplicates \
        INPUT={input.bam} \
        OUTPUT={output.deduped_bam} \
        CREATE_INDEX=true \
        METRICS_FILE={output.mark_duplicate_metrics} \
        TMP_DIR={params.output_dir}/temp_dir \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R08_base_recalibration_report:
    # "$JAVA" -d64 -Xmx8g -jar
    input:
        bam = expand("{output_dir}/07_mark_duplicates/{sample}.deduped.bam",
        sample = config['sample_name'],
        output_dir = config['output_dir']),
        known_1000g_indel_vcf = expand("{input_dir}/refs/{indel_vcf_1000g}",
                           indel_vcf_1000g = config['indel_vcf_1000g'],
                           input_dir = config['input_dir']),
        known_dbsnp_snp_vcf = expand("{input_dir}/refs/{snp_vcf_dbsnp}",
                           snp_vcf_dbsnp = config['snp_vcf_dbsnp'],
                           input_dir = config['input_dir']),
        reference = expand("{input_dir}/refs/{reference}",
                           reference = config['reference_DNA'],
                           input_dir = config['input_dir'])
    output:
        recalibration_report = expand(
        "{output_dir}/08_base_recalibration_report/{sample}_recalibration_report.grp",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    params:
        threads = 8
    log:
        expand("{output_dir}/logs/R08_base_recalibration_report.log",
        output_dir = config['output_dir'])
    shell:
        """
        /opt/gatk-4.0.8.1/gatk BaseRecalibrator \
        --reference {input.reference} \
        --input {input.bam} \
        --known-sites {input.known_1000g_indel_vcf} \
        --known-sites {input.known_dbsnp_snp_vcf} \
        --output {output.recalibration_report} \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R09_analyze_base_quality_covariates:
    # log "GATK AnalyzeCovariates - Generating comparison of base score
    # accuracy before & after recalibration"
    # recalibration_plots="$output/$sample_dna_id.recal_plots.pdf"
    # "$JAVA" -d64 -Xmx8g -jar "$gatk_jar"
    #     -before "$recalibration_report" \
    #     -after "$post_recalibration_report" \
    #     -plots "$recalibration_plots" \
    #     >> $output/log.txt 2>&1 \
    #     || true  # Don't abort if this non-critical step fails
    input:
        recalibration_report = expand(
        "{output_dir}/08_base_recalibration_report/{sample}_recalibration_report.grp",
        sample = config['sample_name'],
        output_dir = config['output_dir']),
        reference = expand("{input_dir}/refs/{reference}",
                    reference = config['reference_DNA'],
                    input_dir = config['input_dir'],
                    output_dir = config['output_dir'])
    output:
        recalibration_plots_pdf = expand(
        "{output_dir}/09_post_recalibration_report/{sample}_recalibration_report.pdf",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R09_post_recalibration_report.log",
                output_dir = config['output_dir'])
    shell:
        """
        /opt/gatk-4.0.8.1/gatk AnalyzeCovariates \
        --before {input.recalibration_report} \
        --plots {output.recalibration_plots_pdf} \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R10_print_recalibrated_bam:
    # "$JAVA" -d64 -Xmx8g -jar "$gatk_jar" -T PrintReads \
    #     -R "$ref_genome_fasta" \
    #     -I "$deduped_bam" \
    #     -BQSR "$recalibration_report" \
    #     -o "$recalibrated_bam" \
    #     -nct $threads \
    #     >> $output/log.txt 2>&1
    input:
        deduped_bam = expand(
            "{output_dir}/07_mark_duplicates/{sample}.deduped.bam",
            sample = config['sample_name'],
            output_dir = config['output_dir']),
        reference = expand("{input_dir}/refs/{reference}",
            reference = config['reference_DNA'],
            input_dir = config['input_dir'],
            output_dir = config['output_dir']),
        recalibration_report = expand(
            "{output_dir}/08_base_recalibration_report/{sample}_recalibration_report.grp",
            sample = config['sample_name'],
            output_dir = config['output_dir'])
    output:
        recalibrated_bam= expand(
            "{output_dir}/10_recalibrated_bam/{sample}_recalibrated.bam",
            sample = config['sample_name'],
            output_dir = config['output_dir'])
    log:
        expand("{output_dir}/logs/R10_post_recalibration_report.log",
        output_dir = config['output_dir'])
    shell:
        """
        /opt/gatk-4.0.8.1/gatk ApplyBQSR \
        --reference {input.reference} \
        --input {input.deduped_bam} \
        --bqsr-recal-file {input.recalibration_report} \
        --output {output.recalibrated_bam} \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
