# RNA mapping and refining alignments based on VaDiR
# Adds read groups, sorts, marks duplicates, splits reads that span splice
# junctions, creates index, realigns around known indels & SNPs, reassigns
# mapping qualities, and recalibrates base quality scores.
# This step is done once for each sample
# Tools used Java
#   SAMtools
#   Picard
#   GATK
# INPUT
#   Config file at <current directory>/CONFIG.sh
#   STAR mapping at <output_dir>/<sample dir>/star-pass2/<sample>.Aligned.out.bam
#   Reference genome at <ref_genome_fasta>
#   Known variant lists at
#     <known_mills_indel_vcf>, <known_1000g_indel_vcf>, <known_dbsnp_snp_vcf>
# pass2_dir="$output_dir/star-pass2"
# log "Output directory:  $pass2_dir"
# mkdir -p "$pass2_dir"
rule R99_create_star_reference:
rule R01_star_map:
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

rule R02_mark_adapters:
    # Recomended java settings: -d64 -Xmx16G
    input:
        ubam = expand("{output_dir}/01_ubam/{sample}.ubam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    params:
        sample_name = config['sample_name']
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
         /home/rstudio/miniconda3/bin/picard MarkIlluminaAdapters \
         I={input.ubam} \
         O={output.sam} \
         M={output.adapter_metrics} \
         TMP_DIR=output/temp_dir \
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
        sample_name = config['sample_name']
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
         /home/rstudio/miniconda3/bin/picard SamToFastq \
         I={input.unaligned_sam} \
         FASTQ={output.interleaved_fastq} \
         CLIPPING_ATTRIBUTE=XT \
         CLIPPING_ACTION=2 \
         INTERLEAVE=true \
         NON_PF=true \
         TMP_DIR=output/temp_dir \
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
        /home/rstudio/miniconda3/bin/bwa mem -M -t 8 {input.reference} \
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
         /home/rstudio/miniconda3/bin/picard MergeBamAlignment \
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
rule R06_add_read_groups_to_bam:
    #  "$JAVA" -d64 -Xmx16g -jar "$picard_jar" AddOrReplaceReadGroups
    input:
        bam = expand(
        "{output_dir}/05_merged_mapped_and_unampped_dna_reads/{sample}.merged.bam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
    params:
        sample_name = config['sample_name']
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
        /home/rstudio/miniconda3/bin/picard AddOrReplaceReadGroups \
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
rule R07_mark_duplicates:
    # log "Picard MarkDuplicates - Marking duplicates & creating index"
    # dedupe_metrics="$output_sample_dir/dedupMetrics.txt"
    # deduped_bam="$output_sample_dir/reads_deduped.bam"
    # deduped_bai="$output_sample_dir/reads_deduped.bai"
    #
    # "$JAVA" \
    #     -d64 -Xmx24g -jar "$picard_jar" MarkDuplicates \
    #     INPUT="$grouped_sorted_bam" \
    #     OUTPUT="$deduped_bam" \
    #     CREATE_INDEX=true \
    #     VALIDATION_STRINGENCY=SILENT \
    #     METRICS_FILE="$dedupe_metrics" \
    #     >> $output_sample_dir/log.txt 2>&1
    #
    # rm -f "$grouped_sorted_bam"
    #  "$JAVA" -d64 -Xmx16g -jar "$picard_jar" MarkDuplicates
    input:
        bam = expand("{output_dir}/06_add_read_groups_to_bam/{sample}.bam",
        sample = config['sample_name'],
        output_dir = config['output_dir'])
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
        /home/rstudio/miniconda3/bin/picard MarkDuplicates \
        INPUT={input.bam} \
        OUTPUT={output.deduped_bam} \
        CREATE_INDEX=true \
        METRICS_FILE={output.mark_duplicate_metrics} \
        TMP_DIR=output/temp_dir \
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
        /home/rstudio/tool_bin/gatk-4.0.8.1/gatk BaseRecalibrator \
        --reference {input.reference} \
        --input {input.bam} \
        --known-sites {input.known_1000g_indel_vcf} \
        --known-sites {input.known_dbsnp_snp_vcf} \
        --output {output.recalibration_report} \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """
rule R09_analyze_base_quality_covariates:
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
        /home/rstudio/tool_bin/gatk-4.0.8.1/gatk AnalyzeCovariates \
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
        /home/rstudio/tool_bin/gatk-4.0.8.1/gatk ApplyBQSR \
        --reference {input.reference} \
        --input {input.deduped_bam} \
        --bqsr-recal-file {input.recalibration_report} \
        --output {output.recalibrated_bam} \
         > {log} 2>&1 || (e=$?; cat {log}; exit $e)
        """





rule R03_add_or_replace_read_groups:
log "Picard AddOrReplaceReadGroups - Adding read groups & sorting"
grouped_sorted_bam="$output_sample_dir/reads_grouped_sorted.bam"

"$JAVA" \
    -d64 -Xmx24g -jar "$picard_jar" AddOrReplaceReadGroups \
    INPUT="$output_dir/star-pass2/$sample_id.Aligned.out.bam" \
    OUTPUT="$grouped_sorted_bam" \
    SORT_ORDER=coordinate \
    RGID="$sample_id" \
    RGLB=stranded \
    RGPL=illumina \
    RGPU="barcode" \
    RGSM="$sample_id" \
    VALIDATION_STRINGENCY=SILENT \
    > $output_sample_dir/log.txt 2>&1


rule R99_generate_ref_genome_dict:
#-------------------------------------------------------------------------------
# Generate reference genome dictionary & index if they don't exist
# http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference
#-------------------------------------------------------------------------------

# To match the file names expected by GATK:
# Use shell parameter expansion to remove fasta extension & replace with .dict.
# Simply append .fai for the index (don't replace fasta extension).
# ref_genome_dict=${ref_genome_fasta%.*}.dict
# ref_genome_index="$ref_genome_fasta.fai"
#
# if [ ! -e "$ref_genome_dict" ]; then
#     log "Picard CreateSequenceDictionary - Generating dictionary file for reference genome FASTA"
#     "$JAVA" \
#         -d64 -Xmx24g -jar "$picard_jar" CreateSequenceDictionary \
#         R="$ref_genome_fasta" \
#         O="$ref_genome_dict" \
# 	>> $output_sample_dir/log.txt 2>&1
# fi
# if [ ! -e "$ref_genome_index" ]; then
#     log "SAMtools faidx - Generating index file for reference genome FASTA"
#     "$samtools" faidx "$ref_genome_fasta" > $output_dir/log.txt 2>&1
# fi


rule R04_gatk_split_n_cigar_reads:
#-------------------------------------------------------------------------------
# Refine with GATK
#-------------------------------------------------------------------------------
# log "GATK SplitNCigarReads - Splitting reads into grouped exon segments & hard-clipping sequences that overhang into intronic regions"
# split_bam="$output_sample_dir/reads_split.bam"
# split_bai="$output_sample_dir/reads_split.bai"
#
# "$JAVA" \
#     -d64 -Xmx24g -jar "$gatk_jar" -T SplitNCigarReads \
#     -R "$ref_genome_fasta" \
#     -I "$deduped_bam" \
#     -o "$split_bam" \
#     -U ALLOW_N_CIGAR_READS \
#     -rf ReassignOneMappingQuality \
#     -RMQF 255 \
#     -RMQT 60 \
#     >> $output_sample_dir/log.txt 2>&1
#
# rm -f "$deduped_bam"
# rm -f "$deduped_bai"
#
# log "GATK RealignerTargetCreator - Preparing to realign around known indels/SNPs"
# interval_file="$output_sample_dir/realigner_target.intervals"
#
# "$JAVA" \
#     -d64 -Xmx24g -jar "$gatk_jar" -T RealignerTargetCreator \
#     -I "$split_bam" \
#     -R "$ref_genome_fasta" \
#     -o "$interval_file" \
#     -known "$known_mills_indel_vcf" \
#     -known "$known_1000g_indel_vcf" \
#     -nt 2 \
#     >> $output_sample_dir/log.txt 2>&1

rule R05_gatk_indel_realigner:
# log "GATK IndelRealigner - Realigning & reassigning STAR mapping quality 255 to GATK-compatible 60"
# realigned_bam="$output_sample_dir/reads_realigned.bam"
# realigned_bai="$output_sample_dir/reads_realigned.bai"
#
# mkdir -p "$java_tmp_dir"
#
# "$JAVA" \
#     -d64 -Xmx24g -Djava.io.tmpdir="$java_tmp_dir" -jar "$gatk_jar" -T IndelRealigner \
#     -I "$split_bam" \
#     -R "$ref_genome_fasta" \
#     -targetIntervals "$interval_file" \
#     -o "$realigned_bam" \
#     -known "$known_mills_indel_vcf" \
#     -known "$known_1000g_indel_vcf" \
#     --consensusDeterminationModel KNOWNS_ONLY \
#     -LOD 0.4 \
#     >> $output_sample_dir/log.txt 2>&1
#
# rm -f "$interval_file"
# rm -f "$split_bam"
# rm -f "$split_bai"



rule R06_gatk_base_recalibrator:
log "GATK BaseRecalibrator - Recalibrating base scores & generating recalibration report"
recalibration_report="$output_sample_dir/recalibration_report.grp"

"$JAVA" \
    -d64 -Xmx24g -jar "$gatk_jar" -T BaseRecalibrator \
    -R "$ref_genome_fasta" \
    -I "$realigned_bam" \
    -knownSites "$known_mills_indel_vcf" \
    -knownSites "$known_1000g_indel_vcf" \
    -knownSites "$known_dbsnp_snp_vcf" \
    -o "$recalibration_report" \
    -nct $threads \
    >> $output_sample_dir/log.txt 2>&1

log "GATK BaseRecalibrator - Generating post-recalibration report"
post_recalibration_report="$output_sample_dir/post_recalibration_report.grp"

"$JAVA" \
    -d64 -Xmx24g -jar "$gatk_jar" -T BaseRecalibrator \
    -R "$ref_genome_fasta" \
    -I "$realigned_bam" \
    -knownSites "$known_mills_indel_vcf" \
    -knownSites "$known_1000g_indel_vcf" \
    -knownSites "$known_dbsnp_snp_vcf" \
    -BQSR "$recalibration_report" \
    -o "$post_recalibration_report" \
    -nct $threads \
    >> $output_sample_dir/log.txt 2>&1

log "GATK AnalyzeCovariates - Generating comparison of base score accuracy before & after recalibration"
recalibration_plots="$output_sample_dir/recalibration_plots.pdf"

# Try to load R so we can generate recalibration plots

"$JAVA" \
    -d64 -Xmx24g -jar "$gatk_jar" -T AnalyzeCovariates \
    -R "$ref_genome_fasta" \
    -before "$recalibration_report" \
    -after "$post_recalibration_report" \
    -plots "$recalibration_plots" \
    >> $output_sample_dir/log.txt 2>&1 \
    || true

log "GATK PrintReads - Writing recalibrated scores to BAM"
recalibrated_bam="$output_sample_dir/$sample_id.refined.bam"

"$JAVA" \
   -d64 -Xmx24g -jar "$gatk_jar" -T PrintReads \
   -R "$ref_genome_fasta" \
   -I "$realigned_bam" \
   -BQSR "$recalibration_report" \
   -o "$recalibrated_bam" \
   -nct $threads \
   >> $output_sample_dir/log.txt 2>&1

rm -f "$realigned_bam"
rm -f "$realigned_bai"
rm -f "$recalibration_report"
rm -f "$post_recalibration_report"
