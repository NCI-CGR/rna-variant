# VaDiR Workflow
# alignment rule

# RNA star basic
# RNA refine mapping
# DNA bwa
# DNA refine mapping
# variant calling rules (use include)

# SNPir
# RVBoost
# Mutect2
# merging and filtering of variants

# Merge VCFs
# compare DNA and RNA (mpileup)
# Endo Annotate (?)
# Filter variants (PERL)

import sys
import os

configfile: "config.yaml"

localrules:all
#localrules means that it is not submitted to cluster but run on the
#current node

include: "rules/alignments_DNA.py"

rule VaDiR_summary_report:
    input:
       recalibrated_bam= expand(
            "{output_dir}/10_recalibrated_bam/{sample}_recalibrated.bam",
            sample = config['sample_name'],
            output_dir = config['output_dir'])




