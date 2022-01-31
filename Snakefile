#!/usr/bin/env python3
import pathlib2

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

###########
# GLOBALS #
###########

#containers
trinity_container = 'docker://trinityrnaseq/trinityrnaseq:2.11.0'

#########
# RULES #
#########

rule target:
    input:
        #expand('output/variants_{analysis}/filtered_output.vcf', analysis=("exposed", "microcosm", "location")),
        expand('output/dtu_{analysis}/dtu_{analysis}_.dexseq.results.dat', analysis=("exposed", "microcosm", "location")),
        "output/dexseq/dtu_exposed/exposure/sig_degs_annots.csv"


## exposure dtu ##

rule dtu_exposure_analysis:
    input:
        dds = "output/dexseq/dtu_exposed/dtu_exposed.dds",
        trinotate = "data/trinotate_annots.csv"
    output:
        dds_res = "output/dexseq/dtu_exposed/results.dds",
        sig_res = "output/dexseq/dtu_exposed/exposure/sig_degs.csv",
        sig_res_annots = "output/dexseq/dtu_exposed/exposure/sig_degs_annots.csv"
    log:
        "output/logs/dexseq/dtu_exposure_analysis.log"
    threads:
        30
    script:
        'src/exposed_dtu/exposure_analysis.R'

rule dtu_exposure_dds:
    input:
        sample_info_file = "data/sample_table_exposed.txt",
        full_info_file = "data/full_sample_table.csv",
        counts_info_file = "output/dtu_exposed/dtu_exposed_.sample_table",
        flatfile = "output/supertranscripts/trinity_genes.gtf.dexseq.gff"
    output:
        dds = "output/dexseq/dtu_exposed/dtu_exposed.dds"
    log:
        "output/logs/dexseq/dtu_exposure_dds.log"
    threads:
        30
    script:
        'src/exposed_dtu/make_dds.R'

############################
## diff. transcript usage ##
############################

#https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage
##pathlib to stop all output in base directory
rule dtu:
    input:
        fasta = str(pathlib2.Path(resolve_path('output/supertranscripts/'), 'trinity_genes.fasta')),
        gtf = str(pathlib2.Path(resolve_path('output/supertranscripts/'), 'trinity_genes.gtf')),
        sample_table = str(pathlib2.Path(resolve_path('data/'), 'sample_table_{analysis}.txt'))
    output:
        'output/dtu_{analysis}/dtu_{analysis}_.dexseq.results.dat'
    params:
        wd = 'output/dtu_{analysis}',
        prefix = 'dtu_{analysis}_'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'), 'dtu_{analysis}.log'))
    singularity:
        trinity_container
    threads:
        30
    shell:
        'cd {params.wd} || exit 1 ; '
        '/usr/local/bin/trinityrnaseq/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl '
        '--genes_fasta {input.fasta} '
        '--genes_gtf {input.gtf} '
        '--samples_file {input.sample_table} '
        '--out_prefix {params.prefix} '
        '--aligner STAR '
        '--CPU {threads} '
        '&> {log}'

#####################
## variant calling ##
#####################

##analyse variants https://murraycadzow.github.io/2021-obss-day3/04-popgen/index.html
##filter variants for depth - https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration
##GATK RNAseq best practices https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-

#https://github.com/trinityrnaseq/trinityrnaseq/wiki/Variant-Calling

##trinity page mentions filtering out homozygous variants BUT
##only applies if using a single sample that was used for assembly
##not when assembly contains multiple samples


##should variant calling be one sample at a time? - pipeline doesn't show multiple samples
##joint calling is a different pipeline - should repeat this cycling through each sample
##GATK - This workflow is designed to operate on a set of samples (uBAM files) one-at-a-time; joint calling RNAseq is not supported.

rule variants:
    input:
        fasta = 'output/supertranscripts/trinity_genes.fasta',
        gtf = 'output/supertranscripts/trinity_genes.gtf',
        sample_table = 'data/sample_table_{analysis}.txt'
    output:
        'output/variants_{analysis}/filtered_output.vcf'
    params:
        wd = 'output/variants_{analysis}'
    log:
        'output/logs/variants_{analysis}.log'
    singularity:
        trinity_container
    threads:
        30
    shell:
        '/usr/local/bin/trinityrnaseq/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py '
        '--st_fa {input.fasta} '
        '--st_gtf {input.gtf} '
        '-S {input.sample_table} '
        '-o {params.wd} '
        '-t {threads} '
        '-m 88193656874 '
        '&> {log}'

##############################
## supertranscript assembly ##
##############################

##pathlib2 to get full paths - have to cd to get output where I want
rule supertranscripts:
    input:
        str(pathlib2.Path(resolve_path('data/asw-transcriptome/output/trinity/'), "Trinity.fasta"))
    output:
        'output/supertranscripts/trinity_genes.fasta'
    params:
        wd = 'output/supertranscripts'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'), 'supertranscripts.log'))
    singularity:
        trinity_container
    threads:
        20
    shell:
        'cd {params.wd} || exit 1 ; '
        '/usr/local/bin/trinityrnaseq/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py '
        '--trinity_fasta {input} '
        '--incl_malign '
        '&> {log}'
