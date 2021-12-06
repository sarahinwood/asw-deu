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
        expand('output/st_variants_{analysis}/filtered_output.vcf', analysis=("exposed", "microcosm", "location")),
        expand('output/dtu_{analysis}/dtu.dexseq.results.dat', analysis=("exposed", "microcosm", "location"))

############################
## diff. transcript usage ##
############################

#https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage

rule dtu:
    input:
        fasta = 'output/supertranscripts/trinity_genes.fasta',
        gtf = 'output/supertranscripts/trinity_genes.gtf',
        sample_table = 'data/sample_table_{analysis}.txt'
    output:
        'output/dtu_{analysis}/dtu.dexseq.results.dat'
    params:
        prefix = 'output/dtu_{analysis}/dtu'
    log:
        'output/logs/dtu_{analysis}.log'
    singularity:
        trinity_container
    threads:
        20
    shell:
        '/usr/local/bin/trinityrnaseq/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl '
        '--genes_fasta {input.fasta} '
        '--genes_gtf {input.gtf} '
        '--samples_file {input.sample_table} '
        '--out_prefix {params.prefix} '
        '--aligner STAR '
        '&> {log}'

#####################
## variant calling ##
#####################

##analyse variants https://murraycadzow.github.io/2021-obss-day3/04-popgen/index.html
##filter variants for depth - https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration

#https://github.com/trinityrnaseq/trinityrnaseq/wiki/Variant-Calling

##trinity page mentions filtering out homozygous variants BUT
##only applies if using a single sample that was used for assembly
##not when assembly contains multiple samples

rule st_variants:
    input:
        fasta = 'output/supertranscripts/trinity_genes.fasta',
        gtf = 'output/supertranscripts/trinity_genes.gtf',
        sample_table = 'data/sample_table_{analysis}.txt'
    output:
        'output/st_variants_{analysis}/filtered_output.vcf'
    params:
        wd = 'output/st_variants_{analysis}'
    log:
        'output/logs/st_variants_{analysis}.log'
    singularity:
        trinity_container
    threads:
        20
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