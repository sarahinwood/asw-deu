#!/usr/bin/env python3
import pathlib2
import peppy

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

#containers
trinity_container = 'docker://trinityrnaseq/trinityrnaseq:2.11.0'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.12:0.0.1'

#########
# RULES #
#########

rule target:
    input:
        expand("output/dexseq/dtu_{analysis}_analysis/asw.dds", analysis=("exposed", "location")), #"microcosm"
       	#expand("output/dexseq/dtu_{analysis}_analysis/asw_{analysis}_res.dds", analysis=("exposed", "location")),
        "output/dexseq/dtu_exposed_location_analysis/asw_exposed_location_res.dds"

##dexseq vignette
    ##https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#7_Visualization

rule dtu_exposed_location_analysis: ##only using half as much bcc as dtu_analysis?
    input:
        trinotate_file = "data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv",
        dds_file = "output/dexseq/dtu_exposed_analysis/asw.dds"
    output:
        dds_res = "output/dexseq/dtu_exposed_location_analysis/asw_exposed_location_res.dds"
    log:
        "output/logs/dexseq/dtu_exposed_location_analysis.log"
    params:
        threads=30
    threads:
        30
    singularity:
        bioconductor_container
    script:
        'src/exposed_dtu/exposed_location_analysis.R'

rule dtu_analysis:
    input:
        trinotate_file = "data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv",
        dds_file = "output/dexseq/dtu_{analysis}_analysis/asw.dds"
    output:
        dds_res = "output/dexseq/dtu_{analysis}_analysis/asw_{analysis}_res.dds"
    log:
        "output/logs/dexseq/dtu_{analysis}_analysis.log"
    threads:
        30
    singularity:
        bioconductor_container
    script:
        'src/{wildcards.analysis}_dtu/{wildcards.analysis}_analysis.R'

rule dtu_make_dual_dds:
    input:
        sample_info_file = "data/sample_table_{analysis}.txt",
        full_info_file = "data/full_sample_table.csv",
        flatfile = "data/asw-mh-transcriptome/supertranscripts.gtf.dexseq.gff",
        counts_info_file = "output/dtu_{analysis}/dtu_{analysis}_.sample_table"
    output:
        dual_dxd = "output/dexseq/dtu_{analysis}_analysis/dual.dds",
        asw_dxd = "output/dexseq/dtu_{analysis}_analysis/asw.dds",
        mh_dxd = "output/dexseq/dtu_{analysis}_analysis/mh.dds"
    params:
        full_path = "output/dtu_{analysis}/"
    log:
        "output/logs/dexseq/dtu_{analysis}_dds.log"
    singularity:
        bioconductor_container
    script:
        'src/make_dds.R'

############################
## diff. transcript usage ##
############################

## this rule should generate flatfile for future stuff - don't need results file as contains mh also

#https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage
##pathlib to stop all output in base directory
rule dtu:
    input:
        fasta = str(pathlib2.Path(resolve_path('data/asw-mh-transcriptome'), 'supertranscripts.fasta')),
        gtf = str(pathlib2.Path(resolve_path('data/asw-mh-transcriptome'), 'supertranscripts.gtf')),
        sample_table = str(pathlib2.Path(resolve_path('data/'), 'sample_table_{analysis}.txt'))
    output:
        sample_table = 'output/dtu_{analysis}/dtu_{analysis}_.sample_table'
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


