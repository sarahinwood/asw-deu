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
       	# location/exposed pairwise
        expand("output/dexseq/dtu_{analysis}_analysis_le/sig_degs_annots.csv", analysis=("location", "exposed")),
        # location-sp exposure analysis
        expand("output/dexseq/dtu_exposed_{location}_analysis/sig_degs_annots.csv", location=("Dunedin", "Ruakura")),
        # attack/para pairwise
        expand("output/dexseq/dtu_{status}_analysis_ap/sig_degs_annots.csv", status=("parasitism", "attack"))

# dexseq vignette
    # https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#7_Visualization

##########################
## attack/para analyses ##
##########################

rule dtu_parasitism_attack_analysis:
    input:
        trinotate_file = "data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv",
        dds_file = "output/dexseq/dtu_location_analysis/asw.dds"
    output:
        dds_res = "output/dexseq/dtu_{status}_analysis_ap/dtu_{status}_res.dds",
        all_res = "output/dexseq/dtu_{status}_analysis_ap/all_res.csv",
        sig_dtu_annots = "output/dexseq/dtu_{status}_analysis_ap/sig_degs_annots.csv",
        dtu_plots = "output/dexseq/dtu_{status}_analysis_ap/dtu_{status}_dexseq.pdf"
    log:
        "output/logs/dexseq/dtu_{status}_analysis.log"
    params:
        bp_threads = 8
    threads:
        8
    singularity:
        bioconductor_container
    script:
        'src/{wildcards.status}_dtu/{wildcards.status}_analysis.R'

#################################
## exposed & location analyses ##
#################################

rule dtu_exposed_location_analysis:
    input:
        trinotate_file = "data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv",
        dds_file = "output/dexseq/dtu_exposed_analysis/asw.dds"
    output:
        dds_res = "output/dexseq/dtu_exposed_{location}_analysis/exposed_{location}_res.dds",
        all_res = "output/dexseq/dtu_exposed_{location}_analysis/all_res.csv",
        sig_dtu_annots = "output/dexseq/dtu_exposed_{location}_analysis/sig_degs_annots.csv",
        dtu_plots = "output/dexseq/dtu_exposed_{location}_analysis/dtu_exposed_{location}_dexseq.pdf"
    log:
        "output/logs/dexseq/dtu_exposed_{location}_analysis.log"
    params:
        location = "{location}",
        bp_threads = 8
    threads:
        8
    singularity:
        bioconductor_container
    script:
        'src/exposed_location_dtu/exposed_location_analysis.R'

rule dtu_exposed_or_location_analysis:
    input:
        trinotate_file = "data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv",
        dds_file = "output/dexseq/dtu_{analysis}_analysis/asw.dds"
    output:
        dds_res = "output/dexseq/dtu_{analysis}_analysis_le/asw_{analysis}_res.dds",
        all_res = "output/dexseq/dtu_{analysis}_analysis_le/all_res.csv",
        sig_dtu_annots = "output/dexseq/dtu_{analysis}_analysis_le/sig_degs_annots.csv",
        dtu_plots = "output/dexseq/dtu_{analysis}_analysis_le/dtu_{analysis}_dexseq.pdf"
    params:
        bp_threads = 8
    threads:
        8
    log:
        "output/logs/dexseq/dtu_{analysis}_analysis.log"
    singularity:
        bioconductor_container
    script:
        'src/{wildcards.analysis}_dtu/{wildcards.analysis}_analysis.R'

####################
## Make DDS files ##
####################

rule dtu_make_dual_dds: # didn't write parasitism/attack status info
    input:
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

#################
## Trinity DTU ##
#################

# this rule should generate flatfile for future stuff - don't need results file as contains mh also
    # https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage
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


