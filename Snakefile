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
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.14:0.0.1' # DEXSeq v
old_bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.12:0.0.1'

#########
# RULES #
#########

rule target:
    input:
        # location/exposed pairwise
        expand("output/dexseq/dtu_{analysis}_analysis/sig_degs_annots.csv", analysis=("location", "exposed")),
        # location-sp exposure analysis
        expand("output/dexseq/dtu_exposed_{location}_analysis/sig_degs_annots.csv", location=("Dunedin", "Ruakura")),
        # attack/para pairwise
        expand("output/dexseq/dtu_{status}_analysis_ap/sig_degs_annots.csv", status=("parasitism", "attack")),
        # location-sp attack analysis
        #expand("output/dexseq/dtu_attack_{location}_analysis/sig_degs_annots.csv", location=("Dunedin", "Ruakura")),
        # GO/Pfam enrichment
        expand('output/term_overrepresentation/dtu_{analysis}/{analysis}_Pfam_overrep.csv', analysis=("location_analysis", "parasitism_analysis_ap", "exposed_Dunedin_analysis", "exposed_analysis")), #"exposed_Ruakura_analysis" - no pfam , "attack_analysis_ap", "attack_Dunedin_analysis", "attack_Ruakura_analysis", "exposed_analysis" - not enough genes
        expand('output/term_overrepresentation/dtu_{analysis}/{analysis}_BlastP_GO_overrep.csv', analysis=("location_analysis", "parasitism_analysis_ap", "attack_analysis_ap", "attack_Dunedin_analysis", "exposed_analysis")), #  "exposed_Dunedin_analysis", "exposed_Ruakura_analysis", "attack_Ruakura_analysis" - no enrichment, "attack_Dunedin_analysis", "exposed_analysis" - not enough genes
        expand('output/term_overrepresentation/dtu_{analysis}/{analysis}_overrep_plot_BlastP_GO_new.pdf', analysis=("location_analysis", "parasitism_analysis_ap", "exposed_analysis"))

# dexseq vignettes
    # https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#7_Visualization

rule dtu_location_design_analyses:
    input:
        trinotate_file = "data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv",
        dds_file = "output/dexseq/dtu_location_analysis/asw.dds"
    output:
        dds_res = "output/dexseq/dtu_location_{design}_analysis/asw_{design}_res.dds",
        all_res = "output/dexseq/dtu_location_{design}_analysis/all_res.csv",
        sig_dtu_annots = "output/dexseq/dtu_location_{design}_analysis/sig_degs_annots.csv",
        dtu_plots = "output/dexseq/dtu_location_{design}_analysis/dtu_{design}_dexseq.pdf"
    params:
        bp_threads = 8
    threads:
        8
    log:
        "output/logs/dexseq/dtu_location_{design}_analysis.log"
    singularity:
        bioconductor_container
    script:
        'src/location_dtu/location_analysis_{wildcards.design}.R'

################################
## GO/Pfam overrepresentation ##
################################

rule plot_overrepresentation:
    input:
        pfam_file = 'output/term_overrepresentation/dtu_{analysis}/{analysis}_Pfam_overrep.csv',
        go_file  = 'output/term_overrepresentation/dtu_{analysis}/{analysis}_BlastP_GO_overrep.csv'
    output:
        enrich_plot = 'output/term_overrepresentation/dtu_{analysis}/{analysis}_overrep_plot_BlastP_GO_new.pdf'
    singularity:
        old_bioconductor_container
    log:
        'output/logs/term_overrepresentation/plot_{analysis}_overrepresentation.log'
    script:
        'src/term_overrepresentation/plot_overrepresentation.R'

rule term_overrepresentation:
    input:
        DEG_file = 'output/dexseq/dtu_{analysis}/sig_degs_annots.csv',
        term_annot_table_file = 'output/term_overrepresentation/{term}_annots/{term}_annots.csv',
        term_to_gene_file = 'output/term_overrepresentation/{term}_annots/{term}_to_gene.csv',
        term_to_name_file = 'output/term_overrepresentation/{term}_annots/{term}_to_name.csv'
    output:
        enrichment_table = 'output/term_overrepresentation/dtu_{analysis}/{analysis}_{term}_overrep.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/term_overrepresentation/{analysis}_analysis_{term}_enrichment.log'
    script:
        'src/term_overrepresentation/{wildcards.term}_overrepresentation.R'

rule term_to_gene:
    input:
        trinotate_file = 'data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_annotation_report.txt'
    output:
        term_annot_table = 'output/term_overrepresentation/{term}_annots/{term}_annots.csv',
        term_to_gene = 'output/term_overrepresentation/{term}_annots/{term}_to_gene.csv',
        term_to_name = 'output/term_overrepresentation/{term}_annots/{term}_to_name.csv',
        term_sizes = 'output/term_overrepresentation/{term}_annots/{term}_sizes.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/term_overrepresentation/{term}_term_to_gene.log'
    script:
        'src/term_overrepresentation/{wildcards.term}_to_geneID.R'

##########################
## attack/para analyses ##
##########################

rule dtu_attack_location_analysis:
    input:
        trinotate_file = "data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv",
        dds_file = "output/dexseq/dtu_location_analysis/asw.dds"
    output:
        dds_res = "output/dexseq/dtu_attack_{location}_analysis/attack_{location}_res.dds",
        all_res = "output/dexseq/dtu_attack_{location}_analysis/all_res.csv",
        sig_dtu_annots = "output/dexseq/dtu_attack_{location}_analysis/sig_degs_annots.csv",
        dtu_plots = "output/dexseq/dtu_attack_{location}_analysis/dtu_attack_{location}_dexseq.pdf"
    log:
        "output/logs/dexseq/dtu_attack_{location}_analysis.log"
    params:
        location = "{location}",
        bp_threads = 8
    threads:
        8
    singularity:
        bioconductor_container
    script:
        'src/attack_location_dtu/attack_location_analysis.R'

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
        dds_res = "output/dexseq/dtu_{analysis}_analysis/asw_{analysis}_res.dds",
        all_res = "output/dexseq/dtu_{analysis}_analysis/all_res.csv",
        sig_dtu_annots = "output/dexseq/dtu_{analysis}_analysis/sig_degs_annots.csv",
        dtu_plots = "output/dexseq/dtu_{analysis}_analysis/dtu_{analysis}_dexseq.pdf"
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
        flatfile = "output/all_supertranscripts/all_supertranscripts.gtf.dexseq.gff",
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

# this rule should generate flatfile and mapped bams for future stuff - don't need DEXseq results file as contains mh also
    # https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage
rule dtu:
    input:
        fasta = str(pathlib2.Path(resolve_path('output/all_supertranscripts'), 'all_supertranscripts.fasta')),
        gtf = str(pathlib2.Path(resolve_path('output/all_supertranscripts'), 'all_supertranscripts.gtf')),
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
