// Import subworkflows/modules
include { INPUT_CHECK   } from '../subworkflows/input_check'
include { ANALYSE_MT    } from '../subworkflows/analyse_MT'

// Create path for dbs
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
ch_genome_fasta             = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
ch_genome_fai               = Channel.fromPath(params.fasta_fai).map {it -> [[id:it[0].simpleName], it]}.collect()
ch_genome_dictionary        = Channel.fromPath(params.dict).map {it -> [[id:it[0].simpleName], it]}.collect()
ch_mt_intervals             = Channel.value([])
ch_mtshift_intervals        = Channel.value([])
ch_genome_bwamem2index      = Channel.fromPath(params.bwamem2).map {it -> [[id:it[0].simpleName], it]}.collect()
ch_mtshift_bwamem2index     = Channel.fromPath(params.bwamem2mt).map {it -> [[id:it[0].simpleName], it]}.collect()
ch_mtshift_fasta            = Channel.fromPath(params.shift_fasta).map {it -> [[id:it[0].simpleName], it]}.collect()
ch_mtshift_dictionary       = Channel.fromPath(params.shift_dict).map {it -> [[id:it[0].simpleName], it]}.collect()
ch_mtshift_fai              = Channel.fromPath(params.shift_fai).map {it -> [[id:it[0].simpleName], it]}.collect()
ch_mtshift_backchain        = Channel.fromPath(params.shift_backchain).map {it -> [[id:it[0].simpleName], it]}.collect()
ch_snpsift_dbs_gnomad       = params.snpsift_gnomad                    ? Channel.fromPath(params.snpsift_gnomad).map { it -> [[id:'snpsift_annotate_gnomad'], it] }.collect(): Channel.value([[],[]])
ch_snpsift_dbs_gnomad_tbi   = params.snpsift_gnomad_tbi                ? Channel.fromPath(params.snpsift_gnomad_tbi).map { it -> [[id:'snpsift_annotate_gnomad'], it] }.collect(): Channel.value([[],[]])
ch_snpsift_mitomap_disease       = params.snpsift_mitomap_disease                ? Channel.fromPath(params.snpsift_mitomap_disease).map { it -> [[id:'snpsift_annotate_gnomad'], it] }.collect(): Channel.value([[],[]])
ch_snpsift_mitomap_disease_tbi   = params.snpsift_mitomap_disease_tbi            ? Channel.fromPath(params.snpsift_mitomap_disease_tbi).map { it -> [[id:'snpsift_annotate_gnomad'], it] }.collect(): Channel.value([[],[]])
ch_snpsift_mitomap_polymorphism       = params.snpsift_mitomap_polymorphism                ? Channel.fromPath(params.snpsift_mitomap_polymorphism).map { it -> [[id:'snpsift_annotate_gnomad'], it] }.collect(): Channel.value([[],[]])
ch_snpsift_mitomap_polymorphism_tbi   = params.snpsift_mitomap_polymorphism_tbi            ? Channel.fromPath(params.snpsift_mitomap_polymorphism_tbi).map { it -> [[id:'snpsift_annotate_gnomad'], it] }.collect(): Channel.value([[],[]])
ch_snpsift_mitomap_tip       = params.snpsift_mitomap_tip                ? Channel.fromPath(params.snpsift_mitomap_tip).map { it -> [[id:'snpsift_annotate_gnomad'], it] }.collect(): Channel.value([[],[]])
ch_snpsift_mitomap_tip_tbi   = params.snpsift_mitomap_tip_tbi            ? Channel.fromPath(params.snpsift_mitomap_tip_tbi).map { it -> [[id:'snpsift_annotate_gnomad'], it] }.collect(): Channel.value([[],[]])
ch_igvreport_ideogram        = params.igvreport_ideogram            ? Channel.fromPath(params.igvreport_ideogram).collect(): Channel.value([])
ch_blacklist                 = params.blacklist ? Channel.fromPath(params.blacklist).map { it -> [[id:'blacklisted_sites'], it] }.collect(): Channel.value([[],[]])
ch_blacklist_idx                 = params.blacklist_idx ? Channel.fromPath(params.blacklist_idx).map { it -> [[id:'blacklisted_sites'], it] }.collect(): Channel.value([[],[]])

workflow GATK_MITOCHONDRIA {

    INPUT_CHECK (
        ch_input
    )

    // REFORMAT CHANNEL
    ch_reads_meta = INPUT_CHECK.out.reads.map{meta, bam_bai->
        [
            meta,bam_bai[0],bam_bai[1]
        ]
    }

    // RUN WORKFLOW
    ANALYSE_MT (
        ch_reads_meta,
        ch_genome_bwamem2index,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dictionary,
        ch_mt_intervals,
        ch_mtshift_bwamem2index,
        ch_mtshift_fasta,
        ch_mtshift_dictionary,
        ch_mtshift_fai,
        ch_mtshift_intervals,
        ch_mtshift_backchain,
        INPUT_CHECK.out.case_info,        
        ch_snpsift_dbs_gnomad,
        ch_snpsift_dbs_gnomad_tbi,
        ch_snpsift_mitomap_disease,
        ch_snpsift_mitomap_disease_tbi,
        ch_snpsift_mitomap_polymorphism,
        ch_snpsift_mitomap_polymorphism_tbi,
        ch_snpsift_mitomap_tip,
        ch_snpsift_mitomap_tip_tbi,
        ch_igvreport_ideogram,
        ch_blacklist,
        ch_blacklist_idx
    )
}