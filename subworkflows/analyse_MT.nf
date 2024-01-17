//
// Analyse MT
//
include { CONVERT_MT_BAM_TO_FASTQ                      } from './mitochondria/convert_mt_bam_to_fastq'
include { ALIGN_AND_CALL_MT                            } from './mitochondria/align_and_call_MT'
include { ALIGN_AND_CALL_MT as ALIGN_AND_CALL_MT_SHIFT } from './mitochondria/align_and_call_MT'
include { PICARD_LIFTOVERVCF                           } from '../modules/nf-core/picard/liftovervcf/main'
include { ANNOTATE_SNPSIFT                             } from './mitochondria/annotate_snpsift'

workflow ANALYSE_MT {
    take:
        ch_cram_crai               // channel: [mandatory] [ val(meta), file(bam), file(bai) ]
        ch_genome_bwamem2_index  // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai            // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_dict           // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mt_intervals          // channel: [mandatory] [ path(interval_list) ]
        ch_mtshift_bwamem2index  // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mtshift_dict          // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_fai           // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mtshift_intervals     // channel: [mandatory] [ path(interval_list) ]
        ch_mtshift_backchain     // channel: [mandatory] [ val(meta), path(back_chain) ]
        ch_case_info             // channel: [mandatory] [ val(case_info) ]
        ch_snpsift_gnomad
        ch_snpsift_gnomad_tbi
        ch_snpsift_mitomap_disease
        ch_snpsift_mitomap_disease_tbi
        ch_snpsift_mitomap_polymorphism
        ch_snpsift_mitomap_polymorphism_tbi
        ch_snpsift_mitomap_tip
        ch_snpsift_mitomap_tip_tbi
        ch_igvreport_ideogram
        ch_blacklist
        ch_blacklist_idx
        
    main:
        ch_versions = Channel.empty()

        // PREPARING READS FOR MT ALIGNMENT
        CONVERT_MT_BAM_TO_FASTQ (
            ch_cram_crai,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dict
        )

        // MT ALIGNMENT  AND VARIANT CALLING
        ALIGN_AND_CALL_MT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.cram,
            ch_genome_bwamem2_index,
            ch_genome_fasta,
            ch_genome_dict,
            ch_genome_fai,
            ch_mt_intervals
        )
        
        
        ALIGN_AND_CALL_MT_SHIFT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.cram,
            ch_mtshift_bwamem2index,
            ch_mtshift_fasta,
            ch_mtshift_dict,
            ch_mtshift_fai,
            ch_mtshift_intervals
        )

        // LIFTOVER VCF FROM REFERENCE MT TO SHIFTED MT
        PICARD_LIFTOVERVCF (
            ALIGN_AND_CALL_MT_SHIFT.out.vcf,
            ch_genome_dict,
            ch_genome_fasta,
            ch_mtshift_backchain
        )


        // MT MERGE AND ANNOTATE VARIANTS
        ANNOTATE_SNPSIFT(            
            ALIGN_AND_CALL_MT.out.vcf,
            PICARD_LIFTOVERVCF.out.vcf_lifted,
            ch_genome_fasta,
            ch_genome_dict,
            ch_genome_fai,
            ch_case_info,
            ch_snpsift_gnomad,
            ch_snpsift_gnomad_tbi,
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
        //
        /*

        ch_versions = ch_versions.mix(CONVERT_MT_BAM_TO_FASTQ.out.versions)
        ch_versions = ch_versions.mix(ALIGN_AND_CALL_MT.out.versions)
        ch_versions = ch_versions.mix(ALIGN_AND_CALL_MT_SHIFT.out.versions)
        ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF.out.versions.first())


    emit:
        vcf            = MERGE_ANNOTATE_MT.out.vcf              // channel: [ val(meta), path(vcf) ]
        tbi            = MERGE_ANNOTATE_MT.out.tbi              // channel: [ val(meta), path(tbi) ]
        stats          = ALIGN_AND_CALL_MT.out.stats            // channel: [ val(meta), path(stats) ]
        filt_stats     = ALIGN_AND_CALL_MT.out.filt_stats       // channel: [ val(meta), path(tsv) ]
        mt_del_result  = ALIGN_AND_CALL_MT.out.mt_del_result    // channel: [ val(meta), path(txt) ]
        stats_sh       = ALIGN_AND_CALL_MT_SHIFT.out.stats      // channel: [ val(meta), path(stats) ]
        filt_stats_sh  = ALIGN_AND_CALL_MT_SHIFT.out.filt_stats // channel: [ val(meta), path(tsv) ]
        eklipse_del    = ALIGN_AND_CALL_MT.out.eklipse_del      // channel: [ val(meta), path(csv) ]
        eklipse_genes  = ALIGN_AND_CALL_MT.out.eklipse_genes    // channel: [ val(meta), path(csv) ]
        eklipse_circos = ALIGN_AND_CALL_MT.out.eklipse_circos   // channel: [ val(meta), path(png) ]
        haplog         = MERGE_ANNOTATE_MT.out.haplog           // channel: [ val(meta), path(txt) ]
        report         = MERGE_ANNOTATE_MT.out.report           // channel: [ path(html) ]
        txt            = ALIGN_AND_CALL_MT.out.txt              // channel: [ val(meta), path(txt) ]
        html           = ALIGN_AND_CALL_MT.out.html             // channel: [ val(meta), path(html) ]
        txt_sh         = ALIGN_AND_CALL_MT_SHIFT.out.txt        // channel: [ val(meta), path(txt) ]
        html_sh        = ALIGN_AND_CALL_MT_SHIFT.out.html       // channel: [ val(meta), path(html) ]
        versions       = ch_versions                            // channel: [ path(versions.yml) ]
    */
}