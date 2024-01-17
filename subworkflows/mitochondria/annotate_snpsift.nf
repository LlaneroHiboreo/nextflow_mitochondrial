include { SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_GNOMAD                   } from '../../modules/nf-core/snpsift/annotate/main'
include { SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_MITOMAP_DISEASE          } from '../../modules/nf-core/snpsift/annotate/main'
include { SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_MITOMAP_POLYMORPHISMS    } from '../../modules/nf-core/snpsift/annotate/main'
include { SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_MITOMAP_TIP              } from '../../modules/nf-core/snpsift/annotate/main'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_LIFT_UNLIFT_MT             } from '../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_MT         } from '../../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_MT                       } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT                                 } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_MT                         } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT2                                } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_MT                           } from '../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX as TABIX_TABIX_MERGE                              } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_MT3                                } from '../../modules/nf-core/tabix/tabix/main'
include { VCFANNO as VCFANNO_MT                                         } from '../../modules/nf-core/vcfanno/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_PASS                       } from '../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_MITOMAP                    } from '../../modules/nf-core/bcftools/filter/main'
//include { VCF2TABULAR                                                   } from '../../modules/local/vcf2tabular'
include { IGVREPORTS as IGVREPORTS_NO_FILTERED                          } from '../../modules/local/igvreports'
include { IGVREPORTS as IGVREPORTS_FILTERED                             } from '../../modules/local/igvreports'

workflow ANNOTATE_SNPSIFT {

    take:
        ch_vcf1                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_vcf2                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_dict         // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_case_info           // channel: [mandatory] [ val(case_info) ]
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

        // JOIN blacklist and index
        ch_blacklist = ch_blacklist.join(ch_blacklist_idx)
        
        ch_vcfs = ch_vcf1
            .join(ch_vcf2, remainder: true)
            .map{ meta, vcf1, vcf2 ->
            [meta, [vcf1, vcf2]]
        }
        GATK4_MERGEVCFS_LIFT_UNLIFT_MT( ch_vcfs, ch_genome_dict)

        // Filtering Variants
        GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.vcf
            .join(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_filt_vcf }
        GATK4_VARIANTFILTRATION_MT (ch_filt_vcf, ch_genome_fasta, ch_genome_fai, ch_genome_dict, ch_blacklist)

        // Spliting multiallelic calls
        GATK4_VARIANTFILTRATION_MT.out.vcf
            .join(GATK4_VARIANTFILTRATION_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_in_split }
        SPLIT_MULTIALLELICS_MT (ch_in_split, ch_genome_fasta)
        TABIX_TABIX_MT(SPLIT_MULTIALLELICS_MT.out.vcf)

        // Removing duplicates and merging if there is more than one sample
        SPLIT_MULTIALLELICS_MT.out.vcf
            .join(TABIX_TABIX_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_in_remdup }
        REMOVE_DUPLICATES_MT(ch_in_remdup, ch_genome_fasta)
        TABIX_TABIX_MT2(REMOVE_DUPLICATES_MT.out.vcf)

        REMOVE_DUPLICATES_MT.out.vcf
            .collect{it[1]}
            .ifEmpty([])
            .toList()
            .set { file_list_vcf }

        TABIX_TABIX_MT2.out.tbi
            .collect{it[1]}
            .ifEmpty([])
            .toList()
            .set { file_list_tbi }

        ch_case_info
            .combine(file_list_vcf)
            .combine(file_list_tbi)
            .set { ch_rem_dup_vcf_tbi }

        ch_rem_dup_vcf_tbi.branch {
            meta, vcf, tbi ->
                single: vcf.size() == 1
                    return [meta, vcf]
                multiple: vcf.size() > 1
                    return [meta, vcf, tbi]
            }.set { ch_case_vcf }

        BCFTOOLS_MERGE_MT( ch_case_vcf.multiple,
            ch_genome_fasta,
            ch_genome_fai,
            []
        )

        BCFTOOLS_MERGE_MT.out.merged_variants
            .mix(ch_case_vcf.single)
            .set { ch_annotation_in }
        
        TABIX_TABIX_MERGE(ch_annotation_in)

        // ANNOTATE
        ch_complete_annotation = ch_annotation_in.join(TABIX_TABIX_MERGE.out.tbi).collect()
        snpsift_dbs_gnomad = ch_snpsift_gnomad.join(ch_snpsift_gnomad_tbi)
        SNPSIFT_ANNOTATE_GNOMAD(ch_complete_annotation, snpsift_dbs_gnomad)
        
        ch_snpsift_mitomap_disease_dbs = ch_snpsift_mitomap_disease.join(ch_snpsift_mitomap_disease_tbi).collect()
        SNPSIFT_ANNOTATE_GNOMAD.out.vcf.map{meta,vcf->
            [meta,vcf,[]]
        }.set{ch_snpsift_mitomap_disease_dbs_in}

        SNPSIFT_ANNOTATE_MITOMAP_DISEASE(ch_snpsift_mitomap_disease_dbs_in,ch_snpsift_mitomap_disease_dbs)

        ch_snpsift_mitomap_polymorphism_dbs = ch_snpsift_mitomap_polymorphism.join(ch_snpsift_mitomap_polymorphism_tbi)
        SNPSIFT_ANNOTATE_MITOMAP_DISEASE.out.vcf.map{meta,vcf->
            [meta, vcf, []]
        }.set{ch_snpsift_mitomap_polymorphism_dbs_in}
        SNPSIFT_ANNOTATE_MITOMAP_POLYMORPHISMS(ch_snpsift_mitomap_polymorphism_dbs_in,ch_snpsift_mitomap_polymorphism_dbs)
        
        ch_snpsift_mitomap_tip_dbs = ch_snpsift_mitomap_tip.join(ch_snpsift_mitomap_tip_tbi).collect()
        SNPSIFT_ANNOTATE_MITOMAP_POLYMORPHISMS.out.vcf.map{meta,vcf->
            [meta, vcf, []]
        }.set{ch_snpsift_mitomap_tip_dbs_in}
        SNPSIFT_ANNOTATE_MITOMAP_TIP(ch_snpsift_mitomap_tip_dbs_in,ch_snpsift_mitomap_tip_dbs)
        
        // VCF2TABULAR
        //VCF2TABULAR(SNPSIFT_ANNOTATE_MITOMAP_TIP.out.vcf)

        // FILTER
        BCFTOOLS_FILTER_PASS(SNPSIFT_ANNOTATE_MITOMAP_TIP.out.vcf)
        BCFTOOLS_FILTER_MITOMAP(BCFTOOLS_FILTER_PASS.out.vcf)

        // GENERATE REPORT
        IGVREPORTS_NO_FILTERED(BCFTOOLS_FILTER_PASS.out.vcf, ch_genome_fasta, ch_igvreport_ideogram)
        IGVREPORTS_FILTERED(BCFTOOLS_FILTER_MITOMAP.out.vcf, ch_genome_fasta, ch_igvreport_ideogram)
    emit:
        versions  = ch_versions                    // channel: [ path(versions.yml) ]
}