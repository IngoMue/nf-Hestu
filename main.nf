#!/usr/bin/env nextflow

log.info """\

         ================================================
         nf-Hestu - A nextflow workflow for Heterozygosity estimation using genotype likelihoods
         https://github.com/IngoMue/nf-Hestu     
         Author: Ingo A. Müller
         ================================================
         |refseq....................: ${params.refseq}
         |chrlist...................: ${params.chrlist}
         |bamlist...................: ${params.bamlist}
         |outdir....................: ${params.outdir}
         |
         |Print population plots?...: ${params.inclPopPlots}
         ================================================
         """
         .stripIndent()

process calcMinMax {
    tag "Calculate minimum and maximum DoCs"
    publishDir "${params.outdir}/00_DoC/", pattern: '*.{list}', mode: 'copy'
    label 'Low_res'

    input:
      file(bam_list)

    output:
      file 'bams_DoCs.list'

    script:
    """
        01_calcMinMax.R $bam_list
    """
}

process indSAF {
    tag "Calculate site allele frequency likelihood for $sample_id "
    publishDir "${params.outdir}/01_SAF/${sample_id}", pattern: '*', mode: 'copy'

    input:
      tuple val(sample_id), file(bam_file), file(bam_index), val(DoC), val(minDoC), val(maxDoC)

    output:
      tuple val(sample_id), file("${sample_id}.saf.gz"), file("${sample_id}.saf.idx"), file("${sample_id}.saf.pos.gz")

    script:
    """
      angsd -i $bam_file -out $sample_id -doSaf 1 -GL 1 -rf ${params.chrlist} \
      -doMajorMinor 1 -ref $params.refseq -anc $params.refseq -doMaf 1 -doCounts 1 \
      -setMinDepth $minDoC -setMaxDepth $maxDoC -minQ 20 -minMapQ 20 -uniqueOnly 1 \
      -only_proper_pairs 1 -remove_bads 1 -baq 1 -C 50 -P ${task.cpus}
    """
}

process indSFS {
    tag "Calculate site frequency spectrum for $sample_id "
    publishDir "${params.outdir}/02_SFS/", pattern: '*', mode: 'copy'

    input:
      tuple val(sample_id), file(saf_file), file(saf_idx), file(saf_pos)

    output:
      file("${sample_id}.sfs")

    script:
    """
      realSFS $saf_idx -P ${task.cpus} -fold 1 > ${sample_id}.sfs
    """
}

process HetEst {
    tag "Calculate Heterozygosity"
    publishDir "${params.outdir}/03_Het/", pattern: '*', mode: 'move'
    label 'Low_res'

    input:
      path(sfs_file_list)

    output:
      file('*')

    script:
    """
      02_Hestu.R $params.bamlist $sfs_file_list
    """
}

process HetEst_pop {
    tag "Calculate Heterozygosity"
    publishDir "${params.outdir}/03_Het/", pattern: '*', mode: 'move'
    label 'Low_res'

    input:
      path(sfs_file_list)

    output:
      file('*')

    script:
    """
      03_Hestu_pop.R $params.bamlist $sfs_file_list
    """
}

workflow {

    calcMinMax(file(params.bamlist))

    inds = calcMinMax.out.splitCsv(header:true, sep:'\t')
		.map { row -> tuple(row.ind, file(row.bam), file("${row.bam}.bai"), row.DoC, row.minDoC, row.maxDoC) }

    indSAF(inds) | indSFS
    
    if (params.inclPopPlots == true) {  

    HetEst_pop(indSFS.out.collect())

    } else {
        
    HetEst(indSFS.out.collect())
    
    }

}
