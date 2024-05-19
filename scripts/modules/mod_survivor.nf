

exclude_list=file(params.exclude_list)

params.enabled=true


params.vsuffix="*.vcf"


process truvari_merge_s {
    cpus 2
    memory '60 GB'
    input:
    path(all_vcfs)
    output:
      tuple path("${out}"), path("${out}.tbi")
    script:
      out = "${params.single}.vcf.gz"
      ref = params.ref
"""
bcftools merge -m none ${params.vsuffix} | bcftools annotate -x INFO/SR,FORMAT/PR,FORMAT/PL,FORMAT/RCR,FORMAT/RCL,FORMAT/DR,FORMAT/GL,FORMAT/DV,FORMAT/DV,FORMAT/RC,FORMAT/RR,FORMAT/RV,FORMAT/RDCN,FORMAT/SR,FORMAT/SQ,FORMAT/RP,FORMAT/QR,FORMAT/AO,FORMAT/RS,FORMAT/DHFC,FORMAT/QA,FORMAT/AS,FORMAT/ASC,FORMAT/DHSP,FORMAT/RO,FORMAT/DHFFC,FORMAT/AP,FORMAT/DP,FORMAT/AB,FORMAT/DHBFC,FORMAT/DHSP,FORMAT/SHQ,FORMAT/FT -Ob -o merge.bcf
bcftools filter -i "ALT='<DEL>'"  -Oz -o dels.vcf.gz merge.bcf
bcftools index -t dels.vcf.gz
truvari collapse -p 0 -P 0.75 -i dels.vcf.gz -o d_mrg.vcf -c d_col.vcf -f $ref >& del.log
bcftools sort -m 10G d_mrg.vcf -Ob -o d_mrg.bcf
bcftools index d_mrg.bcf
bcftools filter -e "ALT='<DEL>'"  -Oz -o others.vcf.gz merge.bcf
bcftools index -t others.vcf.gz
truvari collapse -p 0.9 -P 0.8 -i others.vcf.gz -o oth_mrg.vcf -c oth_col.vcf -f $ref >& oth.log
bcftools sort -m 10G  oth_mrg.vcf -Ob -o oth_mrg.bcf
bcftools index oth_mrg.bcf
bcftools concat -a d_mrg.bcf oth_mrg.bcf -Ob -o combined.bcf
bcftools index  combined.bcf
bcftools sort -m 10G  combined.bcf -Oz -o ${out}
bcftools index -t ${out}
"""
}

process truvari_merge {
    cpus 2
    memory '60 GB'
    input:
    path(all_vcfs)
    output:
      tuple path("${out}"), path("${out}.tbi")
    script:
      out = "${params.single}.vcf.gz"
      ref = params.ref
"""
bcftools merge -m none *vcf.gz | bcftools annotate -x INFO/SR,FORMAT/PR,FORMAT/PL,FORMAT/RCR,FORMAT/RCL,FORMAT/DR,FORMAT/GL,FORMAT/DV,FORMAT/DV,FORMAT/RC,FORMAT/RR,FORMAT/RV,FORMAT/RDCN,FORMAT/SR,FORMAT/SQ,FORMAT/RP,FORMAT/QR,FORMAT/AO,FORMAT/RS,FORMAT/DHFC,FORMAT/QA,FORMAT/AS,FORMAT/ASC,FORMAT/DHSP,FORMAT/RO,FORMAT/DHFFC,FORMAT/AP,FORMAT/DP,FORMAT/AB,FORMAT/DHBFC,FORMAT/DHSP,FORMAT/SHQ,FORMAT/FT -Ob -o merge.bcf
bcftools filter -i "ALT='<DEL>'"  -Oz -o dels.vcf.gz merge.bcf
bcftools index -t dels.vcf.gz
truvari collapse -p 0 -P 0.75 -i dels.vcf.gz -o d_mrg.vcf -c d_col.vcf -f $ref >& del.log
bcftools sort -m 10G d_mrg.vcf -Ob -o d_mrg.bcf
bcftools index d_mrg.bcf
bcftools filter -e "ALT='<DEL>'"  -Oz -o others.vcf.gz merge.bcf
bcftools index -t others.vcf.gz
truvari collapse -p 0.9 -P 0.8 -i others.vcf.gz -o oth_mrg.vcf -c oth_col.vcf -f $ref >& oth.log
bcftools sort -m 10G  oth_mrg.vcf -Ob -o oth_mrg.bcf
bcftools index oth_mrg.bcf
bcftools concat -a d_mrg.bcf oth_mrg.bcf -Ob -o combined.bcf
bcftools index  combined.bcf
bcftools sort -m 10G  combined.bcf -Oz -o ${out}
bcftools index -t ${out}
"""
}


process tmerge_step1_merge {
    cpus 2
    memory '60 GB'
    input:
      path(all_vcfs)
    output:
      tuple val(base), path("merge.bcf")
    script:
    if (params.single=="") 
        out="${all_vcfs[0].simpleName}.vcf.gz"
    else
	out="${params.single}.vcf.gz"
    base=out.replaceAll(".vcf.gz","")
    ref = params.ref
   """
   bcftools ${params.cmd} ${params.vsuffix} | bcftools annotate -x INFO/SR,FORMAT/PR,FORMAT/PL,FORMAT/RCR,FORMAT/RCL,FORMAT/DR,FORMAT/GL,FORMAT/DV,FORMAT/DV,FORMAT/RC,FORMAT/RR,FORMAT/RV,FORMAT/RDCN,FORMAT/SR,FORMAT/SQ,FORMAT/RP,FORMAT/QR,FORMAT/AO,FORMAT/RS,FORMAT/DHFC,FORMAT/QA,FORMAT/AS,FORMAT/ASC,FORMAT/DHSP,FORMAT/RO,FORMAT/DHFFC,FORMAT/AP,FORMAT/DP,FORMAT/AB,FORMAT/DHBFC,FORMAT/DHSP,FORMAT/SHQ,FORMAT/FT -Ob -o merge.bcf
    """
}

process tmerge_step2_filter {
    cpus 2
    input:
     tuple val(base), path("merge.bcf")
    output:
      tuple val(base), path("dels.vcf.gz"), path("dels.vcf.gz.tbi"), emit: del_ch
      tuple val(base), path("others.vcf.gz"), path("others.vcf.gz.tbi"), emit: oth_ch
    script:
    """
       bcftools filter -i "ALT='<DEL>'"  -Oz -o dels.vcf.gz merge.bcf
       bcftools index -t dels.vcf.gz
       bcftools filter -e "ALT='<DEL>'"  -Oz -o others.vcf.gz merge.bcf
       bcftools index -t others.vcf.gz
    """
}


process tmerge_step3a_del_collapse {
    memory '60 GB'
    input:
      tuple val(base), path(dvcf), path(di)
    output:
      tuple val(base), path("d_mrg.bcf"), path("d_mrg.bcf.csi")
    script:
      ref = params.ref
    """
    truvari collapse -p 0 -P 0.75 -i $dvcf -o d_mrg.vcf -c d_col.vcf -f $ref >& del.log
    bcftools sort -m 20G d_mrg.vcf -Ob -o d_mrg.bcf
    bcftools index d_mrg.bcf
    """
}



process tmerge_step3b_oth_collapse {
    memory '60 GB'
    input:
      tuple val(base), path(ovcf), path(oi)
    output:
      tuple val(base),path("oth_mrg.bcf"), path("oth_mrg.bcf.csi")
    script:
      ref = params.ref
    """
    truvari collapse -p 0.9 -P 0.8 -i $ovcf -o oth_mrg.vcf -c oth_col.vcf -f $ref >& oth.lg
    bcftools sort -m 20G  oth_mrg.vcf -Ob -o oth_mrg.bcf
    bcftools index oth_mrg.bcf
    """
}

process tmerge_step4_merge {
    memory '40 GB'
    input:
      tuple val(base), path(files), path(indices)
    output:
      tuple path(out), path("${out}.tbi")
    publishDir params.out_dir, mode:'copy'
    script:
      out = "${base}.vcf.gz"
     """
     bcftools concat -a *bcf -Ob -o combined.bcf
     bcftools index  combined.bcf
     bcftools sort -m 30G  combined.bcf -Oz -o ${out}
     bcftools index -t ${out}
     """
}


// Take VCFs for same person with multiple tools
process survivorByTool {
    tag { base }
    memory '5 GB'
    maxForks 80
    errorStrategy 'finish'
    input:
       tuple val(base), path('sample?.vcf') 
       each support 
    output:
       tuple val(support), file(out) 
    clusterOptions = '--constraint=avx2'
  script:
    out = "${base}.vcf"
    """
    /bin/hostname
    ls *vcf > merge_list
    SURVIVOR merge merge_list 1000 $support 0 0 0 50 $out
    """
}



process mergeSamples {
    tag { "overall" }
    memory '80 GB'
    clusterOptions = '--constraint=avx2'
    input:
       tuple val(support), file(vcfs) 
    output:
       tuple val(support), path("${out}"), path("${out}.tbi")
    publishDir "${params.out_dir}/", mode: 'copy', enabled: params.enabled
  script:
    out = "res${params.sub_tag}_${support}.vcf.gz"
    """
    ls *vcf > vcf_list
    SURVIVOR merge vcf_list 1000 1 1 1 0 50 all.vcf
    decolonise.py all.vcf | bcftools sort -Oz -o ${out}
    tabix ${out}
    """
}


process filterOut {  // legacy code -- this now done in clean up -- remove at a later stage
    memory "8 GB"
    input:
      tuple val(support), path(vcf), path(index) 
    output:
    tuple val(support), path("${out}"), path("${out}.tbi") , emit: filtered
    publishDir "${params.out_dir}/", mode: 'copy', enabled: params.enabled
    script:
       out = "res${params.sub_tag}_${support}.vcf.gz"
       """
         bcftools view -T ^${exclude_list} $vcf -Oz -o $out
         tabix ${out}
       """
}


workflow tmerge {
    take: vcfs
    main:
      tmerge_step1_merge(vcfs) | tmerge_step2_filter
      tmerge_step3a_del_collapse(tmerge_step2_filter.out.del_ch)
      tmerge_step3b_oth_collapse(tmerge_step2_filter.out.oth_ch)
      tmerge_step4_merge((tmerge_step3a_del_collapse.out.mix(tmerge_step3b_oth_collapse.out)).groupTuple())
    emit:
      tmerge_step4_merge.out
}

workflow doSurvivor {
    take: grouped
    take: support_vals
    main:
      all_vcfs = grouped.map {  it -> it[1] }.\
          flatten().map { file -> [file.simpleName.replace("-"+params.sub_tag,""), file] }.\
          groupTuple()
      survivorByTool(all_vcfs, support_vals) | \
	  groupTuple | \
	  mergeSamples 
    emit:
       mergeSamples.out
}
