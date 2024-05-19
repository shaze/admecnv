#!/usr/bin/env nextflow


nextflow.enable.dsl=2




ref  =file(params.ref, type:'file')

map_file = file(params.map_file, type:'file')


out_dir=params.out_dir
max_forks=params.max_forks




key_fnames = file(params.gluster_suggestion)




include {  initNodeSuggestion } from "../utils/gluster.nf"





bams = Channel.fromFilePairs(params.input) { it -> it.simpleName }.map { [it[0],it[1][0], it[1][1]] }



vcfs = Channel.fromFilePairs("${params.indiv_vcf_dir}/*.vcf.gz*")  { it -> it.simpleName }
              .map { [it[0],it[1][0], it[1][1]] }




node_suggestion = initNodeSuggestion(params.gluster_suggestion)


process cnvIndiv {
  errorStrategy 'finish'
  cpus 4
  clusterOptions { node_suggestion[bam.getName()] }
  input:
    tuple val(sample), path(bam), path(bai), path(sv_vcf), path(sv_index)
  output:
    tuple path("${sample}_cnv.bcf"), path("${sample}_cnv.bcf.csi")
  """
     hostname 
     /usr/bin/time -f "%e %M" delly cnv -o ${sample}_cnv.bcf -g $ref -m $map_file  -l $sv_vcf  $bam
  """

}


process mergeCNVSites {
  maxForks params.max_forks
  memory '24GB'
  cpus 2
  input:
  path(bcfs) 
  output:
    tuple path("merged_cnv.bcf"),  path("merged_cnv.bcf.csi")
script: 
  """

  export omp_num_threads=2
  delly merge -e -p -o merged_cnv.bcf -m 1000 -n 100000 *bcf
  
  export omp_num_threads=2

  
  """
}


process genotypeCNV {
  memory '9GB'
  maxForks 40
  clusterOptions { node_suggestion[bam.getName()] }  
  cpus 4
  input:
    tuple  path(sites), path(siteidx), val(sample), path(bam), path(ind)
  output:
    path("${sample}_geno_cnv.bcf*") 
  publishDir  "${params.out_cnv_dir}/indivs"
  script:
  """
     export omp_num_threads=2
     delly cnv -u -v $sites -g $ref -m $map_file -o ${sample}_geno_cnv.bcf $bam
  """
}


process mergeCNV {
  maxForks params.max_forks
  memory '24GB'
  cpus 2
  input:
  path(bcfs) 
  output:
    tuple path("merged_cnv.bcf"), path("merged_cnv.bcf.csi") 
  """
  export omp_num_threads=2
  bcftools merge -m id -O b -o merged_cnv.bcf *.bcf
  bcftools index merged_cnv.bcf
  """
}



process filterCNV {
   memory '16GB'
   cpus 4
  input:
    tuple path(unfiltered), path(unfiltered_index) 
  output:
    tuple  path("delly_cnv.bcf"), path("delly_cnv.bcf.csi")
  publishDir params.out_cnv_dir, mode: "copy"
  """
  delly classify -f germline -o delly_cnv.bcf $unfiltered
  """
}



process splitCNVs {
    input:
      tuple path(bcf), path(index)
    output:
      path("splits/*vcf.gz")
    publishDir "${params.out_cnv_dir}/cnv_refined_indivs"
    script:
    """
      mkdir splits
      bcftools +split -Oz -o splits $bcf
    """
}


process stitchCNVs {
    maxForks 20
    input:
      path(vcfgz)
    output:
      path(out)
    publishDir "${params.out_cnv_dir}/stitched_cnvs"
    script:
    out="${vcfgz.simpleName}.vcf"
    """
       stitch.py $vcfgz $out
    """
}


workflow {
    main:
    inputCh = bams.join(vcfs).randomSample(2000)
    cnvIndiv(inputCh) | collect |  \
	mergeCNVSites | combine(bams) | randomSample(10000) | \
	genotypeCNV | collect | \
	mergeCNV | \
	filterCNV  | \
	splitCNVs | flatten | \
	stitchCNVs
}
    
   




