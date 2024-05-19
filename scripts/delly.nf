#!/usr/bin/env nextflow


nextflow.enable.dsl=2




ref  =file(params.ref, type:'file')

map_file = file(params.ref, type:'file')


out_dir=params.out_dir
max_forks=params.max_forks




key_fnames = file(params.gluster_suggestion)




include {  initNodeSuggestion } from "../utils/gluster.nf"





bams = Channel.fromFilePairs(params.input) { it -> it.simpleName }.map { [it[0],it[1][0], it[1][1]] }.randomSample(2000)



node_suggestion = initNodeSuggestion(params.gluster_suggestion)


process delly {
    tag { sample }
    memory '9 GB'
    maxForks max_forks
    clusterOptions { node_suggestion[bam.getName()] }
    cpus 1
    input:
       tuple  val(sample), path(bam), path(bai) 
    output:
      tuple path("${bam.baseName}.bcf"), path("${bam.baseName}.bcf.csi") 
  script:
    """
    export omp_num_threads=1
    delly call \
        --outfile ${bam.baseName}.bcf \
        --genome ${ref} \
        ${bam}
    """
}


process merge_sites {
  memory '32 GB'
  cpus 16
  input:
    path(bcfs) 
  output:
    path("sites.bcf") 
  script:
     """
     export omp_num_threads=16
     delly merge -o sites *bcf
     mv sites sites.bcf
     """
}



process genotype {
  memory '9GB'
  clusterOptions { node_suggestion[bam.getName()] }  
  cpus 1
  input:
    tuple  path(sites), val(sample), path(bam), path(ind)
  output:
    path("${sample}_geno.bcf*") 
  publishDir "${params.out_dir}/indivs"
  script:
  """
     export omp_num_threads=1
     delly call  -v $sites  -g ${ref}  -o ${sample}_geno.bcf $bam
  """
}


process merge_genos {
  maxForks params.max_forks
  memory '24GB'
  cpus 2
  input:
  path(bcfs) 
  output:
    tuple path("merged.bcf"), path("merged.bcf.csi") 
  """
  export omp_num_threads=2
  bcftools merge -m id -Ob -o merged *bcf
  mv merged merged.bcf
  bcftools index merged.bcf
  """
}



process filterSV {
   memory '16GB'
   cpus 4
  input:
    tuple path(unfiltered), path(unfiltered_index) 
  output:
    tuple  path("delly.bcf"), path("delly.bcf.csi")
  publishDir params.out_dir, mode: "copy"
  """
  delly filter -f germline  -o delly.bcf $unfiltered
  """
}


process splitSV {
   memory '16GB'
   cpus 4
  input:
    tuple path(unfiltered), path(unfiltered_index) 
  output:
     path("indivs/*")
  publishDir params.indiv_out_dir, mode: "copy"
  """
  mkdir indivs
  bcftools +split $delly -Oz  -o indivs
  """
}


process cnvIndiv {
  input:
    tuple val(sample), path(bcf_sv), path(sv_index), path(bam), path(bai)
  output:
    tuple path("${sample}_cnv.bcf"), path("${sample}_cnv.bcf.csi")   
  """
     delly cnv -o ${sample}_cnv.bcf -g $ref -m $map_file  -l $bcf_csv  $bam
  """

}


process mergeCNVSites {
  maxForks params.max_forks
  memory '24GB'
  cpus 2
  input:
  path(bcfs) 
  output:
    tuple path("merged_cnv.bcf"), path("merged_cnv.bcf.csi") 
  """
  export omp_num_threads=2
  delly merge -e -p -o merged_cnv.bcf -m 1000 -n 100000 *bcf
  bcftools index merged_cnv.bcf
  """
}


process genotypeCNV {
  memory '9GB'
  clusterOptions { node_suggestion[bam.getName()] }  
  cpus 1
  input:
    tuple  path(sites), val(sample), path(bam), path(ind)
  output:
    path("${sample}_geno_cnv.bcf*") 
  publishDir "${params.out_dir}/indivs"
  script:
  """
     export omp_num_threads=1
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
    tuple path("merged.bcf"), path("merged.bcf.csi") 
  """
  export omp_num_threads=2
  bcftools merge -m id -o b -o merged_cnv.bcf $bcfs
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
  publishDir params.out_dir, mode: "copy"
  """
  delly classify -f germline -o delly_cnv.bcf $unfiltered
  bcftools index delly.bcf
  """
}

  

workflow {
    main:
    delly(bams)  | flatten | collect | \
	merge_sites  | combine(bams) | \
	genotype | collect | \
	merge_genos | \
	filterSV | \
	splitSV

}
    
   




