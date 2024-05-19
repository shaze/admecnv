

nextflow.enable.dsl=2
ref=file(params.ref)
excludes = file (params.excludes)
out = params.out

bams = Channel.fromFilePairs(params.input) { it -> it.simpleName }.map { [it[0],it[1][0], it[1][1]] }.randomSample(2000)




include {  initNodeSuggestion } from "../utils/gluster.nf"





process smoove {
    cpus 2
    maxForks 15
    clusterOptions = '--constraint=avx2'
    input:
       tuple  val(sample), path(bam), path(bai) 
    output:
      path("*vcf.gz*")
    publishDir params.output
    script:
      """
      smoove call --outdir ./ --exclude $excludes --name $sample\
	    --fasta $ref -p 1 --genotype *bam
      """
}


process mergeSites {
  label 'singularity'
  input:
     path(vcfs)
  output:
     path("res/*vcf.gz*")
  script:
     """
     hostname
     mkdir res
     smoove merge --name merged -f $ref --outdir res *vcf.gz 
     """
}

process smooveGT {
    label 'singularity'
    maxForks 40
    cpus 2
    input:
       tuple path(merged), val(sample), path(bam), path(bai),\
             path(the_ref), path(the_index)
     output:
       path("res/*")
     script:
      
      """
      mkdir res
      smoove genotype -d -x -p 1 --name $out --outdir   res  --fasta $the_ref --vcf $merged $bam
      """
} 


process paste {
    label 'singularity'
    input:
      path("sample????.vcf.gz")
    output:
       path(vcfgz)
    script:
       vcfgz="combine.vcf.gz"       
       index = "${vcfgz}.tbi"
       """
       for x in *.vcf.gz; do 
           bcftools index \$x;
       done
       smoove paste --name combine *.vcf.gz
       """
}


process annotate {
    cpus 2
    input:
      tuple path(gff), path(vcfgz)
    output:
      tuple path(out), path("${out}.tbi")
    publishDir params.out_dir, mode: 'copy'
    script:
      out = "smoove.vcf.gz"
    """
     smoove annotate --gff ${params.out_dir} $vcfgz| \
     bgzip -c > $out
     tabix $out
    """
}


process split {
    input:
      tuple path(vcfgz), path(index)
    output:
      path("raw/*vcf.gz")
    script:
    """
      mkdir raw
      bcftools +split $vcfgz -Ob -o raw
    """
}


process checkQuality {
    maxForks 12
    input:
      path(raw)
    output:
      path(good)
    publishDir "${params.out_dir}/filtered"
    script:
     good = raw.simpleName+".vcf.gz"
    
    """
       bcftools f$raw  -i "FORMAT/SHQ==4 | FORMAT/GQ==1/1"" -o test.vcf.gz -Oz 
       stitch.py -a test.vcf.gz - | bgzip -c > $good
    """
}


workflow {
    refs= Channel.fromFilePairs("${params.ref}{,.fai}").map { it[1]}
    smoove(bams) | collect | mergeSites | combine(bams) | combine(refs) |
        smooveGT  | collect |
	paste |
	annotate |
	split |
	flatten |
	checkQuality 
}
