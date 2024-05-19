
nextflow.enable.dsl=2

params.gluster_suggestion = ""
params.max_forks = 60

include {  initNodeSuggestion } from "../utils/gluster.nf"


params.ref = "Homo_sapiens_assembly38.fasta"

ref        = file(params.ref)
the_ref    = Channel.fromPath(params.ref)


print params.gluster_suggestion

bam = Channel.fromFilePairs(params.input, type:'any').randomSample(5000)
// .randomSample(5000) can be delteted. The 5000 is arbitrary and should be bigger than the
// number of BAMs. Our BAMs are in our gluster file systyem and are not always randomly distributed
// across bricks. Hence we randomise the inputs


// Following line can be omitted -- useful on our cluster to take work to node where
// BAM is stored rather than arbitrary node on the gluster system
// Need to delete clusterOptions if you delete this
node_suggestion = initNodeSuggestion(params.gluster_suggestion)


process setup {
    cpus 8
    memory  '6GB'
    storeDir params.store_dir
    input:
      path(ref)
    output:
      path("${ref}.{fasta,bwt,amb,ann,pac,sa,dict}"), includeInputs: true
      """
        gridss -s setupreference -r $ref
      """
}

process preprocess {
    maxForks 24
    cpus 10
    memory '12.GB'
    publishDir "preprocess"
    clusterOptions { node_suggestion[bam.getName()] }
    input:
      tuple val(base), path(bam), path(bai), path(dict) 
    output:
      tuple val(base), path(bam), path(bai), path("${bam}.gridss.working") 
    """
    hostname
    gridss -r $ref -t 10  -s preprocess $bam 
    """

}



process assemble {
    maxForks params.max_forks
    memory '32.GB'
    cpus 8
    input:
      tuple val(base), path(bam), path(bai), path(working)
    output:
      tuple val(sname), path(outbam), path(bam), path(bai),\
            path(working), path("assemble.bam.gridss.working")
    script:
    sname =  bam.simpleName
    outbam = "assemble.bam"
    """
      hostname
      gridss -s assemble -a assemble.bam  -t 8   -r $ref $bam
    """
}


process vcfcall {
    cpus 8
    maxForks 20
    memory '16.GB'
    errorStrategy 'finish'
    input:
      tuple val(base), path(abam), path(ibam), path(bai), path(wdir), path(adir)
    output:
      path "${vcfout}.gz"
    publishDir "gridss_vcf2"
    script:
    vcfout = "${base}.vcf"
    """
    gridss -t 8 -s call  -a $abam -r $ref $ibam --output $vcfout
    bgzip $vcfout
    """
} 

workflow {
    bam = Channel.fromFilePairs(params.input, type:'any')\
	   {it -> it.simpleName}.\
           map { it -> [it[0], it[1][0], it[1][1]] }
    setup(the_ref)
    dict=setup.out.flatten().filter( ~/.*.dict$/ )
    preprocess(bam.combine(dict)) |  assemble | vcfcall
}
