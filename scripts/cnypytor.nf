
nextflow.enable.dsl=2


include {  initNodeSuggestion } from "../utils/gluster.nf"


tag = params.tag

bam = Channel.fromFilePairs(params.input, type:'any') { it -> it.simpleName }.randomSample(5000)


// .randomSample(5000) can be delteted. The 5000 is arbitrary and should be bigger than the
// number of BAMs. Our BAMs are in our gluster file systyem and are not always randomly distributed
// across bricks. Hence we randomise the inputs


// Following line can be omitted -- useful on our cluster to take work to node where
// BAM is stored rather than arbitrary node on the gluster system
// Need to delete clusterOptions if you delete this
node_suggestion = initNodeSuggestion(params.gluster_suggestion)







def memOfBam(bam) {
  b=bam
  if (b.contains("X160") || b.contains("X1509") || b.contains("BMGF")) {
    return '90.GB'
  } else {
    return params.init_mem
  }
}




// This is the main process and us used to provide output to survivor

process cnvpytorProcess {
  scratch params.scratch_dir
  maxForks params.max_forks
  errorStrategy 'finish'
  cpus    params.init_cores
  memory  { memOfBam(base) }
  clusterOptions { node_suggestion[bam.getName()] }
  time '1d'
  input:
    tuple val(b), file(bam), file(bai) 
  output:
     tuple val(base), file(pytor), file(call), path(vcf)
  publishDir  "${params.out_dir}/calls"
  script:
     base=bam.simpleName
     pytor = "${base}.pytor"
     call  = "${base}.call"
     vcf = "${base}${tag}.vcf"
    """
     cnvpytor -j ${params.init_cores} -root $pytor -rd $bam
     cnvpytor -j ${params.init_cores} -root $pytor -his ${params.histos}
     cnvpytor -j ${params.init_cores} -root $pytor -partition ${params.histos}
     cnvpytor -j ${params.init_cores} -root $pytor -genotype ${params.histos}
     cnvpytor -j ${params.init_cores} -root $pytor -call ${params.histos} > $call         
     cnvpytor -root $pytor -view 1000 <<ENDL
      set Q0_range 0 0.5
      set size_range 500 inf
      set print_filename $vcf
      print calls
ENDL
    """
}

// This is used to create a relatively high quality merged data set that represents what would be produced
// if you only used cnvpyor

process mergeCalls {
    cpus 8
    memory "15.G"
    input:
     path(pytors)
    output:
     path(vcf)
    script:
      vcf = "merged.vcf"
    """
     cnvpytor -root $pytor -view ${params.histos} <<ENDL
        set Q0_range 0 0.5
        set p_range 0 0.0001
        set p_N 0 0.5
        set size_range 1000 inf
        set print_filename $vcf
        print calls
ENDL
    """
}



// cnvpytor doesn't produce compliant VCF. We ned to reheader the VCF so that all the chromosomes
// and contigs are in the header of the file (.vcf.gz produced to standard output) and second
// the CNs are sometimes fractional and so we need to allow for that 
process reheader {
    input:
    path(merged)
    output:
    path(out)
    script:
    out="reheadered.vcf.gz"
    """
      bcftools reheader  -f $ref   $merged |\
      zcat | sd "##FORMAT=<ID=CN,Number=1,Type=Integer," "##FORMAT=<ID=CN,Number=1,Type=Float," | bgzip > $out
    """
}




// The results need to be sorted
process sort {
    input:
    path(reheader)
    output:
    path(out)
    script:
    out="cnvpytor_merged.vcf.gz"
    """
      bcftools sort -Oz -o $out $reheader
    """
}



bam_pairs = Channel.fromFilePairs("${params.input}/*{.bam,.bai}", size:2) \
                { it.getName().replaceAll(".bam.*","")} \
		. map { it -> [it[0], it[1][0], it[1][1]] }


workflow  {
  main:
    cnvpytorProcess(bam_pairs) // | map { it -> it[1] } | collect  | mergeCalls | reheader | sort 
	      
}
