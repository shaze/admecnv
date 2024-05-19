#!/usr/bin/env nextflow


nextflow.enable.dsl=1

ref = file(params.ref, type:'file')
refi = file(params.fai, type:'file')


src = params.src




params.variant = 'diploidSV'
params.out     = "manta_out"



params.first_chrom=1
out = params.out




def getNodesOfBricks(fname) {
  cmd = "getfattr -n glusterfs.pathinfo -e text ${fname}";
  msg=cmd.execute().text;
  def matcher = msg =~ /(<POSIX.*)/;
  def bricks = matcher[0][0].tokenize(" ")
  nodes = []
  for (b : bricks ) {
    if (b =~ /.*arbiter.*/) continue
    matcher  = b =~ /.*:(.*):.*/; 
    node = matcher[0][1]
    matcher = node  =~ /(.*?)\..*/;
    if (matcher)
      node=matcher[0][1]
    nodes << node
  }
  return nodes
}


possible_states = ['idle','alloc','mix' ]
free_states = ['idle','mix']

def getStatus(nodes) {
  node_states ='sinfo -p batch -O NodeHost,StateCompact'.execute().text.split("\n")
  state_map = [:]
  possible  = []
  num_free  = 0
  for (n : node_states) {
    line=n.split()
    the_node=line[0]
    the_state=line[1]
    state_map[the_node]=the_state
    if  (the_state in possible_states) possible << the_node
    if  ( !(the_node in nodes)) continue;
    if  (the_state in free_states) num_free++;
  }
  return [num_free,possible]
}


def nodeOption(fname,aggression=1,other="") {
  nodes = getNodesOfBricks(fname)
  state = getStatus(nodes)
  possible=state[1]
  if ((possible.intersect(nodes)).size()<aggression)
    return "${other}"
  else {
    possible=possible - nodes;
    options="--exclude="+possible.join(',')+" ${other}"
    return options
  }
}



pairs = Channel.fromFilePairs("$src*{.bam,.bam.bai}", size:2) { it -> it.getName().replaceAll(/.bam.*/,"") }
	.map { [it[0],it[1][0], it[1][1], nodeOption(it[1][0])] }
        .randomSample(2000)


max_forks = params.max_forks


process runManta {
	cpus params.manta_parallel
        memory '5GB'
	tag {idSample}
	maxForks max_forks
	input:
	  file(ref) 
	  file(refi)
	  tuple idSample, file(bam), file(bai), val(cluster_opt) from pairs
        publishDir "manta_vcf", mode:'copy'
	output:
	  set idSample, file("${idSample}.vcf.gz"),
      	                file("${idSample}.vcf.gz.tbi") \
            into mantaOutput
       clusterOptions { cluster_opt }
	"""
	/opt/exp_soft/bioinf/manta/bin/configManta.py \
	    --bam ${idSample}.bam \
	    --referenceFasta ${ref} \
	    --runDir .
	./runWorkflow.py -j ${params.manta_parallel} -g 5
        cp results/variants/${params.variant}.vcf.gz ${idSample}.vcf.gz
        cp results/variants/${params.variant}.vcf.gz.tbi ${idSample}.vcf.gz.tbi
	"""
}

