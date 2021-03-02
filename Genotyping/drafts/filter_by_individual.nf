#!/usr/bin/env nextflow

// Input files
params.input_vcf = "input/vcf_by_sample/*.vcf.gz"
params.shuffle   = "tmp/shuffle"

// Freebayes parameters
params.mincov = 2

// Output files
params.out_vcf = "results/filtered_vcf_by_sample"

// Environment
params.env = "env.sh"

Channel
  .fromPath(params.input_vcf)
  .map { file -> tuple(file.getSimpleName(), file) }
  .set { vcf_channel }

process run_filter {

  publishDir params.out_vcf, mode: 'move'

  errorStrategy 'retry'

  maxRetries 2

  scratch true

  stageOutMode 'move'

  executor 'sge'
  penv 'smp'
  maxForks 90
  cpus 1
  time {12.hour * ((2 * task.attempt) -  1)}
  clusterOptions "-l h_vmem=${ [30, 60, 120][ task.attempt - 1] }G"

  input:
  file env from file(params.env)
  file shuffle from file(params.shuffle)
  set sample_id, file(vcf) from vcf_channel
  val mincov from params.mincov

  output:
  file "results_${sample_id}" into final_vcf

  """

  source ${env}
  tabix -p vcf $vcf
  $baseDir/SNP-genotype-filtering.v2.HTZ2Random.MT.zsh \
    $vcf $mincov $sample_id ${shuffle} results_${sample_id}

  """
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
