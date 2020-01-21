#configfile: "config.yaml"
#proj = ["SCOP_9"]
#run_name =  ["190707_A00642_0027_AH7MWNDRXX"]
sample_name = ["ZA10","ZA12"]
read = ["R1","R2"]
lane = ["1","2"]
#FC_ID = ["H7MWNDRXX"]

rule all:
 input:
  expand("/projects/dylan/count_mats/SCOP_9/{sample_name}/counts_filtered/spliced.mtx", sample_name=sample_name),
  expand("/projects/dylan/gen_fastq_output/fastqc/SCOP_9/{sample_name}_{lane}/fastqc/results.html", sample_name=sample_name, lane=lane)
  #expand("/projects/dylan/gen_fastq_output/fastqc/SCOP_9/{sample_name}_{lane}/fastqc/results.html", sample_name=sample_name, lane=lane)
  #expand("/projects/dylan/gen_fastq_output/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S1_L00{lane}_{read}_001.fastq.gz",
  #sample_name=sample_name, read=read, lane=lane)

rule mkfastq:
 input:
  "/data/sc-seq/SCOP_9/190707_A00642_0027_AH7MWNDRXX" # SCOP_9 can be saved as project ID #190707... can be saved as run ID
 output:
  "/projects/dylan/gen_fastq_output/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S1_L00{lane}_{read}_001.fastq.gz"
 shell:
  """
  rm -rf /projects/dylan/gen_fastq_output/SCOP_9
  mkdir /projects/dylan/gen_fastq_output/SCOP_9
  cd /projects/dylan/gen_fastq_output/SCOP_9
  cellranger mkfastq --id 190707_A00642_0027_AH7MWNDRXX --csv {input}/cellranger_samplesheet.csv --run {input}
  """

rule cat_fastq:
 input:
  "/projects/dylan/gen_fastq_output/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S1_L00{lane}_{read}_001.fastq.gz"
 output:

 shell:
  """

  """

rule fastqc:
 input:
  read1="/projects/dylan/gen_fastq_output/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S1_L00{lane}_R1_001.fastq.gz",
  read2="/projects/dylan/gen_fastq_output/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S1_L00{lane}_R2_001.fastq.gz"
 output:
  "/projects/dylan/gen_fastq_output/fastqc/SCOP_9/{sample_name}_{lane}/fastqc/results.html"
 shell:
  """
  mkdir -p /projects/dylan/gen_fastq_output/fastqc/SCOP_9
  fastqc --outdir="/projects/dylan/gen_fastq_output/fastqc/SCOP_9" {input.read1} {input.read2}
  """

rule kallisto_count:
 input:
  expand("/projects/dylan/gen_fastq_output/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S1_L00{lane}_{read}_001.fastq.gz",
  lane=lane, read=read, sample_name=sample_name)
 output:
  "/projects/dylan/count_mats/SCOP_9/{sample_name}/counts_filtered/spliced.mtx"
 run:
  "kb count --h5ad -i index.idx -g t2g.txt -x run_chemistry -o {outfolder} \
  -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 2 \
  {input}"

rule kite_count:
 input:
 output:
 run:


# #
# # input:
# #  expand("/projects/dylan/gen_fastq_output/{proj}/{run_name}/outs/fastq_path/{FC_ID}/{sample_name}/{sample_name}_S1_LOO{lane}_{read}_001.fastq.gz", read=read)
# # output:
# #  expand("/projects/dylan/count_mats/{proj}/{sample_name}", sample_name=sample_name)
# # run:
# #   "kb count --h5ad -i {index}.idx -g {t2g}.txt -x {run_chemistry} -o {outfolder} \
# #   -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 2 \
# #  SRR6470906_S1_L001_R1_001.fastq.gz \
# #   SRR6470906_S1_L001_R2_001.fastq.gz \
# #   SRR6470906_S1_L002_R1_001.fastq.gz \
# #   SRR6470906_S1_L002_R2_001.fastq.gz"
