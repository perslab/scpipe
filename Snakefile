from glob import glob

"""
Dependencies:
  
  - Cellranger 
    Cmd: 
      `module use /tools/modules/
      module load cellranger/3.0.0`

  - fastqc

  - Kallisto, Alevin, Starsolo etc. 
"""

#configfile: "config.yaml"
OUT_DIR = "/home/cbmr/xbq246/20200117-dylan_pipeline/development/outs/gen_fastq_output/"
project = ["SCOP_9"]
run_name =  "190707_A00642_0027_AH7MWNDRXX"
DATA_DIR = ["/data/sc-seq"]
sample_name = ["ZA10","ZA12"]
read = ["R1","R2"]
lane = ["1","2"]
#FC_ID = ["H7MWNDRXX"]

dict_sample_name_s = {"ZA10": "S1", "ZA12": "S2"}
dict_sample_name_lane = {"ZA10": "L002", "ZA12": "L002"}

list_target_files = []

for s_name in sample_name:
  tmp = expand("/home/cbmr/xbq246/20200117-dylan_pipeline/development/outs/gen_fastq_output/{project}/{run_name}/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_{S_index}_L00{lane}_{read}_001.fastq.gz",
        project = project,
        run_name = run_name,
        sample_name = s_name,
        S_index = dict_sample_name_s[s_name],
        read = read,
        lane = 2) #should be lane = lane, but apparently lane = 1 doesn't get produced.
  list_target_files.extend(tmp)


rule all:
  input: 
    list_target_files,
    #expand("/home/cbmr/xbq246/20200117-dylan_pipeline/development/outs/gen_fastq_output/{project}/{run_name}/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S1_L00{lane}_{read}_001.fastq.gz",
    #  project = project,
    #  run_name = run_name,
    #  sample_name = sample_name,
    #  read = read,
    #  lane = lane), #rule mkfastq
    #expand("/home/cbmr/xbq246/20200117-dylan_pipeline/development/outs/gen_fastq_output/fastqc/{project}/{sample_name}_{lane}/fastqc/results.html",
    expand("/home/cbmr/xbq246/20200117-dylan_pipeline/development/outs/gen_fastq_output/fastqc/{project}/{run_name}/{sample_name}",
      project = project,
      run_name = run_name,
      sample_name = sample_name,
      lane = 2) #rule fastqc
  #Final outputs of Kallisto
  #expand("/projects/dylan/count_mats/SCOP_9/{sample_name}/counts_filtered/spliced.mtx", sample_name=sample_name),
  #expand("/projects/dylan/gen_fastq_output/fastqc/SCOP_9/{sample_name}_{lane}/fastqc/results.html", sample_name=sample_name, lane=lane)
  
  #expand("/projects/dylan/gen_fastq_output/fastqc/SCOP_9/{sample_name}_{lane}/fastqc/results.html", sample_name=sample_name, lane=lane)
  #expand("/projects/dylan/gen_fastq_output/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S1_L00{lane}_{read}_001.fastq.gz",
  #sample_name=sample_name, read=read, lane=lane)


"""
  📌 Dev notes 22/1 2 pm
  issue with Cellranger:
  The issue is that if a dir doesn't exist, then snakemake will create it
  cellranger needs to make its own output folder, otherwise it gets angry

  Reference thread: https://twitter.com/tangming2005/status/1056283068732981248
  Solution suggested by Johannes (not implemented yet): https://twitter.com/johanneskoester/status/1059532816537522176?s=20
    
    Implemented w.  a `rm -rf {params.OUT_DIR}{params.project}{params.run_name}` as per Johannes suggestion
    
  Alt solution: touch {output}
  The latter might confuse even more.
"""
rule mkfastq:
  input:
    #"/data/sc-seq/SCOP_9/190707_A00642_0027_AH7MWNDRXX" # SCOP_9 can be saved as project ID #190707... can be saved as run ID
    expand("{DATA_DIR}/{project}/{run_name}", 
      DATA_DIR = DATA_DIR,
      project = project,
      run_name = run_name)
  output:
    list_target_files
  #"""
  #
  #  📌 Dev note 23/1 1 am
  #
  #  Issue with file name matching:
  #  Apparently for ZA10 it's S1, while ZA12 it's S2
  #  
  #  ZA10_S1_L002_I1_001.fastq.gz  ZA10_S1_L002_R1_001.fastq.gz  ZA10_S1_L002_R2_001.fastq.gz
  #  ZA12_S2_L002_I1_001.fastq.gz  ZA12_S2_L002_R1_001.fastq.gz  ZA12_S2_L002_R2_001.fastq.gz
  #
  #  Also, not entirely sure how lane 1 (L001) comes into this picture. 
  #  Dylan defines a lane 1 and lane 2, but really the output doesn't have a lane 1
  # 
  #  This will probably call for some conditional string matching
  #
  #
  #  📌 Dev note 23/1 7:30 pm
  #
  #  I suspect that a lot of hassle with. straing matching can be avoided if we simply opt for using
  #  the dir that Cellranger makes locally (something something ./{run_name}/ and __{run_name}.mro). 
  #  Talk with Dylan about this 💬
  #
  #"""
  params:
    OUT_DIR = OUT_DIR,
    run_name = run_name,
    project = project #will be put in a lambda w. card
  shell:
    "rm -rf {params.OUT_DIR}{params.project}/{params.run_name} &&\
    module load cellranger/3.0.0 &&\
    cellranger mkfastq \
    --id {params.run_name} \
    --csv {input}/cellranger_samplesheet.csv \
    --run {input} \
    --output-dir {params.OUT_DIR}{params.project}/{params.run_name}/outs/fastq_path/ &&\
    mkdir {params.OUT_DIR}{params.project}/{params.run_name}/tmp &&\
    mv {params.run_name} {params.OUT_DIR}{params.project}/{params.run_name}/tmp && mv __{params.run_name}.mro {params.OUT_DIR}{params.project}/{params.run_name}/tmp"
 
#
#  📌 Dev note 22/1, 1pm: I wanted to label these files as temp() files in the output, but then that apparently messes up the wildcards in rule mkfastq
#  input:
#    expand("{run_name}/", run_name = run_name),
#    expand("__{run_name}.mro", run_name = run_name)
#  shell:
#    "rm -rf {input} | echo 'deleted tmp files"
#
#  📌📌 Dev note 23/1, 8pm: These files might actually be useful! Maybe design might be much simpler. 


rule cat_fastq:
  """ 📝 Concatenates multiple .fastq.gz files into one big fastq file.

  Globs all files w. same sample_name and same read id, but that might have different lanes.
  Combined using simple `cat` in shell.
  
  Input: 
    read1 = glob(.../{sample_name}_S[0-9]*_L[0-9]*_R1_001.fastq.gz)
    read2 = glob(.../{sample_name}_S[0-9]*_L[0-9]*_R2_001.fastq.gz)

  Output:
    {sample_name}__concat__R1_001.fastq.gz
    {sample_name}__concat__R2_001.fastq.gz

  Example:
    Input: 
      read1 = ['.../ZA10_S1_L1_R1_001.fastq.gz', '.../ZA10_S1_L2_R1_001.fastq.gz', '.../ZA10_S1_L3_R1_001.fastq.gz']
      read2 = ['.../ZA10_S1_L1_R2_001.fastq.gz', '.../ZA10_S1_L2_R2_001.fastq.gz']
  
    Output:
      .../ZA10__concat__R1_001.fastq.gz
      .../ZA10__concat__R2_001.fastq.gz
  """
  input:
    read1 = lambda wildcards: glob(OUT_DIR + "{project}/{run_name}/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S[0-9]*_L[0-9]*_R1_001.fastq.gz".format(\
      project = wildcards.project,
      run_name = str(wildcards.run_name),
      sample_name = wildcards.sample_name)),
    read2 = lambda wildcards: glob(OUT_DIR + "{project}/{run_name}/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}_S[0-9]*_L[0-9]*_R2_001.fastq.gz".format(\
      project = wildcards.project,
      run_name = str(wildcards.run_name), 
      sample_name = wildcards.sample_name))
  output:
    cat_R1 = OUT_DIR + "{project}/{run_name}/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}__concat__R1_001.fastq.gz",
    cat_R2 = OUT_DIR + "{project}/{run_name}/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}__concat__R2_001.fastq.gz"
  run:
    if len(input.read1) > 0 and len(input.read2) > 0: 
      r1 = " ".join(input.read1)
      r2 = " ".join(input.read2)
      shell("touch {output.cat_R1}")
      shell("cat " + r1 + " > {output.cat_R1}")
      shell("touch {output.cat_R2}")
      shell("cat " + r2 + " > {output.cat_R2}")



rule fastqc:
  """ 📝 Runs fastqc w. the concatenated .fastqc.gz files
   
  Input: 
   read1 = .../*__concat__R1_001.fastq.gz file from rule cat_fastq
   read2 = .../*__concat__R2_001.fastq.gz file from rule cat_fastq
   
   wildcards: {project}, 
               {run_name},
               {sample_name},
   
   lambda is mostly for ensureing that we stick to 1 {sample_name} during pr. rule execution. 
   E.g. so that we don't do QC where we match different samples such as ZA10__concat__R1_001.fastq.gz and ZA12__concat__R2_001.fastq.gz
  
  Output: 
   fastqc directory
  
  Example:
   Input:
     read1 = <BASE_OUT_DIR> + '/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/ZA10/ZA10__concat__R1_001.fastq.gz'
     read2 = <BASE_OUT_DIR> + '/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/ZA10/ZA10__concat__R2_001.fastq.gz'
  
   Output:
     <BASE_OUT_DIR> + '/fastqc/SCOP_9/190707_A00642_0027_AH7MWNDRXX/ZA10/''
  
   shell:
     "module load fastqc/0.11.5 &&\
     mkdir <OUT_DIR>/fastqc/SCOP_9/190707_A00642_0027_AH7MWNDRXX/ZA10/ &&\
     fastqc --outdir=<BASE_OUT_DIR>/fastqc/SCOP_9 \
     '<BASE_OUT_DIR>/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/ZA10/ZA10__concat__R1_001.fastq.gz' \
     '<BASE_OUT_DIR>/SCOP_9/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/ZA10/ZA10__concat__R2_001.fastq.gz'
  """
  input:
    read1 = lambda wildcards: OUT_DIR + "{project}/{run_name}/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}__concat__R1_001.fastq.gz".format(\
      project = wildcards.project,
      run_name = str(wildcards.run_name),
      sample_name = wildcards.sample_name) ,
    read2 = lambda wildcards: OUT_DIR + "{project}/{run_name}/outs/fastq_path/H7MWNDRXX/{sample_name}/{sample_name}__concat__R2_001.fastq.gz".format(\
      project = wildcards.project,
      run_name = str(wildcards.run_name),
      sample_name = wildcards.sample_name)
  output:
    directory(OUT_DIR + "fastqc/{project}/{run_name}/{sample_name}")
  #"/projects/dylan/gen_fastq_output/fastqc/{project}/{sample_name}_{lane}/fastqc/results.html"
  params:
    project = project,
    run_name = run_name,
    OUT_DIR = OUT_DIR,
    sample_name = lambda wildcards: wildcards.sample_name
  shell: 
    "module load fastqc/0.11.5 &&\
    mkdir {params.OUT_DIR}/fastqc/{params.project}/{params.run_name}/{params.sample_name} &&\
    fastqc --outdir={params.OUT_DIR}/fastqc/{params.project}/{params.run_name}/{params.sample_name} {input.read1} {input.read2}"


rule kallisto:
  input:
    "placeholder"
  output:
    "/projects/dylan/rausch_abay_vlmc/data/kallisto/ZA10/"
  shell:
    "{Kallisto} bus \
    -i /home/cbmr/jph712/projects/serup_velocytoLoompy/velocity_getting_started/index.idx \
    -o /projects/dylan/rausch_abay_vlmc/data/kallisto/ZA10/ \
    -x 10xv3 \
    -t 50 \
    /data/sc-10x/data-mkfastq/190707_A00642_0027_AH7MWNDRXX/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/ZA10/ZA10_S1_L002_R1_001.fastq.gz /data/sc-10x/data-mkfastq/190707_A00642_0027_AH7MWNDRXX/190707_A00642_0027_AH7MWNDRXX/outs/fastq_path/H7MWNDRXX/ZA10/ZA10_S1_L002_R2_001.fastq.gz"
    """
    📌 Dev Notes:
      Lots of things that can be put into a config
      doc: https://pachterlab.github.io/kallisto/manual#bus
      
      {Kallisto} = /tools/anaconda/envs/lhv464/kallisto/lib/python3.7/site-packages/kb_python/bins/linux/kallisto/kallisto
      bus           Generate BUS files for single-cell data
      -i, --index=STRING          Filename for the kallisto index to be constructed
      -o, --output-dir=STRING       Directory to write output to
      -x, --technology=STRING       Single-cell technology used 
      -t, --threads=INT             Number of threads to use (default: 1)
      <FASTQ files>                 ex. ZA10_S1_L002_R1_001.fastq.gz  ZA10_S1_L002_R2_001.fastq.gz
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

#rule kite_count:
# input:
# output:
# run:


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

