import sys

t_file = sys.argv[1]
n_file = sys.argv[2]
out_path = sys.argv[3]

header = "group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath\n"
group_id = t_file.split("/")[-4]
subject_id = group_id.replace("P-","")

sample_id_tumor="SARC-"+group_id+"-T-"+subject_id
sample_id_normal="SARC-"+group_id+"-N-"+n_file.split("/")[-3].split("-")[10]
n_type=",normal,dna,bam"
t_type=",tumor,dna,bam"

onco_file = out_path+group_id+".csv"

oncoanalyzer_csv=open(onco_file,'w')
oncoanalyzer_csv.write(header)
content_n=group_id+","+subject_id+","+sample_id_normal+n_type+","+n_file+"\n"
content_t=group_id+","+subject_id+","+sample_id_tumor+t_type+","+t_file+"\n"

oncoanalyzer_csv.write(content_n)
oncoanalyzer_csv.write(content_t)
oncoanalyzer_csv.close()

slurm_file_name=out_path+group_id+".sh"
slurm_file=open(slurm_file_name,'w')
slurm_content1="#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --tasks-per-node=25\n#SBATCH --job-name=oncoanalyzer_"+group_id+"\n#SBATCH --time=5-15:00:00\n"
slurm_content2="nextflow run /nfs/PIPELINE/oncoanalyser/main.nf -profile singularity -config /nfs/PIPELINE/nf_references/GRCh37/reference.config --mode wgts --genome CustomGenome --genome_type no_alt --genome_version 37 --force_genome --ref_data_genome_gtf /nfs/PIPELINE/nf_references/GRCh37/2.7.3a/gencode.v38.annotation_corrected.gtf --input "+onco_file+" --outdir "+out_path
slurm_file.write(slurm_content1)
slurm_file.write(slurm_content2)
slurm_file.close()
