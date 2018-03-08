#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=50G
#SBATCH --time=0-8:00:00     # 2 minutes
#SBATCH --output=/scratch/pooja/runSLAMdunk
#SBATCH --job-name=runSLAMdunk
#SBATCH --array=1-87
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at


module load cutadapt/1.9.1-foss-2017a-python-2.7.13

OUTDIR=/scratch/pooja/mapping_dr10_backgroundCW/
mkdir -p "$OUTDIR"

arrayfile=`ls /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/allData/quantseq/ | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

INDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/allData/quantseq/
module load cutadapt/1.9.1-foss-2017a-python-2.7.13

#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -o "$OUTDIR"/"$arrayfile"_adapterTrimmed.fastq -m 18 --trim-n "$INDIR"/"$arrayfile"

module load python/2.7.13-foss-2017a
#pip install slamdunk
#export PATH=$PATH:$HOME/.local/bin/"

OUTDIR=/scratch/pooja/mapping/
mkdir -p "$OUTDIR"



#/groups/ameres/bioinformatics/tools/slamdunk_v.0.2.4/slamdunk/bin/slamdunk  all -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -o "$OUTDIR" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countingWindoes_ensemblAnnotation.bed  -fb /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/proteinCoding_annotatedUTRs.bed -a 5 -5 12 -n 100 -t 15 -mq 0 -mi 0.95 -m -rl 88 /scratch/pooja/mapping_dr10_backgroundCW/combinedFile_TTTATG.fastq.gz_adapterTrimmed.fastq

OUTDIR=/scratch/pooja/mapping_lowExpression/
mkdir -p "$OUTDIR"

/groups/ameres/bioinformatics/tools/slamdunk_v.0.2.4/slamdunk/bin/slamdunk  all -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -o "$OUTDIR" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/countingWindows_lowexpression.bed  -fb /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/proteinCoding_annotatedUTRs.bed -a 5 -5 12 -n 100 -t 15 -mq 0 -mi 0.95 -m -rl 88 /scratch/pooja/mapping_dr10_backgroundCW/combinedFile_TTTATG.fastq.gz_adapterTrimmed.fastq
