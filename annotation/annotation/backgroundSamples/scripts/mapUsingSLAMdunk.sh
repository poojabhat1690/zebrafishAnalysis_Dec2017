#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=100G
#SBATCH --time=0-8:00:00     # 2 minutes
#SBATCH --output=/scratch/pooja/runSLAMdunk
#SBATCH --job-name=runSLAMdunk
#SBATCH --array=1-87
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at






OUTDIR=/scratch/pooja/mapping/
mkdir -p "$OUTDIR"

arrayfile=`ls /scratch/pooja/mapping_dr10_backgroundCW/ | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`


#/groups/ameres/bioinformatics/tools/slamdunk_v.0.2.4/slamdunk/bin/slamdunk  all -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -o "$OUTDIR" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/output/final90percent/allAnnotations.bed  -fb /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/output/final90percent/countingWindows_transcriptionalOutput.bed -a 5 -5 12 -n 100 -t 15 -mq 0 -mi 0.95 -m -rl 88 /scratch/pooja/mapping_dr10_backgroundCW/"$arrayfile"


module load python/2.7.13-foss-2017a
OUTDIR=/scratch/pooja/mapping_lowExpression/
mkdir -p "$OUTDIR"

/groups/ameres/bioinformatics/tools/slamdunk_v.0.2.4/slamdunk/bin/slamdunk  all -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -o "$OUTDIR" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/countingWindows_lowexpression.bed  -fb /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/proteinCoding_annotatedUTRs.bed -a 5 -5 12 -n 100 -t 15 -mq 0 -mi 0.95 -m -rl 88 /scratch/pooja/mapping_dr10_backgroundCW/"$arrayfile"
