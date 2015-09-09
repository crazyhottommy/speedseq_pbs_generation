__author__ = 'Ming Tang'

import csv
import re

datadir = '/scratch/bcb/fpbarthel/TCGA/WGS/GBM/'
ref = '/Users/Tammy/projects/verhaak_lab/GBM_structural_variants/data/2015-06-23-Floris_bamfiles_info/GBM_bamfiles.csv'
results_dir = '/scratch/genomic_med/mtang1/scratch/results/'

# the bamfiles.csv file from Floris contains the tumor and the paired normal bam file info
# the UUID in column 5 and column 8 are the folder names containing the bam file in HPC cluster

# keep the normal, tumor barcode in two lists
normal_barcode = []
tumor_barcode = []

with open(ref, "r") as f:
    next(f) # skip the header
    reader = csv.reader(f, delimiter=',')
    ## need a backslash to escape a literal backslash \\
    for row in reader:
        normal_barcode.append(row[6])
        tumor_barcode.append(row[0])
    # the normal bam has duplicates, because there are primary and recurrent tumor matching it
    # return a unique and sorted list
    normal_barcode = sorted(set(normal_barcode))
    # just in case
    tumor_barcode = sorted(set(tumor_barcode))


# make a dictionary with key is the normal, value is a list of [primary, recurrent]
normal_tumor_pairs = {}

pattern = r'(TCGA-\d+-\d+)(.+)'

for normal in normal_barcode:
    result = re.search(pattern, normal)
    case = result.group(1)
    primary_pattern = r'(%s-01\w)-.+' % case   # 01 is the primary tumor
    recurrent_pattern = r'(%s-02\w)-.+' % case # 02 is the recurrent tumor
    tumor_list = []
    for tumor in tumor_barcode:
        primary = re.search(primary_pattern, tumor)
        recurrent = re.search(recurrent_pattern, tumor)
        if primary:
            tumor_list.append(primary.group(0))
        if recurrent:
            tumor_list.append(recurrent.group(0))
    print normal, tumor_list
    normal_tumor_pairs[normal] = tumor_list

########################### How many cases for normal, primary, recurrent triplet
n=0
for normal, tumor in normal_tumor_pairs.items():
    if len(tumor) ==2:
        print tumor
        n+=1
print "there are %d triplets" %n
###########################


##### generate the pbs files

for normal, tumor in normal_tumor_pairs.items():
    if len(tumor) == 1:
        job_name = normal[0:12]
        normal_bam = normal
        primary_bam = tumor[0]
        job_string = """
#PBS -N %s
#PBS -l nodes=1:ppn=6,walltime=12:00:00
#PBS -l mem=18g
#PBS -M mtang1@mdanderson.org
#PBS -m abe
#PBS -d /scratch/genomic_med/mtang1/scratch/results/%s
#PBS -o /scratch/genomic_med/mtang1/scratch/results/logs
#PBS -e /scratch/genomic_med/mtang1/scratch/results/logs
#PBS -V

speedseq sv \\
    -o %s \\
    -x /scratch/genomic_med/mtang1/scratch/LOWPASS_WGS_LGG1/ref_genome/ceph18.b37.lumpy.exclude.2014-01-15.bed \\
    -B %s.bam,%s.bam \\
    -D %s.discordants.bam,%s.discordants.bam \\
    -S %s.splitters.bam,%s.splitters.bam \\
    -R /scratch/genomic_med/mtang1/scratch/LOWPASS_WGS_LGG1/ref_genome/human_g1k_v37.fasta \\
    -g \\
    -d \\
    -t 6 \\
    -T temp
        """ % (job_name, job_name, job_name,
               results_dir + normal_bam + "/" + normal_bam, results_dir + primary_bam + "/" + primary_bam,
               results_dir + normal_bam + "/" + normal_bam, results_dir + primary_bam + "/" + primary_bam,
               results_dir + normal_bam + "/" + normal_bam, results_dir + primary_bam + "/" + primary_bam)
        ofile= open("pbs_dir/%s_sv.pbs" % job_name, "w")
        ofile.write(job_string)
        ofile.close()
    elif len(tumor) == 2:
        job_name = normal[0:12]
        normal_bam = normal
        primary_bam = tumor[0]
        recurrent_bam = tumor[1]
        job_string = """
#PBS -N %s
#PBS -l nodes=1:ppn=6,walltime=18:00:00
#PBS -l mem=18g
#PBS -M mtang1@mdanderson.org
#PBS -m abe
#PBS -d /scratch/genomic_med/mtang1/scratch/results/%s
#PBS -o /scratch/genomic_med/mtang1/scratch/results/logs
#PBS -e /scratch/genomic_med/mtang1/scratch/results/logs
#PBS -V

speedseq sv \\
    -o %s \\
    -x /scratch/genomic_med/mtang1/scratch/LOWPASS_WGS_LGG1/ref_genome/ceph18.b37.lumpy.exclude.2014-01-15.bed \\
    -B %s.bam,%s.bam,%s.bam \\
    -D %s.discordants.bam,%s.discordants.bam,%s.discordant.bam \\
    -S %s.splitters.bam,%s.splitters.bam,%s.splitters.bam \\
    -R /scratch/genomic_med/mtang1/scratch/LOWPASS_WGS_LGG1/ref_genome/human_g1k_v37.fasta \\
    -g \\
    -d \\
    -t 6 \\
    -T temp
        """ % (job_name, job_name, job_name,
               results_dir + normal_bam + "/" + normal_bam, results_dir + primary_bam + "/" + primary_bam, results_dir + recurrent_bam + "/" + recurrent_bam,
               results_dir + normal_bam + "/" + normal_bam, results_dir + primary_bam + "/" + primary_bam, results_dir + recurrent_bam + "/" + recurrent_bam,
               results_dir + normal_bam + "/" + normal_bam, results_dir + primary_bam + "/" + primary_bam, results_dir + recurrent_bam + "/" + recurrent_bam)
        ofile= open("pbs_dir/%s_sv.pbs" % job_name, "w")
        ofile.write(job_string)
        ofile.close()
