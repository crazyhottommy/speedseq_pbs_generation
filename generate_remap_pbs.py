
__author__ = 'Ming Tang'

import csv

datadir = '/scratch/bcb/fpbarthel/TCGA/WGS/GBM/'
ref = '/Users/Tammy/projects/verhaak_lab/GBM_structural_variants/data/2015-06-23-Floris_bamfiles_info/GBM_bamfiles.csv'


# the bamfiles.csv file from Floris contains the tumor and the paired normal bam file info
# the UUID in column 5 and column 8 are the folder names containing the bam file in HPC cluster

# create pbs submission files for tumor samples
with open(ref, "r") as f:
    next(f) # skip the header
    reader = csv.reader(f, delimiter=',')
    ## need a backslash to escape a literal backslash \\
    for row in reader:
        barcode = row[0]
        job_dir = row[4]
        bam = row[5]
        print barcode, job_dir, bam
        job_string = """
#PBS -N %s
#PBS -l nodes=1:ppn=12,walltime=24:00:00
#PBS -l mem=64g
#PBS -M mtang1@mdanderson.org
#PBS -m abe
#PBS -d /scratch/genomic_med/mtang1/scratch/results/%s
#PBS -o /scratch/genomic_med/mtang1/scratch/results/logs
#PBS -e /scratch/genomic_med/mtang1/scratch/results/logs
#PBS -V

speedseq realign -o %s -v -t 12 -T temp -M 40 /scratch/genomic_med/mtang1/scratch/LOWPASS_WGS_LGG1/ref_genome/human_g1k_v37.fasta \\
%s
        """ % (barcode, barcode, barcode, datadir + job_dir + '/' + bam)
        print job_string
        ofile= open("%s.pbs" % barcode, "w")
        ofile.write(job_string)
        ofile.close()

# create pbs submission files for normal paired samples, some normal pairs will be duplicated
# because there are recurrent tumors that paired with the same normal. It is fine here, the same
# normal pbs file will be over-written.

with open(ref, "r") as f:
    next(f) # skip the header
    reader = csv.reader(f, delimiter=',')
    ## need a backslash to escape a literal backslash \\
    for row in reader:
        barcode = row[6]
        job_dir = row[7]
        bam = row[8]
        print barcode, job_dir, bam
        job_string = """
#PBS -N %s
#PBS -l nodes=1:ppn=12,walltime=24:00:00
#PBS -l mem=64g
#PBS -M mtang1@mdanderson.org
#PBS -m abe
#PBS -d /scratch/genomic_med/mtang1/scratch/results/%s
#PBS -o /scratch/genomic_med/mtang1/scratch/results/logs
#PBS -e /scratch/genomic_med/mtang1/scratch/results/logs
#PBS -V

speedseq realign -o %s -v -t 12 -T temp -M 40 /scratch/genomic_med/mtang1/scratch/LOWPASS_WGS_LGG1/ref_genome/human_g1k_v37.fasta \\
%s
        """ % (barcode, barcode, barcode, datadir + job_dir + '/' + bam)
        print job_string
        # write the pbs file into the folder
        ofile = open("%s.pbs" % barcode, "w")
        ofile.write(job_string)
        ofile.close()
