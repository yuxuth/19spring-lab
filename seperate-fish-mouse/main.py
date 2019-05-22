import HTSeq
import gzip
import subprocess

"""
This script is used to seperate fish reads from a fastq file containing reads of
both fish and mouse.
(1) Creating a .bam file by mapping the original fastq file to the fish reference genome.
(2) Finding all matched reads in the bam file using HTseq is_matched methods.
(3) Using the dictionary produced by (2), iterate through R1 and R2 fastq files to
    seperate fish reads from mouse.
"""

BAM_FILE = 'result.bam'
FASTQ_R1 = 'RNA_10_USPD16095262-AK399-AK12890_HY2VKCCXY_L8_1.fq.gz'
FASTQ_R2 = 'RNA_10_USPD16095262-AK399-AK12890_HY2VKCCXY_L8_2.fq.gz'
# FASTQ_R2 = 'test.fq'


def find_matched_reads_in_bam(filename):
    """
    This method takes the name of a bam file as input, returning a list of names
    of the reads that is matched to the fish reference genome.
    """
    matched = [] # Names of all matched reads.

    all_reads = 0
    aligned_reads = 0

    bam_reader = HTSeq.BAM_Reader(filename)
    for align in bam_reader:
        all_reads+=1
        if align.aligned:
            # if(align.read.name in matched):
            #     print(align.read.name)
            matched.append(align.read.name)
            aligned_reads+=1
            # print(align.read.name)

    print("Total reads:",all_reads)
    print("Aligned reads:", aligned_reads)
    print("Unique aligned:",len(set(matched)))
    return set(matched)

def catagorize_fastq(matched,filename):
    """
    Catagoriz each read of a fastq file based on whether it is aligned to
    the reference genome (therefore name appears in matched) or not.
    Output two files - "filename_aligned" & "filename_unaligned".
    """
    fastq_reader = HTSeq.FastqReader(filename)
    counter = 0

    aligned_output = open(filename[:-6]+"_aligned.fq", "w" )
    unaligned_output = open(filename[:-6]+"_unaligned.fq", "w" )

    for read in fastq_reader:
        if read.name.split(" ")[0] in matched:
            counter+=1
            read.write_to_fasta_file(aligned_output)
        else:
            read.write_to_fasta_file(unaligned_output)

    aligned_output.close()
    unaligned_output.close()

    # Compress as .gz
    subprocess.call(["gzip", filename[:-6]+"_aligned.fq"])
    subprocess.call(["gzip", filename[:-6]+"_unaligned.fq"])

    print("Find %d aligns in fastq" % counter)
    return 0;

def main():
    matched = find_matched_reads_in_bam(BAM_FILE)
    catagorize_fastq(matched,FASTQ_R1)
    catagorize_fastq(matched,FASTQ_R2)

if __name__== "__main__":
    main()
