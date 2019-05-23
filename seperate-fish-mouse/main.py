import HTSeq
import gzip
import subprocess
from argparse import ArgumentParser

"""
seperate fq based on align statues
"""

# BAM_FILE = 'result.bam'
# FASTQ_R1 = 'RNA_10_USPD16095262-AK399-AK12890_HY2VKCCXY_L8_1.fq.gz'
# FASTQ_R2 = 'RNA_10_USPD16095262-AK399-AK12890_HY2VKCCXY_L8_2.fq.gz'
# FASTQ_R2 = 'test.fq'


def find_matched_reads_in_bam(filename):
    """
    return mapped reads name as set
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

def catagorize_fastq(matched,r1,r2): ## process two files
    """
    spit reads by name based on mapped reads name set
    """
    fastq_reader1 = HTSeq.FastqReader(r1)
    fastq_reader2 = HTSeq.FastqReader(r2)
    filename1 = r1.split('.')[0]
    filename2 = r2.split('.')[0]
    counter = 0
    aligned_output1 = open(filename1+"_aligned.fq", "w" )
    unaligned_output1 = open(filename1+"_unaligned.fq", "w" )
    
    aligned_output2 = open(filename2+"_aligned.fq", "w" )
    unaligned_output2 = open(filename2+"_unaligned.fq", "w" )
    
    for read1,read2 in zip(fastq_reader1,fastq_reader2 ):
        if read1.name.split(" ")[0] in matched:
            counter += 1
            read1.write_to_fasta_file(aligned_output1)
            read2.write_to_fasta_file(aligned_output2)
        else:
            read1.write_to_fasta_file(unaligned_output1)
            read2.write_to_fasta_file(unaligned_output2)

    aligned_output1.close()
    unaligned_output1.close()
    aligned_output2.close()
    unaligned_output2.close()

    # Compress as .gz
    subprocess.call(["gzip", filename1+"_aligned.fq"])
    subprocess.call(["gzip", filename1+"_unaligned.fq"])
    subprocess.call(["gzip", filename2+"_aligned.fq"])
    subprocess.call(["gzip", filename2+"_unaligned.fq"])

    print("Find %d aligns in fastq" % counter)
    return 0;

def main():
    parser = ArgumentParser(description='provide fastq and bam files') ## add parser
    parser.add_argument('-r1', help='r1', required=True)
    parser.add_argument('-r2',  help='r2', required=True)
    parser.add_argument('-bam',  help='aligned bam file', required=True)
    options = parser.parse_args()
    FASTQ_R1 = options.r1
    FASTQ_R2 = options.r2
    BAM_FILE = options.bam
    matched = find_matched_reads_in_bam(BAM_FILE)
    catagorize_fastq(matched,FASTQ_R1, FASTQ_R2)

if __name__== "__main__":
    main()
