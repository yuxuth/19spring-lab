import os
import sys
input_file = sys.argv[1]
output_file = sys.argv[2]
# input_file = '200.vcf'
# output_file = '200_test.txt'
input = open(input_file,'r')
output = open(output_file,'w')

format = ''

while True:
    line = input.readline()
    fields = line.strip('\n').split('\t')
    if not str(fields[0]).startswith('#'):
        fields = line.strip('\n').split('\t')
        format = fields[len(fields)-1]
        # print (format)
        break
format = format.split(':')
input.close()

format = ['.'] * len(format)

with open(input_file) as input:
    for line in input:
        fields = line.strip('\n').split('\t')
        if str(fields[0]).startswith('##') :
            output.writelines('\t'.join(fields) + '\n')
            continue
        if str(fields[0]).startswith('#CHROM') :
            fields.append('B6')
            output.writelines('\t'.join(fields) +'' +'\n')
            continue
        fields.append(':'.join(format))
        output.writelines('\t'.join(fields) +'' +'\n')
input.close()
output.close()
    