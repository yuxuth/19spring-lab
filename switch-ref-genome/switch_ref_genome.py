# input_file = sys.argv[1]
# output_file = sys.argv[2]
input_file = 'balb_b6j_exon.vcf'
# input_file = 'test_case1.vcf'
output_file = 'balb_b6_exon_compared_to_b6j.vcf'
input = open(input_file,'r')
output = open(output_file,'w')

def set_first_entry_to_11(line):
    """
    If the input line is a dot, "1/1" will be returned. Otherwise, Set the first
    entry to 1/1 while keeping the rest the same. For example,
    1/0:1:59:127:0:255,178,0:290,192,0:60:2:59:0,0,28,31:0:.:-0.693147
    will be transformed to
    1/1:1:59:127:0:255,178,0:290,192,0:60:2:59:0,0,28,31:0:.:-0.693147
    """
    if line==".":
        return "1/1:1:44:127:0:255,129,0:287,142,0:60:2:44:0,0,25,19:0:.:-0.693146"
    all_entries = line.strip('\n').split(':')
    all_entries[0] = "1/1"
    return ':'.join(all_entries)

for line in input:
    fields = line.strip('\n').split('\t')

    # If current line is a part of header, keep same.
    if str(fields[0]).startswith('##') :
        output.writelines('\t'.join(fields) + '\n')
        continue
    if str(fields[0]).startswith('#CHROM') :
        fields[-1] = "B6" # Change name of the last column
        output.writelines('\t'.join(fields) + '\n')
        continue

    # Check if the line is split correctly.
    if len(fields)!=11:
        print("Error: length of fields not equal to 11.")
        print("Current line split: ", fields)
        raise SyntaxError(line)

    # If last column is a dot, keep everything same.
    if fields[-1]=="." :
        output.writelines('\t'.join(fields) + '\n')
        continue

    # If the 2nd last column is a dot, we assume the last column is always 1/1, then
    # (1) switch alt/ref;
    # (2) set the 2nd last colum to 1/1;
    # (3) set the last colum to 1/1.
    if fields[-2]=="." :
        fields[3],fields[4] = fields[4][0], fields[3] # Switch ref and alt alleles.
        fields[-2]= set_first_entry_to_11(fields[-2])
        fields[-1]= set_first_entry_to_11(fields[-1])
        output.writelines('\t'.join(fields) + '\n')
        continue

    # If neither the two last columns are dot, we assume both of them are 1/1, then
    # (1) switch alt/ref;
    # (2) set the 2nd last column to ".";
    # (3) set the last colum to 1/1.
    else:
        fields[3],fields[4] = fields[4][0], fields[3] # Switch ref and alt alleles.
        fields[-2] = "."
        fields[-1] = set_first_entry_to_11(fields[-1])
        output.writelines('\t'.join(fields) + '\n')

input.close()
output.close()
