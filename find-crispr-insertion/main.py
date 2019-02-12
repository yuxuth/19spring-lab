

P1_SEQUENCE = "PPPP"
P2_SEQUENCE = "pppp"
INSERT_SEQUENCE = "--------"

# 4 types of reads are stored in four lists
perfectly_inserted = []
not_inserted = [] # HMR
not_cleaved = []
noise = []
exception = [] # Can be handled manually,
# Exception List Includes:
# (1) Flawed/Incomplete insertion: the sequence is inserted but there are more.
# (2) Index of p1 is bigger than index of p2.
# (3) Only p1/p2 is found.
# (4) Only target sequence is found.


def read_from_sample_data(filename):
    f = open(filename, "r")
    reads = []
    while True:
        read = f.readline()
        if not read:
            break
        reads.append(read)
    f.close()
    return reads


def read_from_fastq(filename):
    reads = []
    with open(filename, 'r') as f:
        while True:
            f.readline()
            read = f.readline()
            if not read:
                break
            else:
                reads.append(read)
            f.readline()
            f.readline()
    return reads
def catagorize(reads):
    """
    Catagorize each read stored and add to the corresponding list.
    """
    for each in reads:
        if(P1_SEQUENCE in each and P2_SEQUENCE in each):
            p1_index = each.find(P1_SEQUENCE)
            p2_index = each.find(P2_SEQUENCE)
            if(p1_index<p2_index):
                if(p1_index+len(P1_SEQUENCE)==p2_index):
                    not_cleaved.append(each)
                else:
                    # The p1/p2 site is cleaved.
                    fragment = each[p1_index+len(P1_SEQUENCE):p2_index]
                    if(fragment==INSERT_SEQUENCE):
                        perfectly_inserted.append(each)
                    elif(INSERT_SEQUENCE in fragment):
                        # The sequence is inserted but there are more --> exception.
                        exception.append((each,1))
                    else:
                        not_inserted.append(each)
            else:
                # Index of p1 is bigger than index of p2.
                exception.append((each,2))
        elif(P1_SEQUENCE in each or P2_SEQUENCE in each):
            # Only p1/p2 is found --> exception.
            exception.append((each,3))
        elif(INSERT_SEQUENCE in each):
            # Only target sequence is found --> exception.
            exception.append((each,4))
        else:
            # Nothing is found --> noise.
            noise.append(each)
    return

def display(total):
    """
    Display result.
    """

    print("Perfectly inserted:  %.1f%%" % (len(perfectly_inserted)*100.0/total))
    print("Not inserted:        %.1f%%" % (len(not_inserted)*100.0/total))
    print("Not cleaved:         %.1f%%" % (len(not_cleaved)*100.0/total))
    print("Noise:               %.1f%%" % (len(noise)*100.0/total))
    print("Exception:           %.1f%%" % (len(exception)*100.0/total))

    print()
    for each in exception:
        if(each[1]==1):
            print("The sequence is inserted but there are more:")
        elif(each[1]==2):
            print("Index of p1 is bigger than index of p2:")
        elif(each[1]==3):
            print("Only p1/p2 is found:")
        elif(each[1]==4):
            print("Only target sequence is found:")
        print(each[0])
        print()
    return


def main():

    # Read data from data source.
    # filename = "sample_data.txt"
    # reads = read_from_sample_data(filename)
    filename = "data/trimmed-mClta_F.fastq"
    reads = read_from_fastq(filename)
    for each in reads:
        print(each)
    # total = len(reads)

    # catagorize(reads)

    # display(total)


if __name__ == '__main__':
    main()
