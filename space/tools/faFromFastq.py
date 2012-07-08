def main(fastq, fasta):
    fa = open(fasta, 'w')
    fq = open(fastq, 'r')
    lines = fq.readlines()
    fq.close()

    for i in range(len(lines),):
        if i % 4 == 0:
            print >> fa, '>' + lines[i][1:],
        elif i % 4 == 1:
            print >> fa, lines[i],
        elif i % 4 == 2 or i % 4 == 3:
            continue
    fa.close()
    fq.close()

if __name__ == '__main__':

    import sys

    infile = sys.argv[1]
    outfile = sys.argv[2]
    
    main(infile, outfile)
