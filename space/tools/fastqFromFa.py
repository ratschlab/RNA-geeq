def main(fasta, fastq):
    fq = open(fastq, 'w')
    fa = open(fasta, 'r')
    lines = fa.readlines()
    fa.close()

    name = ''

    for i in range(len(lines),):
        if i % 2 == 0:
            name = lines[i][1:]
            print >> fq, '@' + name,
        else:
            print >> fq, lines[i].upper(),
            print >> fq, '+' + name,
            print >> fq, (len(lines[i]) - 1) * 'g'
    fq.close()

if __name__ == '__main__':

    import sys

    infile = sys.argv[1]
    outfile = sys.argv[2]
    
    main(infile, outfile)
