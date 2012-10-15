#!/usr/bin/python
"""

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2009-2011 Andre Kahles
  Copyright (C) 2009-2011 by Andre Kahles

  This program generates overviews of different statistical parameters 
  of a given alignment files.
  
  For detailed usage information type:

    python genAlignmentStatistic.py

"""


import sys
import re
import subprocess

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--alignment', dest='align', metavar='FILE', help='alignment file in BAM  (default) or SAM (-S) format', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-R', '--ignore_multimappers', dest='ignore_multi', action='store_true', help='ignores multimappers, if marked in SAM file [off]', default=False)
    optional.add_option('-g', '--genome', dest='genome_fasta', metavar='FILE1[,FILE2, ...]', help='Genome sequence as fasta file(s) ', default='-')
    optional.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    optional.add_option('-X', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismathes [unlimited]', default=-1)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default='0')
    optional.add_option('-o', '--outfile_base', dest='outfile_base', metavar='PATH', help='basedir for outfiles written', default='-')
    optional.add_option('-s', '--samtools_path', dest='samtools', metavar='PATH', help='absolute path to samtools, if not in PATH', default='samtools')
    optional.add_option('-S', '--SAM', dest='sam', action='store_true', help='input alignment is in SAM format [off]', default=False)
    optional.add_option('-M', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [unlimited]', default=-1)
    optional.add_option('-x', '--use-x', dest='x', action='store_true', help='Use graphical output [off]', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def read_fasta(options):

    """Parses genome information from infiles in fasta-format"""

    infiles = options.genome_fasta.split(',')

    genome = dict()
    curr_chr = ''
    for infile in infiles:
        for line in open(infile, 'r'):
            if line[0] == '>':
                chr_name = line.strip()[1:]
                genome[chr_name] = ''
                curr_chr = chr_name
                if options.verbose:
                    print >> sys.stderr, "parsing %s ..." % curr_chr
                continue

            if curr_chr != '':
                genome[curr_chr] += line.strip().lower()
            else:
                print >> sys.stderr, "File in %s has no valid Fasta format!\n" % infile
                exit(-1)
                
    return genome

def main():
    """Main function generating the alignment statistics."""

    ### parse options from argument vector
    options = parse_options(sys.argv)

    ### initializations
    number_of_exons = []
    intron_pos = []
    counter = 0
    filter_counter = 0
    unspliced = 0
    readlen = 0

    if options.outfile_base != '-':
        outfile_base = (options.outfile_base + '/' + options.align.split('/')[-1])
    else:
        outfile_base = options.align

    if options.genome_fasta != '-':
        genome = read_fasta(options)

    ### set filter tags for filenames
    if options.min_exon_len > 0:
        outfile_base += '_me%s' % options.min_exon_len
    if options.max_mismatches > -1:
        outfile_base += '_mm%s' % options.max_mismatches
    if options.max_intron_len < 100000000:
        outfile_base += '_mi%s' % options.max_intron_len
 
    mismatches = []
    deletions = []
    insertions = []
    qualities = []

    if options.align == '-':    
        infile = sys.stdin
    else:
        if options.sam:
            infile = open(options.align, 'r')
        else:
            file_handle = subprocess.Popen([options.samtools, 'view', options.align], stdout=subprocess.PIPE) 
            infile = file_handle.stdout

    for line in infile:
        if line[0] in ['@', '#' ] or line[:2] == 'SQ':
            continue
        if options.verbose and counter % 10000 == 0:
            print 'lines read: [ %s (taken: %s / filtered: %s)]' % (counter, counter - filter_counter, filter_counter)
        counter += 1
        sl = line.strip().split('\t')
        
        if len(sl) < 11 or (options.ignore_multi and int(sl[1]) & 256 == 256):
            filter_counter += 1
            continue

        if sl[9] != '*':
            readlen = len(sl[9])
            read = sl[9].upper()
        else:
            print >> sys.stderr, 'Error: No read sequence given in alignment file!'
            sys.exit(-1)

        if options.max_mismatches > -1:
            cont_flag = False
            try:
                for opt in sl[11:]:
                    if (opt[:3] == 'NM:' and int(opt[5:]) > options.max_mismatches):
                        cont_flag = True
                        break
 
                if cont_flag:
                    filter_counter += 1
                    continue
            except KeyError:
                print >> sys.stderr, 'Error: Alignment file lacks mismatch informaition!'
                sys.exit(1)

        if options.min_exon_len > 0:
            cont_flag = False
            __cig = sl[5]
            __cig = re.sub('[0-9]*I', '', __cig) 
            for _cig in __cig.strip().split('N'):
                if sum([int('0'+i) for i in re.split('[^0-9]', '0' + _cig + 'Z0')][:-2]) < options.min_exon_len:
                    cont_flag = True
                    break
            if cont_flag:
                filter_counter += 1
                continue

        if options.max_intron_len > -1:
            cont_flag = False
            (op, size) = (re.split('[0-9]*', sl[5])[1:], re.split('[^0-9]', sl[5])[:-1])
            size = [int(i) for i in size]

            for o in range(len(op)):
                if op[o] == 'N' and (size[o] > options.max_intron_len):
                    cont_flag = True
                    break
            if cont_flag:
                filter_counter += 1
                continue

        ### count exons / segments in read
        idx = sl[5].count('N') + 1
        try:
            number_of_exons[idx] += 1
        except IndexError:
            number_of_exons.extend([0] * (idx - len(number_of_exons) - 1))
            number_of_exons.append(1)

        ### check, if read is reversed -> must change coordinates
        if (int(sl[1]) & 16) == 16:
            _reversed = readlen - 1
        else:
            _reversed = 0

        ### count intron distribution for spliced reads
        ### the intron position is measured as the length of the first exon/segment (0-based position counting)
        if sl[5].find('N') == -1:
            unspliced += 1
        else:
            ### handle deletions - they do not affect block length
            rl = sl[5]
            rl = re.sub('[0-9]*D', '', rl)
            rl = re.sub('[MISH]', '$', rl) ### for this analysis softclips and hardclips are counted as positions in the original read
            exon = rl.split('N')[0]
            exon_list = exon.split('$')

            ### determine intron position (always position of the FIRST intron)
            ### in case of alignment to minus strand position is reversed
            exon_len = sum([int(i) for i in exon_list[:-1]])
            idx = abs(_reversed - exon_len)
            try:
                intron_pos[idx] += 1
            except IndexError:
                intron_pos.extend([0] * (idx - len(intron_pos) - 1))
                intron_pos.append(1)
    
        ### build up mismatch-statistics 
        if options.genome_fasta != '-':

            if not genome.has_key(sl[2]):
                print >> sys.stderr, 'Chromosome name %s could not be found in %s' % (sl[2], options.gio_file)
                sys.exit(-1)

            (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
            size = [int(i) for i in size]
            chrm_pos = 0    # position in chrm
            read_pos = 0    # position in the actual read
            hardclipped_read_pos = 0
            
            broken = False

            for pos in range(len(size)):
                if op[pos] == 'M':
                    gen_start = int(sl[3]) - 1
                    gen = genome[sl[2]][gen_start + chrm_pos : gen_start + chrm_pos + size[pos]].upper()
                    for p in range(size[pos]):
                        if not (read_pos + p) < len(read):
                            print '\naborted analysis of following line, due to block size inconsistencies:'
                            print (line + '\n')
                            broken = True
                            break
                        if gen[p] != read[read_pos + p]:
                            idx = abs(_reversed - (hardclipped_read_pos + read_pos + p))
                            try:
                                mismatches[idx] += 1
                            except IndexError:
                                mismatches.extend([0] * (idx - len(mismatches) - 1))
                                mismatches.append(1)
                    if broken:
                        break
                    chrm_pos += size[pos]
                    read_pos += size[pos]
                elif op[pos] == 'I': # insertions
                    for _p in range(size[pos]):
                        idx = abs(_reversed - (read_pos + _p + hardclipped_read_pos))
                        try:
                            insertions[idx] += 1
                        except KeyError:
                            insertions.extend([0] * (idx - len(insertions) - 1))
                            insertions.append(1)
                    read_pos += size[pos]
                elif op[pos] == 'D': # deletions
                    idx = abs(_reversed - read_pos - hardclipped_read_pos)
                    try:
                        deletions[idx] += 1  # count only one deletion, not depending on number of positions deleted. ...size[pos]
                    except KeyError:
                        deletions.extend([0] * (idx - len(deletions) - 1))
                        deletions.append(1) # size [pos]
                    chrm_pos += size[pos]
                elif op[pos] == 'N': # introns
                    chrm_pos += size[pos]
                elif op[pos] == 'S': # softclips
                    read_pos += size[pos]
                    chrm_pos += size[pos]
                elif op[pos] == 'H': # hardclips
                    hardclipped_read_pos += size[pos]

        ### build up quality distribution
        if len(sl) > 10 and sl[10] != '*':
            for _p in sl[10]:
                idx = ord(_p)
                try:
                    qualities[idx] += 1
                except IndexError:
                    qualities.extend([0] * (idx - len(qualities) - 1))
                    qualities.append(1)

    import matplotlib
    if not options.x:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(intron_pos)
    plt.title('Intron position distribution')
    plt.xlabel('read position')
    plt.ylabel('occurrences')
    plt.show()
    plt.figure(2)
    plt.plot(number_of_exons) ## --> make box plot
    plt.title('Number of exons')
    plt.xlabel('number of exons')
    plt.ylabel('occurrences')
    plt.show()
    plt.figure(3)
    plt.plot(mismatches)
    plt.title('Positional mismatch distribution')
    plt.xlabel('read position')
    plt.ylabel('occurrences')
    plt.show()
    plt.figure(4)
    plt.plot(insertions)
    plt.title('Positional insertion distribution')
    plt.xlabel('read position')
    plt.ylabel('occurrences')
    plt.show()
    plt.figure(5)
    plt.plot(deletions)
    plt.title('Positional deletion distribution')
    plt.xlabel('read position')
    plt.ylabel('occurrences')
    plt.show()
    plt.figure(6)
    plt.plot(qualities)
    plt.title('Quality distribution')
    plt.xlabel('phred score')
    plt.ylabel('occurrences')
    plt.show()
    print 'number of exons: \n%s' % str(number_of_exons)
    print '%s reads were unspliced' % unspliced

if __name__ == '__main__':
    main()
