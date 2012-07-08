"""This script generates statistical overviews for a given alignment. """
import sys
import re
import genome_utils
import cPickle
import plot

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--alignment', dest='align', metavar='FILE', help='alignment file in sam format', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-H', '--histogram', dest='histogram', metavar='PATH', help='destination PATH of evaluation plots', default='-')
    optional.add_option('-R', '--ignore_multireads', dest='multireads', metavar='FILE', help='file containing the multireads to ignore', default='-')
    optional.add_option('-r', '--readfile', dest='readfile', metavar='FILE', help='file containing all reads in FASTQ format, sorted by read_id', default='-')
    optional.add_option('-O', '--organism', dest='organism', metavar='STRING', help='one of [human|drosophila|elegans], default elegans', default='elegans')
    optional.add_option('-n', '--number_of_reads_to_align', dest='given_reads', metavar='INT', type='int', \
        help='number of reads that had to be aligned - needed to calculate amount of unmapped reads', default=0)
    optional.add_option('-g', '--genome', dest='gio_file', metavar='FILE', help='Genome information object', default='-')
    optional.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    optional.add_option('-X', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismathes [10000]', default=10000)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default='0')
    optional.add_option('-o', '--outfile_base', dest='outfile_base', metavar='PATH', help='basedir for outfiles written', default='-')
    optional.add_option('-M', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [100000000]', default='100000000')
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def main():
    """Main function generating the alignment statistics."""

    counter = 0
    filter_counter = 0
    f_line_counter = 0
    options = parse_options(sys.argv)

    ### initializations
    number_of_exons = dict()
    intron_pos = dict()
    unspliced = 0
    readlen = 0
    plotlist = plot.Plotlist()

    if options.outfile_base != '-':
        outfile_base = (options.outfile_base + '/' + options.align.split('/')[-1])
    else:
        outfile_base = options.align

    if options.readfile != '-':
        fastq = open(options.readfile, 'r')
        fastq_sl = fastq.readline().strip().split('\t')


    ### set filter tags for filenames
    if options.min_exon_len > 0:
        outfile_base += '_me%s' % options.min_exon_len
    if options.max_mismatches < 10000:
        outfile_base += '_mm%s' % options.max_mismatches
    if options.max_intron_len < 100000000:
        outfile_base += '_mi%s' % options.max_intron_len
 
    multireads = set()
    if options.multireads != '-':
        print '\nParsing multireads from file %s' % options.multireads
        print '-----------------------------------------'
        for line in open(options.multireads, 'r'):
            _l = line.strip().split('\t')
            multireads.add((_l[0], int(_l[1])))

    mismatches = dict()
    deletions = dict()
    insertions = dict()
    qualities = dict()

    if options.gio_file != '-':
        gio = genome_utils.GenomeInfo(options.gio_file)
        gio.contig_names.sort()
        genome = dict()
       
    for line in open(options.align, 'r'):
        if line[0] in ['@', '#' ] or line[:2] == 'SQ':
            continue
        if options.verbose and counter % 10000 == 0:
            print 'lines read: [ %s (taken: %s / filtered: %s)]' % (counter, counter - filter_counter, filter_counter)
        counter += 1
        sl = line.strip().split('\t')
        
        if len(sl) < 11 or (sl[0], int(sl[1]) & 128) in multireads:
            filter_counter += 1
            continue

        if options.readfile != '-':    
            while (fastq_sl[0] != sl[0] and len(fastq_sl) > 0 and fastq_sl[0] != ''):
                fastq_sl = fastq.readline().strip().split('\t')
                f_line_counter += 1
                if options.verbose and (f_line_counter % 10000) == 0:
                    print '[line in pseudofastq: %s]' % (f_line_counter)
            if (fastq_sl[0] != sl[0]):
                print >> sys.stderr, 'pseudofastq not sorted or not matching the current infile!'
                sys.exit(1)
            if int(sl[1]) & 64 == 64 or sl[0][-2:] == '/1':
                read = fastq_sl[1]
            elif int(sl[1]) & 128 == 128 or sl[0][-2:] == '/2':
                read = fastq_sl[2]
            readlen = len(read)
            if int(sl[1]) & 16 == 16:
                read = genome_utils.reverse_complement(read)
        elif sl[9] != '*':
            readlen = len(sl[9])
            read = sl[9].upper()
        else:
            print >> sys.stderr, 'No read sequence given in SAM, nor in additional readfile!'
            sys.exit(-1)

        if options.max_mismatches < 10000:
            cont_flag = False
            try:
                for opt in sl[11:]:
                    if (opt[:3] == 'NM:' and int(opt[5:]) > options.max_mismatches):
                        cont_flag = True
                        break
 
                if cont_flag:
                    filter_counter += 1
                    continue
            except:
                print >> sys.stderr, 'No mismatch information available or read string missing in %s' % options.align
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

        if options.max_intron_len < 100000000:
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
        try:
            number_of_exons[sl[5].count('N') + 1] += 1
        except KeyError:
            number_of_exons[sl[5].count('N') + 1] = 1

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
            try:
                intron_pos[abs(_reversed - exon_len)] += 1
            except KeyError:
                intron_pos[abs(_reversed - exon_len)] = 1
    
        ### build up mismatch-statistics 
        if options.gio_file != '-':
            sl[2] = re.sub('Cel_', '', sl[2])
            sl[2] = re.sub('.1\]', '', sl[2])
            sl[2] = re.sub('chr', '', sl[2])
            if sl[2] == 'M':
                sl[2] = 'MT'
            if options.organism == 'drosophila':
                sl[2] = re.sub('MT', 'dmel_mitochondrion_genome', sl[2])

            if not genome.has_key(sl[2]):
                print 'Reading chromosome %s' % sl[2]
                try:
                    fgen = open('%s/genome/%s.flat' % (gio.basedir, sl[2]))
                    genome[sl[2]] = fgen.readline()
                    fgen.close()
                except:
                    print >> sys.stderr, 'Chromosome name %s could not be found in %s' % (sl[2], options.gio_file)
                    sys.exit(1)

            (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
            size = [int(i) for i in size]
            chrm_pos = 0        # position in chrm
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
                            try:
                                mismatches[abs(_reversed - (hardclipped_read_pos + read_pos + p))] += 1
                            except KeyError:
                                mismatches[abs(_reversed - (hardclipped_read_pos + read_pos + p))] = 1
                    if broken:
                        break
                    chrm_pos += size[pos]
                    read_pos += size[pos]
                elif op[pos] == 'I': # insertions
                    for _p in range(size[pos]):
                        try:
                            insertions[abs(_reversed - (read_pos + _p + hardclipped_read_pos))] += 1
                        except KeyError:
                            insertions[abs(_reversed - (read_pos + _p + hardclipped_read_pos))] = 1
                    read_pos += size[pos]
                elif op[pos] == 'D': # deletions
                    try:
                        deletions[abs(_reversed - read_pos - hardclipped_read_pos)] += 1 # count only one deletion, not depending on number of positions deleted. ...size[pos]
                    except KeyError:
                        deletions[abs(_reversed - read_pos - hardclipped_read_pos)] = 1 #size[pos]
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
                try:
                    qualities[ord(_p)] += 1
                except KeyError:
                    qualities[ord(_p)] = 1

    max_len = max([76, len(mismatches.keys()), len(insertions.keys()), len(deletions.keys()), len(qualities.keys())])

    if len(mismatches.keys()) == 0:
        mismatches = dict()
        for i in range(max_len):
            mismatches[i] = 0
    elif len(mismatches.keys()) < max_len:
        for i in range(1, max_len + 1):
            if not mismatches.has_key(i):
                mismatches[i] = 0
    if len(insertions.keys()) == 0:
        insertions = dict()
        for i in range(max_len):
            insertions[i] = 0
    elif len(insertions.keys()) < max_len:
        for i in range(1, max_len + 1):
            if not insertions.has_key(i):
                insertions[i] = 0
    if len(deletions.keys()) == 0:
        deletions = dict()
        for i in range(max_len):
            deletions[i] = 0
    elif len(deletions.keys()) < max_len:
        for i in range(1, max_len + 1):
            if not deletions.has_key(i):
                deletions[i] = 0
    if len(intron_pos.keys()) == 0:
        intron_pos = dict()
        for i in range(max_len):
            intron_pos[i] = 0
    elif len(intron_pos.keys()) < max_len:
        for i in range(1, max_len + 1):
            if not intron_pos.has_key(i):
                intron_pos[i] = 0
    if len(qualities.keys()) == 0:
        qualities = dict()
        for i in range(max_len):
            qualities[i] = 0
    elif len(qualities.keys()) < max_len:
        for i in range(1, max_len + 1):
            if not qualities.has_key(i):
                qualities[i] = 0

    plotlist.append(plot.Plot('Intron position distribution', intron_pos, xlabel='read position', ylabel='occurrences', xrange=(0, 77)))
    plot_data(intron_pos, outfile_base, 'intron_position')
    plotlist.append(plot.Plot('Number of exons', number_of_exons, xlabel='number of exons', ylabel='occurrences', plottype='boxes'))
    plot_data(number_of_exons, outfile_base, 'num_of_exons')
    plotlist.append(plot.Plot('Positional mismatch distribution', mismatches, xlabel='read position', ylabel='occurrences', xrange=(0, 77)))
    plot_data(mismatches, outfile_base, 'mismatches')
    plotlist.append(plot.Plot('Positional insertion distribution', insertions, xlabel='read position', ylabel='occurrences', xrange=(0, 77)))
    plot_data(insertions, outfile_base, 'insertions')
    plotlist.append(plot.Plot('Positional deletion distribution', deletions, xlabel='read position', ylabel='occurrences', xrange=(0, 77)))
    plot_data(deletions, outfile_base, 'deletions')
    plotlist.append(plot.Plot('Quality distribution', qualities, xlabel='phred score', ylabel='occurrences', xrange=(0, 77)))
    plot_data(qualities, outfile_base, 'qualities')
    print 'number of exons: \n%s' % str(number_of_exons)
    print '%s reads were unspliced' % unspliced

    ### plot data
    if options.histogram != '-':
        plotlist.plot(options.histogram)

    if options.readfile != '-':
        fastq.close()

    ### write plotlists to file, for later evaluation
    plot_out = (outfile_base + '.statistics.plotlist')
    cPickle.dump(plotlist, open(plot_out, 'w'))
 
def plot_data(data, outfile_base, tag):
    """Plots data to outfile"""

    outfile = open(outfile_base + '.' + tag, 'w')

    _line = ''
    if type(data) == list:
        for d in data:
            _line += (str(d) + '\t')
    elif type(data) == dict:
        idx = data.keys()
        idx.sort()
        for d in idx:
            _line += (str(data[d]) + '\t')
    print >> outfile, _line[:-1] 
    outfile.close()


if __name__ == '__main__':
    main()
