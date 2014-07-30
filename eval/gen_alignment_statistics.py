"""This script generates statistical overviews for a given alignment. """
import sys
import re
import cPickle
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy as sp
import pdb

from utils import *

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--alignment', dest='align', metavar='FILE', help='alignment file in sam format', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-R', '--ignore_multireads', dest='multireads', metavar='FILE', help='file containing the multireads to ignore', default='-')
    optional.add_option('-g', '--genome', dest='genome', metavar='FILE', help='genome in fasta or hdf5 format (needs ending .hdf5 for latter)', default='-')
    optional.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    optional.add_option('-X', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismathes [-]', default=None)
    optional.add_option('-M', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [100000000]', default='100000000')
    optional.add_option('-I', '--ignore_missing_chr', dest='ignore_missing_chr', action='store_true', help='ignore chromosomes missing in the annotation', default=False)
    optional.add_option('-s', '--shift_start', dest='shift_start', action='store_false', help='turn shifting start of softclips to accomodate for old bug OFF - it is usually ON!', default=True)
    optional.add_option('-b', '--bam_input', dest='bam_input', action='store_true', help='input has BAM format - does not work for STDIN', default=False)
    optional.add_option('-S', '--samtools', dest='samtools', metavar='PATH', help='if SAMtools is not in your PATH, provide the right path here (only neccessary for BAM input)', default='samtools')
    optional.add_option('-o', '--outfile_base', dest='outfile_base', metavar='PATH', help='basedir for outfiles written', default='-')
    optional.add_option('-l', '--lines', dest='lines', metavar='INT', type='int', help='maximal number of alignment lines to read [-]', default=None)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    optional.add_option('-d', '--debug', dest='debug', action='store_true', help='print debugging output', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def main():
    """Main function generating the alignment statistics."""

    ### get command line arguments
    options = parse_options(sys.argv)

    if options.outfile_base != '-':
        if options.align != '-':
            outfile_base = os.path.join(options.outfile_base, os.path.basename(options.align))
        else:
            outfile_base = os.path.join(options.outfile_base, 'results')
    else:
        if options.align != '-':
            outfile_base = options.align
        else:
            print >> sys.stderr, "Please provide an output file basedir when reading from STDIN!"
            sys.exit(-1)

    ### set filter tags for filenames
    if options.min_exon_len > 0:
        outfile_base += '_me%s' % options.min_exon_len
    if options.max_mismatches is not None:
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

    ### load genome
    if options.genome != '-':
        if options.genome.split('.')[-1] == 'hdf5':
            genome = hdf52dict(options.genome)
            for g in genome:
                genome[g] = str(genome[g])
        else:
            genome = read_fasta(options.genome)
       
    ### set infile streams
    if options.align != '-':
        infiles = []
        for fi in options.align.strip(',').split(','):
            if options.bam_input:
                file_handle = subprocess.Popen([options.samtools, 'view', options.align], stdout=subprocess.PIPE) 
                infiles.append(file_handle.stdout)
            else:
                infiles.append(open(options.align, 'r'))
    else:
        infiles = [sys.stdin]

    ### initializations
    filter_counter = 0
    unspliced = 0
    readlen = 0
    max_readlen = 30
    last_id = ('', 0)
    mismatches = sp.zeros((len(infiles), 1000), dtype='int')
    deletions = sp.zeros((len(infiles), 1000), dtype='int')
    insertions = sp.zeros((len(infiles), 1000), dtype='int')
    qualities_per_pos = sp.zeros((len(infiles), 1000), dtype='int')
    qualities = sp.zeros((len(infiles), 200), dtype='int')
    number_of_exons = sp.zeros((len(infile), 100), dtype='int')
    intron_pos = sp.zeros((len(infiles), 1000), dtype='int')

    ### iterate over infiles
    for f, infile in enumerate(infiles):
        for counter, line in enumerate(infile):
            if line[0] in ['@', '#' ] or line[:2] == 'SQ':
                continue
            if options.lines is not None and counter > options.lines:
                break
            if options.verbose and counter > 0 and counter % 10000 == 0:
                print 'lines read: [ %s (taken: %s / filtered: %s)]' % (counter, counter - filter_counter, filter_counter)
            sl = line.strip().split('\t')
            
            if len(sl) < 11 or (sl[0], int(sl[1]) & 128) in multireads:
                filter_counter += 1
                continue

            quality_string = ''
            if sl[9] != '*':
                readlen = len(sl[9])
                read = sl[9].upper()
                max_readlen = max(readlen, max_readlen)
            else:
                print >> sys.stderr, 'No read sequence given in SAM'
                sys.exit(-1)

            if options.max_mismatches is not None:
                cont_flag = False
                mm = get_mm(sl)
                if mm == -1:
                    print >> sys.stderr, 'No mismatch information available or read string missing in %s' % options.align
                    sys.exit(1)
                elif mm > options.max_mismatches:
                    filter_counter += 1
                    continue

            if options.min_exon_len > 0:
                cont_flag = False
                __cig = sl[5]
                __cig = re.sub('[0-9]*[IHS]', '', __cig) 
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

            curr_id = (sl[0], int(sl[1]) & 196)

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
            if options.genome != '-':

                (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
                size = [int(i) for i in size]
                chrm_pos = 0    # position in chrm
                read_pos = 0    # actual position in the read
                clipped_read_pos = 0
                
                for pos in range(len(size)):
                    if op[pos] == 'M':
                        gen_start = int(sl[3]) - 1
                        try:
                            gen = genome[sl[2]][gen_start + chrm_pos : gen_start + chrm_pos + size[pos]].upper()
                        except:
                            if options.ignore_missing_chr:
                                continue
                            else:
                                print >> sys.stderr, 'Chromosome name %s could not be found in %s' % (sl[2], options.genome)
                                sys.exit(1)
                        for p in range(size[pos]):
                            try:
                                if gen[p] != read[read_pos + p]:
                                    try:
                                        mismatches[abs(_reversed - (clipped_read_pos + read_pos + p))] += 1
                                    except KeyError:
                                        mismatches[abs(_reversed - (clipped_read_pos + read_pos + p))] = 1
                            except IndexError:
                                if options.debug:
                                    print >> sys.stderr, 'gen: %s' % gen
                                    print >> sys.stderr, 'read: %s' % read
                                    print >> sys.stderr, 'pos in gen: %i' % p
                                    print >> sys.stderr, 'pos in read: %i' % (read_pos + p)
                                    pdb.set_trace()
                                else:
                                    print >> sys.stderr, 'Index Error in line:\n %s' % line
                                    sys.exit(1)
                        chrm_pos += size[pos]
                        read_pos += size[pos]
                    elif op[pos] == 'I': # insertions
                        for _p in range(size[pos]):
                            try:
                                insertions[abs(_reversed - (read_pos + _p + clipped_read_pos))] += 1
                            except KeyError:
                                insertions[abs(_reversed - (read_pos + _p + clipped_read_pos))] = 1
                        read_pos += size[pos]
                    elif op[pos] == 'D': # deletions
                        try:
                            deletions[abs(_reversed - read_pos - clipped_read_pos)] += 1 # count only one deletion, not depending on number of positions deleted. ...size[pos]
                        except KeyError:
                            deletions[abs(_reversed - read_pos - clipped_read_pos)] = 1 #size[pos]
                        chrm_pos += size[pos]
                    elif op[pos] == 'N': # introns
                        chrm_pos += size[pos]
                    elif op[pos] == 'S': # softclips
                        read_pos += size[pos]
                        if options.shift_start:
                            chrm_pos += size[pos]
                    elif op[pos] == 'H': # hardclips
                        clipped_read_pos += size[pos]

            ### build up quality distribution
            if curr_id != last_id:
                if quality_string == '' and len(sl) > 10 and sl[10] != '*':
                    quality_string = sl[10]
                    if (int(sl[1])) & 16 == 16:
                        quality_string = quality_string[::-1]

                if quality_string != '':
                    for _pidx, _p in enumerate(quality_string):
                        qualities[ord(_p)] += 1
                        qualities_per_pos[_pidx] += ord(_p)

            last_id = curr_id

    ### plotting
    fig = plt.figure(figsize=(10, 15), dpi=300)
    gs = gridspec.GridSpec(4, 2)
    axes = []
    ### intron positions
    axes.append(plt.subplot(gs[0, 0]))
    axes[-1].plot(sp.arange(max_readlen), intron_pos[:max_readlen])
    axes[-1].set_xlabel('read position')
    axes[-1].set_ylabel('frequency')
    axes[-1].set_title('Split Position Distribution')
    ### number of exons
    axes.append(plt.subplot(gs[0, 1]))
    axes[-1].bar(sp.arange(number_of_exons.shape[0]), number_of_exons)
    axes[-1].set_xlabel('nunber of segments')
    axes[-1].set_ylabel('frequency')
    axes[-1].set_title('Number of Segments')
    ### mismatch distribution
    axes.append(plt.subplot(gs[1, 0]))
    axes[-1].plot(sp.arange(max_readlen), mismatches[:max_readlen])
    axes[-1].set_xlabel('read position')
    axes[-1].set_ylabel('mismatches')
    axes[-1].set_title('Mismatch Distribution')
    ### insertion distribution
    axes.append(plt.subplot(gs[1, 1]))
    axes[-1].plot(sp.arange(max_readlen), insertions[:max_readlen])
    axes[-1].set_xlabel('read position')
    axes[-1].set_ylabel('insertions')
    axes[-1].set_title('Insertion Distribution')
    ### deletion distribution
    axes.append(plt.subplot(gs[2, 0]))
    axes[-1].plot(sp.arange(max_readlen), deletions[:max_readlen])
    axes[-1].set_xlabel('read position')
    axes[-1].set_ylabel('deletions')
    axes[-1].set_title('Deletion Distribution')
    ### quality distribution
    axes.append(plt.subplot(gs[2, 1]))
    axes[-1].plot(sp.arange(qualities.shape[0]), qualities)
    axes[-1].set_xlabel('phred score')
    axes[-1].set_ylabel('frequency')
    axes[-1].set_title('Quality Distribution')
    ### deletion distribution
    axes.append(plt.subplot(gs[3, 0]))
    axes[-1].plot(sp.arange(max_readlen), qualities_per_pos[:max_readlen])
    axes[-1].set_xlabel('read position')
    axes[-1].set_ylabel('avg. quality')
    axes[-1].set_title('Positional Quality Distribution')

    if options.verbose:
        print 'number of exons: \n%s' % str(number_of_exons)
        print '%s reads were unspliced' % unspliced

    plt.tight_layout()

    ### plot data
    plt.savefig('overview.pdf', format='pdf')

    if options.align != '-':
        infile.close()

    ### write plotlists to file, for later evaluation
    #plot_out = (outfile_base + '.statistics.plotlist')
    #cPickle.dump(plotlist, open(plot_out, 'w'))
 
if __name__ == '__main__':
    main()
