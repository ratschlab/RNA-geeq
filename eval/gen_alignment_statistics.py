"""This script generates statistical overviews for a given alignment. """
import sys
import os
import re
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy as sp
import numpy.random as npr
import h5py
import time
import pdb

from modules.utils import *
from modules.plotting import *
from optparse import OptionParser, OptionGroup

def parse_options(argv, parser):

    """Parses options from the command line """

    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-g', '--genome', dest='genome', metavar='FILE', help='genome in fasta or hdf5 format (needs ending .hdf5 for latter)', default='-')
    optional.add_option('-I', '--ignore_missing_chr', dest='ignore_missing_chr', action='store_true', help='ignore chromosomes missing in the annotation', default=False)
    optional.add_option('-s', '--shift_start', dest='shift_start', action='store_false', help='turn shifting start of softclips to accomodate for old bug OFF - it is usually ON!', default=True)
    optional.add_option('-b', '--bam_input', dest='bam_input', action='store_true', help='input has BAM format - does not work for STDIN', default=False)
    optional.add_option('-S', '--samtools', dest='samtools', metavar='PATH', help='if SAMtools is not in your PATH, provide the right path here (only neccessary for BAM input)', default='samtools')
    optional.add_option('-o', '--outfile_base', dest='outfile_base', metavar='PATH', help='basename for outfiles written [align_stats]', default='align_stats')
    optional.add_option('-L', '--legend', dest='legend', action='store_true', help='put legend into plots [off]', default=False)
    optional.add_option('-l', '--lines', dest='lines', metavar='INT', type='int', help='maximal number of alignment lines to read [-]', default=None)
    optional.add_option('-r', '--random', dest='random', metavar='FLOAT', type='float', help='probability to accept an input line -- effective subsampling [1.0]', default=1.0)
    optional.add_option('-m', '--max_readlength', dest='max_readlen', metavar='INT', type='int', help='maximal read length to be considered [200]', default=200)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    optional.add_option('-d', '--debug', dest='debug', action='store_true', help='print debugging output', default=False)
    parser.add_option_group(optional)

    return parser.parse_args()

def get_tags(sl):
    """Extract tags from SAM line and return as dict"""

    #return  dict(z for z in [(x[0], int(x[2])) if x[1] == 'i' else (x[0], float(x[2])) if x[1] == 'f' else (x[0], x[2]) for x in [y.split(':') for y in sl]])

    tags = dict()
    for s in sl:
        ssl = s.split(':')
        #if ssl[1] == 'i':
        #    tags[ssl[0]] = int(ssl[2])
        #elif ssl[1] == 'f':
        #    tags[ssl[0]] = float(ssl[2])
        #else:
        tags[ssl[0]] = ssl[2]
    return tags

def main():
    """Main function generating the alignment statistics."""

    ### get command line arguments
    parser = OptionParser(usage="%prog [options] LIST OF ALIGNMENT FILES")
    (options, args) = parse_options(sys.argv, parser)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    ### load genome
    if options.genome != '-':
        if options.genome.split('.')[-1] == 'hdf5':
            genome = hdf52dict(options.genome)
            for g in genome:
                genome[g] = str(genome[g])
        else:
            genome = read_fasta(options.genome)
       
    infiles = args

    ### check, if infile is hdf5, in this case only do the plotting
    if infiles[0].endswith('hdf5'):
        for i, fname in enumerate(infiles):
            print >> sys.stdout, 'Loading counts from hdf5 %s' % fname
            h5_in = h5py.File(fname)
            if i == 0:
                plot_info = h5_in['plot_info'][:]
                counts = dict()
                filelist = h5_in['files'][:]
                for key in h5_in:
                    if key in ['files', 'plot_info']:
                        continue
                    counts[key] = h5_in[key][:]
            else:
                filelist = sp.r_[filelist, h5_in['files'][:]]
                for key in h5_in:
                    if key in ['files', 'plot_info']:
                        continue
                    if len(h5_in[key].shape) > 1 and h5_in[key].shape[1] > counts[key].shape[1]:
                        counts[key] = sp.c_[counts[key], sp.zeros((counts[key].shape[0], h5_in[key].shape[1] - counts[key].shape[1]))]
                        counts[key] = sp.r_[counts[key], h5_in[key][:]]
                    elif len(h5_in[key].shape) > 1 and h5_in[key].shape[1] < counts[key].shape[1]:
                        tmp = h5_in[key][:]
                        tmp = sp.c_[tmp, sp.zeros((tmp.shape[0], counts[key].shape[1] - h5_in[key].shape[1]))]
                        counts[key] = sp.r_[counts[key], tmp]
                    else:
                        counts[key] = sp.r_[counts[key], h5_in[key][:]]
            h5_in.close()
    else:
        ### initializations
        filter_counter = 0
        unspliced = 0
        readlen = 0
        max_readlen = 30
        counts = dict()
        for category in ['mismatches', 'deletions', 'insertions', 'qualities_per_pos', 'intron_pos', 'min_seg_len']:
            counts[category] = sp.zeros((len(infiles), options.max_readlen), dtype='int')
        counts['qualities'] = sp.zeros((len(infiles), 80), dtype='int')
        counts['number_of_segments'] = sp.zeros((len(infiles), 10), dtype='int')
        counts['deletion_lens'] = sp.zeros((len(infiles), 500), dtype='int')
        counts['insertion_lens'] = sp.zeros((len(infiles), 500), dtype='int')
        counts['multimappers'] = sp.zeros((len(infiles), 1000), dtype='int')
        for category in ['unaligned_reads', 'primary_alignments', 'secondary_alignments', 'unique_alignments', 'non_unique_alignments']:
            counts[category] = sp.zeros((len(infiles), ), dtype='int')

        t0 = time.time()

        ### iterate over infiles
        for f, fname in enumerate(infiles):
            ### open infile handle
            if fname == '-':
                infile = sys.stdin
            elif options.bam_input:
                fh = subprocess.Popen([options.samtools, 'view', fname], stdout=subprocess.PIPE)
                infile = fh.stdout
            else:
                infile = open(fname, 'r')

            taken_ids = set()

            if options.verbose:
                print >> sys.stdout, 'Parsing alignments from %s' % fname

            for counter, line in enumerate(infile):
                if line[0] in ['@', '#' ] or line[:2] == 'SQ':
                    continue
                if options.lines is not None and counter > options.lines:
                    break
                  
                if options.verbose and counter > 0 and counter % 100000 == 0:
                    t1 = time.time()
                    print 'lines read: [ %s (taken: %s / filtered: %s)] ... took %i sec' % (counter, counter - filter_counter, filter_counter, t1 - t0)
                    t0 = t1
                sl = line.strip().split('\t')

                if options.random < 1.0:
                    if npr.rand() > options.random and not sl[0] in taken_ids:
                        continue
                    else:
                        taken_ids.add(sl[0])
                
                if len(sl) < 11:
                    filter_counter += 1
                    continue

                ### check if unmapped
                if ((int(sl[1]) & 4) == 4):
                    counts['unaligned_reads'][f] +=1
                    continue

                if sl[9] != '*':
                    readlen = len(sl[9])
                    read = sl[9].upper()
                    max_readlen = max(readlen, max_readlen)
                else:
                    print >> sys.stderr, 'No read sequence given in SAM'
                    sys.exit(-1)

                is_secondary = ((int(sl[1]) & 256) == 256)
                if is_secondary:
                    counts['secondary_alignments'][f] += 1
                else:
                    counts['primary_alignments'][f] += 1

                tags = get_tags(sl[11:])
                if 'NH' in tags:
                    if int(tags['NH']) == 1:
                        counts['unique_alignments'][f] += 1
                    else:
                        counts['non_unique_alignments'][f] += 1
                    counts['multimappers'][f, int(tags['NH'])] += 1

                is_reversed = ((int(sl[1]) & 16) == 16)

                ### check, if read is reversed -> must change coordinates
                if is_reversed:
                    _reversed = readlen - 1
                else:
                    _reversed = 0

                ### record min segment length for spliced alignments
                if 'N' in sl[5]: 
                    __cig = sl[5]
                    __cig = re.sub('[0-9]*[IHS]', '', __cig) 
                    min_sl = min([sum([int('0'+i) for i in re.split('[^0-9]', '0' + _cig + 'Z0')][:-2]) for _cig in __cig.strip().split('N')])
                    counts['min_seg_len'][f, min_sl] += 1
                    ### count exons / segments in read
                    counts['number_of_segments'][f, sl[5].count('N') + 1] += 1

                    ### count intron distribution for spliced reads
                    ### the intron position is measured as the length of the first exon/segment (0-based position counting)
                    ### handle deletions - they do not affect block length
                    rl = sl[5]
                    rl = re.sub('[0-9]*D', '', rl)
                    rl = re.sub('[MISH]', 'M', rl) ### for this analysis softclips and hardclips are counted as positions in the original read
                    segm_len = sp.cumsum([sp.array(x.split('M')[:-1], dtype='int').sum() for x in ('%s0' % rl).split('N')])

                    ### in case of alignment to minus strand position is reversed
                    for s in segm_len[:-1]:
                        counts['intron_pos'][f, abs(_reversed - s)] += 1
                else:
                    unspliced += 1
                    ### count exons / segments in read
                    counts['number_of_segments'][f, 1] += 1

                ### build up mismatch-statistics from genome if MD tag is not available 
                (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
                size = [int(i) for i in size]
                chrm_pos = 0    # position in chrm
                read_pos = 0    # actual position in the read
                clipped_read_pos = 0
                
                for pos in range(len(size)):
                    if op[pos] == 'M' and options.genome != '-':
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
                                    counts['mismatches'][f, abs(_reversed - (clipped_read_pos + read_pos + p))] += 1
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
                        counts['insertion_lens'][f, size[pos]] += 1
                        _p = abs(_reversed - (read_pos + clipped_read_pos))
                        counts['insertions'][f, _p:_p + size[pos]] += 1
                       # for _p in range(size[pos]):
                       #     counts['insertions'][f, abs(_reversed - (read_pos + _p + clipped_read_pos))] += 1
                        read_pos += size[pos]
                    elif op[pos] == 'D': # deletions
                        counts['deletion_lens'][f, size[pos]] += 1
                        counts['deletions'][f, abs(_reversed - read_pos - clipped_read_pos)] += 1 # count only one deletion, not depending on number of positions deleted. ...size[pos]
                        chrm_pos += size[pos]
                    elif op[pos] == 'N': # introns
                        chrm_pos += size[pos]
                    elif op[pos] == 'S': # softclips
                        read_pos += size[pos]
                        if options.shift_start:
                            chrm_pos += size[pos]
                    elif op[pos] == 'H': # hardclips
                        clipped_read_pos += size[pos]

                ### build up quality distribution (only for primary alignments as this is a property of the key)
                ### do it only for 1% of the reads as it is too costly otherwise
                if not is_secondary and npr.random < 0.01:
                    if len(sl) > 10 and sl[10] != '*':
                        if is_reversed:
                            quality_string = sl[10][::-1]
                        else:
                            quality_string = sl[10]

                        for _pidx, _p in enumerate(quality_string):
                            counts['qualities'][f, ord(_p)] += 1
                            counts['qualities_per_pos'][f, _pidx] += ord(_p)

            ### clean up
            if fname != '-':
                infile.close()
            del taken_ids


        ### truncate counts to max non-zero x
        for c in counts:
            if len(counts[c].shape) > 1:
                max_idx = 0
                for i in range(counts[c].shape[0]):
                    idx = sp.where(counts[c][i, :] > 0)[0]
                    if idx.shape[0] > 0:
                        max_idx = max(max_idx, min(idx[-1] + 1, counts[c].shape[1]))
                    else:
                        max_idx = counts[c].shape[1]
                counts[c] = counts[c][:, :max_idx]
            else:
                idx = sp.where(counts[c] > 0)[0]
                if idx.shape[0] > 0:
                    max_idx = min(idx[-1] + 1, counts[c].shape[0])
                    counts[c] = counts[c][:max_idx]

        ### collect plot_info
        ### [data_field, plot_type, transformation, x-label, y-label, title']
        plot_info = [
            ['intron_pos', 'plot', '', 'read position', 'frequency', 'Split Position Distribution'],
            ['number_of_segments', 'bar', 'log10', 'number of segments', 'frequency', 'Number of Segments'],
            ['mismatches', 'plot', '', 'read position', 'mismatches', 'Mismatch Distribution'],
            ['insertions', 'plot', '', 'read position', 'insertions', 'Insertion Distribution'],
            ['deletions', 'plot', '', 'read position', 'deletions', 'Deletion Distribution'],
            ['qualities', 'plot', '', 'phred score', 'fequency', 'Quality Value Distribution'],
            ['qualities_per_pos', 'plot', '', 'read position', 'avg. quality', 'Position-wise Quality Distribution'],
            ['deletion_lens', 'plot', '', 'deletion length', 'frequency', 'Deletion Length Distribution'],
            ['insertion_lens', 'plot', '', 'deletion length', 'frequency', 'Insertion Length Distribution'],
            ['min_seg_len', 'plot', '', 'shortest segment length', 'frequency', 'Shortest Segment Length Distribution'],
            ['multimappers', 'plot', '', 'number of hits', 'frequency', 'Distribution of Alignment Ambiguity'],
            ['primary_alignments', 'bar', '', 'sample', 'number of alignments', 'Number of Primary Alignments'],
            ['secondary_alignments', 'bar', '', 'sample', 'number of alignments', 'Number of Secondary Alignments'],
            ['unaligned_reads', 'bar', '', 'sample', 'number of unaligned reads', 'Number of Unaligned Reads'],
            ['unique_alignments', 'bar', '', 'sample', 'number of unique alignments', 'Number of Unique Alignments'],
            ['non_unique_alignments', 'bar', '', 'sample', 'number of non-unique alignments', 'Number of Non-unique Alignments'],
            ]
        plot_info = sp.array(plot_info, dtype='str')

        ### store output as HDF5 file
        h5_out = h5py.File('%s.hdf5' % options.outfile_base, 'w')
        h5_out.create_dataset(name='files', data=sp.array(infiles, dtype='str'))
        h5_out.create_dataset(name='plot_info', data=plot_info)
        for key in counts:
            h5_out.create_dataset(name=key, data=counts[key], dtype='int')

        h5_out.close()

        filelist = infiles
    
    ### plotting
    fig = plt.figure(figsize=(15, 2*plot_info.shape[0]), dpi=300)
    gs = gridspec.GridSpec((plot_info.shape[0] + 1) / 2, 2)
    cmap = plt.get_cmap('jet')
    norm = plt.Normalize(0, len(infiles))  
    axes = []
    label_list = ['...' + x[-12:] if len(x) > 12 else x for x in filelist]
    for i in range(plot_info.shape[0]):
        axes.append(plt.subplot(gs[i / 2, i % 2]))
        if options.legend:
            plot(counts[plot_info[i, 0]], plot_info[i, :], ax=axes[-1], labels=label_list)
        else:
            plot(counts[plot_info[i, 0]], plot_info[i, :], ax=axes[-1])

    plt.tight_layout()

    ### plot data
    plt.savefig(options.outfile_base + '.overview.pdf', format='pdf')

if __name__ == '__main__':
    main()
