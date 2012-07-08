"""This script evaluates the error distribution in a given error model"""

import sys
import genome_utils
import re

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--alignment', dest='alignment', metavar='FILE', help='alignment file', default='-')
    required.add_option('-H', '--histogram', dest='histogram', metavar='PATH', help='print error histograms to PATH ', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-f', '--format', dest='format', metavar='STRING', help='alignment file format sam or bed [default: sam]', default='sam')
    optional.add_option('-g', '--genome', dest='gio_file', metavar='FILE', help='Genome information object', default='-')
#    optional.add_option('-O', '--output', dest='outputfile', metavar='FILE', help='file to store the evaluation output', default='-')
    
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()

    if len(argv) < 3:
        parser.print_help()
        sys.exit(2)

    return options

def read_alignment(options):
    """Extracts quality und substitutions information from an alignment given in bed-format """
    qualities = dict()
    substitutions = dict()
    quality_per_pos = dict()
    avg_quality_per_pos = dict()
    line_counter = 0
    
    if options.format == 'sam' and options.gio_file == '-':
        print >> sys.stderr, 'Please specify a genome information object to complete information missing in the SAM file - option -g'
        sys.exit(-1)


    if options.gio_file != '-':
        gio = genome_utils.GenomeInfo(options.gio_file)
        gio.contig_names.sort()
        genome = dict()

    for line in open(options.alignment, 'r'):
        if line[0] == '#':
            continue
        line_counter += 1
        sl = line.strip().split('\t')
        if options.format == 'bed':
            if len(sl) < 9:
                continue
            read = sl[12]
            quality = sl[13]
        elif options.format == 'sam':
            if len(sl) < 9:
                continue

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
            
            read = sl[9]
            _read = ''
            gen_start = int(sl[3]) - 1

            for pos in range(len(size)):
                if op[pos] == 'M':
                    gen = genome[sl[2]][gen_start + chrm_pos : gen_start + chrm_pos + size[pos]].upper()
                    for p in range(size[pos]):
                        if gen[p] != read[read_pos + p]:
                           _read += '[%s%s]' % (gen[p], read[read_pos + p])
                        else:
                            _read += read[read_pos + p]
                    chrm_pos += size[pos]
                    read_pos += size[pos]
                elif op[pos] == 'I': # insertions
                    for _p in range(size[pos]):
                        _read += '[-%s]' % read[read_pos + _p]
                    read_pos += size[pos]
                elif op[pos] == 'D': # deletions
                    for _p in range(size[pos]):
                        _read += '[%s-]' % genome[sl[2]][gen_start + chrm_pos + _p]
                    chrm_pos += size[pos]
                elif op[pos] == 'N': # introns
                    chrm_pos += size[pos]
                elif op[pos] == 'S': # softclips
                    read_pos += size[pos]
                    chrm_pos += size[pos]

            read = _read
            quality = sl[10]
            
        start = read.find('[', 0)
        
        ### check for substitutions
        while start > -1:
            ### skip indels
            if read[start + 1] == '-' or read[start + 2] == '-':
                start = read.find('[', start + 1)
                continue
            if substitutions.has_key((read[start + 1], read[start + 2])):
                substitutions[(read[start + 1], read[start + 2])] += 1
            else:
                substitutions[(read[start + 1], read[start + 2])] = 1
            start = read.find('[', start + 1)

        ### read quality values
        offset = 0
        subst = False

        for q in range(len(quality)):
            if read[q + offset] == '[':
                ### skip insertions
                if read[q + offset + 1] == '-':
                    offset += 3
                    continue 
                c = read[q + offset + 1]
                offset += 3
                subst = True
            else:
                c = read[q + offset]
            if not qualities.has_key((c, quality[q])):
                qualities[(c, quality[q])] = [0, 0]
            if subst:
                qualities[(c, quality[q])][0] += 1
                subst = False
            qualities[(c, quality[q])][1] += 1
 
            ### build up list of quality distribution per read length
            try:
                avg_quality_per_pos[q] += ord(quality[q])
            except KeyError:
                avg_quality_per_pos[q] = ord(quality[q])
            try:
                quality_per_pos[q][quality[q]] += 1
            except KeyError:
                try:
                    quality_per_pos[q][quality[q]] = 1
                except KeyError:
                    quality_per_pos[q] = {quality[q]:1}

    ### average quality values
    for key in avg_quality_per_pos.keys():
        avg_quality_per_pos[key] /= line_counter

 
    return (qualities, substitutions, quality_per_pos, avg_quality_per_pos, line_counter) 

def plot_histograms(options, qualities, substitutions, quality_per_pos, avg_quality_per_pos, line_counter):
    """Plots the error distribution relative to intron positions using gnuplot. """

    import tempfile
    import os
    import subprocess

    plot_tmp_name = tempfile.mkstemp()[1]
    plot_tmp = open(plot_tmp_name, 'w')
    plotpath = options.histogram

    ### build up list of error probabilities per position and avg. error probabilities per quality value
    ### we assume, that the nucleotides are equally distributed 
    error_prob_per_pos = dict()
    error_prob_per_quality = dict()
    for pos in quality_per_pos.keys():
        if not error_prob_per_pos.has_key(pos):
            error_prob_per_pos[pos] = 0
        for qual in quality_per_pos[pos].keys():
            avg_substitution_rate = 0
            for base in ['A', 'C', 'G', 'T']:
                try:
                    avg_substitution_rate +=  (0.25 * float(qualities[(base, qual)][0]) / float(qualities[(base, qual)][1]))
                except KeyError:
                    pass
            error_prob_per_pos[pos] += (float(quality_per_pos[pos][qual]) / float(line_counter) * avg_substitution_rate)
            if not error_prob_per_quality.has_key(qual):    
                error_prob_per_quality[qual] = avg_substitution_rate
            #error_prob_per_quality[qual] += (float(quality_per_pos[pos][qual]) / float(line_counter) * avg_substitution_rate / len(quality_per_pos.keys()))
            #error_prob_per_quality[qual] += (float(quality_per_pos[pos][qual]) / float(line_counter) * avg_substitution_rate)
    print >> plot_tmp, 'set terminal png size 1024,1024 enhanced'
    print >> plot_tmp, 'set output \"' + plotpath + '\"'
    print >> plot_tmp, 'set multiplot'

    print >> plot_tmp, 'reset'
    print >> plot_tmp, 'set title \"Substitions\"'
    print >> plot_tmp, 'set size 0.5, 0.5'
    print >> plot_tmp, 'set origin 0.0, 0.0'
    print >> plot_tmp, 'set xrange [-0.5:3.5]'
    print >> plot_tmp, 'set xlabel \"to base\"' 
    print >> plot_tmp, 'set ylabel \"rel. number of substitutions\"' 
    print >> plot_tmp, 'set xtics (\"A\" 0, \"C\" 1, \"G\" 2, \"T\" 3)'

    labelset = ['A', 'C', 'G', 'T']
    sums = [0, 0, 0, 0]
    for idx in range(len(labelset)):
        for idx2 in range(len(labelset)):
            if substitutions.has_key((labelset[idx], labelset[idx2])):
                sums[idx] += substitutions[(labelset[idx], labelset[idx2])]

    plot = 'plot '
    data = '\n'
    cum = [sums[i] for i in range(len(sums))]
    for idx2 in range(len(labelset)):
        plot += '\"-\" with boxes fs solid 0.7 title \"from %s\",' % labelset[idx2]
        for idx in range(len(labelset)):
            if substitutions.has_key((labelset[idx], labelset[idx2])):
                data += (str(float(cum[idx]) / max(sums[idx], 1)) + '\n')
                cum[idx] -= substitutions[(labelset[idx], labelset[idx2])]
                #data += (str(float(substitutions[(labelset[idx], labelset[idx2])]) / max(sums[idx], 1)) + '\n')
            else:
                data += '0\n'
        data += 'e\n'
    print >> plot_tmp, plot[:-1], data
    print >> plot_tmp, 'set key'
 
    print >> plot_tmp, 'reset'
    print >> plot_tmp, 'set title \"Avg. quality per pos\"'
    print >> plot_tmp, 'set size 0.5, 0.5'
    print >> plot_tmp, 'set origin 0.5, 0.0'
    print >> plot_tmp, 'set xlabel \"position\"' 
    print >> plot_tmp, 'set ylabel \"avg. quality\"' 
    print >> plot_tmp, 'set key'

    plot = 'plot \"-\" with boxes title \"avg. quality\",'
    data = '\n'
    index = avg_quality_per_pos.keys()
    index.sort()
    for pos in index:
        data += (str(avg_quality_per_pos[pos]) + '\n')
    data += 'e\n'
    print >> plot_tmp, plot[:-1], data

    print >> plot_tmp, 'reset'
    print >> plot_tmp, 'set title \"Error probability per pos\"'
    print >> plot_tmp, 'set size 0.5, 0.5'
    print >> plot_tmp, 'set origin 0.0, 0.5'
    print >> plot_tmp, 'set xlabel \"position\"' 
    print >> plot_tmp, 'set ylabel \"error prob.\"' 
    print >> plot_tmp, 'set key'

    plot = 'plot \"-\" with boxes title \"error prob.\",'
    data = '\n'
    index = error_prob_per_pos.keys()
    index.sort()
    for pos in index:
        data += (str(error_prob_per_pos[pos]) + '\n')
    data += 'e\n'
    print >> plot_tmp, plot[:-1], data

    print >> plot_tmp, 'reset'
    print >> plot_tmp, 'set title \"Average error probability per quality\"'
    print >> plot_tmp, 'set size 0.5, 0.5'
    print >> plot_tmp, 'set origin 0.5, 0.5'
    print >> plot_tmp, 'set xlabel \"quality value\"' 
    print >> plot_tmp, 'set ylabel \"error prob.\"' 
    print >> plot_tmp, 'set key'

    plot = 'plot \"-\" with boxes title \"avg. error prob.\",'
    data = '\n'
    index = error_prob_per_quality.keys()
    index.sort()
    tics = 'set xtics ('
    for tic in range(len(index)):
        if not index[tic] in set(['\"', '%']):
            tics += '\"%s\" %i,' % (index[tic], tic)
#        if index[tic] == '\"':
#            tics += '\"\\"\" %i,' % tic
#        elif index[tic] == '\%':
#            tics += '\"-\" %i,' % tic
        else:
            tics += '\"??\" %i,' % tic

    tics = tics[:-1] + ')'
    print >> plot_tmp, tics
    for qual in index:
        data += (str(error_prob_per_quality[qual]) + '\n')
    data += 'e\n'
    print >> plot_tmp, plot[:-1], data


    print >> plot_tmp, 'unset multiplot'

    plot_tmp.close()
    subprocess.call(['gnuplot', plot_tmp_name])
    #print plot_tmp_name
    os.remove(plot_tmp_name)


def main():
    """Function controlling the program flow ... """
    options = parse_options(sys.argv)

    (qualities, substitutions, quality_per_pos, avg_quality_per_pos, line_counter) = read_alignment(options)

    plot_histograms(options, qualities, substitutions, quality_per_pos, avg_quality_per_pos, line_counter)


if __name__ == '__main__':
    main()

