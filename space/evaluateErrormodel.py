"""This script evaluates the error distribution in a given error model"""

import sys

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-e', '--errormodel', dest='error', metavar='FILE', help='error model', default='-')
    required.add_option('-q', '--quality', dest='quality', metavar='FILE', help='quality file', default='-')
    required.add_option('-H', '--histogram', dest='histogram', metavar='PATH', help='print error histograms to PATH ', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-o', '--offset', dest='offset', metavar='INT', type='int', help='offset for quality range', default=0)
    
    parser.add_option_group(required)

    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    if options.quality == '-':
        options.quality = options.error + '.quality'
    
    return options

def plot_histograms(options, qualities, substitutions):
    """Plots the error distribution relative to intron positions using gnuplot. """

    import tempfile
    import os
    import subprocess

    plot_tmp_name = tempfile.mkstemp()[1]
    plot_tmp = open(plot_tmp_name, 'w')
    plotpath = options.histogram

    ### build up list of quality distribution per read length
    quality_per_pos = dict()
    avg_quality_per_pos = dict()
    line_counter = 0
    for line in open(options.quality, 'r'):
        line_counter += 1
        line.strip()
        for c in range(len(line)):
            try:
                avg_quality_per_pos[c] += max(ord(line[c]) - options.offset, 0)
            except KeyError:
                avg_quality_per_pos[c] = max(ord(line[c]) - options.offset, 0)
            
            try:
                quality_per_pos[c][line[c]] += 1
            except KeyError:
                try:
                    quality_per_pos[c][line[c]] = 1
                except KeyError:
                    quality_per_pos[c] = {line[c]:1}


    ### average quality values
    for key in avg_quality_per_pos.keys():
        avg_quality_per_pos[key] /= float(line_counter)

    ### build up list of error probabilities per position and error probabilities per quality value
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
                    avg_substitution_rate +=  (0.25 * float(qualities[(base, qual)][0]) / float(max(qualities[(base, qual)][1], 1)))
                except KeyError:
                    pass
            error_prob_per_pos[pos] += (float(quality_per_pos[pos][qual]) / float(line_counter) * avg_substitution_rate)
            if not error_prob_per_quality.has_key(qual):    
            #    error_prob_per_quality[qual] = 0
                error_prob_per_quality[qual] = avg_substitution_rate 
            #error_prob_per_quality[qual] += (float(quality_per_pos[pos][qual]) / float(line_counter) * avg_substitution_rate / len(quality_per_pos.keys()))
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
        #plot += '\"-\" with boxes title \"from %s\",' % labelset[idx]
        plot += '\"-\" with boxes fs solid 0.7 title \"from %s\",' % labelset[idx2]
        for idx in range(len(labelset)):
            if substitutions.has_key((labelset[idx], labelset[idx2])):
                data += (str(float(cum[idx]) / max(sums[idx], 1)) + '\n')
                cum[idx] -= substitutions[(labelset[idx], labelset[idx2])]
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

    plot = 'plot \"-\" with boxes title \"error prob.\",'
    data = '\n'
    index = error_prob_per_quality.keys()
    index.sort()
    tics = 'set xtics ('
    for tic in range(len(index)):
        if index[tic] == '\n':
            continue
        if not index[tic] in set(['\"', '%', ']', '[', '\\', '`']):
            tics += '\"%s\" %i,' % (index[tic], tic)
        else:
            tics += '\"%s\" %i,' % ('l', tic)
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
    import pickle

    options = parse_options(sys.argv)

    # dummy contains only statistics, how many lines have been processed
    ((qualities, substitutions), dummy) = pickle.load(file(options.error))

    plot_histograms(options, qualities, substitutions)


if __name__ == '__main__':
    main()

