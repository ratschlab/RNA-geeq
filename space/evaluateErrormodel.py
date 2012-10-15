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

    ### build up list of quality distribution per read length
    quality_per_pos = dict()
    avg_quality_per_pos = []
    line_counter = 0
    for line in open(options.quality, 'r'):
        line_counter += 1
        line.strip()
        for c in range(len(line)):
            try:
                avg_quality_per_pos[c] += max(ord(line[c]) - options.offset, 0)
            except IndexError:
                avg_quality_per_pos.extend([0] * (c - len(avg_quality_per_pos) - 1))
                avg_quality_per_pos.append(max(ord(line[c]) - options.offset, 0))
            
            try:
                quality_per_pos[c][line[c]] += 1
            except KeyError:
                try:
                    quality_per_pos[c][line[c]] = 1
                except KeyError:
                    quality_per_pos[c] = {line[c]:1}


    ### average quality values
    for idx in range(len(avg_quality_per_pos)):
        avg_quality_per_pos[idx] /= float(line_counter)

    ### build up list of error probabilities per position and error probabilities per quality value
    ### we assume, that the nucleotides are equally distributed 
    error_prob_per_pos = []
    error_prob_per_quality = []
    for pos in sorted(quality_per_pos.keys()):
        
        for qual in quality_per_pos[pos].keys():
            avg_substitution_rate = 0
            for base in ['A', 'C', 'G', 'T']:
                try:
                    avg_substitution_rate +=  (0.25 * float(qualities[(base, qual)][0]) / float(max(qualities[(base, qual)][1], 1)))
                except KeyError:
                    pass
            try:
                error_prob_per_pos[pos] += (float(quality_per_pos[pos][qual]) / float(line_counter) * avg_substitution_rate)
            except IndexError:
                error_prob_per_pos.extend([0] * (pos - len(avg_quality_per_pos) - 1))
                error_prob_per_pos.append(float(quality_per_pos[pos][qual]) / float(line_counter) * avg_substitution_rate)

            try:
                error_prob_per_quality[qual] = avg_substitution_rate 
            except IndexError:
                error_prob_per_quality.extend([0] * (qual - len(error_prob_per_quality) - 1))
                error_prob_per_quality[qual] = avg_substitution_rate 

    import matplotlib
    if not options.x:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ### pre-processing of substitution matrix
    labelset = ['A', 'C', 'G', 'T']
    sub_sum = dict()
    for idx in labelset:
        for idx2 in labelset:
            if not substitutions.has_key((idx, idx2)):
                substitutions[(idx, idx2)] = 0
        sub_sum[idx] = sum([substitutions[idx, idx2] for idx2 in labelset]) 
    for idx in labelset:
        for idx2 in labelset:
            substitutions[idx, idx2] /= float(sub_sum[idx])
 
    plt.subplot(2, 2, 1)
    colors = {'A':'red', 'C':'green', 'G':'blue', 'T':'yellow'}
    cum_sum = [0, 0, 0, 0]
    p = []
    for label2 in labelset:
        p.append(plt.bar(range(4), [substitutions[label1, label2] for label1 in labelset], bottom=cum_sum, color=colors[label2]))
        cum_sum = [cum_sum[i] + substitutions[labelset[i], label2] for i in range(len(labelset))]
    plt.title('Substitutions')
    plt.xlabel('to base')
    plt.ylabel('rel. number of substitutions')
    plt.xlim((0, 5.5))
    plt.xticks([0.5, 1.5, 2.5, 3.5], ('A', 'C', 'G', 'T'))
    plt.legend((p[0], p[1], p[2], p[3]), ('A', 'C', 'G', 'T'), loc=1, borderaxespad=0.)

    plt.subplot(2, 2, 2)
    plt.bar(range(len(avg_quality_per_pos)), avg_quality_per_pos, color='w')
    plt.title('Avg. quality per pos')
    plt.xlabel('position')
    plt.ylabel('avg. quality')
    
    plt.subplot(2, 2, 3)
    plt.bar(range(len(error_prob_per_pos)), error_prob_per_pos)
    plt.title('Error probability per pos')
    plt.xlabel('position')
    plt.ylabel('error prob.')

    plt.subplot(2, 2, 4)
    plt.bar(range(len(error_prob_per_quality)), error_prob_per_quality)
    plt.title("Average error probability per quality") 
    plt.xlabel("quality value")
    plt.ylabel("error prob.")

    plt.subplots_adjust(hspace=0.4)
    plt.savefig(options.histogram, format='pdf')
    plt.show()

def main():
    """ Main calling routine"""
    import pickle

    options = parse_options(sys.argv)

    ((qualities, substitutions), dummy) = pickle.load(file(options.error))

    plot_histograms(options, qualities, substitutions)


if __name__ == '__main__':
    main()

