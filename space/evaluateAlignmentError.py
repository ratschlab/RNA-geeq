#!/usr/bin/python
"""

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2009-2011 Andre Kahles
  Copyright (C) 2009-2011 by Andre Kahles
  This program evaluates the error distribution in a given error model
  
  For detailed usage information type:

    python evaluateAlignmentError.py 

"""

import sys
import re
import subprocess

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--alignment', dest='alignment', metavar='FILE', help='alignment file', default='-')
    required.add_option('-H', '--histograms', dest='histogram', metavar='FILE', help='print error histograms in pdf format into FILE', default='-')
    required.add_option('-g', '--genome', dest='genome', metavar='FILE', help='Genome sequence in FASTA format', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-s', '--samtools', dest='samtools', metavar='STRING', help='location of samtools on the system', default='samtools')
    optional.add_option('-f', '--format', dest='format', metavar='STRING', help='alignment file format sam or bed [default: sam]', default='sam')
    optional.add_option('-x', '--use-x', dest='x', action='store_true', help='Use graphical output [off]', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='Set verbose output [off]', default=False)
    
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()

    if len(argv) < 4:
        parser.print_help()
        sys.exit(2)

    return options

def read_fasta(options):

    """Parses genome information from infiles in fasta-format"""

    if options.genome.count(',') > 0:
        infiles = options.genome.split(',')
    else:
        infiles = [options.genome]

    genome = dict()
    curr_chr = ''
    for infile in infiles:
        if options.verbose:
            print >> sys.stdout, 'Processing %s' % infile
        for line in open(infile, 'r'):
            if line[0] == '>':
                chr_name = line.strip()[1:]
                genome[chr_name] = []
                curr_chr = chr_name
                if options.verbose:
                    print >> sys.stdout, '... reading %s' % curr_chr
                continue

            if curr_chr != '':
                genome[curr_chr].append(line.strip().lower())
            else:
                print >> sys.stderr, "File in %s has no valid Fasta format!\n" % infile
                exit(-1)
                
    ### doing it that way is apparently the fastest solution
    for chrm in genome.keys():
        genome[chrm] = ''.join(genome[chrm])
    return genome


def parse_alignment(options):
    """Extracts quality und substitution information from an alignment given in different alignment formats"""

    qualities = dict()
    substitutions = dict()
    quality_per_pos = dict()
    avg_quality_per_pos = dict()
    line_counter = 0
    
    if options.alignment != '-':
        if options.alignment[-3:] == 'bam' or options.format == 'bam':
            infile_handle = subprocess.Popen([options.samtools, 'view', options.alignment], stdout=subprocess.PIPE)
            infile = infile_handle.stdout
        else:
            infile = open(options.alignment, 'r')
    else:
        infile = sys.stdin

    genome = read_fasta(options)

    for line in infile:

        ### skip comments and header sections
        if line[0] in ['#', '@'] or line[:2] == 'SQ':
            continue

        if options.verbose and line_counter % 1000 == 0:
            print >> sys.stdout, 'read %i lines' % line_counter

        line_counter += 1
        sl = line.strip().split('\t')

        if options.format == 'bed':
            if len(sl) < 9:
                continue
            read = sl[12]
            quality = sl[13]
        elif options.format in ['sam', 'bam']:
            if len(sl) < 9:
                print >> sys.stderr, line
                print >> sys.stderr, 'Line incomplete: %i - skipping!\n' % line_counter
                continue

            (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
            size = [int(i) for i in size]
            chrm_pos = 0    # position in chrm
            read_pos = 0    # position in the actual read
            
            read = sl[9]
            _read = ''
            gen_start = int(sl[3]) - 1

            ### build up augmented read string containing all substitution information
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
                    _read += read[read_pos:read_pos+size[pos]]
                    read_pos += size[pos]
                    #chrm_pos += size[pos] #--> new SAM standard describes start-pos as first MATCHING position

            read = _read
            quality = sl[10]
            
        ### check for substitutions
        start = read.find('[', 0)
        while start > -1:
            ### skip indels
            if read[start + 1] == '-' or read[start + 2] == '-':
                start = read.find('[', start + 1)
                continue
            try:
                substitutions[(read[start + 1], read[start + 2])] += 1
            except KeyError:
                substitutions[(read[start + 1], read[start + 2])] = 1
            start = read.find('[', start + 1)

        ### read quality values
        if int(sl[1]) & 16 == 16:
            strand_offset = len(quality) - 1
        else:
            strand_offset = 0
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
            q_ = abs(strand_offset - q)
            try:
                avg_quality_per_pos[q_] += ord(quality[q])
            except KeyError:
                avg_quality_per_pos[q_] = ord(quality[q])
            try:
                quality_per_pos[q_][quality[q]] += 1
            except KeyError:
                try:
                    quality_per_pos[q_][quality[q]] = 1
                except KeyError:
                    quality_per_pos[q_] = {quality[q]:1}

    if options.format == 'bam':
        infile_handle.kill()

    ### average quality values
    for key in avg_quality_per_pos.keys():
        avg_quality_per_pos[key] /= line_counter

 
    return (qualities, substitutions, quality_per_pos, avg_quality_per_pos, line_counter) 

def plot_histograms(options, qualities, substitutions, quality_per_pos, avg_quality_per_pos, line_counter):
    """Plots the error distribution relative to intron positions using gnuplot. """

    import matplotlib

    if not options.x:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

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
    plt.bar(range(len(avg_quality_per_pos)), [avg_quality_per_pos[pos] for pos in sorted(avg_quality_per_pos.keys())], color='w')
    plt.title('Avg. quality per pos')
    plt.xlabel('position')
    plt.ylabel('avg. quality')

    plt.subplot(2, 2, 3)
    plt.bar(range(len(error_prob_per_pos)), [error_prob_per_pos[pos] for pos in sorted(error_prob_per_pos.keys())])
    plt.title('Error probability per pos')
    plt.xlabel('position')
    plt.ylabel('error prob.')

    plt.subplot(2, 2, 4)
    plt.bar(range(len(error_prob_per_quality)), [error_prob_per_quality[qual] for qual in sorted(error_prob_per_quality.keys())])
    plt.title("Average error probability per quality") 
    plt.xlabel("quality value")
    plt.ylabel("error prob.")

    plt.subplots_adjust(hspace=0.4)
    plt.savefig(options.histogram, format='pdf')
    plt.show()


def main():
    """Function controlling the program flow ... """
    options = parse_options(sys.argv)
    if options.histogram == '-':
        print >> sys.stderr, "Please provide an output filename for the plots! - Use option -H <file>"
        sys.exit(-1)

    (qualities, substitutions, quality_per_pos, avg_quality_per_pos, line_counter) = parse_alignment(options)

    plot_histograms(options, qualities, substitutions, quality_per_pos, avg_quality_per_pos, line_counter)


if __name__ == '__main__':
    main()

