"""This script compares an alignment junction list to annotated junctions. """

import cPickle
import sys
import re
import subprocess
import plot

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-I', '--annotation_introns', dest='anno_int', metavar='FILE', help='annotation intron list', default='-')
    required.add_option('-C', '--annotation_coverage', dest='anno_cov', metavar='FILE', help='annotation coverage', default='-')
    required.add_option('-i', '--alignment_introns', dest='align_int', metavar='FILE', help='alignment intron list', default='-')
    required.add_option('-c', '--alignment_coverage', dest='align_cov', metavar='FILE', help='alignment coverage', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-a', '--alignment', dest='align', metavar='FILE', help='alignment file in sam format - can be given instead of alignment intron file', default='-')
    optional.add_option('-S', '--store', dest='store', action='store_true', help='stores generated intron lists and coverage maps', default=False)
    optional.add_option('-s', '--sam', dest='sam', action='store_true', help='input alignment is in sam format not in bam', default=False)
    optional.add_option('-F', '--full_analysis', dest='full_analysis', action='store_true', help='generates overlap and coverage plots', default=False)
    optional.add_option('-s', '--strand_specific', dest='strand_specific', action='store_true', help='data is strand specific', default=False)
    optional.add_option('-R', '--ignore_multireads', dest='multireads', metavar='FILE', help='file containing the multireads to ignore', default='-')
    optional.add_option('-p', '--performance_log', dest='performance_log', action='store_true', help='store the intron recovery performance in extra log [off]', default=False)
    optional.add_option('-g', '--coverage_gff', dest='cov_gff', action='store_true', help='store coverage info as GFF [off]', default=False)
    optional.add_option('-M', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [100000000]', default='100000000')
    optional.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    optional.add_option('-X', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismathes [10000]', default=10000)
    optional.add_option('-m', '--min_read_coverage', dest='min_coverage', metavar='INT', type='int', help='minimal coverage to count an intron [1]', default='1')
    optional.add_option('-O', '--max_overlap', dest='max_overlap', metavar='INT', type='int', help='maximal evalutated exon overlap [5]', default='5')
    optional.add_option('-H', '--histogram', dest='histogram', metavar='PATH', help='destination PATH of evaluation plots', default='-')
    optional.add_option('-x', '--exclude_chrms', dest='exclude_chrm', metavar='STRINGLIST', help='list of comma separated chromosomes to exclude from evaluation', default='-')
    #optional.add_option('-A', '--align_out', dest='align_out', metavar='PATH', help='destination for generated alignment intron list', default='-')
    optional.add_option('-o', '--outfile_base', dest='outfile_base', metavar='PATH', help='basedir for outfiles written', default='-')
    optional.add_option('-d', '--details', dest='details', action='store_true', help='destination for generated alignment intron list', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 3 or '-' in [options.anno_int, options.anno_cov] :
        parser.print_help()
        sys.exit(2)

    return options

def print_to_datafile(data, file):
    """Prints data to file"""

    outf = open(file, 'w')
    _line = ''
    for d in data:
        _line = (str(d) + '\t')
    print >> outf, _line[:-1]
    outf.close()

def build_intron_list(options, multireads):
    """Builds up an intron list from the given alignment and writes it to a file. Returns a map of all covered intron positions """

    ### parse SAM/BAM:
    ### 0      1     2       3             4     5      6        7        8       9    10    11
    ### QUERY  FLAG  REFSEQ  POS(1-based)  MAPQ  CIGAR  MATEREF  MATEPOS  INSIZE  SEQ  QUAL  OPT
    intron_lists = dict()
    cov_map = dict()
    counter = 0
    filter_counter = 0

    if options.align == '-':    
        infile = sys.stdin
    else:
        if options.sam:
            infile = open(options.align, 'r')
        else:
            file_handle = subprocess.Popen([options.samtools, 'view', options.align], stdout=subprocess.PIPE) 
            infile = file_handle.stdout

    for line in infile:
        if line[0] in ['@', '#'] or line[:2] == 'SQ':
            continue
        if counter % 10000 == 0 and options.verbose:
            print 'lines read: [ %s (taken: %s / filtered: %s)]' % (counter, counter - filter_counter, filter_counter)
        counter += 1

        sl = line.strip().split('\t')
    
        if (sl[0], (int(sl[1]) & 128 + int(sl[1]) & 64)) in multireads:
            filter_counter += 1
            continue

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

        if not cov_map.has_key(sl[2]):
            cov_map[sl[2]] = dict()

        ### process CIGAR string

        ### handle Insertions - they do not affect intron length
        rl = sl[5]
        rl = re.sub('[0-9]*I', '', rl)
        
        (op, size) = (re.split('[0-9]*', rl)[1:], re.split('[^0-9]', rl)[:-1])
        size = [int(i) for i in size]
        offset = int(sl[3]) - 1
        covered = 0
        for p in range(len(op)):
            if op[p] in ['M', 'D']:
                ### handle starting position
                try:
                    ### cov_map[chromosome][position] = [number of all alignments, number of starting alignments]
                    cov_map[sl[2]][offset + covered][0] += 1
                    cov_map[sl[2]][offset + covered][1] += 1
                except KeyError:
                    cov_map[sl[2]][offset + covered] = [1, 1]
                covered += 1
                ### handle following positions
                for pp in range(1, size[p]):                
                    try:
                        cov_map[sl[2]][offset + covered][0] += 1
                    except KeyError:
                        cov_map[sl[2]][offset + covered] = [1, 0]
                    covered += 1
            elif op[p] == 'N':
                ### update intron lists
                ### intron = (strand [-1, 1], start, end)  ... end in the pythonic way
                if options.strand_specific:
                    if (int(sl[1]) & 128) == 128:
                        intron = ((((int(sl[1]) & 16) / -8) + 1) * -1, offset + covered, offset + covered + size[p])
                    else:
                        intron = (((int(sl[1]) & 16) / -8) + 1, offset + covered, offset + covered + size[p])
                else:
                    intron = (0, offset + covered, offset + covered + size[p])
                ### intron_lists[chromosome][intron] = [read1, read2, ...]
                try:
                    intron_lists[sl[2]][intron].add(sl[0])
                except KeyError:
                    try:
                        intron_lists[sl[2]][intron] = set([sl[0]])
                    except KeyError:
                        intron_lists[sl[2]] = {intron:set([sl[0]])}
 
                covered += size[p]
            ### handle softclips and hardclips --> means ignoring them :)
            elif op[p] == 'S':
                covered += size[p]
            elif op[p] == 'H':
                continue
    
    print '\n parsed file %s' % options.align
    print 'read %s lines' % counter
    if options.max_mismatches < 10000 or options.min_exon_len > 0:
        print 'filter criteria:'
        print '    max mismatches: %s' % options.max_mismatches
        print '    min exon length: %s\n' % options.max_mismatches
        print 'filtered %s lines' % filter_counter

    return (cov_map, intron_lists)


def evaluate_intron_overlap(options, annotation, not_matched, plotlist):
    """This part evaluates the overlap between annotation and alignment at the intron boundaries."""

    total_overlap = 0
    total_enclosing = 0
    total_introns = 0
    total_offset = 0

    not_matched_idx = not_matched.keys()
    not_matched_idx.sort()

    distance_5p = dict()
    distance_3p = dict()

    for chrm in not_matched_idx:
        ### load introns
        align_idx = list(not_matched[chrm])
        anno_idx = annotation[chrm].keys()

        ### sort both intron lists by starting position and intron length
        align_idx.sort(lambda u, v: v[2] - u[2])
        align_idx.sort(lambda u, v: u[1] - v[1])
        anno_idx.sort(lambda u, v: v[2] - u[2])
        anno_idx.sort(lambda u, v: u[1] - v[1])
        
        ### initialize counters
        counter_start = 0
        enclosing_introns = 0 #set()
        overlap_introns = 0
        completely_offset = 0
        
        min_anno = 0

        ### we search for each intron not matching the annotation exactly, for wich of the annotated introns it laps least into an exon
        for align in align_idx:
            ### go to next annotated intron overlapping the current one
            while counter_start < len(anno_idx) and anno_idx[counter_start][2] < align[1] :
                counter_start += 1
           
            anno_counter = counter_start
            min_overlap = 10000000
            if anno_counter < len(anno_idx):
                min_5p_overlap = abs(anno_idx[anno_counter][1] - align[1])
                min_3p_overlap = abs(anno_idx[anno_counter][2] - align[2])
            else: # take overlap to last annotated intron
                min_5p_overlap = abs(anno_idx[-1][1] - align[1])
                min_3p_overlap = abs(anno_idx[-1][2] - align[2])

            no_match = True

            ### iterate over all overlapping introns and search for the one with the smallest overlap
            ### (overlap means in this context the amount of the read lapping over the boundaries of the intron)
            while (anno_counter < len(anno_idx)) and (anno_idx[anno_counter][1] < align[2]) :
                no_match = False 
                if abs(anno_idx[anno_counter][1] - align[1]) + abs(anno_idx[anno_counter][2] - align[2]) < min_overlap:
                    min_anno = anno_counter
                    min_overlap = abs(anno_idx[anno_counter][1] - align[1]) + abs(anno_idx[anno_counter][2] - align[2]) 
                ### min 5p overlap
                if abs(anno_idx[anno_counter][1] - align[1]) < min_5p_overlap:
                    min_5p_overlap = abs(anno_idx[anno_counter][1] - align[1])
                ### min 3p overlap
                if abs(anno_idx[anno_counter][2] - align[2]) < min_3p_overlap:
                    min_3p_overlap = abs(anno_idx[anno_counter][2] - align[2])

                anno_counter += 1
            
            ### handle case, where no annotated intron overlaps and check, if next intron would be nearer than the last
            if no_match and (anno_counter + 1) < len(anno_idx):
                if abs(anno_idx[anno_counter + 1][1] - align[1]) + abs(anno_idx[anno_counter + 1][2] - align[2]) < min_overlap:
                    min_anno = anno_counter + 1
                    min_overlap = abs(anno_idx[min_anno][1] - align[1]) + abs(anno_idx[min_anno][2] - align[2]) 
                min_5p_overlap = min(min_5p_overlap, abs(anno_idx[anno_counter + 1][1] - align[1]))
                min_3p_overlap = min(min_3p_overlap, abs(anno_idx[anno_counter + 1][2] - align[2]))
    

            ### update distance dict
            try:
                distance_5p[min_5p_overlap] += 1
            except KeyError:
                distance_5p[min_5p_overlap] = 1
            try:
                distance_3p[min_3p_overlap] += 1
            except KeyError:
                distance_3p[min_3p_overlap] = 1

            ### test for inclusion
            if (anno_idx[min_anno][1] <= align[1] and anno_idx[min_anno][2] >= align[2]):
                enclosing_introns += 1
            ### test for overlap
            elif (anno_idx[min_anno][1] <= align[1] and anno_idx[min_anno][2] > align[1]) or \
                (anno_idx[min_anno][2] >= align[2]):
                overlap_introns += 1
            else:
                completely_offset += 1
                continue

            if options.details:
                print '---------------------------'
                print 'in chrm %s: %s' % (chrm, annotation[chrm][anno_idx[min_anno]])
                print 'overlap: %s %s' % (anno_idx[min_anno][1] - align[1], anno_idx[min_anno][2] - align[2]) 
                print 'length difference: %s' % ((anno_idx[min_anno][2] - anno_idx[min_anno][1]) - (align[2] - align[1])) 
                print 'read set: %s' % str(align)
                print 'anno set: %s' % str(anno_idx[min_anno])

        print '---------------------------'
        print 'in chrm %s: found %s introns overlapping an annotated intron (for %s introns not matching the annotation)' % (chrm, overlap_introns, len(not_matched[chrm]))
        print 'in chrm %s: found %s introns enclosing an annotated intron (for %s introns not matching the annotation)' % (chrm, enclosing_introns, len(not_matched[chrm]))
        print 'in chrm %s: found %s introns not overlapping to any annotated intron (for %s introns not matching the annotation)' % (chrm, completely_offset, len(not_matched[chrm]))
        total_overlap += overlap_introns
        total_enclosing += enclosing_introns
        total_offset += completely_offset
        total_introns += len(not_matched[chrm])

    print '---------------------------\n---------------------------'
    print 'TOTAL: found %s introns overlapping an annotated intron (for %s introns not matching the annotation)' % (total_overlap, total_introns)
    print '       found %s introns enclosing an annotated intron (for %s introns not matching the annotation)' % (total_enclosing, total_introns)
    print '       found %s introns not overlapping to any annotated intron (for %s introns not matching the annotation)' % (total_offset, total_introns)
    print '---------------------------'

    
    ### write plotlists
    if len(distance_5p.keys()) == 0:
        distance_5p[0] = 0
    if len(distance_3p.keys()) == 0:
        distance_3p[0] = 0
    
    plotlist.append(plot.Plot('Splice site distance distribution 3 prime', distance_3p, xlabel='distance', ylabel='occurrences', xlog=True, ylog=True))
    plotlist.append(plot.Plot('Splice site distance distribution 5 prime', distance_5p, xlabel='distance', ylabel='occurrences', xlog=True, ylog=True))

def evaluate_exon_overlap(options, plotlist, anno_exon_map, align_cov_map, total_alignments):
    """This part evaluates the 3 prime and 5 prime overlap to the annotated exons. """

    ### count total positions
    ### build up all covered positions
    unique_exon_map = dict()
    anno_total_position_set = dict()

    for chrm in anno_exon_map.keys():
        unique_exon_map[chrm] = dict()
        anno_total_position_set[chrm] = set()
        ### iterate over all annotated transcripts
        for trans in anno_exon_map[chrm].keys():
            ### iterate over contained exons
            for exon in anno_exon_map[chrm][trans].keys():
                try:
                    if unique_exon_map[chrm][exon] == 3:
                        continue
                    else:
                        unique_exon_map[chrm][exon] |= anno_exon_map[chrm][trans][exon]
                except KeyError:
                    unique_exon_map[chrm][exon] = anno_exon_map[chrm][trans][exon]
                anno_total_position_set[chrm].union(set(range(exon[0], exon[1])))

    max_overlap = options.max_overlap

    p5_overlap = [0 for i in range(max_overlap + 1)]
    p3_overlap = [0 for i in range(max_overlap + 1)]
    intronic_overlap = [0 for i in range(max_overlap + 1)]    
    p5_intronic_overlap = [0 for i in range(max_overlap + 1)]
    p3_intronic_overlap = [0 for i in range(max_overlap + 1)]


    for chrm in align_cov_map.keys():
        if not unique_exon_map.has_key(chrm):
            continue
        for exon in unique_exon_map[chrm].keys():
            try:
                last_cov = align_cov_map[chrm][exon[1] + max_overlap]
            except KeyError:
                last_cov = [0, 0]

            for overlap in range(max_overlap, 0, -1):
                ### handle 5 prime
                ### count exons starting at position
                if not (exon[0] - overlap) in anno_total_position_set[chrm]:
                    try:
                        add_cov = align_cov_map[chrm][exon[0] - overlap][1]
                        p5_overlap[overlap] += add_cov
                        if unique_exon_map[chrm][exon] in [2, 3]:
                            p5_intronic_overlap[overlap] += add_cov
                    except:
                        pass

                ### handle 3 prime
                ### count exons ending at position
                if not (exon[1] - 1 + overlap) in anno_total_position_set[chrm]:
                    try:
                        curr_cov = align_cov_map[chrm][exon[1] - 1 + overlap]
                        add_cov = curr_cov[0] + (last_cov[1] - last_cov[0])
                        last_cov = curr_cov
                        p3_overlap[overlap] += add_cov
                        if unique_exon_map[chrm][exon] in [1, 3]:
                            p3_intronic_overlap[overlap] += add_cov
                    except KeyError:
                        last_cov = [0, 0]
                else:
                    try:
                        last_cov = align_cov_map[chrm][exon[1] - 1 + overlap]
                    except KeyError:
                        last_cov = [0, 0]

    intronic_overlap = []
    for overlap in range(max_overlap):
        intronic_overlap.append(p5_intronic_overlap[overlap] + p3_intronic_overlap[overlap])

    plotlist.append(plot.Plot('5 prime overlap to exon boundary (abs)', p5_overlap, xlabel='overlap', ylabel='occurrence', plottype='boxes'))
    plotlist.append(plot.Plot('3 prime overlap to exon boundary (abs)', p3_overlap, xlabel='overlap', ylabel='occurrence', plottype='boxes'))
    plotlist.append(plot.Plot('5 prime overlap to exon boundary into intron (abs)', p5_intronic_overlap, xlabel='overlap', ylabel='occurrence', plottype='boxes'))
    plotlist.append(plot.Plot('3 prime overlap to exon boundary into intron (abs)', p3_intronic_overlap, xlabel='overlap', ylabel='occurrence', plottype='boxes'))
    plotlist.append(plot.Plot('total intronic overlap to exon boundary (abs)', intronic_overlap, xlabel='overlap', ylabel='occurrence', plottype='boxes'))

    ### build up relative measurement (division by total_covered_positions)
    p5_overlap_rel = []
    for ol in p5_overlap:
        p5_overlap_rel.append(float(ol) / max(float(total_alignments), 1.0))
    p3_overlap_rel = []
    for ol in p3_overlap:
        p3_overlap_rel.append(float(ol) / max(float(total_alignments), 1.0))
    p5_intronic_overlap_rel = []
    for ol in p5_intronic_overlap:
        p5_intronic_overlap_rel.append(float(ol) / max(float(total_alignments), 1.0))
    p3_intronic_overlap_rel = []
    for ol in p3_intronic_overlap:
        p3_intronic_overlap_rel.append(float(ol) / max(float(total_alignments), 1.0))
    intronic_overlap_rel = []
    for ol in intronic_overlap:
        intronic_overlap_rel.append(float(ol) / max(float(total_alignments), 1.0))

    plotlist.append(plot.Plot('5 prime overlap to exon boundary (rel)', p5_overlap_rel, xlabel='overlap', ylabel='occurrence', plottype='boxes'))
    plotlist.append(plot.Plot('3 prime overlap to exon boundary (rel)', p3_overlap_rel, xlabel='overlap', ylabel='occurrence', plottype='boxes'))
    plotlist.append(plot.Plot('5 prime overlap to exon boundary into intron (rel)', p5_intronic_overlap_rel, xlabel='overlap', ylabel='occurrence', plottype='boxes'))
    plotlist.append(plot.Plot('3 prime overlap to exon boundary into intron (rel)', p3_intronic_overlap_rel, xlabel='overlap', ylabel='occurrence', plottype='boxes'))
    plotlist.append(plot.Plot('total intronic overlap to exon boundary (rel)', intronic_overlap_rel, xlabel='overlap', ylabel='occurrence', plottype='boxes'))



def evaluate_coverage(options, align_cov_map, plotlist):
    """Function for evaluation of positional and transcript coverage. """
    
    anno_exon_map = cPickle.load(open(options.anno_cov, 'r'))

    ### initialize counting dictionaries
    covered_transcripts = dict([(i * 10, 0) for i in range(11)]) 

    transcript_count = 0

    for chrm in align_cov_map.keys():
        ### iterate over annotated transcripts
        if not anno_exon_map.has_key(chrm):
            continue
        transcript_count += len(anno_exon_map[chrm].keys())
        ### align_cov_set = set(covered position in chromosome 'chrm')
        align_cov_set = set(align_cov_map[chrm].keys())
        for trans in anno_exon_map[chrm].keys():
            anno_trans = set()
            ### iterate over positions covered by annotated transcript
            for exon in anno_exon_map[chrm][trans].keys():
                anno_trans = anno_trans.union(set(range(exon[0], exon[1])))

            contig_len = len(anno_trans)
            intersect_len = len(anno_trans.intersection(align_cov_set))

            ### add to transcript-wise coverage map
            for i in range(100, -10, -10):
                if float(intersect_len) / float(contig_len) * 100 >= i:
                    covered_transcripts[i] += 1
                    break

    ### counter for all covered positions
    total_alignments = 0

    ### coverage distribution
    cov_distr = dict()
    ### align_cov_map[chromosom][position] = [number of all alignments, number of starting alignments]
    for chrm in align_cov_map.keys():
        for key in align_cov_map[chrm].keys():
            try:
                cov_distr[align_cov_map[chrm][key][0]] += 1
            except KeyError:
                cov_distr[align_cov_map[chrm][key][0]] = 1
            total_alignments += align_cov_map[chrm][key][1]

    ### analysis of covered  but not annotated positions
    #non_covered_pos = set(align_cov_map.keys()).difference(anno_total_positions_set)
    
    plotlist.append(plot.Plot('Coverage distribution', cov_distr, xlabel='coverage', ylabel='occurence', ylog=True))

    print 'Coverage analysis: '
    print '--------------------------'

    print '%s of %s transcripts were covered to 100 percent' % (covered_transcripts[100], transcript_count)
    for c in range(90, -10, -10):
        print '%s of %s transcripts were covered between %i and %i percent' % (covered_transcripts[c], transcript_count, c, c + 10)

    plotlist.append(plot.Plot('Transcript coverage', covered_transcripts, xlabel='percent covered', ylabel='number of transcripts'))

    ### analyze exon overlaps
    evaluate_exon_overlap(options, plotlist, anno_exon_map, align_cov_map, total_alignments)

def build_intron_support_list(alignment_list, plotlist):
    """Builds up a dictionairy with number of supp. reads as keys and number of introns as values. """
    intron_support_list = dict()

    max_len = 0
    for chrm in alignment_list.keys():
        for intron in alignment_list[chrm]:
            try:
                intron_support_list[len(alignment_list[chrm][intron])] += 1
            except KeyError:
                intron_support_list[len(alignment_list[chrm][intron])] = 1
            max_len = max(len(alignment_list[chrm][intron]), max_len)

    ### sanitize for plotting
    diff = set(range(1, max_len)).difference(set(intron_support_list.keys()))
    for d in diff:
        intron_support_list[d] = 0

    plotlist.append(plot.Plot('Reads supporting an intron', intron_support_list, xlabel='number of supp. reads', ylabel='number of introns', xlog=True, ylog=True))

def store_intron_lists(options, align_cov_map, alignment_list):
    """Stores intron information for later usage."""
    
    # TODO STUB
    pass

    #outfile = open(options.align + ".introns", "w")

    
    #outfile.close()

def build_coverage_gff(align_cov_map, outfile_base):

    gff_out = open(outfile_base + '_coverage.gff', 'w')
    for chrm in align_cov_map.keys():
        counter = 0
        last_pos = -2
        last_cov = -1
        region = None
        subregion = None
        for pos in sorted(align_cov_map[chrm].keys()):
            cov = align_cov_map[chrm][pos][0]
            ### extend region
            if pos == last_pos + 1:
                region[1] += 1
                ### extend subregion
                if cov == last_cov:
                    subregion[-1][2] += 1
                ### start new subregion
                else:
                    subregion.append([pos, cov, 1])
                    last_cov = cov
            ### start new region
            else:
                counter += 1
                last_cov = -1
                if region != None and subregion != None:
                    _mean_cov = sum([i[1]*i[2] for i in subregion])
                    _mean_cov /= region[1]
                    print >> gff_out, '%s\tFML_evalTools\tCoverage_region\t%s\t%s\t%s\t.\t.\tid "region_%s_%s"' % (chrm, str(region[0] + 1), str(region[0] + region[1]), _mean_cov, chrm, counter)
                    for i in range(len(subregion)):
                        print >> gff_out, '%s\tFML_evalTools\tSubregion\t%s\t%s\t%s\t.\t.\tid "subregion_%s_%s_%s"; parent_id "region_%s_%s"' % (chrm, str(subregion[i][0] + 1), str(subregion[i][0] + subregion[i][2]), min(subregion[i][1], 1000), chrm, counter, str(i + 1), chrm, counter)
                        
                region = [pos, 1]
                subregion = [[pos, cov, 1]]
                last_cov = cov
            last_pos = pos

    gff_out.close()

def main():
    """Main function for comparing two intron lists """    
    
    options = parse_options(sys.argv)

    plotlist = plot.Plotlist()


    if options.outfile_base != '-':
        if options.align == '-':
            outfile_base = (options.outfile_base + '/' + options.align_int.split('/')[-1])
        else:
            outfile_base = (options.outfile_base + '/' + options.align.split('/')[-1])
    else:
        if options.align == '-':
            outfile_base = options.align_int
        else:
            outfile_base = options.align

    ### set filter tags for filenames
    if options.min_exon_len > 0:
        outfile_base += '_me%s' % options.min_exon_len
    if options.max_mismatches < 10000:
        outfile_base += '_mm%s' % options.max_mismatches
    if options.max_intron_len < 100000000:
        outfile_base += '_mi%s' % options.max_intron_len
    if options.min_coverage > 1:
        outfile_base += '_mc%s' % options.min_coverage

    if options.align_int != '-':
        alignment_list = cPickle.load(open(options.align_int, 'r'))

        ### load alignment coverage vector
        align_cov_vec = cPickle.load(open(options.align_cov, 'r'))

        ### gen align_cov_map from vector
        align_cov_map = dict()
        for chrm in align_cov_vec.keys():
            for idx in range(len(align_cov_vec[chrm])):
                align_cov_map[chrm][align_cov_vec[chrm][0]] = [align_cov_vec[chrm][1], align_cov_vec[chrm][2]]
    else:
        ### check for multireads to handle
        multireads = set()
        if options.multireads != '-':
            print '\nParsing multireads from file %s' % options.multireads
            print '-----------------------------------------'
            for line in open(options.multireads, 'r'):
                _l = line.strip().split('\t')
                multireads.add((_l[0], int(_l[1])))
        ### parse alignmentfile for creating intron_list and coverage_map
        (align_cov_map, alignment_list) = build_intron_list(options, multireads)
        
        if options.store:
            store_intron_lists(options, align_cov_map, alignment_list)
    
    annotation_list = cPickle.load(open(options.anno_int, 'r'))

    ### filter annotated and predicted introns for excluded chromosomes
    if options.exclude_chrm != '-':
        _ex_chrm = options.exclude_chrm.strip().split(',')
        ### handle leading or trailing commas
        if _ex_chrm[0] == '':
            _ex_chrm = _ex_chrm[1:]
        if _ex_chrm[-1] == '':
            _ex_chrm = _ex_chrm[:-1]
        for chrm in _ex_chrm:
            if annotation_list.has_key(chrm):
                del annotation_list[chrm]

            if alignment_list.has_key(chrm):
                del alignment_list[chrm]

    ### filter intron lists for max intron length
    print '\nFiltering intron list for max intron len'
    print '-----------------------------------------'
    skipped = 0
    for chrm in annotation_list.keys():
        skiplist = set()
        for intron in annotation_list[chrm].keys():
            if (intron[2] - intron[1]) > options.max_intron_len:
                skiplist.add(intron)
        for intron in skiplist:
            del annotation_list[chrm][intron]
        skipped += len(skiplist)
    print '%s introns removed from annotation' % skipped

    skipped = 0
    for chrm in alignment_list.keys():
        skiplist = set()
        for intron in alignment_list[chrm].keys():
            if (intron[2] - intron[1]) > options.max_intron_len:
                skiplist.add(intron)
        for intron in skiplist:
            del alignment_list[chrm][intron]
        skipped += len(skiplist)
    print '%s introns removed from alignment' % skipped
    del skiplist

    ### filter intron lists for min coverage
    if options.min_coverage > 1:
        print '\nFiltering intron list for min support '
        print '-----------------------------------------'
        skipped = 0
        for chrm in alignment_list.keys():
            skiplist = set()
            for intron in alignment_list[chrm].keys():
                if len(alignment_list[chrm][intron]) < options.min_coverage:
                    skiplist.add(intron)
            for intron in skiplist:
                del alignment_list[chrm][intron]
                if annotation_list[chrm].has_key(intron):
                    del annotation_list[chrm][intron]
            skipped += len(skiplist)
        print '%s introns removed from alignment\n' % skipped
        del skiplist


    ### match intron lists
    non_matched = dict()
    total_precision = float(0)
    total_recall = float(0)
    key_count = 0
    for chrm in annotation_list.keys():
        if alignment_list.has_key(chrm):
            matches = len(set(annotation_list[chrm].keys()).intersection(set(alignment_list[chrm].keys())))
            #if options.details:
            #    print (set(annotation_list[chrm].keys()) - set(annotation_list[chrm].keys()).intersection(set(alignment_list[chrm].keys())))
            print '-----------------------------'
            print ' in Chromosome %s ' % chrm 
            print '-----------------------------'
            print 'recall: %s' % (float(matches) / float(max(1, len(annotation_list[chrm].keys()))))
            print 'precision: %s' % (float(matches) / float(max(1, len(alignment_list[chrm].keys()))))
            total_precision += (float(matches) / float(max(1, len(alignment_list[chrm].keys()))))
            total_recall += (float(matches) / float(max(1, len(annotation_list[chrm].keys()))))
            non_matched[chrm] = set(alignment_list[chrm].keys()).difference(set(annotation_list[chrm].keys()))
            ### do not include chromosomes with zero values into average
            if matches > 0:
                key_count += 1
   
    total_precision /= max(1.0, float(key_count))
    total_recall /= max(1.0, float(key_count))
    print '-----------------------------'
    print ' average over all chromosomes '
    print '-----------------------------'
    print 'recall: %s' % total_recall
    print 'precision: %s' % total_precision

    if options.performance_log:
        outf = open(outfile_base + '_performance.log', 'w')
        non_matched = dict()
        total_precision = float(0)
        total_recall = float(0)
        key_count = 0
        _recall_line = ''
        _precision_line = ''
        _header_line = ''
        for chrm in annotation_list.keys():
            if alignment_list.has_key(chrm):
                matches = len(set(annotation_list[chrm].keys()).intersection(set(alignment_list[chrm].keys())))
                #print >> outf, '-----------------------------'
                #print >> outf, ' in Chromosome %s ' % chrm 
                #print >> outf, '-----------------------------'
                _header_line += (str(chrm) + '\t')
                _recall_line += (str(float(matches) / float(max(1, len(annotation_list[chrm].keys())))) + '\t')
                _precision_line += (str(float(matches) / float(max(1, len(alignment_list[chrm].keys())))) + '\t')
                total_precision += (float(matches) / float(max(1, len(alignment_list[chrm].keys()))))
                total_recall += (float(matches) / float(max(1, len(annotation_list[chrm].keys()))))
                non_matched[chrm] = set(alignment_list[chrm].keys()).difference(set(annotation_list[chrm].keys()))
                ### do not include chromosomes with zero values into average
                if matches > 0:
                    key_count += 1
        total_precision /= max(1.0, float(key_count))
        total_recall /= max(1.0, float(key_count))
       # print >> outf, '-----------------------------'
       # print >> outf, ' average over all chromosomes '
       # print >> outf, '-----------------------------'
       # print >> outf, 'recall: %s' % total_recall
       # print >> outf, 'precision: %s' % total_precision
        _header_line += 'average'
        _recall_line += str(total_recall)
        _precision_line += str(total_precision)
        print >> outf, _header_line
        print >> outf, _recall_line
        print >> outf, _precision_line
        outf.close() 
    
    ### get coverage gff
    if options.cov_gff:
        build_coverage_gff(align_cov_map, outfile_base)

    if options.full_analysis:
        ### evaluate overlap only for those reads not exactly matching the annotation introns
        evaluate_intron_overlap(options, annotation_list, non_matched, plotlist)

        ### coverage analysis of annotated introns
        evaluate_coverage(options, align_cov_map, plotlist)

        ### build list of intron supporting reads
        build_intron_support_list(alignment_list, plotlist)

        ### plot histograms
        if options.histogram != '-':
            plotlist.plot(options.histogram)

        ### write plotlists to file, for later evaluation
        plot_out = (outfile_base + '.introns.plotlist')
        cPickle.dump(plotlist, open(plot_out, 'w'))
    
if __name__ == '__main__':
    main()
