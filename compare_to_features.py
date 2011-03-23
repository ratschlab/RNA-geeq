"""

Authors: Andre Kahles
Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

This script compares an alignment to the annotation. """

import cPickle
import sys
import plot
import pdb

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-i', '--annotation_introns', dest='anno_int', metavar='FILE', help='annotation intron list', default='-')
    required.add_option('-f', '--features', dest='features', metavar='FILE', help='alignment intron features', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-s', '--strand_specific', dest='strand_specific', action='store_true', help='data is strand specific', default=False)
    optional.add_option('-S', '--splice_consensus', dest='splice_consensus', action='store_true', help='enforce introns to have splice consensus', default=False)
    optional.add_option('-p', '--performance_log', dest='performance_log', action='store_true', help='store the intron recovery performance in extra log [off]', default=False)
    optional.add_option('-I', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [100000000]', default='100000000')
    optional.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    optional.add_option('-X', '--max_feat_mismatches', dest='max_feat_mismatches', metavar='INT', type='int', help='max number of mismatches for feat generation[80]', default=80)
    optional.add_option('-M', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismacthes [1000000]', default=1000000)
    optional.add_option('-c', '--min_coverage', dest='min_coverage', metavar='INT', type='int', help='minimal coverage to count an intron [0]', default='0')
    optional.add_option('-C', '--min_coverage_2', dest='min_coverage_2', metavar='INT', type='int', help='minimal coverage after applying all other filters [1]', default='1')
    optional.add_option('-H', '--histogram', dest='histogram', metavar='PATH', help='destination PATH of evaluation plots', default='-')
    optional.add_option('-x', '--exclude_chrms', dest='exclude_chrm', metavar='STRINGLIST', help='list of comma separated chromosomes to exclude from evaluation', default='-')
    optional.add_option('-E', '--exclude_introns', dest='exclude_introns', metavar='STRINGLIST', help='list of comma separated intron files to exclude from submitted features', default='-')
    optional.add_option('-F', '--full_analysis', dest='full_analysis', action='store_true', help='generates overlap and coverage plots', default=False)
    optional.add_option('-o', '--outfile_base', dest='outfile_base', metavar='PATH', help='basedir for outfiles written', default='-')
    optional.add_option('-d', '--details', dest='details', action='store_true', help='destination for generated alignment intron list', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 3 or '-' in [options.features, options.anno_int] :
        parser.print_help()
        sys.exit(2)

    return options

def build_intron_list(options):
    """Builds up an intron list from the given alignment features."""

    counter = 0
    filter_counter = 0
    intron_lists = dict()

    for line in open(options.features, 'r'):
        
        if counter % 10000 == 0 and options.verbose:
            print 'lines read: [ %s (taken: %s / filtered: %s)]' % (counter, counter - filter_counter, filter_counter)
        counter += 1

        sl = line.strip().split('\t')
        chrm = sl[0]

        ### find valid intron
        if int(sl[3]) < options.min_coverage:
            filter_counter += 1
            continue 

        if options.splice_consensus:
            ### load genome if necessary
            if not options.genome_dict.has_key(chrm):
                if chrm in options.gio.contig_names:
                    if options.verbose:
                        print 'Reading chromosome %s' % chrm
                    fgen = open('%s/genome/%s.flat' % (options.gio.basedir, chrm))
                    options.genome_dict[chrm] = fgen.readline()
                    fgen.close()
                elif options.ignore_missing_chrm:
                    continue
                else:       
                    print >> sys.stderr, 'Chromosome Names do not match genome information file. Chromosome %s can not be found!' % chrm
                    exit(2)



        for _sl in sl[5:]:
            (val, count) = _sl.split(':')
            ex = int(val) / (options.max_feat_mismatches + 1)
            mm = int(val) % (options.max_feat_mismatches + 1)
            if ex >= options.min_exon_len and mm <= options.max_mismatches:
                ### intron = (strand [-1, 1], start, end)  ... end in the pythonic way
                intron = (0, int(sl[1]), int(sl[2]))
                try:
                    intron_lists[chrm][intron] += int(count)
                except KeyError:
                    try:
                        intron_lists[chrm][intron] = int(count)
                    except KeyError:
                        intron_lists[chrm] = {intron:int(count)}
            else:
                filter_counter += 1
            

    print '\n parsed file %s' % options.features
    print 'read %s lines' % counter
    if options.max_mismatches < 1000000 or options.min_exon_len > 0 or options.min_coverage > 1:
        print 'filter criteria:'
        print '    max mismatches: %s' % options.max_mismatches
        print '    min exon length: %s' % options.min_exon_len
        print '    min coverage: %s\n' % options.min_coverage
        print 'filtered %s lines' % filter_counter

    return intron_lists


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

def main():
    """Main function for comparing two intron lists """    
    
    options = parse_options(sys.argv)

    plotlist = plot.Plotlist()


    if options.outfile_base != '-':
        outfile_base = (options.outfile_base + '/' + options.features.split('/')[-1])
    else:
        outfile_base = options.features

    ### set filter tags for filenames
    if options.min_exon_len > 0:
        outfile_base += '_me%s' % options.min_exon_len
    if options.max_mismatches < 1000000:
        outfile_base += '_mm%s' % options.max_mismatches
    if options.max_intron_len < 100000000:
        outfile_base += '_mi%s' % options.max_intron_len
    if options.min_coverage > 0:
        outfile_base += '_mc%s' % options.min_coverage

    ### load genome, if necessary
    if options.splice_consensus:
        import genome_utils
        options.gio = genome_utils.GenomeInfo(options.genome)
        options.gio.contig_names.sort()
        options.genome_dict = dict()

    ### parse feature file for creating intron_list and coverage_map
    alignment_list = build_intron_list(options)
    
    ### get list of annotated introns
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

    ### filter predicted introns for excluded introns
    if options.exclude_introns != '-':
        _ex_introns = options.exclude_introns.strip().split(',')
        ### handle leading or trailing commas
        if _ex_introns[0] == '':
            _ex_introns = _ex_introns[1:]
        if _ex_introns[-1] == '':
            _ex_introns = _ex_introns[:-1]
        for _infile in _ex_introns:
            _ex_intron = cPickle.load(open(_infile, 'r'))
            for chrm in _ex_intron.keys():
                if alignment_list.has_key(chrm):
                    for _intron in _ex_intron[chrm].keys():
                        try:
                            del alignment_list[chrm][_intron]
                        except KeyError:
                            continue

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
    if options.min_coverage_2 > 1:
        print '\nFiltering intron list for min support after filtering'
        print '-----------------------------------------------------'
        skipped = 0
        for chrm in alignment_list.keys():
            skiplist = set()
            for intron in alignment_list[chrm].keys():
                if alignment_list[chrm][intron] < options.min_coverage_2:
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
            _curr_recall = float(matches) / float(max(1, len(annotation_list[chrm].keys())))
            _curr_precision = float(matches) / float(max(1, len(alignment_list[chrm].keys())))
            _curr_fscore = (2 * _curr_precision * _curr_recall) / float(max(1, _curr_precision + _curr_recall))
            print '-----------------------------'
            print ' in Chromosome %s ' % chrm 
            print '-----------------------------'
            print 'recall: %s' % _curr_recall
            print 'precision: %s' % _curr_precision
            print 'F-Score: %s' % _curr_fscore
            total_precision += _curr_precision
            total_recall += _curr_recall
            non_matched[chrm] = set(alignment_list[chrm].keys()).difference(set(annotation_list[chrm].keys()))
            ### do not include chromosomes with zero values into average
            if matches > 0:
                key_count += 1
   
    total_precision /= max(1.0, float(key_count))
    total_recall /= max(1.0, float(key_count))
    total_fscore = (2 * total_precision * total_recall) / float(max(1, total_precision + total_recall))
    print '-----------------------------'
    print ' average over all chromosomes '
    print '-----------------------------'
    print 'recall: %s' % total_recall
    print 'precision: %s' % total_precision
    print 'F-Score: %s' % total_fscore

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
        _header_line += 'average'
        _recall_line += str(total_recall)
        _precision_line += str(total_precision)
        print >> outf, _header_line
        print >> outf, _recall_line
        print >> outf, _precision_line
        outf.close() 
    
    if options.full_analysis:
        ### evaluate overlap only for those reads not exactly matching the annotation introns
        evaluate_intron_overlap(options, annotation_list, non_matched, plotlist)

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
