"""This script compares an alignment to the annotation. """

import cPickle
import sys
import os
import re
import scipy as sp
from parser import *

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--alignment', dest='alignment', metavar='FILE', help='alignment file in sam/bam format or alignment intron summary in pickle format', default='-')
    required.add_option('-A', '--annotation', dest='anno', metavar='FILE', help='annotation file in gtf or gff format', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-F', '--full_analysis', dest='full_analysis', action='store_true', help='generates overlap and coverage plots', default=False)
    optional.add_option('-s', '--strand_specific', dest='strand_specific', action='store_true', help='data is strand specific', default=False)
    optional.add_option('-R', '--ignore_reads', dest='exclude_set', metavar='FILE', help='file containing the reads to ignore', default='-')
    optional.add_option('-p', '--performance_log', dest='performance_log', action='store_true', help='store the intron recovery performance in extra log [off]', default=False)
    optional.add_option('-M', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [100000000]', default='100000000')
    optional.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    optional.add_option('-X', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismathes [10000]', default=10000)
    optional.add_option('-b', '--bam_input', dest='bam_input', action='store_true', help='input has BAM format - does not work for STDIN', default=False)
    optional.add_option('-t', '--samtools', dest='samtools', metavar='PATH', help='if SAMtools is not in your PATH, provide the right path here (only neccessary for BAM input)', default='samtools')
    optional.add_option('-m', '--min_read_coverage', dest='min_coverage', metavar='INT', type='int', help='minimal coverage to count an intron [1]', default='1')
    optional.add_option('-O', '--max_overlap', dest='max_overlap', metavar='INT', type='int', help='maximal evalutated exon overlap [5]', default='5')
    optional.add_option('-H', '--histogram', dest='histogram', metavar='PATH', help='destination PATH of evaluation plots', default='-')
    optional.add_option('-x', '--exclude_chrms', dest='exclude_chrm', metavar='STRINGLIST', help='list of comma separated chromosomes to exclude from evaluation', default='-')
    optional.add_option('-S', '--show_sources', dest='show_sources', action='store_true', help='only show available sources of gff file', default=False)
    optional.add_option('-L', '--sources', dest='sources', help='list of comma-separated soources to use from annotation', default='-')
    optional.add_option('-o', '--outfile_base', dest='outfile_base', metavar='PATH', help='basedir for outfiles written', default='-')
    optional.add_option('-l', '--lines', dest='lines', metavar='INT', type='int', help='maximal number of alignment lines to read [-]', default=None)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 3 or '-' in [options.anno, options.alignment] :
        parser.print_help()
        sys.exit(2)

    return options

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

    if options.outfile_base != '-':
        outfile_base = (options.outfile_base + '/' + options.alignment.split('/')[-1])
    else:
        outfile_base = options.alignment

    ### set filter tags for filenames
    if options.min_exon_len > 0:
        outfile_base += '_me%s' % options.min_exon_len
    if options.max_mismatches < 10000:
        outfile_base += '_mm%s' % options.max_mismatches
    if options.max_intron_len < 100000000:
        outfile_base += '_mi%s' % options.max_intron_len
    if options.min_coverage > 1:
        outfile_base += '_mc%s' % options.min_coverage

    ### load or compute annotated intron list
    if options.anno.endswith('pickle'):
        (annotation_list, exon_map) = cPickle.load(open(options.anno, 'r'))
    elif os.path.exists(options.anno + '.pickle'):
        (annotation_list, exon_map) = cPickle.load(open(options.anno + '.pickle', 'r'))
    else:
        (annotation_list, exon_map) = intron_list_from_annotation(options)
        cPickle.dump((annotation_list, exon_map), open(options.anno + '.pickle', 'wb'), -1)

    ### load or compute alignment intron list
    if options.alignment.endswith('pickle'):
        (alignment_list, align_cov_vec) = cPickle.load(open(options.alignment, 'r'))
    elif os.path.exists(options.alignment + '.pickle'):
        (alignment_list, align_cov_vec) = cPickle.load(open(options.alignment + '.pickle', 'r'))
    else:
        ### check for blacklist to handle
        if options.exclude_set != '-':
            exclude_set = set()
            print '\nParsing exclude list from file %s' % options.exclude_set
            print '-----------------------------------------'
            for line in open(options.exclude_set, 'r'):
                _l = line.strip().split('\t')
                exclude_set.add((_l[0], int(_l[1])))
            (align_cov_vec, alignment_list, alignment_support) = intron_list_from_alignment(options, exclude_set)
        else:
            ### parse alignmentfile for creating intron_list and coverage_map
            (align_cov_vec, alignment_list, alignment_support) = intron_list_from_alignment(options)
            cPickle.dump((align_cov_vec, alignment_list, alignment_support), open(options.alignment + '.pickle', 'wb'), -1)
        
    ### filter annotated and predicted introns for excluded chromosomes
    if options.exclude_chrm != '-':
        k_idx = sp.where(~sp.in1d(aligment_list[:, 0], options.exclude_chrm.strip().strip(',').split(',')))[0]
        alignment_list = alignment_list[k_idx, :]

    ### filter intron lists for max intron length
    print '\nFiltering intron lists for max intron len'
    print '-----------------------------------------'
    tmp = annotation_list[:, 2:].astype('int')
    k_idx = sp.where(tmp[:, 1] - tmp[:, 0] <= options.max_intron_len)[0]
    print '%s introns removed from annotation\n' % (tmp.shape[0] - k_idx.shape[0])
    annotation_list = annotation_list[k_idx, :]

    tmp = alignment_list[:, 2:].astype('int')
    k_idx = sp.where(tmp[:, 1] - tmp[:, 0] <= options.max_intron_len)[0]
    print '%s introns removed from alignment\n' % (tmp.shape[0] - k_idx.shape[0])
    alignment_list = alignment_list[k_idx, :]
    alignment_support = alignment_support[k_idx, :]

    ### filter intron lists for min coverage
    if options.min_coverage > 1:
        print '\nFiltering intron list for min support '
        print '-----------------------------------------'
        k_idx = sp.where(alignment_support < options.min_coverage)[0]
        print '%s introns removed from alignment\n' % (alignment_list.shape[0] - k_idx.shape[0])
        alignment_list = alignment_list[k_idx, :]
        alignment_support = alignment_support[k_idx, :]

    ### remove strand for non-strand-specific comparison
    if not options.strand_specific:
        alignment_list = alignment_list[:, [0,2,3]]
        annotation_list = annotation_list[:, [0, 2, 3]]

    ### match intron lists
    precision = []
    recall = []
    pdb.set_trace()
    chrms = sp.unique(annotation_list[:, 0])
    for chrm in chrms:
        c_idx = sp.where(alignment_list[:, 0] == chrm)[0]
        if c_idx.shape[0] > 0:
            a_idx = sp.where(annotation_list[:, 0] == chrm)[0]
            matches = sp.sum(sp.in1d(row_strings(alignment_list[c_idx, :]), row_strings(annotation_list[a_idx, :])))
            precision.append(float(matches) / float(max(1, c_idx.shape[0])))
            recall.append(float(matches) / float(max(1, a_idx.shape[0])))
            print '-----------------------------'
            print ' in Chromosome %s ' % chrm 
            print '-----------------------------'
            print 'recall: %.6f' % recall[-1]
            print 'precision: %.6f' % precision[-1]
            #non_matched[chrm] = set(alignment_list[chrm].keys()).difference(set(annotation_list[chrm].keys()))
   
    ### compute total precision
    matches = sp.sum(sp.in1d(row_strings(alignment_list), row_strings(annotation_list)))
    total_recall = float(matches) / max(1, float(sp.sum(sp.in1d(annotation_list[:, 0], sp.unique(alignment_list[:, 0])))))
    total_precision = float(matches) / max(1, float(alignment_list.shape[0]))
    print '-----------------------------'
    print ' average over all chromosomes '
    print '-----------------------------'
    print 'recall: %.6f' % total_recall
    print 'precision: %.6f' % total_precision

    if options.performance_log:
        outf = open(outfile_base + '_performance.log', 'w')
        print >> outf, '\t'.joun(chrms) + '\taverage'
        print >> outf, '\t'.join(['%.6f' % x for x in recall]) + '\t%.6f' % total_recall 
        print >> outf, '\t'.join(['%.6f' % x for x in precision]) + '\t%.6f' % total_precision
        outf.close() 
    
    ### get coverage gff
    if False: #options.cov_gff:
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
