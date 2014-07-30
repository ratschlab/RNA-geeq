"""This tool parses several alignments in SAM format and scans them for multireads. The multireads are stored in the given outfile. The input files MUST be SORTED by read_id! """
import sys
import re
import pdb

from utils import *

def parse_options(argv):
    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-o', '--output', dest='output', metavar='FILE', help='outfile for multireads or list of files [combined,align1,align2,...]', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-M', '--max_mismatch', dest='max_mismatches', metavar='INT', type='int', help='maximal number of mismatches in the score [30]', default=30)
    optional.add_option('-N', '--max_reads', dest='max_reads', metavar='INT', type='int', help='maximal number of readse to be evaluated [1000000000]', default=1000000000)
    optional.add_option('-S', '--score', dest='score', metavar='FILE', help='outfile for the scoring line', default='-')
    optional.add_option('-L', '--labels', dest='labels', metavar='FILE', help='label list', default='-')
    optional.add_option('-x', '--skip_mt', dest='skip_mt', metavar='BOOL', action='store_true', help='skip mtdna for analysis', default=False)
    optional.add_option('-s', '--samtools', dest='samtools', metavar='PATH', help='if SAMtools is not in your PATH, provide the right path here (only neccessary for BAM input)', default='samtools')
    optional.add_option('-E', '--eval_strategy', dest='eval_strategy', metavar='INT',  type='int', help='evaluation strategy (0-penalty score, 1-stratum by mm, 2-stratum)', default=0)
    optional.add_option('-i', '--ignore_clipped', dest='ignore_clipped', metavar='BOOL', action='store_true', help='ignore clipped alignments', default=False)
    optional.add_option('-e', '--eval_only', dest='eval_only', metavar='BOOL', action='store_true', help='suppresses writing multimappers to files in -o', default=False)
    optional.add_option('-O', '--outfile_only', dest='outfile_only', metavar='BOOL', action='store_true', help='suppresses evaluation only writing multimappers to files in -o', default=False)
    optional.add_option('-b', '--bam_input', dest='bam_input', metavar='BOOL', action='store_true', help='input files have BAM format', default=False)
    optional.add_option('-d', '--debug', dest='debug', metavar='BOOL', action='store_true', help='print debug info to screen', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', metavar='BOOL', action='store_true', help='be verbose', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) <  2 :
        parser.print_help()
        sys.exit(2)

    return (options, args)


def process_line(sl, options):
    """Processes evaluation of current read line (sl)."""
    
	### split CIGAR
    (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
    try:
        size = [int(i) for i in size]
    except ValueError:
        print >> sys.stderr, 'Integer conversion error in CIGAR of %s' % str(sl)
        sys.exit(-1)
	
	### initialization
    length = 0
    start_pos = int(sl[3])
    exon_range = set()
    chrm = sl[2]
    chrm = chrm.replace('chr', '')
    if chrm == 'M':
        chrm = 'MT'

	### extract exonic positions from CIGAR
    exon_sum = 0
    for pp in range(len(op)):
        if op[pp] == 'H' and pp == 0:
            start_pos -= size[pp]
        elif op[pp] in ['N', 'H', 'S', 'M', 'D']:
            length += size[pp]

        if not op[pp] in ['N', 'I']:
            seg_start = start_pos + exon_sum
            exon_range = exon_range.union(set(range(seg_start, seg_start + size[pp])))
        if op[pp] != 'I':
            exon_sum += size[pp]
    end_pos = (start_pos + length)

    ### get mismatches
    mismatches = get_mm(sl)
    if mismatches == -1:
        print >> sys.stderr, "One of given files lacks mismatch information! - Bailing out.\n%s" % (str(sl))
        sys.exit(-1)

	### compute score for later stratum definition
    if options.eval_strategy == 0:
        score = float(1) - (float(min(mismatches, options.max_mismatches)) / float(options.max_mismatches))
        assert(score <= 1)
    elif options.eval_strategy in [1, 2, 3]:
        score = min(mismatches, options.max_mismatches)

    if options.debug:
        print 'process line: \n %s\n' % str(sl)

    return ((chrm, int(sl[1]) & 128, start_pos, end_pos), (score, sl[0], sum(list(exon_range))), exon_range)

def cmp4(a, b): 
    """Compares fourth element of a and b """
    return cmp(b[4], a[4])

def process_multimappers_by_mm(scored_submissions, score, options, alignments, alignments2, exons, max_score):
    """Processes the evaluation of a multi mapping alignment set in strata defined by mismatches"""
    
    ### choose strategy (0 - ; 1 - )
    strategy = 1

    ### bundle submission strata by mm
    mm_list_all = dict()
    mm_list_sub = dict()
    exon_list_sub = dict()
    exon_list_all = dict()

    ### merge alignments into strata
    ### mm_list_all -> merged list of all submissions
    ### mm_list_sub -> all strata of a single submission
    ### exon lists, respectively
    for subm in scored_submissions.keys():
        mm_list_sub[subm] = dict()
        exon_list_sub[subm] = dict()
        if len(exons[subm]) != len(scored_submissions[subm]):
            print >> sys.stderr, 'Number of alignments and exonic position lists is different!'
            print >> sys.stderr, str(exons)
            print >> sys.stderr, str(scored_submissions)
            sys.exit(-1)

        idx = 0
        for s in scored_submissions[subm]:
            try:
                mm_list_all[s[4]].add(s)
            except KeyError:
                mm_list_all[s[4]] = set([s])
            try:
                mm_list_sub[subm][s[4]].add(s)
            except KeyError:
                mm_list_sub[subm][s[4]] = set([s])
            exon_list_sub[subm][s] = exons[subm][idx]
            exon_list_all[s] = exons[subm][idx]
            idx += 1

    ### determine max possible score
	### that is, a submission can explain at most all strata in the list to 100 percent
    for mm in mm_list_all:
        max_score[mm] += 1

    ### strategy 0 merges the alignments by mismatches into strata
    ### it augments each submission with alignments overlapping the own alignments to at least 90 percent
    if strategy == 0:
        ### which stratum counts for sensitivity?
        ###     a stratum counts for sensitivity, if itself or a worse 
        ###     (more mm) stratum is  present in the submission
        for subm in scored_submissions.keys():
            for mm in sorted(mm_list_all.keys()):
                if mm_list_sub[subm].has_key(mm):
                    alignments[subm][mm] += 1
                    alignments2[subm][mm] += len(mm_list_sub[subm][mm])
                elif len(mm_list_sub[subm].keys()) > 0 and mm < max(mm_list_sub[subm].keys()):
                    alignments[subm][mm] += 1
                else:
                    break
        
        ### check for overlapping alignments, and count them as positives
        ### augment alignments of each submission with alignments overlapping
        ### the merged stratum list with at least 90 percent exonic positions
        for subm in scored_submissions.keys():
            for mm in mm_list_sub[subm].keys():
                add_set = set()
                for align1 in (mm_list_all[mm] - mm_list_sub[subm][mm]):
                    for align2 in mm_list_sub[subm][mm]:
                        if float(len(exon_list_all[align1].intersection(exon_list_sub[subm][align2]))) / float(len(exon_list_all[align1])) > 0.9:
                            add_set.add(align1)
                mm_list_sub[subm][mm] = mm_list_sub[subm][mm].union(add_set)

        ### compute score for each stratum as fraction of 'explained' alignments
        ### in the current stratum of the mergd list
        for subm in scored_submissions.keys():
            for mm in mm_list_sub[subm].keys():
                score[subm][mm] += (float(len(mm_list_sub[subm][mm].intersection(mm_list_all[mm])))/float(len(mm_list_all[mm])))

    ### strategy 1 counts how many alignments of a stratum of the merged list
    ### overlap to the alignments of a submission by more than 90 percent
    elif strategy == 1:
        for subm in scored_submissions.keys():
            if len(scored_submissions[subm]) == 0:
                continue
            max_mm = max([v[4] for v in scored_submissions[subm]])

            for mm in sorted(mm_list_all.keys()):
                if mm >  max_mm:
                    break
                score_count = 0
                for align1 in (mm_list_all[mm]):
                    for align2 in (scored_submissions[subm]):
                        if align1[0] == align2[0] and float(len(exon_list_all[align1].intersection(exon_list_sub[subm][align2]))) / float(max(len(exon_list_all[align1]), 1)) > 0.9:
                            score_count  += 1
                            break
                score[subm][mm] += (float(score_count)/float(len(mm_list_all[mm])))
                alignments[subm][mm] += 1
                alignments2[subm][mm] += 1 

    if options.debug:
        print '_________________'
        print score
        print alignments
        print '+++++++++++++++++'

    return (score, alignments, alignments2, max_score)
        
def process_multimappers_by_aligns(scored_submissions, score, options, alignments, exons, max_score):
    """Processes the evaluation of a multi mapping alignment set separated by aligns"""
    
    ### form strata by score
    strata_list_all = dict()
    strata_list_sub = dict()
    exon_list_sub = dict()
    exon_list_all = dict()
    
    ### merge alignments into strata
    ### strata_list_all -> merged list of all submissions
    ### strata_list_sub -> all strata of a submission
    ### exon lists, respectively
    for subm in scored_submissions.keys():
        strata_list_sub[subm] = dict()
        exon_list_sub[subm] = dict()
        idx = 0
        for s in scored_submissions[subm]:
            try:
                strata_list_all[s[4]].add(s)
            except KeyError:
                strata_list_all[s[4]] = set([s])
            try:
                strata_list_sub[subm][s[4]].add(s)
            except KeyError:
                strata_list_sub[subm][s[4]] = set([s])
            exon_list_sub[subm][s] = exons[subm][idx]
            exon_list_all[s] = exons[subm][idx]
            idx += 1

    ### compute scores
    for subm in scored_submissions.keys():
		### initializations
        length = 0
        skipped = 0
        pos = 0
        merge_list = []

        for stratum in sorted(strata_list_all.keys()):
            if strata_list_sub[subm].has_key(stratum):
                length += len(strata_list_sub[subm][stratum])

                _strata_all = strata_list_all[stratum]
                if len(merge_list) > 0:
                    for i in merge_list:
                        _strata_all = _strata_all.union(strata_list_all[i])
                    merge_list = []

                ### check for overlapping alignments and count them as true positives
                add_set = set()
                for align1 in (_strata_all - strata_list_sub[subm][stratum]):
                    for align2 in (strata_list_sub[subm][stratum]):
                        if align1[0] == align2[0] and float(len(exon_list_all[align1].intersection(exon_list_sub[subm][align2]))) / float(len(exon_list_all[align1])) > 0.9:
                            add_set.add(align1)
                _strata_sub = strata_list_sub[subm][stratum].union(add_set)

                try:
                    score[subm][length] += (float(len(_strata_sub.intersection(_strata_all)))/float(len(_strata_all)))
                except KeyError:
                    score[subm][length] = (float(len(_strata_sub.intersection(_strata_all)))/float(len(_strata_all)))
                try:
                    alignments[subm][length] += 1
                except KeyError:
                    alignments[subm][length] = 1
            else:
                if len(strata_list_sub[subm]) > (pos - skipped):
                    ### merge current stratum with next
                    merge_list.append(stratum)
                skipped += 1
            pos += 1

    if options.debug:
        print '_________________'
        print score
        print alignments
        print '+++++++++++++++++'

    return (score, alignments, max_score)

def process_multimappers_stratum_score(scored_submissions, score, max_score, alignments2, options):
    """This function scores each submission with an additive score of how complete the list of all merged strata has been predicted."""
    
    ### form strata by score
    strata_list_all = dict()
    strata_list_sub = dict()

    ### merge alignments into strata
    ### strata_list_all -> merged list of all submissions
    ### strata_list_sub -> all strata of a single submission
    ### exon lists, respectively
    for subm in scored_submissions.keys():
        strata_list_sub[subm] = dict()
        for s in scored_submissions[subm]:
            try:
                strata_list_all[s[4]].add(s)
            except KeyError:
                strata_list_all[s[4]] = set([s])
            try:
                strata_list_sub[subm][s[4]].add(s)
            except KeyError:
                strata_list_sub[subm][s[4]] = set([s])

    ### determine max possible score
	### that is, a submission can explain at most all strata in the list to 100 percent
    ### and thus add 1/stratum-nr to the score
    
    for stratum in strata_list_all:
    #for stratum in range(len(strata_list_all)):
        max_score[stratum] += (1.0/float(stratum + 1))

    ### the score for each stratum is weighed by 1/stratum-nr and each alignment contributes 1/stratum-size
    ### only exact stratum matches are counted
    for subm in scored_submissions.keys():
        for stratum in strata_list_sub[subm].keys():
			### stratum score as explained fraction weighed with the stratum number
            score[subm][stratum] += (float(len(strata_list_sub[subm][stratum])) / float(max(1, len(strata_list_all[stratum]))) / float(stratum + 1))
            alignments2[subm][stratum] += 1
            #score[subm][stratum] += (float(len(strata_list_sub[subm][stratum])) / float(max(1, len(strata_list_all[stratum]))) / float(sorted(strata_list_all.keys()).index(stratum) + 1))
            #alignments2[subm][sorted(strata_list_all.keys()).index(stratum)] += 1 
    return (score, alignments2, max_score)

def process_multimappers(scored_submissions, score, options):
    """Processes the evaluation of a multi mapping alignment set """

    ranking = set()

    for s in scored_submissions.keys():
        ### unite all submission's multimapping alignments
        ranking = ranking.union(set(scored_submissions[s]))
        if options.debug:
            print scored_submissions[s]

    ### rank list by score (highest score first)
    ranking = list(ranking)
    ranking.sort(cmp4)

    if options.debug:
        print '\n%s\n\n' % str(ranking)

    for s in scored_submissions.keys():
        if len(scored_submissions[s]) == 0:
            continue

        ### sort by score (highest first)
        scored_submissions[s].sort(cmp4)
        ### choose end position as highest index in ranking with same score as worst alignment in submission
        ranking_end = 0
        for i in range(ranking.index(scored_submissions[s][-1]), len(ranking)):
            if ranking[i][4] < scored_submissions[s][-1][4]:
                break
            ranking_end = i

        for i in range(ranking_end + 1):
            if ranking[i] in scored_submissions[s]:
                score[s] += ranking[i][4]
            else:
                ### do not penalize undetected alignments, if submission has an overlapping alignment
                overlap = False
                score_start = ranking[i][2]
                score_end = ranking[i][3]
                for j in range(len(scored_submissions[s])):
                    if scored_submissions[s][:2] != ranking[i][:2]:
                        continue
                    sub_start = scored_submissions[s][j][2]
                    sub_end = scored_submissions[s][j][3]
                    if (sub_start >= score_start and sub_start < score_end) or \
                        (sub_end > score_start and sub_end <= score_end):
                        overlap =  True
                        break
                if not overlap:
                    score[s] -= ranking[i][4]

    if options.debug:
        print 'score: %s\n' % score
    return score

def main():
    """Main function ... """
   
    (options, args) = parse_options(sys.argv)

    import locale
    locale.setlocale(locale.LC_ALL, 'C')

    ### check given labels    
    if options.labels != '-':
        label_list = options.labels.strip().strip(',').split(',')
        if len(label_list) > 1 and len(label_list) != len(args):
            print >> sys.stderr, "ERROR: Labels for all single files must be given"
            sys.exit(-1)

    ### initializations
    closed_files = 0
    read_counter = 0
    multireads = 0
    min_read_id = '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    file_list = dict()
    stored_lines = dict()
    scores = dict()
    if options.eval_strategy in [1, 2, 3]:
        alignments = dict()
        alignments2 = dict()
        max_score = dict()

    ### create input and output streams
    ### initialize lists
    if not options.eval_only:
        ### parse output file list
        output_list = options.output.strip().strip(',').split(',')
        if len(output_list) > 1 and len(output_list) != (len(args) + 1):
            print >> sys.stderr, "ERROR: Output files for all single files must be given"
            sys.exit(-1)
        out_files = {-1 : open(output_list[0], 'w')}

    ### check input files
    skipped_labels = []
    skipped_aligns = []
    if options.bam_input:
        import subprocess
        process_handles = dict()
    for idx in range(len(args)):
        if options.bam_input:
            process_handle = subprocess.Popen([options.samtools, 'view', '-h', args[idx]], stdout=subprocess.PIPE)
            _file = process_handle.stdout
        else:
            _file = open(args[idx], 'r')
        _line = _file.readline()
        ### seek for first valid SAM line - skip header lines
        while len(_line) > 0 and _line[0] == '@':
            _line = _file.readline()
        if len(_line) == 0:
            print >> sys.stderr, 'file %s has no valid lines' % args[idx]
            _file.close()
            sys.exit(-1)
        ### check number of available fields
        sl = _line.strip().split()
        if len(sl) < 9:
            print >> sys.stderr, 'file %s has too few fields' % args[idx]
            print >> sys.stderr, _line
            _file.close()
            sys.exit(2)
        ### check files for mismatch info
        mm = get_mm(sl)
        if mm == -1:
            print >> sys.stderr, '\nskip alignment %s - lacks mismatch information\n' % args[idx]
            skipped_aligns.append(args[idx])
            skipped_labels.append(label_list[idx])
            _file.close()
            continue

        ### check for minimal read ID
        if sl[0] < min_read_id:
            min_read_id = sl[0]

        ### close file handles
        _file.close()
        if options.bam_input:
            process_handle.kill()
            del process_handle

        ### init data structures depending on eval strategy
        stored_lines[idx] = []
        if options.eval_strategy == 0:
            scores[idx] = 0
        elif options.eval_strategy in [1, 3]:
            scores[idx] = dict()
            alignments[idx] = dict()
            alignments2[idx] = dict()
            max_score[idx] = dict()
            for mm in range(options.max_mismatches + 1):
                scores[idx][mm] = 0.0
                alignments[idx][mm] = 0
                alignments2[idx][mm] = 0
                max_score[mm] = 0
        elif options.eval_strategy == 2:
            scores[idx] = dict()
            alignments[idx] = dict()
            alignments2[idx] = dict()
            max_score[idx] = dict()
            for mm in range(1, options.max_mismatches + 1):
                scores[idx][mm] = 0.0
                alignments[idx][mm] = 0
                alignments2[idx][mm] = 0
                max_score[mm] = 0
 
        ### open alignment files
        if options.bam_input:
            process_handles[idx] = subprocess.Popen([options.samtools, 'view', '-h', args[idx]], stdout=subprocess.PIPE)
            file_list[idx] = process_handles[idx].stdout
        else:
            file_list[idx] = open(args[idx], 'r')
        if not options.eval_only:
            if len(output_list) > 1:
                out_files[idx] = open(output_list[idx + 1], 'w')
            else:
                out_files[idx] = open(args[idx] + '.multireads', 'w')
   
    ### remove skipped alignments 
    for alignment in skipped_aligns:
        args.remove(alignment)
    if options.verbose:
        print 'arguments:\n%s' % str(args)
    ### remove labels of skipped alignments
    if options.labels != '-':
        for label in skipped_labels:
            label_list.remove(label)
        if options.verbose:
            print 'labels:\n%s' % str(label_list)

    assert(options.eval_only or (len(file_list) + 1) == len(out_files))

    ### parse alignments as long as there is at least one file stream open
    ### or the loop has been broken because of other reasons
    while closed_files < len(file_list):
        if read_counter > options.max_reads:
            break
        if options.verbose and (read_counter % 1000 in range(len(file_list))):
            print 'processed %s reads - found %s multimappers [ %s percent ]' % (read_counter, multireads, float(multireads) * float(100) / float(max(read_counter, 1)))  

        curr_aligns_1 = []
        curr_aligns_2 = []
        scored_aligns_1 = dict()
        scored_aligns_2 = dict()
        exons_1 = dict()
        exons_2 = dict()
        curr_exons_1 = dict()
        curr_exons_2 = dict()
        curr_id = ''

        ### iterate over submission files
        for subm in file_list.keys():
            if file_list[subm].closed:
                continue

            scored_aligns_1[subm] = []
            scored_aligns_2[subm] = []
            exons_1[subm] = []
            exons_2[subm] = []

            ### parse for valid lines - otherwise close file stream
            ### parse only new files, iff stored_lines have been processed
            ### read alignments from stored next lines
            if len(stored_lines[subm]) >= 1: 
                if stored_lines[subm][0] <= min_read_id:
                    # (_a, _b, _c) = ((chrm, strand, start, stop), (score, read_id, checksum), exon_range)
                    (_a, _b, _c) = process_line(stored_lines[subm], options)

                    ### skip analysis for mtdna
                    if not (options.skip_mt and (_a[0] in ['M', 'MT'])):
                        if (int(stored_lines[subm][1]) & 128 == 0):
                            curr_aligns_1.append(_a)
                            scored_aligns_1[subm].append(_a + _b)
                            exons_1[subm].append(_c)
                            try:
                                curr_exons_1[_a[0]].append(_c)
                            except KeyError:
                                curr_exons_1[_a[0]] = [_c]
                        else:
                            curr_aligns_2.append(_a)
                            scored_aligns_2[subm].append(_a + _b)
                            exons_2[subm].append(_c)
                            try:
                                curr_exons_2[_a[0]].append(_c)
                            except KeyError:
                                curr_exons_2[_a[0]] = [_c]

                        if stored_lines[subm][0].strip() != '':
                            curr_id = stored_lines[subm][0].strip()
                    stored_lines[subm] = []

            ### if no line stored, read new ones        
            if len(stored_lines[subm]) == 0: 
                line = file_list[subm].readline()

                while (len(line) > 0 and (line[0] == '@' or line[0] == '#')):
                    line = file_list[subm].readline()

                if line == '':
                    file_list[subm].close()
                    if options.bam_input:
                        process_handles[subm].kill()
                    closed_files += 1
                    continue

                sl = line.strip().split()

                if len(sl) < 9:
                    file_list[subm].close()
                    if options.bam_input:
                        process_handles[subm].kill()
                    closed_files += 1
                    continue

                ### read alignments from submission as long as read id stays the same
                while sl[0] <= min_read_id:
                    # (_a, _b, _c) = ((chrm, strand, start, stop), (score, read_id, checksum), exon_range)
                    (_a, _b, _c) = process_line(sl, options)
                    ### skip analysis for mtdna
                    if options.skip_mt and (_a[0] in ['M', 'MT']):
                        sl = file_list[subm].readline().strip().split()
                        if len(sl) < 9:
                            file_list[subm].close()
                            if options.bam_input: 
                                process_handles[subm].kill()
                            closed_files += 1
                            break
                        continue

                    ### ignore clipped alignments
                    if options.ignore_clipped and (sl[5].find('H') > -1 or sl[5].find('S') > -1):
                        sl = file_list[subm].readline().strip().split()
                        if len(sl) < 9:
                            file_list[subm].close()
                            if options.bam_input: 
                                process_handles[subm].kill()
                            closed_files += 1
                            break
                        continue

                    if (int(sl[1]) & 128 == 0):
                        curr_aligns_1.append(_a)
                        scored_aligns_1[subm].append(_a + _b)
                        exons_1[subm].append(_c)
                        try:
                            curr_exons_1[_a[0]].append(_c)
                        except KeyError:
                            curr_exons_1[_a[0]] = [_c]
                    else:
                        curr_aligns_2.append(_a)
                        scored_aligns_2[subm].append(_a +_b)
                        exons_2[subm].append(_c)
                        try:
                            curr_exons_2[_a[0]].append(_c)
                        except KeyError:
                            curr_exons_2[_a[0]] = [_c]

                    if sl[0].strip() != '':
                        curr_id = sl[0].strip()

                    sl = file_list[subm].readline().strip().split()

                    if len(sl) < 9:
                        file_list[subm].close()
                        if options.bam_input: 
                            process_handles[subm].kill()
                        closed_files += 1
                        break

                #if options.debug:
                #    print 'submission %s' % subm
                #    print scored_aligns_1
                #    print scored_aligns_2

                ### store next line only if further stored line has been processed
                if len(sl) > 1:
                    stored_lines[subm] = sl
                ### close, if no next line
                else:
                    file_list[subm].close()
                    if options.bam_input: 
                        process_handles[subm].kill()
                    closed_files += 1

                ### print multireads on single submission level
                if not options.eval_only:
                    ### if we found a multimapper (in a single submission) --> log multimapper
                    if len(scored_aligns_1[subm]) > 1:
                        broken = False
                        for k in range(len(scored_aligns_1[subm]) - 1):
                            for l in range(k + 1, len(scored_aligns_1[subm])):
                                # the mappings share no single exonic position or are located on different chromosomes 
                                if (scored_aligns_1[subm][k][0] != scored_aligns_1[subm][l][0]) or len(set(exons_1[subm][k]).intersection(exons_1[subm][l])) == 0:
                                    broken = True
                                    break
                            if broken:
                                break
                        if broken:
                            assert (curr_id != '')
                            if not options.eval_only:
                                print >> out_files[subm], '%s\t64' % curr_id

                    if len(scored_aligns_2[subm]) > 1:
                        broken = False
                        for k in range(len(scored_aligns_2[subm]) - 1):
                            for l in range(k + 1, len(scored_aligns_2[subm])):
                                # the mappings share no single exonic position or are located on different chromosomes 
                                if (scored_aligns_2[subm][k][0] != scored_aligns_2[subm][l][0]) or len(set(exons_2[subm][k]).intersection(exons_2[subm][l])) == 0:
                                    broken = True
                                    break
                            if broken:
                                break
                        if broken:
                            assert (curr_id != '')
                            if not options.eval_only:
                                print >> out_files[subm], '%s\t128' % curr_id
           
        min_read_id_list = [stored_lines[v][0] for v in stored_lines.keys() if len(stored_lines[v]) > 0]
        if min_read_id_list != None and len(min_read_id_list) > 0:
            min_read_id = min(min_read_id_list)

        ### if we found a multimapper (over the union of all submissions) --> log multimapper
        if len(set(curr_aligns_1)) > 1:
            ## check for mutually overlapping alignments
            ## curr_aligns = (chr, pair, start, end, score)
            broken = False
            #for k in range(len(curr_aligns_1) - 1):
            #    for l in range(k + 1, len(curr_aligns_1)):
            #        if not ((curr_aligns_1[k][2] >= curr_aligns_1[l][2] and curr_aligns_1[k][2] < curr_aligns_1[l][3]) or \
            #                (curr_aligns_1[k][3] > curr_aligns_1[l][2] and curr_aligns_1[k][3] <= curr_aligns_1[l][3]) or \
            #                (curr_aligns_1[k][3] >= curr_aligns_1[l][3] and curr_aligns_1[k][2] <= curr_aligns_1[l][2])):
            if len(curr_exons_1.keys()) > 1:
                broken = True
            else:
                _chrm = curr_exons_1.keys()[0]
                for k in range(len(curr_exons_1[_chrm]) - 1):
                    for l in range(k + 1, len(curr_exons_1[_chrm])):
                        if len(set(curr_exons_1[_chrm][k]).intersection(curr_exons_1[_chrm][l])) == 0:
                            broken = True
                            break
                    if broken:
                        break
            if broken:
                assert (curr_id != '')
                if not options.eval_only:
                    print >> out_files[-1], '%s\t64' % curr_id
                multireads += 1
        
                if not options.outfile_only:
                    ### run further analysis
                    if options.eval_strategy == 0:
                        scores = process_multimappers(scored_aligns_1, scores, options)
                    elif options.eval_strategy == 1:
                        (scores, alignments, alignments2, max_score) = process_multimappers_by_mm(scored_aligns_1, scores, options, alignments, alignments2, exons_1, max_score)
                    elif options.eval_strategy == 2:
                        (scores, alignments, max_score) = process_multimappers_by_aligns(scored_aligns_1, scores, options, alignments, exons_1, max_score)
                    elif options.eval_strategy == 3:
                        (scores, alignments2, max_score) = process_multimappers_stratum_score(scored_aligns_1, scores, max_score, alignments2, options)

        if len(set(curr_aligns_2)) > 1:
            ## check for mutually overlapping alignments
            ## curr_aligns = (chr, pair, start, end, score)
            broken = False
            #for k in range(len(curr_aligns_2) - 1):
            #    for l in range(k + 1, len(curr_aligns_2)):
            #        if not ((curr_aligns_2[k][2] >= curr_aligns_2[l][2] and curr_aligns_2[k][2] < curr_aligns_2[l][3]) or \
            #                (curr_aligns_2[k][3] > curr_aligns_2[l][2] and curr_aligns_2[k][3] <= curr_aligns_2[l][3]) or \
            #                (curr_aligns_2[k][3] >= curr_aligns_2[l][3] and curr_aligns_2[k][2] <= curr_aligns_2[l][2])):
            if len(curr_exons_2.keys()) > 1:
                broken = True
            else:
                _chrm = curr_exons_2.keys()[0]
                for k in range(len(curr_exons_2[_chrm]) - 1):
                    for l in range(k + 1, len(curr_exons_2[_chrm])):
                        if len(set(curr_exons_2[_chrm][k]).intersection(curr_exons_2[_chrm][l])) > 0:
                            broken = True
                            break
                    if broken:
                        break
            if broken:
                assert (curr_id != '')
                if not options.eval_only:
                    print >> out_files[-1], '%s\t128' % curr_id
                multireads += 1

                if not options.outfile_only:
                    ### run further analysis
                    if options.eval_strategy == 0:
                        scores = process_multimappers(scored_aligns_2, scores, options)
                    elif options.eval_strategy == 1:
                        (scores, alignments, alignments2, max_score) = process_multimappers_by_mm(scored_aligns_2, scores, options, alignments, alignments2, exons_2, max_score)
                    elif options.eval_strategy == 2:
                        (scores, alignments, max_score) = process_multimappers_by_aligns(scored_aligns_2, scores, options, alignments, exons_2, max_score)
                    elif options.eval_strategy == 3:
                        (scores, alignments2, max_score) = process_multimappers_stratum_score(scored_aligns_1, scores, max_score, alignments2, options)

        if len(set(curr_aligns_1)) > 0:
            read_counter += 1
        if len(set(curr_aligns_2)) > 0:
            read_counter += 1


    ### close everything which is not closed yet
    for subm in file_list:
        if not file_list[subm].closed:
            file_list[subm].close()
            if options.bam_input: 
                process_handles[subm].kill()
    if not options.eval_only:
        for subm in out_files:
            if not out_files[subm].closed:
                out_files[subm].close()

    ### print score to outfile
    if options.score != '-':
        score_file = options.score.replace('.dat', '')
        score_file += ('.strat_' + str(options.eval_strategy) + '.dat')
        score_out = open(score_file, 'w')

        if options.eval_strategy == 0:
            _line = ''
            for sc in scores:
                _line += (str(sc) + '\t')
            print >> score_out, _line[:-1]
        elif options.eval_strategy in [1, 2, 3]:
            ### check equal length of scores
            idx_list = set()
            for subm in scores.keys():
                idx_list = idx_list.union(set(scores[subm].keys()))

            for subm in scores.keys():
                _line = ''
                for sc in idx_list:
                    try:
                        #_line += (str(float(scores[subm][sc]) / float(max(1, alignments[subm][sc]))) + '\t')
                        _line += (str(float(scores[subm][sc])) + '\t')
                    except KeyError:
                        _line += ('0.0\t')
                print >> score_out, _line[:-1]
        #elif options.eval_strategy == 3:
        #    for subm in scores.keys():
        #        _line = ''
        #        for stratum in scores[subm]:
        #            _line += (str(scores[subm][stratum]) + '\t')
        #        print >> score_out, _line[:-1]

        score_out.close()
    
        
        if options.eval_strategy in [1, 2, 3]:
            support_file = options.score.replace('.dat','') 
            support_file += ('.support.strat_' + str(options.eval_strategy) + '.dat')
            support_out = open(support_file, 'w')
            idx_list = set()
            for subm in scores.keys():
                idx_list = idx_list.union(set(scores[subm].keys()))

            for subm in scores.keys():
                _line = ''
                for sc in idx_list:
                    try:
                        _line += (str(alignments2[subm][sc]) + '\t')
                    except KeyError:
                        _line += ('0\t')
                print >> support_out, _line[:-1]
            support_out.close()

            max_outname = options.score.replace('.dat','') 
            max_outname += ('.max_score.strat_' + str(options.eval_strategy) + '.dat')
            max_out = open(max_outname, 'w')
            idx_list = set()
            for subm in scores.keys():
                idx_list = idx_list.union(set(scores[subm].keys()))

            _line = ''
            for sc in idx_list:
                try:
                    _line += (str(max_score[sc]) + '\t')
                except KeyError:
                    _line += ('0\t')
            print >> max_out, _line[:-1]
            max_out.close()

        if options.labels != '-':
            out_name = options.score
            out_name = out_name.replace('.dat','')
            out_name += ('.strat_' + str(options.eval_strategy) + '.labels') 
            label_out = open(out_name, 'w')
            _line = ''
            for subm in scores.keys():
                _line += (label_list[subm] + '\t')
            print >> label_out, _line[:-1]
            label_out.close()

    ### print summary information
    print 'processed %s reads - found %s multimappers [ %s percent ]' % (read_counter, multireads, float(multireads) * float(100) / float(max(read_counter, 1)))  

if __name__ == '__main__':
    main()
