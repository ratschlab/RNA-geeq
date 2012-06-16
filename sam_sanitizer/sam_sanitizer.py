"""This script checks a given file in SAM format for correctness
and sanitizes format issues as far as possible.

SAM format according to:
   http://samtools.sourceforge.net/SAM1.pdf

Version: 0.0.1

Written by: Andre Kahles (Friedrich Miescher Laboratory, Tuebingen, Germany)
            andre.kahles'at'tuebingen.mpg.de
            2010-08-01
"""

import sys
import re
import genome_utils
import cPickle

class Profile(object):
    """Represents an error profile containing all error information"""

    def __init__(self):
        """Data type initialisation"""
        self.store_dict = dict()
        ### general errors
        self.store_dict['type_error'] = dict()    # dict of fields with type error, containing the lines with the error
        self.store_dict['range_error'] = dict()    # dict of fields with range error, containing the lines with the error
        self.store_dict['read_len'] = []          # list of lines with too long/short read strings
        self.store_dict['quality_len'] = []          # list of lines with too long/short quality strings
        self.store_dict['read_vs_qual_len'] = []    # list of lines with not matching read and quality lengths
        self.store_dict['wrong_mate_info'] = []     # list of pairs of lines where mate information does not fit (start and id)
        self.store_dict['quality_score_1'] = []    # list of lines where quality score is neither sanger, nor illumina
        self.store_dict['quality_score_2'] = False    # different quality scores are used in the same file (e.g., Sanger and Illumina)
        self.store_dict['chrm_names'] = []         # list of lines where chromosom/reference names do not fit scheme / regexp / given list file
        self.store_dict['pairinfo_in_id'] = []    # list of lines where pair information is stored in read id
        self.store_dict['insert_size'] = []        # list of pairs of lines where inferred insert size seems to be wrong
        self.store_dict['start_one_based'] = []     # list of lines where start is not one based
        self.store_dict['too_many_mismatches'] = []   # list of lines, where number of re-aligned mismatches exceeds given maximum
        ### header errors
        self.store_dict['header_missing'] = False # header information is completely missing
        self.store_dict['header_missing_file_format'] = False    # file format info in header is missing
        self.store_dict['header_missing_seq_name'] = False    # sequence name (chromosome) is missing
        self.store_dict['header_missing_seq_len'] = False     # sequence (chromosome) length is missing
        self.store_dict['header_missing_rg_id'] = False         # read group id is missing
        self.store_dict['header_missing_rg_pool'] = False        # read group pool information is missing
        self.store_dict['header_missing_pg_name'] = False        # program name is missing
        ### flag errors
        self.store_dict['flag_paired_nonpair'] = [] # list of lines where unpaired reads have been mapped as pair
        self.store_dict['flag_same_pair_info'] = [] # list of pairs of linenumbers where mapped pairs have same pair info (both 64 or both 128)
        self.store_dict['flag_both_pair_info'] = [] # list of lines where one read has both pair information (64 and 128)
        self.store_dict['flag_missing_pair_info'] = [] # pair information is missing for paired reads
        self.store_dict['flag_mate_mapping'] = [] # list of pairs of linenumbers where mate is indicated a unmapped but is mapped (vice versa)
        self.store_dict['flag_mate_orientation'] = [] # list of pairs of linenumbers where orientation and mate orientation do not fit
        self.store_dict['flag_multimapper'] = [] # list of lines, where flag for multimappers is missing
        ### cigar errors
        self.store_dict['cigar_code'] = dict()     # dict over non standard letters used, containing the lines of occurance
        self.store_dict['cigar_length'] = []        # list of lines with not matching cigar lengths
        self.store_dict['cigar_orientation'] = [] # list of lines with wrong CIGAR orientation
        ### tag errors
        self.store_dict['tag_type'] = dict()    # dict of tags where tag has non valid type, contains list of affected lines
        self.store_dict['tag_name'] = dict()    # dict of tags where tag name is not standrad conform, contains list of affected lines
        self.store_dict['tag_value'] = dict()    # dict of tags where tag value is not valid, contains list of affected lines
        self.store_dict['tag_wrong_mismatch'] = []    # list of lines, where given mismatch info differs from re-aligned mismatch info

        
    def add_type_error(self, key, line):
        """Adds type error to profile struct"""
        try:
            self.store_dict['type_error'][key].append(line)
        except KeyError:
            self.store_dict['type_error'][key] = [line]
            
    def add_range_error(self, key, line):
        """Adds range error to profile struct"""
        try:
            self.store_dict['range_error'][key].append(line)
        except KeyError:
            self.store_dict['range_error'][key] = [line]
            
    def add_cigar_code_error(self, key, line):
        """Adds cigar code error to profile struct"""
        try:
            self.store_dict['cigar_code'][key].append(line)
        except KeyError:
            self.store_dict['cigar_code'][key] = [line]
            
    def add_tag_value_error(self, key, line):
        """Adds cigar code error to profile struct"""
        try:
            self.store_dict['tag_value'][key].append(line)
        except KeyError:
            self.store_dict['tag_value'][key] = [line]
            
    def add_tag_name_error(self, key, line):
        """Adds cigar code error to profile struct"""
        try:
            self.store_dict['tag_name'][key].append(line)
        except KeyError:
            self.store_dict['tag_name'][key] = [line]
            
    def add_tag_type_error(self, key, line):
        """Adds cigar code error to profile struct"""
        try:
            self.store_dict['tag_type'][key].append(line)
        except KeyError:
            self.store_dict['tag_type'][key] = [line]
            
    def print_overview(self, options):
        """This function prints an overview of all non empty error fields in the profile and gives
        additional information to the found errors and advice how to correct them."""

        if options.phase == 1:
            ### general errors
            if len(self.store_dict['type_error']) > 0:
                print "Wrong data type:"
                for key in self.store_dict['type_error']:
                    print '\terror in %s - %i lines' % (key, len(self.store_dict['type_error'][key])) 
                print "\n"
            if len(self.store_dict['range_error']) > 0:
                print "Wrong data range:"
                for key in self.store_dict['range_error']:
                    print '\terror in %s - %i lines' % (key, len(self.store_dict['range_error'][key])) 
                print "\n"
            if len(self.store_dict['read_len']) > 0:
                print "Lines with too long or too short read strings:"
                print '\t affected lines: %i\n' % len(self.store_dict['read_len']) 
            if len(self.store_dict['quality_len']) > 0:
                print "Lines with too long or too short quality strings:"
                print '\t affected lines: %i\n' % len(self.store_dict['quality_len']) 
            if len(self.store_dict['read_vs_qual_len']) > 0:
                print "Lines where lengths of read and quality string do not match:"
                print '\t affected lines: %i\n' % len(self.store_dict['read_vs_qual_len']) 
            if len(self.store_dict['quality_score_1']) > 0:
                print "Quality score seems to be neither Sanger nor Illumina:"
                print '\t affected lines: %i\n' % len(self.store_dict['quality_score_1']) 
            if self.store_dict['quality_score_2']:
                print "Different quality scores seem to be used (Sanger and Illumina):"
                print '\t YES\n' 
            if len(self.store_dict['chrm_names']) > 0:
                print "Chromosome names do not fit given pattern:"
                print '\t affected lines: %i\n' % len(self.store_dict['chrm_names']) 
            if len(self.store_dict['pairinfo_in_id']) > 0:
                print "Pairing information is stored in read ID:"
                print '\t affected lines: %i\n' % len(self.store_dict['pairinfo_in_id']) 
            if len(self.store_dict['wrong_mate_info']) > 0:
                print 'Mate information ofalignment pair is not consistent:'
                print '\t affected lines: %i\n' % len(self.store_dict['wrong_mate_info'])
            if len(self.store_dict['insert_size']) > 0:
                print 'Inferred insert siz seems to be wrong:'
                print '\t affected lines: %i\n' % len(self.store_dict['insert_size'])
            if len(self.store_dict['start_one_based']) > 0:
                print 'Alignment start position is not 1-based:'
                print '\t affected lines: %i\n' % len(self.store_dict['start_one_based'])
            ### header errors
            if self.store_dict['header_missing']:
                print 'Header information is missing\n'
            if self.store_dict['header_missing_file_format']:
                print 'File format information in header is missing\n'
            if self.store_dict['header_missing_seq_name']:
                print 'Sequence name (chromosome) is missing in header\n' 
            if self.store_dict['header_missing_seq_len']:
                print 'Sequence (chromosome) length is missing in header\n'
            if self.store_dict['header_missing_rg_id']:
                print 'Read group id is missing in header\n'
            if self.store_dict['header_missing_rg_pool']:
                print 'Read group pool information is missing in header\n'
            if self.store_dict['header_missing_pg_name']:
                print 'Program name is missing in header\n'
            ### flag errors
            if len(self.store_dict['flag_paired_nonpair']) > 0:
                print 'Unpaired reads have been mapped as pair:'
                print '\t affected lines: %i\n' % len(self.store_dict['flag_paired_nonpair'])
            if len(self.store_dict['flag_same_pair_info']) > 0:
                print 'Pair information in flag of mapped pair is the same:'
                print '\taffected lines: %i\n' % len(self.store_dict['flag_same_pair_info'])
            if len(self.store_dict['flag_both_pair_info']) > 0:
                print 'Alignment flag has both pair information:'
                print '\t affected lines: %i\n' % len(self.store_dict['flag_both_pair_info'])
            if len(self.store_dict['flag_missing_pair_info']) > 0:
                print 'Pairing information is missing in flag of paired reads:' 
                print '\t affected lines: %i\n' % len(self.store_dict['flag_missing_pair_info'])
            if len(self.store_dict['flag_mate_mapping']) > 0:
                print 'Flag indicates mate as unmapped, but mate is mapped:'
                print '\t affected lines: %i' % len(self.store_dict['flag_mate_mapping'])
            if len(self.store_dict['flag_mate_orientation']) > 0:
                print 'Orientation of read and its mate do not fit:'
                print '\t affected lines: %i\n' % len(self.store_dict['flag_mate_orientation'])
            if len(self.store_dict['flag_multimapped']) > 0:
                print 'Multimapper flag is missing'
                print '\t affected lines: %i\n' % len(self.store_dict['flag_multimapper'])
            ### cigar errors
            if len(self.store_dict['cigar_code']) > 0:
                print 'Following non standard letters have been used in the CIGAR:'
                for key in self.store_dict['cigar_code']:
                    print '\t key %s used in %i lines' % (key, len(self.cigar_code))
            if len(self.store_dict['cigar_length']) > 0:
                print 'Length of CIGAR does not match given read length'
                print '\t affected lines: %i' % len(self.store_dict['cigar_length'])
            if len(self.store_dict['cigar_orientation']) > 0:
                print 'Orientation of CIGAR is wrong'
                print 'affected lines: %i' % len(self.store_dict['cigar_orientation'])
            ### tag errors
            if len(self.store_dict['tag_type']) > 0:
                print 'Following tags do not have a valid type:'
                for key in self.store_dict['tag_type']:
                    print '\t tag with type %s are not valid in %i lines' % (key, self.store_dict['tag_type'])
            if len(self.store_dict['tag_name']) > 0:
                print 'Following tags do not have a valid name:'
                for key in self.store_dict['tag_name']:
                    print '\t tag with name %s are not valid in %i lines' % (key, self.store_dict['tag_name'])
            if len(self.store_dict['tag_value']) > 0:
                print 'Following tags do not have a valid value:'
                for key in self.store_dict['tag_value']:
                    print '\t tag with type %s are not valid in %i lines' % (key, self.store_dict['tag_value'])

### types of errors:
## only some errortypes are a dict with SAM fields as keys, thus the error can occurr in each field
# - flag is contradictory
# -- pair info in mapped pair is the same
# -- mate is unmapped, but indicated as mapped (vice versa)
# -- in mapped pairs orientation of read and orientation of mate's mate are different
# -- flag for multiple alignment missing
# - mate information for mate is wrong / missing
# - id of paired reads is the same (pair info must not be stored in read id)
# - check inferred insert size
# - try re-alignment to verify:
# -- CIGAR string (orientation, usage of letter code)
# -- start is one-based
# -- tags (mismatches, perfect positions)
# - verify header information
# count:
#  - update number of reads, number of aligns and multimappers

# -check correctnes of read seq, orientation and alignment (mismatches) and CIGAR using the given fastq and genome file
# -re-alignment using needleman wunsch


    def store(self, options):
        """Stored profile information for later usage"""

        outfile = open(options.infile, 'w')
        cPickle.dump(self.store_dict, outfile)
        outfile.close()



def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-i', '--infile', dest='infile', metavar='FILE', help='alignment file in sam format', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-o', '--outfile', dest='outfile', metavar='PATH', help='outfile', default='-')
    optional.add_option('-P', '--phase', dest='phase', metavar='INT', type='int', help='phase: 0 - summary, 1 - check, 2 - sanitize [0]', default=0)
    optional.add_option('-g', '--genome_file', dest='genome_file', metavar='FILE', help='genome config file', default='-')
    optional.add_option('-c', '--chrm_cond', dest='chrm_cond', metavar='FILE', help='file with valid chromosome names or regexp per line', default='-')
    optional.add_option('-r', '--readfile', dest='readfile', metavar='FILE', help='file containing all reads in FASTQ format', default='-')
    optional.add_option('-l', '--readlength', dest='readlength', metavar='INT', type='int', help='length of given read [75]', default=75)
    optional.add_option('-M', '--max_mismatch', dest='max_mismatch', metavar='INT', type='int', help='maximal number of mismatches allowed per alignment [10]', default=10)
    optional.add_option('-k', '--k_last_pair', dest='pair_id', metavar='INT', type='int', help='number of id positions at the end discriminating read pairs [2]', default=2)
    optional.add_option('-K', '--k_last_pair_pattern', dest='pair_id_pattern', metavar='STRING', help='pattern of pair id identification [/1|:1|/A|:A,/2|:2|/B|:B]', default='/1|:1|/A|:A,/2|:2|/B|:B')
    optional.add_option('-L', '--lines', dest='lines', metavar='INT', type='int', help='number of sam lines for sanity checking (header excluded) [100]', default=100)
    optional.add_option('-p', '--paired', dest='paired', action='store_true', help='reads were paired in sequencing [not set]', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity [not set]', default=False)
    optional.add_option('-d', '--debug', dest='debug', action='store_true', help='display debug info [not set]', default=False)
    optional.add_option('-H', '--hard', dest='hard', action='store_true', help='try hard to minimize mismatches [not set]', default=False)
    optional.add_option('-N', '--ignore_Ns', dest='ignore_Ns', action='store_false', help='if set, N in a sequence does not count for the distance [set]', default=True)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def ed(str1, str2):
    """Edit distance for two equally long strings. """

    if (len(str1) == len(str2)):
        return sum([(1 - int(str1.upper()[i] == str2.upper()[i] or str1.upper()[i] == 'N' or str2.upper()[i] == 'N')) for i in range(len(str1))])
    else:
        print >> sys.stderr, "strings to compare have not same length!: \n%s \n%s" % (str1, str2)
        sys.exit(2)

def check_bundle(bundle, profile, options):
    """Checks correctness of alignments in one bundle"""
   
    ### check read id
    idx = -1
    for (sl, line_counter) in bundle:
        idx += 1
        # look for pair information
        try:
            flag = int(sl[1])
            if options.paired:
                pair_bundle = [[], []]
                if flag & 192 in [64, 128]:
                    if flag & 64 == 64:
                        pair_bundle[0].append(idx)
                    else:
                        pair_bundle[1].append(idx)
                else:
                    if flag & 192 == 192:
                        profile.flag_both_pair_info.append(line_counter)
                    else:
                        profile.flag_missing_pair_info.append(line_counter)
                    ### try to find pair info in read id
                    if re.match(options.pair_id_pattern.split(',')[0], sl[0][-options.pair_id]) != None:
                        pair_bundle[0].append(idx)
                        profile.pairinfo_in_id.append(line_counter)
                    if re.match(options.pair_id_pattern.split(',')[1], sl[0][-options.pair_id]) != None:
                        pair_bundle[1].append(idx)
                        profile.pairinfo_in_id.append(line_counter)
            else:
                if flag & 192 in [64, 128]:
                    profile.flag_paired_nonpair.append(line_counter)
        except ValueError:
            continue

    #    if options.paired:
    #        for align1 in pair_bundle[0]:
    #            for align2 in pair_bundle[1]:
                    
        profile = check_sam_line(sl, profile, line_counter, options)
        
    return profile

def check_sam_line(sl, profile, line_counter, options):
    """Checking correctness of one SAM line, storing error information in profile"""
    
    # check flag
    try:
        flag = int(sl[1])
        if flag < 0:
            profile.add_range_error('flag', line_counter)
        if flag & 2 == 2 and flag & 1 == 0:
            profile.flag_paired_nonpair.append(line_counter)

    except ValueError:
        profile.add_type_error('flag', line_counter)


    # check chrm / contig
    chrm = sl[2]
    if options.chrm_cond != '-' and re.match(options.chromosomes, chrm) == None:
        profile.chrm_names.append(line_counter)

    # check start pos
    try:
        start = int(sl[3])
        if start < 0:
            profile.add_range_error('start', line_counter)
            start = None
    except ValueError:
        profile.add_type_error('start', line_counter)   
        start = None

    # check mapping quality
    try:
        map_qual = int(sl[4])
        if map_qual < 0 or map_qual > 255:
            profile.add_range_error('map_qual', line_counter)
    except ValueError:
        profile.add_type_error('map_qual', line_counter)

    # check CIGAR
    cigar = sl[5]
    (size, op) = (re.split('[^0-9]', cigar)[:-1], re.split('[0-9]*', cigar)[1:])
    size = [int(i) for i in size]
    cig_len = sum(size)
    if len(set(op).difference(options.cigar_letters)) > 0:
        for letter in set(op).difference(options.cigar_letters):
            profile.add_cigar_code_error(letter, line_counter)
            
    # mate_id
    mate_id = sl[6]
    # mate start
    try:
        mate_start = int(sl[7])
        if mate_start < 0:
            profile.add_range_error('mate_start', line_counter)
    except ValueError:
        profile.add_type_error('mate_start', line_counter)
    # insert size
    try:
        isize = int(sl[8])
    except ValueError:
        profile.add_type_error('isize', line_counter)

    # check sequence
    seq = sl[9]
    if seq != '*' and len(seq) != options.readlength:
        profile.read_len.append(line_counter)

    if cig_len != len(seq):
        profile.cigar_length.append(line_counter)
        
    # check quality
    quality = sl[10]
    if quality != '*' and len(quality) != options.readlength:
        profile.quality_len.append(line_counter)
    if quality != '*' and seq != '*' and len(quality) != len(seq):
        profile.read_vs_quality_len.append(line_counter)

    ### sanger 0+33 - 93+33 (33-126) 
    ### illumina1.3+ 0+64 - 62+64 (64-126)
    ### illumina1.5+ 2+64 - 62+64 (66-126) / special usage for 66
    for q in quality:
        if ord(q) < 33 or ord(q) > 126:
            profile.quality_score_1.append(line_counter)
            break
    min_q = min([q for q in quality])
    if ord(min_q) < options.qual_min:
        options.qual_min = ord(min_q)
    max_q = max([q for q in quality])
    if ord(max_q) < options.qual_max:
        options.qual_max = ord(max_q)
    if options.qual_min < 62 and options.qual_max > 93:
        profile.quality_score_2.append(line_counter)

    # check tags
    tags = sl[11:]
    tag_dict = dict()
    for tag in tags:
        _tag = tag.strip().split(':')   
        cont = False
        if re.match(options.tag_names, _tag[0]) == None:
            profile.add_tag_name_error(_tag[0], line_counter)
            cont = True
        if not _tag[1] in options.tag_types:
            profile.add_tag_type_error(_tag[1], line_counter)
        if (_tag[1] == 'i' and not _tag[1].isdigit()) or \
           (_tag[1] == 'A' and len(_tag[2]) > 1):
            profile.add_tag_value_error(_tag[2], line_counter)
            cont = True
        if _tag[1] == 'H':
            try:
                int(_tag[2], 16)
            except ValueError:
                profile.add_tag_value_error(_tag[2], line_counter)
                cont = True
        if _tag[1] == 'f':
            try:
                float(_tag[2])
            except ValueError:
                profile.add_tag_value_error(_tag[2], line_counter)
                cont = True
        if cont:
            continue
        tag_dict[_tag[0]] = _tag[2]
                
    # run additional checks, if necessary information is available
    if options.genome_file != '-' and seq != '*' and start != None:
        ### load genome if necessary
        if not options.genome.has_key(sl[2]):
            if sl[2] in options.gio.contig_names:
                if options.verbose:
                    print 'Reading chromosome %s' % sl[2]
                fgen = open('%s/genome/%s.flat' % (options.gio.basedir, sl[2]))
                options.genome[sl[2]] = fgen.readline()
                fgen.close()
            else:       
                print >> sys.stderr, 'Chromosome Names do not match genome information file. Chromosome %s can not be found!' % sl[2]
                exit(2)

        ### reconstruct genome_part from CIGAR
        gen = ''
        covered = 0
        skipped = 0
        deletions = 0
        insertions = 0
        for o in range(len(op)):
            if op[o] in ['M', 'S']:
                gen += options.genome[sl[2]][start + covered + skipped : start + covered + skipped + size[o]].upper()
                covered += size[o]
            elif op[o] == 'I':
                gen += (size[o] * 'X') ### insertions go into distance
                insertions += size[o]
            elif op[o] == 'D':
                skipped += size[o]
                deletions += size[o]
            elif op[o] == 'N':
                skipped += size[o]
        
        rev_gen = genome_utils.reverse_complement(gen)

        if (covered + insertions) != len(gen):
            print >> sys.stderr, "WARNING: CIGAR length could net be resembled - maybe reached end of chromosome - ignore line" 
            return profile
            
        if len(gen) > options.readlength:
            print >> sys.stderr, "WARNING: CIGAR was longer then given readlength - ignore line" 
            return profile

        d1 = 0
        d2 = 0
        sc1 = 0
        sc2 = 0
        seq_len = len(seq)

        if op[0] == 'H':
            d1 = size[0]
        if op[-1] == 'H':
            d2 = size[-1]
        if op[0] == 'S':
            d1 = size[0]
            sc1 = size[0]
        if op[-1] == 'S':
            d2 = size[-1]
            sc2 = size[-1]

        dist = [options.readlength for i in range(5)]

        if options.hard:
            dist[0] = ed(seq[d1:seq_len - d2], gen[sc1:seq_len - sc2])
            dist[1] = ed(seq[d1:seq_len - d2], rev_gen[sc1:seq_len - sc2]) 
            dist[2] = ed(seq[d2:seq_len - d1], gen[sc2:seq_len - sc1])
            dist[3] = ed(seq[d2:seq_len - d1], rev_gen[sc2:seq_len - sc1]) 
            min_d_ind = dist.index(min(dist))
            min_d = dist[min_d_ind]
        else:
            if flag & 16 == 0:
                min_d = ed(seq[d1:seq_len - d2], gen[sc1:seq_len - sc2])
            else:
                min_d = ed(seq[d2:seq_len - d1], rev_gen[sc2:seq_len - sc1])

        if tag_dict.has_key('NM') and int(tag_dict['NM']) != min_d:
            profile.tag_wrong_mismatch.append(line_counter)
        if min_d > options.max_mismatche:
            profile.too_many_mismatches.append(line_counter)

    return profile

def main():
    """This function handles the main program flow"""

    options = parse_options(sys.argv)

    line_counter = 0
    header_counter = 0
    skipped_counter = 0
    bundle_dict = dict()

    profile = Profile()

    options.cigar_letters = set(['M', 'D', 'I', 'S', 'H', 'P', 'N'])
    options.qual_min = 256
    options.qual_max = 0
    options.tag_names = 'X.|Y.|Z.|PU|PG|AS|SQ|OQ|E2|U2|MQ|NM|H0|H1|H2|UQ|PQ|NH|IH|HI|MD|CS|CQ|CM|R2|Q2|S2|CC|CP|SM|AM|MF'
    options.tag_types = ['A', 'i', 'H', 'z', 'f']
    options.gio = genome_utils.GenomeInfo(options.genome_file)
    options.gio.contig_names.sort()
    options.genome = dict()
    if options.chrm_cond != '-':
        options.chromosomes = ''
        for line in open(options.chrm_cond):
            options.chromosomes.append(line + '|')
        options.chromosomes = options.chromosomes[:-1]
            

    if options.verbose:
        print "Parsing %i lines from alignment file %s (excluding header lines)" % (options.lines, options.infile)
        print "--------------------------------------------------"

    for line in open(options.infile, 'r'):

        ### check header information
        if line[0] == '@':
            # do some stuff TODO
            header_counter += 1
            continue

        if line_counter > options.lines:
            break
        line_counter += 1

        sl = line.strip().split('\t')

        if len(sl) < 9:
            skipped_counter += 1
            continue
        else:
            ### one bundle constists of all lines having the same read_id (excluding the last k <default: 2> positions)
            if options.paired > 0:
                bundle_id = sl[0][:-options.pair_id]
            else:
                bundle_id = sl[0]
            try:
                bundle_dict[bundle_id].append((sl, line_counter))
            except KeyError:
                bundle_dict[bundle_id] = [(sl, line_counter)]

    if options.verbose:
        print "Parsed %i alignment lines" % line_counter
        print "Parsed %i header lines" % header_counter
        print "Skipped %i lines" % skipped_counter

    for bundle in bundle_dict:
        profile = check_bundle(bundle_dict[bundle], profile, options)

    if int(options.phase) in [1, 2]:
        profile.print_overview(options)
    else:
        profile.store(options)

    #### store all lines in a dict (key = read_id) -> iterate over the dict, important for sorting and analysing paired reads and multimappers
    


if __name__ == '__main__':
    main()
