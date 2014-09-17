import sys
import scipy as sp
import pdb
import subprocess
import re
import time

from sputils import *

def get_tags_gff(tagline):
    """Extract tags from given tagline in GFF format"""

    tags = dict()
    for t in tagline.split(';'):
        tt = t.split('=')
        tags[tt[0]] = tt[1]
    return tags

def get_tags_gtf(tagline):
    """Extract tags from given tagline in GTF format"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags


def parse_anno_from_file(options):
    """This function reads the gff3 input file and returns the information in an
       internal data structure"""

    if options.verbose:
        print >> sys.stderr, "Parsing annotation from %s ..." % options.anno

    ### perform inital parse of the gff file to gather gene - transcript - exon relations
    trans2gene = dict() ### dict with: keys = transcript IDs, values = gene IDs
    gene2trans = dict() ### dict with: keys = gene IDs, values = list of transcript IDs
    trans2exons = dict() ### dict with: keys = transcript IDs, values = list of exon coords
    gene_pos = dict()

    if options.anno.lower().endswith('gff') or options.anno.lower().endswith('gff3'):
        format='gff'
    elif options.anno.lower().endswith('gtf'):
        format='gtf'

    ### only show available sources - if neccessary
    if options.show_sources:
        sources = []
        for line in open(options.anno, 'r'):
            sources.append(line.strip().split('\t')[1])
        sources = sp.unique(sources)
        print 'Parsed file %s\n' % options.anno
        print 'Following sources are available:\n'
        for source in sources:
            print source    
        print '\nUse option -s to specify a comma-separated list of sources (-s source1,source2,source3), otherwise all sources are taken'
        sys.exit(0)

    ### assume all coordinates as 1 based an half-open intervals
    any_source_taken = False
    for line in open(options.anno, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        if options.sources[0] != '-' and not sl[1] in options.sources:
            continue
        any_source_taken = True
        if format == 'gtf':
            tags = get_tags_gtf(sl[8])
        else:
            tags = get_tags_gff(sl[8])
        if sl[2].lower() in ['transcript', 'mrna', 'mirna', 'trna', 'snrna', 'snorna', 'ncrna', 'mrna_te_gene', 'rrna', 'pseudogenic_transcript', 'transposon_fragment']:
            if format == 'gtf':
                try:
                    trans2gene[tags['transcript_id']] = tags['gene_id']
                except:
                    pdb.set_trace()
                try:
                    gene2trans[tags['gene_id']].append(tags['transcript_id'])
                except KeyError:
                    gene2trans[tags['gene_id']] = [tags['transcript_id']]
            else:
                trans2gene[tags['ID']] = tags['Parent']
                try:
                    gene2trans[tags['Parent']].append(tags['ID'])
                except KeyError:
                    gene2trans[tags['Parent']] = [tags['ID']]
        elif sl[2].lower() == 'exon':
            ### assume positions are 1 based and in half open intervals for plus strand
            ### assume positions are 0 based and in half open intervals for minus strand --> this is an artifact of our anno pipeline ...
            if sl[6] == '-':
                strand_offset = 0
            else:
                strand_offset = -1
            if format == 'gtf':
                try:
                    #trans2exons[tags['transcript_id']].append((sl[0], int(sl[3]) + strand_offset, int(sl[4]) + strand_offset - 1)) ### TODO account for half-open intervals
                    trans2exons[tags['transcript_id']].append((sl[0], int(sl[3]), int(sl[4]))) ### TODO account for half-open intervals
                except KeyError:
                    #trans2exons[tags['transcript_id']] = [(sl[0], int(sl[3]) + strand_offset, int(sl[4]) + strand_offset - 1)]
                    trans2exons[tags['transcript_id']] = [(sl[0], int(sl[3]), int(sl[4]))]
            else:
                try:
                    #trans2exons[tags['Parent']].append((sl[0], int(sl[3]) + strand_offset, int(sl[4]) + strand_offset - 1)) ### TODO account for half-open intervals
                    trans2exons[tags['Parent']].append((sl[0], int(sl[3]), int(sl[4]))) ### TODO account for half-open intervals
                except KeyError:
                    #trans2exons[tags['Parent']] = [(sl[0], int(sl[3]) + strand_offset, int(sl[4]) + strand_offset - 1)]
                    trans2exons[tags['Parent']] = [(sl[0], int(sl[3]), int(sl[4]))]
        elif sl[2].lower() == 'gene':
            ### assume positions are 1 based and in half open intervals for plus strand
            ### assume positions are 0 based and in half open intervals for minus strand --> this is an artifact of our anno pipeline ...
            if sl[6] == '-':
                strand_offset = 0
            else:
                strand_offset = -1
            if format == 'gtf':
                try:
                    #gene_pos[tags['gene_id']] = (sl[0], sl[6], int(sl[3]) + strand_offset, int(sl[4]) + strand_offset)
                    gene_pos[tags['gene_id']] = (sl[0], sl[6], int(sl[3]), int(sl[4]))
                except ValueError:
                    print >> sys.stderr, 'Could not read coordinates of gene %s (%s to %s), set to (0 to 1)' % (tags['gene_id'], sl[3], sl[4])
                    gene_pos[tags['gene_id']] = (sl[0], '+', 0, 1)
            else:
                try:
                    #gene_pos[tags['ID']] = (sl[0], sl[6], int(sl[3]) + strand_offset, int(sl[4]) + strand_offset)
                    gene_pos[tags['ID']] = (sl[0], sl[6], int(sl[3]), int(sl[4]))
                except ValueError:
                    print >> sys.stderr, 'Could not read coordinates of gene %s (%s to %s), set to (0 to 1)' % (tags['ID'], sl[3], sl[4])
                    gene_pos[tags['ID']] = (sl[0], '+', 0, 1)

    if not any_source_taken:
        print >> sys.stderr, 'ERROR: The specified sources do not match any of the available sources - Please use option -S to get a list of available sources'
        sys.exit(-1)

    if options.verbose:
        print >> sys.stderr, '...done.'
            
    return {'genes' : gene_pos, 'gene2trans' : gene2trans, 'trans2exons' : trans2exons}


def get_mm_tag(sl, options):
    """Gets a list of tags and returns number of mismatches"""

    for o in sl:
        if o[:3] == 'NM:':
            return int(o[5:])

    print >> sys.stderr, 'No mismatch information available or read string missing in %s' % options.alignment



def intron_list_from_alignment(options, exclude_set=None):
    """Builds up an intron list from the given alignment and writes it to a file. Returns a map of all covered intron positions """

    ### parse sam:
    ### 0      1     2       3             4     5      6        7        8       9    10    11
    ### QUERY  FLAG  REFSEQ  POS(1-based)  MAPQ  CIGAR  MATEREF  MATEPOS  INSIZE  SEQ  QUAL  OPT
    #intron_lists = sp.zeros((0, 4), dtype='str')
    intron_lists = []
    cov_map_start = dict()
    cov_map_cont = dict()
    filter_counter = 0

    if options.bam_input:
        file_handle = subprocess.Popen([options.samtools, 'view', options.alignment], stdout=subprocess.PIPE) 
        infile = file_handle.stdout
    else:
        infile = open(options.alignment, 'r')

    t0 = time.time()
    for counter, line in enumerate(infile):
        if line[0] in ['@', '#'] or line[:2] == 'SQ':
            continue
        if options.lines > 0 and counter > options.lines:
            break
        if counter % 100000 == 0 and options.verbose:
            t1 = time.time()
            print 'lines read: [ %s (taken: %s / filtered: %s)] ... took %i sec' % (counter, counter - filter_counter, filter_counter, t1 - t0)
            t0 = t1

        sl = line.strip().split('\t')
    
        ### unaligned reads
        if (int(sl[1]) & 4) == 4:
            continue

        ### secondary reads
        if (int(sl[1]) & 256) == 256 and not options.include_secondary:
            continue

        if exclude_set is not None and (sl[0], (int(sl[1]) & 128 + int(sl[1]) & 64)) in exclude_set:
            filter_counter += 1
            continue

        if options.max_mismatches is not None:
            cont_flag = False
            mm = get_mm_tag(sl[11:], options)
            if mm > options.max_mismatches:
                filter_counte += 1
                continue

        if options.min_exon_len > 0:
            cont_flag = False
            __cig = sl[5]
            __cig = re.sub('[0-9]*[IHS]', '', __cig) 
            if min([sum([int(i) for i in re.split('[^0-9]', _cig + 'Z0')[:-2]]) for _cig in __cig.strip().split('N')]) < options.min_exon_len:
                filter_counter += 1
                continue

        ### handle Insertions - they do not affect intron length
        rl = re.sub('[0-9]*I', '', sl[5])
        
        (op, size) = (re.split('[0-9]*', rl)[1:], re.split('[^0-9]', rl)[:-1])
        size = [int(i) for i in size]
        offset = int(sl[3])# - 1
        covered = 0
        for p in range(len(op)):
            if op[p] in ['M', 'D']:
                ### handle starting position
                if options.full_analysis:
                    ### cov_map_start[chromosome] = [list if positions where alignments start]
                    ### cov_map_cont[chromosome] = [list if positions covered by alignments]
                    try:
                        cov_map_start[sl[2]].append(offset + covered)
                        cov_map_cont[sl[2]].append(offset + covered)
                    except KeyError:
                        cov_map_start[sl[2]] = [offset + covered]
                        cov_map_cont[sl[2]] = [offset + covered]
                        
                covered += 1
                ### handle following positions
                for pp in range(1, size[p]):                
                    if options.full_analysis:
                        cov_map_cont[sl[2]].append(offset + covered)
                    covered += 1
            elif op[p] == 'N':
                ### update intron lists
                ### intron_lists = [[chrm, strand, start, stop], [chrm, strand, start, stop], ...]
                if options.strand_specific:
                    if (int(sl[1]) & 128) == 128:
                        #intron_lists = sp.r_[intron_lists, [[sl[2], str((((int(sl[1]) & 16) / -8) + 1) * -1), str(offset + covered), str(offset + covered + size[p])]]]
                        intron_lists.append([sl[2], str((((int(sl[1]) & 16) / -8) + 1) * -1), str(offset + covered), str(offset + covered + size[p])])
                    else:
                        #intron_lists = sp.r_[intron_lists, [[sl[2], str(((int(sl[1]) & 16) / -8) + 1), str(offset + covered), str(offset + covered + size[p])]]]
                        intron_lists.append([sl[2], str(((int(sl[1]) & 16) / -8) + 1), str(offset + covered), str(offset + covered + size[p])])
                else:
                    #intron_lists = sp.r_[intron_lists, [[sl[2], '.', str(offset + covered), str(offset + covered + size[p])]]]
                    intron_lists.append([sl[2], '.', str(offset + covered), str(offset + covered + size[p])])
                covered += size[p]
            ### handle softclips and hardclips --> means ignoring them :)
            elif op[p] == 'S':
                covered += size[p]
            elif op[p] == 'H':
                continue
    infile.close()
    
    intron_lists = sp.array(intron_lists, dtype='str')

    ### make introns unique and count
    intron_lists = sort_rows(intron_lists)
    intron_lists_u, f_idx = unique_rows(intron_lists, index=True)
    l_idx = sp.r_[f_idx[1:], intron_lists.shape[0]]
    intron_support = l_idx - f_idx

    ### summarize coverage map positions
    cov_vect = dict()
    for key in cov_map_start:
        tmp = sp.sort(sp.array(cov_map_start[key]))
        starts, f_idx = sp.unique(tmp, return_index=True)
        l_idx = sp.r_[f_idx[1:], tmp.shape[0]]
        starts_cnt = l_idx - f_idx
        
        tmp = sp.sort(sp.array(cov_map_cont[key]))
        cont, f_idx = sp.unique(tmp, return_index=True)
        l_idx = sp.r_[f_idx[1:], tmp.shape[0]]
        cont_cnt = l_idx - f_idx
        
        cov_vect[key] = sp.zeros((cont.shape[0], 3), dtype='int')
        cov_vect[key][:, 0] = cont
        cov_vect[key][:, 1] = cont_cnt
        c_idx = sp.where(sp.in1d(cont, starts))[0]
        cov_vect[key][c_idx, 2] = starts_cnt

    if options.verbose:
        print '\n parsed file %s' % options.alignment
        print 'read %s lines' % counter
        if options.max_mismatches is not None or options.min_exon_len > 0:
            print 'filter criteria:'
            if options.max_mismatches is None:
                print '    max mismatches: not set'
            else:
                print '    max mismatches: %s' % options.max_mismatches
            print '    min exon length: %s\n' % options.max_mismatches
            print 'filtered %s lines' % filter_counter

    return (cov_vect, intron_lists_u, intron_support)

def intron_list_from_annotation(options):

    ### parse gff-file
    anno_dict = parse_anno_from_file(options)
    #intron_lists = sp.zeros((0, 4), dtype='str')
    intron_lists = []
    exon_map = dict()

    for g, gene in enumerate(anno_dict['genes']):
        if options.verbose and g > 0 and g % 100 == 0:
            print '.',
            if g % 1000 == 0:
                print '%i/%i' % (g, len(anno_dict['genes']))
        chrm = anno_dict['genes'][gene][0]
        strand = anno_dict['genes'][gene][1]
        if not chrm in exon_map:
            exon_map[chrm] = dict()
        for trans in anno_dict['gene2trans'][gene]:
            exons = []
            for exon in anno_dict['trans2exons'][trans]:
                exons.append(exon[1:])
            exons = sp.array(exons)
            s_idx = sp.argsort(exons[:, 0])
            exons = exons[s_idx, :]

            for e in range(exons.shape[0] - 1):
                ### intron_lists = [[chrm, strand, start, stop], [chrm, strand, start, stop], ...]
                #intron_lists = sp.r_[intron_lists, [[chrm, strand, str(exons[e, 1] + 1), str(exons[e + 1, 0])]]]
                intron_lists.append([chrm, strand, str(exons[e, 1] + 1), str(exons[e + 1, 0])])
            
            ### we assume, that an exon cannot occurr twice in the same transcript!
            ### the value in the dict is a binary encoding, if the left/right end is intronic 10 = 2 means, 5' end is intronic
            #if len(contig_list) == 1:
            #    exon_map[chrm][trans.id][contig_list[exon]] = 0 ### 00 -> should never occurr
            #elif exon == 0:
            #    exon_map[chrm][trans.id][contig_list[exon]] = 2 ### 10
            #elif exon == len(contig_list) - 1:
            #    exon_map[chrm][trans.id][contig_list[exon]] = 1 ### 01
            #else:
            #    exon_map[chrm][trans.id][contig_list[exon]] = 3 ### 11 

    intron_lists = sp.array(intron_lists, dtype='str')
    ### make intron_list unique
    tmp, u_idx = sp.unique(row_strings(intron_lists), return_index=True)
    intron_lists = intron_lists[u_idx, :]

    return (intron_lists, exon_map)

