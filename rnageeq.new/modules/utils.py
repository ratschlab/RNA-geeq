"""The module utils contains functions of general applicability."""

import subprocess
import sys 
import re

def get_infile_handles(config):
    """Reads the infile strings from the config structure and returns an array with alignment file handles"""

    ### check what kind of input format we have
    infile_handles = []
    if config.infile == '-':
        infiles = [sys.stdin]
    else:
        infiles = []
        for inf in config.alignment.strip().split(','):
            if inf == '':
                continue
            if config.input_format == 'bam':
                infile_handles.append(subprocess.Popen([config.samtools, 'view', '-h', config.infile], stdout=subprocess.PIPE))
                infiles.append(infile_handles[-1].stdout)
            elif config.input_format == 'sam':
                infiles.append(open(config.infile, 'r'))
            else:
                print >> sys.stderr, "\nERROR: Unknown alignment file format: %s" % config.input_format

    return (infiles, infile_handles)

def get_outfile_handles(config):
    """Reads the outfile strings from the config structure and returns an array with output file handles"""

    ### check what kind of output format we have
    if config.outfile == '-':
        outfile = sys.stdout
    else:
        if config.output_format == 'bam':
            outfile_handle = subprocess.Popen([config.samtools, 'view', '-h', config.outfile], stdout=subprocess.PIPE)
            outfile = outfile_handle.stdout
        elif config.output_format == 'sam':
            outfile = open(config.outfile, 'r')
        else:
            print >> sys.stderr, "\nERROR: Unknown output file format: %s" % config.output_format

    return outfile

def get_mismatch_info(config, tag_list):
    """This function extracts the appropriate mismatch information specified in config from the given tag_list"""

    found = False
    mm = 0
    for opt in tag_list:
        if config.variants and opt[:3] in ['XG:', 'XM:']:
            mm += int(opt[5:])
            found = True
        if not config.variants and opt[:3] == 'NM:':
            mm = int(opt[5:])
            found = True
            break
    if not found:
        print >> sys.stderr, 'ERROR: No mismatch information available or read string missing in %s' % config.infile
        if config.variants:
            print >> sys.stderr, 'Alignment has to contain NM tag!'
        else:
            print >> sys.stderr, 'Alignment has to contain tags XM and XG!'
        sys.exit(1)
        ### TODO get the currently processed file id in here somehow
    
    return mm
 
def get_qpalma_score(tag_list):
    """This function extracts the QPALMA score from the given tag_list"""

    found = False
    qp = 0.0
    for opt in tag_list:
        if opt[:5] == 'ZS:f:':
            qp = float(opt[5:])
            found = True
            break

    if not found:
        print >> sys.stderr, 'ERROR: No QPalma Score for sorting available! Bailing out! - This option is only available for PALMapper alignments'
        sys.exit(1)
    
    return qp

def get_min_segment_len(cigar): 
    """This function takes a cigar string and returns the minimal segment length"""

    cigar = re.sub('[0-9]*[IHS]', '', cigar) 
    min_ex_len = sys.maxint
    for _cig in cigar.strip().split('N'):
        min_ex_len = min(sum([int('0'+i) for i in re.split('[^0-9]', '0' + _cig + 'Z0')][:-2]), min_ex_len)

    return min_ex_len
 
def get_filter_from_file(infile, filterset=None):
    """Takes a file containing filter information and returns a filterset object"""

    ### filterset contains filter-ID/comparator pairs as keys and an integer value
    if not filterset:
        filterset = dict()

    ### lines in the filter file need to have following format
    ### tab or space separated:
    ###    <filtername> <integer> <comparator>
    ### comparator is one of: lt, le, gt, ge, eq
    ### commentary lines must start with '#'

    ### iterate over lines and fill struct
    for line in open(infile, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split()

        filterset[(sl[0], sl[2])] = int(sl[1])

    return filterset

def get_filter_from_config(config, filterset=None):
    """Takes the config object, extracts the filter information and returns a filterset object"""

    ### filterset contains filter-ID/comparator pairs as keys and an integer value
    if not filterset:
        filterset = dict()

    if config.max_mismatch > -1:
        filterset[('mismatch', 'le')] = config.max_mismatch
    if config.min_segment_len > -1:
        filterset[('exon_len', 'ge')] = config.min_segment_len
    if config.min_ss_cov > -1:
        filterset[('ss_cov', 'ge')] = config.min_ss_cov
    if config.max_ss_cov > -1:
        filterset[('ss_cov', 'le')] = config.min_ss_cov

    return filterset

def get_segment_features(config, infiles, infile_handles=list()):
    """Takes infiles to read segment information and stores them in a feature dictionary"""

    import cPickle

    ### init variables
    line_counter = 0
    skipped_lines = 0
    features = dict()

    for infile in infiles:
        for line in infile:
            sl = line.strip().split()

            if config.verbose and line_counter % 10000 == 0:
                print '[ lines read: %i / skipped: %i ]\r' % (line_counter, skipped_lines),
            line_counter += 1

            if len(sl) < 9:
                skipped_lines += 1
                continue
            elif sl[5].find('N') == -1:
                skipped_lines += 1
                continue
            elif len(sl) < 12:
                print >> sys.stderr, 'Mismatch information is missing in %s - No tags found!' % config.alignment
                sys.exit(-1)
        
            min_ex_len = get_min_segment_len(sl[5])
            mm = get_mismatch_info(config, sl[11:])

            offset = 0
            (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
            size = [int(i) for i in size]
            
            start = int(sl[3]) - 1
            chrm = sl[2]
            ### TODO: offer a way to substite chromosome names

            for o in range(len(op)):
                if op[o] == 'N':
                    istart = start + offset
                    iend = istart + size[o]
                    try:
                        features[(chrm, istart, iend)][set([min_ex_len, mm])] += 1
                    except KeyError:
                        features[(chrm, istart, iend)][set([min_ex_len, mm])] = 1
                if not op[o] in ['I', 'H', 'S']:
                    offset += size[o]

        infile.close()
        if config.bam_input:
            infile_handles[infiles.index(infile)].kill()

    ### store features if requested
    if config.outfile != '-':
        outfile = open(config.outfile , 'w')
        cPickle.dump(features, outfile, cPickle.HIGHEST_PROTOCOL)
        outfile.close()


