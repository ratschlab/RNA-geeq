""" This module contains all functionality to filter alignments"""
import sys

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


def filter_alignment_sam(infile, filterset, outfile, config):
    """This procedure parses alignment lines from infile and writes passing lines to outfile"""

    counter = 0
    filter_counter = 0

    for line in infile:
        if counter % 1000000 == 0 and config.verbose:
            print >> sys.stderr, "lines read: [ %s (taken: %s / filtered: %s)]" % (counter, counter - filter_counter, filter_counter)
        counter += 1
        if line[0] == '@':
            print >> outfile, line, 
            continue
        sl = line.strip().split()
        if len(sl) < 9:
            print >> sys.stderr, "\nERROR: Alignment line has not expected format: SAM\n%s" % line
            sys.exit(2)

