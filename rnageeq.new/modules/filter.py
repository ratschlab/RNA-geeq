""" This module contains all functionality to filter alignments"""
import sys

def filter_alignment_sam(infile, filterset, outfile, config):
    """This procedure parses alignment lines from infile and writes passing lines to outfile"""

    counter = 0
    filter_counter = 0

    for line in infile:
        if counter % 1000000 == 0 and config.verbose:
            print >> sys.stderr, 'lines read: [ %s (taken: %s / filtered: %s)]' % (counter, counter - filter_counter, filter_counter)
        counter += 1
        if line[0] == '@':
            print >> outfile, line, 
            continue
        sl = line.strip().split()
        if len(sl) < 9:
            print >> sys.stderr, "\nERROR: Alignment line has not expected format: SAM\n%s" % line
            sys.exit(2)

