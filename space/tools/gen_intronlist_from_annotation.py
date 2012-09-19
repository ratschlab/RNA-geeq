
import GFFParser
import cPickle
import sys
import pdb
from Bio.SeqFeature import SeqFeature

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--annotation', dest='anno', metavar='FILE', help='annotation file in gff3 format', default='-')
    required.add_option('-o', '--output', dest='outfile', metavar='FILE', help='annotation intron list', default='-')
    parser.add_option_group(required)

    (options, args) = parser.parse_args()
    
    if len(argv) < 3 or '-' in [options.anno, options.outfile] :
        parser.print_help()
        sys.exit(2)

    return (options, args)


def fix_structure(trans_dict):
    """Adapts the structure of the trans_dict to different gff-file structures"""

    for idx in range(len(trans_dict.features)):
        old_trans_dict = trans_dict.features[idx]
        if trans_dict.features[idx].type == 'exon':
            exon = trans_dict.features[idx]
            trans_dict.features[idx] = SeqFeature(exon.location, type = 'gene', strand = exon.strand, id = exon.id)
            trans_dict.features[idx].sub_features = [SeqFeature(exon.location, type = 'Transcript', strand = exon.strand, id = exon.id)] 
            trans_dict.features[idx].sub_features[0].sub_features = [exon]
        elif len(trans_dict.features[idx].sub_features) > 0 and trans_dict.features[idx].sub_features[0].type == 'exon':
            exon = trans_dict.features[idx]
            trans_dict.features[idx] = SeqFeature(exon.location, type = 'gene', strand = exon.strand, id = exon.id)
            trans_dict.features[idx].sub_features = [exon]

def main():
    """Main function ... """

    (options, args) = parse_options(sys.argv)

    iterator = GFFParser.GFFAddingIterator()    
    examiner = GFFParser.GFFExaminer()

    exon_map = dict()

    id_dict = examiner.available_limits(options.anno)['gff_id']
    intron_lists = dict()
    ### collect all available sources from gff-file
    source_dict = examiner.available_limits(options.anno)['gff_source_type']
    taken_sources = dict()
    ### parse only for exons and let the GFFparser 
    ### infer the respecting parents (otherwise doubled entries occured)
    ### we sanitize the structure later on anyways
    types = ['exon']
    for key in source_dict.keys():
        if key[1] in set(types):
            if taken_sources.has_key(key[0]):
                taken_sources[key[0]] += 1
            else:
                taken_sources[key[0]] = 1
    
    ### print taken_sources
    if len(taken_sources.keys()) == 0:
        print >> sys.stderr, 'No suitable sources found!'
        sys.exit(-1)
    else:
        source_strings = taken_sources.keys()
        print "take sources %s" % source_strings

    ### build up gff-parsing filter
    gff_sources = []
    for source in source_strings:
        gff_sources.extend(zip([source] * len(types), types))

    ### parse gff-file
    for idx in id_dict.keys():
        print 'parsing chromosome %s' % idx
        if len(gff_sources) > 0:
            trans_dict = iterator.get_all_features(options.anno, {'gff_source_type':gff_sources, 'gff_id':idx})
        else:
            trans_dict = iterator.get_all_features(options.anno, {'gff_id':idx})
        ### since we parse only one chromosome, this loop is evaluated only once
        for chrm in trans_dict.keys():
            ### verify/sanitize the created dictionairy
            fix_structure(trans_dict[chrm])
            intron_lists[chrm] = dict()
            for gene in trans_dict[chrm].features:
                for trans in gene.sub_features:
                    if trans.type == 'exon':
                        print "WARNING: Exon on transcript level:"
                        print trans
                        print 'will continue\n'
                        continue
                    elif len(trans.sub_features) > 1: ### at least two exons for one intron ...
                        strand = trans.sub_features[0].strand
                        contig_list = [(trans.sub_features[i].location.nofuzzy_start, trans.sub_features[i].location.nofuzzy_end) for i in range(len(trans.sub_features))]
                        contig_list.sort(lambda u, v:u[0]-v[0])
                        for exon in range(len(contig_list) - 1):
                            ### update intron lists
                            if contig_list[exon][1] - contig_list[exon + 1][0] == 0:
                                continue
                            assert(contig_list[exon][1] < contig_list[exon + 1][0])
                            ### for now strand information is only dummy
                            intron_lists[chrm][(0, contig_list[exon][1], contig_list[exon + 1][0])] = strand
                     
                        ### update exon map
                        for exon in range(len(contig_list)):
                            if not exon_map.has_key(chrm):
                                exon_map[chrm] = dict()

                            if not exon_map[chrm].has_key(trans.id):
                                exon_map[chrm][trans.id] = dict()
                            ### we assume, that an exon cannot occurr twice in the same transcript!
                            ### the value in the dict is a binary encoding, if the left/right end is intronic 10 = 2 means, 5' end is intronic
                            if len(contig_list) == 1:
                                exon_map[chrm][trans.id][contig_list[exon]] = 0 ### 00 -> should never occurr
                            elif exon == 0:
                                exon_map[chrm][trans.id][contig_list[exon]] = 2 ### 10
                            elif exon == len(contig_list) - 1:
                                exon_map[chrm][trans.id][contig_list[exon]] = 1 ### 01
                            else:
                                exon_map[chrm][trans.id][contig_list[exon]] = 3 ### 11 

    outfile = open(options.outfile, 'w')
    cPickle.dump(intron_lists, outfile)
    outfile.close()
    
    outfile = open(options.outfile + '.' + 'cov', 'w')
    cPickle.dump(exon_map, outfile)
    outfile.close()

if __name__ == '__main__':
    main()
