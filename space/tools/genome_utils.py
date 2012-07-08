# copied from palma version of genomic.py

import os
import time
import numpy
from string import maketrans


"""A gene info structure which is a copy of the matlab version.

Written by: Regina Bohnert
"""
class GenomeInfo:
    """ Holds information on the configuration for a particular genome
    @attribute basedir: base directory
    @type basedir: string
    @attribute contig_names: name
    @type contig_names: list of strings
    @attribute flat_fnames: flat file names of the sequence
    @type flat_fnames: list of strings
    @attribute fasta_fnames: fasta file names of the sequence
    @type fasta_fnames: list of strings
    @attribute alphabet: sequence alphabet (DNA, amino acids)
    @type alphabet: string
    @attribute est_fnames: file names of the ESTs
    @type est_fnames: list of strings
    @attribute prot_fnames: file names of the protein sequence (?)
    @type prot_fnames: list of strings
    @attribute cdna_fnames: file names of the cDNAs
    @type cdna_fnames: list of strings
    @attribute annotation_fnames: file names of the annotation
    @type annotation_fnames: list of strings
    @attribute flcdna_fnames: file names of full length cDNAs
    @type flcdna_fnames: list of strings
    @attribute min_intron_len: minimal intron length
    @type min_intron_len: int
    @attribute max_intron_len: maximal intron length
    @type max_intron_len: int
    @attribute merge_est_transcripts: denotes whether EST transcripts are merged
    @type merge_est_transcripts: boolean
    """
    def __init__(self, config_file):
        est_fnames = []
        prot_fnames = []
        cdna_fnames = []
        flcdna_fnames = []
        contig_names = []
        flat_fnames = []
        fasta_fnames = []
        annotation_fnames = []
        alphabet = 'acgt'
        basedir = ""
        max_intron_len = 100000
        min_intron_len = 10
        merge_est_transcripts = 1

        try:
            fsock = open(config_file, mode="r")
        except IOError:
            print "error opening file: " + config_file
            exit(1)
        
        while 1:    
            line = fsock.readline()
            if line=="":
                break
            if not isinstance(line, str):
                break
            if line==" " or line=="\n":
                continue
            if line[0]=="#":
                continue
            if len(line)>8 and line[0:7]=="BASEDIR":
                basedir = line[8:].strip()
            elif len(line)>21 and line[0:21]=="MERGE_EST_TRANSCRIPTS":
                merge_est_transcripts = int(line[21:])
            elif len(line)>15 and line[0:15]=="MAX_DOUBLE_QGAP":
                pass
            elif len(line)>15 and line[0:15]=="MAX_DOUBLE_TGAP":
                pass
            elif len(line)>17 and line[0:17]=="MAX_CDNA_DELETION":
                pass
            elif len(line)>18 and line[0:18]=="MAX_CDNA_INSERTION":
                pass
            elif len(line)>21 and line[0:21]=="TERMINAL_EXON_END_TOL":
                pass
            elif len(line)>14 and line[0:14]=="MIN_INTRON_LEN":
                min_intron_len = int(line[14:]) 
                assert(min_intron_len>0 and min_intron_len<1000) ;
            elif len(line)>14 and line[0:14]=="MAX_INTRON_LEN":
                max_intron_len = int(line[14:]) ; 
                assert(max_intron_len>1000) ;
            elif len(line)>18 and line[0:18]=="MIN_EST_COVER_FRAC":
                pass
            elif len(line)>19 and line[0:19]=="MIN_PROT_COVER_FRAC":
                pass
            elif len(line)>19 and line[0:19]=="MIN_CDNA_COVER_FRAC":
                pass
            elif len(line)>18 and line[0:18]=="BLAT_BEST_HIT_ONLY":
                pass
            elif len(line)>20 and line[0:20]=="BLAT_BEST_HIT_MARGIN":
                pass
            elif len(line)>8 and line[0:8]=="ALPHABET":
                alphabet = line[9:].strip()
            elif len(line)>7 and line[0:7]=="CONTIGS":
                line_nr = int(line[7:])
                for i in range(0,line_nr):
                    line = fsock.readline()
                    elems = line.split()
                    assert(len(elems)==3)
                    contig_names.append(elems[0].strip())
                    flat_fnames.append(os.path.join(basedir,elems[1].strip()))
                    fasta_fnames.append(os.path.join(basedir,elems[2].strip()))
            elif len(line)>8 and line[0:8]=="ESTFILES":
                line_nr = int(line[8:])
                for i in range (0,line_nr):
                    line = fsock.readline()
                    est_fnames.append(os.path.join(basedir,line.strip())) 
            elif len(line)>9 and line[0:9]=="PROTFILES":
                line_nr = int(line[9:])
                for i in range (0,line_nr):
                    line = fsock.readline()
                    prot_fnames.append(os.path.join(basedir,line.strip()))
            elif len(line)>9 and line[0:9]=="CDNAFILES":
                line_nr = int(line[9:])
                for i in range (0,line_nr):
                    line = fsock.readline()
                    cdna_fnames.append(os.path.join(basedir,line.strip()))
            elif len(line)>11 and line[0:11]=="FLCDNAFILES":
                line_nr = int(line[11:])
                for i in range (0,line_nr):
                    line = fsock.readline()
                    flcdna_fnames.append(os.path.join(basedir,line.strip())) 
            elif len(line)>15 and line[0:15]=="ANNOTATIONFILES":
                line_nr = int(line[15:])
                for i in range (0,line_nr):
                    line = fsock.readline()
                    annotation_fnames.append(os.path.join(basedir,line.strip()))
            elif len(line)>0:
                print "cannot understand " + line + " in config file"
                
        fsock.close()
        
        self.est_fnames = est_fnames
        self.prot_fnames = prot_fnames
        self.cdna_fnames = cdna_fnames
        self.flcdna_fnames = flcdna_fnames
        self.contig_names = contig_names
        self.flat_fnames = flat_fnames
        self.fasta_fnames = fasta_fnames
        self.annotation_fnames = annotation_fnames
        self.alphabet = alphabet
        self.basedir = basedir
        self.max_intron_len = max_intron_len
        self.min_intron_len = min_intron_len
        self.merge_est_transcripts = merge_est_transcripts

    def __repr__(self):
        return 'basedir: ' + repr(self.basedir) + ' contig_names: ' + repr(self.contig_names)\
               + ' flat_fnames: ' + repr(self.flat_fnames) + ' fasta_fnames: ' + repr(self.fasta_fnames)\
               + ' alphabet: ' + repr(self.alphabet) + ' est_fnames: ' + repr(self.est_fnames)\
               + ' prot_fnames: ' + repr(self.prot_fnames) + ' cdna_fnames: ' + repr(self.cdna_fnames)\
               + ' annotation_fnames: ' + repr(self.annotation_fnames) + ' flcdna_fnames: ' + repr(self.flcdna_fnames)\
               + ' min_intron_len: ' + repr(self.min_intron_len) + ' max_intron_len: ' + repr(self.max_intron_len)\
               + ' merge_est_transcripts: ' + repr(self.merge_est_transcripts)

    def __str__(self):
        return 'basedir:\t' + str(self.basedir) + '\ncontig_names:\t' + str(self.contig_names)\
               + '\nflat_fnames:\t' + str(self.flat_fnames) + '\nfasta_fnames:\t' + str(self.fasta_fnames)\
               + '\nalphabet:\t' + str(self.alphabet) + '\nest_fnames:\t' + str(self.est_fnames)\
               + '\nprot_fnames:\t' + str(self.prot_fnames) + '\ncdna_fnames:\t' + str(self.cdna_fnames)\
               + '\nannotation_fnames:\t' + str(self.annotation_fnames) + '\nflcdna_fnames:\t' + str(self.flcdna_fnames)\
               + '\nmin_intron_len:\t' + str(self.min_intron_len) + '\nmax_intron_len:\t' + str(self.max_intron_len)\
               + '\nmerge_est_transcripts:\t' + str(self.merge_est_transcripts)



""" this function is 100% compatible to the matlab function, thus it is one based (!)
	use one_based=False if needed, then however the interval is [start,stop) (excluding stop)
"""
def load_genomic(chromosome, strand, start, stop, genome, one_based=True):
    if (type(start)==numpy.ndarray) or (type(stop)==numpy.ndarray):
        assert(len(start) == len(stop))
        assert((start[1:]-stop[:-1]>0).all())
        if strand == '+':
            idx = xrange(len(start))
        else:
            idx = xrange(len(start)-1,-1,-1)

        seq = ''.join([load_genomic(chromosome, strand, start[ix], stop[ix], genome)\
                       for ix in idx])
        return seq

    #fname = '/fml/ag-raetsch/share/databases/genomes/' + genome + '/' + chromosome + '.dna.flat'
    #fname = genome + '/' + chromosome + '.dna.flat'
    fname = genome[chromosome]
    f = open(fname, 'r')
    if one_based:
        f.seek(start-1)
        str = f.read(stop-start+1)
    else:
        f.seek(start)
        str = f.read(stop-start)
    f.close()
        
    if strand=='-':
        return reverse_complement(str)
    elif strand=='+':
        return str
    else:
        print 'strand must be + or -'
        raise KeyError

""" read a table browser ascii output file (http://genome.ucsc.edu/cgi-bin/hgTables) """
def read_table_browser(f):
    table=dict();
    for l in f.readlines():
        if not l.startswith('#'):
            (name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,proteinID,alignID)=l.split('\t')
            exonStarts=[ int(i) for i in exonStarts.split(',')[:-1] ]
            exonEnds=[ int(i) for i in exonEnds.split(',')[:-1] ]

            table[name]={ 'chrom': chrom, 'strand': strand, 'txStart': int(txStart), 'txEnd': int(txEnd), 
            'cdsStart': int(cdsStart), 'cdsEnd': int(cdsEnd), 'exonCount': int(exonCount), 'exonStarts': exonStarts,
            'exonEnds': exonEnds, 'proteinID': proteinID, 'alignID': alignID[:-1] }

    return table

""" get promoter region """
def get_promoter_region(chromosome, strand, gene_start, gene_end, genome, length):

    if strand == '+':
        return load_genomic(chromosome, strand, gene_start, gene_start+length, genome, one_based=False)
    elif strand == '-':
        return load_genomic(chromosome, strand, gene_end, gene_end+length, genome, one_based=False)
    else:
        print 'unknown strand'
        return None

""" reverse + complement a DNA sequence (only letters ACGT are translated!)
    FIXME won't work with all the rest like y... """
def reverse_complement(str):
    t=maketrans('acgtACGT','tgcaTGCA')
    return str[len(str)::-1].translate(t)

""" works only with .fa files that contain a single entry """
def read_single_fasta(fname):
    str=file(fname).read()
    str=str[str.index('\n')+1:].replace('\n','')
    return str

""" writes only single enty .fa files """
def write_single_fasta(fname, name, str, linelen=60):
    header= '>' + name + '\n'
    f=file(fname,'a')
    f.write(header)
    for i in xrange(0,len(str),linelen):
        f.write(str[i:i+linelen]+'\n')
    f.close()

""" read fasta as dictionary """
def read_fasta(f, fasta=dict(), id_only=True):
    import numpy
    str=f.read()
    idx = str.find('>')
    #if idx==-1:
    if 0: # no support for contig list files -> would need extra blastn workaround
        # file name list?
        files = str.split('\n') ;
        for s in files:
            if len(s)==0: continue
            print s.strip() + '\n'
            fasta = read_fasta(file(s.strip()), fasta) 
    else:
        # real fasta file?
        sequences = str.split('>') 
        for s in sequences:
            if len(s)==0: continue
            header = s[0:s.index('\n')]
            if id_only:
                header = header.split(' ')[0] ;
                header = header.split('\t')[0] ;
            seq = s[s.index('\n')+1:].replace('\n','').upper()
            #print 'has_key', fasta.has_key(header),header
            assert(not fasta.has_key(header))
            fasta[header]=seq ;

    return fasta

""" write dictionary fasta """
def write_fasta(f, d, linelen=60):
    for k in sorted(d):
        f.write('>%s\n' % k);
        s = d[k]
        for i in xrange(0, len(s), linelen):
            f.write(s[i:i+linelen] + '\n')

def write_gff(f, (source, version), (seqtype, seqname), descrlist, skipheader=False):
    """ writes a gff version 2 file
        descrlist is a list of dictionaries, each of which contain these fields:
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
    """

    if not skipheader:
        f.write('##gff-version 2\n')
        f.write('##source-version %s %s\n' % (source, version) )

        t=time.localtime()
        f.write("##date %d-%d-%d %d:%d:%d\n" % t[0:6])

    f.write('##Type %s %s\n' % (seqtype, seqname) )

    for d in descrlist:
        f.write('%s\t%s\t%s\t%d\t%d\t%f\t%s\t%d' % (d['seqname'], d['source'], 
                                                                                d['feature'], d['start'], d['end'], 
                                                                                d['score'], d['strand'], d['frame']))
        if d.has_key('attributes'):
            f.write('\t' + d['attributes'])
            if d.has_key('comments'):
                f.write('\t' + d['comments'])
        f.write('\n')


if __name__ == '__main__':
    import sys,os

    table=read_table_browser(file('/fml/ag-raetsch/home/sonne/addnet/tfbs/share/data/wt1_bibliosphere_table_browser_hg17.txt'))
    print table.keys()
    print table[table.keys()[0]]
    d = { 'ahoernchen' : 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',
              'bhoernchen' : 'GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA' }

    write_fasta(sys.stdout, d)
    write_fasta(file('/tmp/test.fa','w'), d)

    d2 = read_fasta(file('/tmp/test.fa'))
    os.unlink('/tmp/test.fa')

    print d
    print d2
    print d == d2

    p=load_genomic('chr5', '+', 100000, 100100,'hg17')
    n=load_genomic('chr1', '-', 3000000, 3001000,'mm7')
    write_single_fasta('bla.fa','bla', 'ACGT')
    n2=read_single_fasta('bla.fa')
