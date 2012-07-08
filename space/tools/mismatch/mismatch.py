import re
import sys
import pdb
import pickle

def compute_qualities(qualities, read, quality):
    #determine read qualities
    i = 0
    bracket = False
    gotbase = False
    base = None
    for c in read:
        if c == '[':
            bracket = True
            continue

        if bracket:
            if c == ']':
                q = (base, quality[i])
                if not qualities.has_key(q):
                    qualities[q] = [0, 0]
                qualities[q][0] += 1
                qualities[q][1] += 1
                bracket = False
                gotbase = False
                base = None
                i += 1
            elif c in ('A', 'C', 'G', 'T', 'N') and not gotbase:
                base = c
                gotbase = True
        else:
            q = (c, quality[i])
            if not qualities.has_key(q):
                qualities[q] = [0, 0]
            qualities[q][1] += 1
            i += 1

def compute_substitutions(substitutions, read):
    idx = read.find('[')
    while idx != -1:
        read = read[idx+1:]
        p = (read[0], read[1])
        if not substitutions.has_key(p):
            substitutions[p] = 0
        substitutions[p] += 1
        idx = read.find('[')

def process_line(qualities, substitutions, l, trunc):
    try:
        read = l.split('read=')[1].split(';')[0]
        quality = l.split('quality=')[1].split()[0]
       
        #### delete deletions, since they have no associated quality
        #while re.search('\[.-\]', read) != None:
        #    st = re.search('\[.-\]', read)
        #    read = read.replace[read[start, start + 4], '']
        #
        #assert len(read) - (3 * len(re.findall('\[', read))) == len(quality)
            ## -> indels skipped anyway for quality sampling

        ### truncate
        if trunc > 0:
            trunc_front = trunc
            trunc_back = trunc
            while len(read[:trunc_front]) - (3 * read[:trunc_front].count('\[')) < trunc:
                trunc_front += 1
            while len(read[-trunc_back:]) - (3 * read[-trunc_back:].count('\[')) < trunc:
                trunc_back += 1
            
            if len(read) < (trunc_front + trunc_back):
                print >> sys.stderr, 'truncation length (%d) MUST NOT exceed half read length (%d)!' % (trunc, (len(read) - (3 * len(re.findall('\[', read)))) / 2)
                
            read = read[trunc_front:-trunc_back]
            quality = quality[trunc:-trunc]
    except:
        print >> sys.stderr
        print >> sys.stderr, l
        #print len(read)
        #print len(quality)
        return #sys.exit(1)

    try:
        compute_substitutions(substitutions, read)
        if read.find('-') == -1:
            compute_qualities(qualities, read, quality)
    except:
        pass

def print_qualities(qualities):
    acgt = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    for q in qualities.iterkeys():
        print acgt[q[0]], ord(q[1]), qualities[q][0], qualities[q][1]

def print_substitutions(substitutions):
    acgt = ('A', 'C', 'G', 'T', 'N', '-')
    for i in acgt:
        for j in acgt:
            p = (i, j)
            r = 0;
            try:
                r = substitutions[p]
            except KeyError:
                pass
            print r,
        print '\n',

def get_quality_matrix(qualities):
    import numpy
    x = numpy.zeros((5,127))
    acgt = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    for q in qualities.iterkeys():
        x[acgt[q[0]], ord(q[1])] = float(qualities[q][0]) / float(qualities[q][1] + 100)

    start =- 1
    stop =- 1
    for i in xrange(127):
        if start == -1 and any(x[:,i]):
            start = i

        if stop == -1 and any(x[:,127-i-1]):
            stop = 127 - i - 1
    
    print start, stop
    x = x[:,start:stop]
    return x

if __name__ == "__main__":
    regex = re.compile('.*read=.*read=.*')
    #regex = re.compile('.*(read=.*read=|read=.*[N-].*;).*')
    qualities = dict()
    substitutions = dict()
    linecount = 0
    processed = 0

    try:
        src = sys.argv[1]
        dst = sys.argv[2]
        trunc = 0
        if len(sys.argv) > 3:
            trunc = int(sys.argv[3])
          
        #src='/fml/ag-raetsch/nobackup/projects/rgasp.tmp/reads/drosophila/CME_W1_CI.8-60-PE/read_sample.result10'
        #src='/fml/ag-raetsch/nobackup/projects/rgasp.tmp/reads/drosophila/CME_W1_CI.8-60-PE/read_sample.result'
        #src='/tmp/y'
        #pdb.set_trace()
        print >> sys.stdout, "processing", src, "writing outputs to", dst
        srcfile = file(src)
        dstfile = file(dst, 'w')
    except:
        print >> sys.stderr, \
    """Usage:
        mismatch.py read_sample.result mismatch.pickle

    where
       read_sample.result are matches of the read sample
       mismatch.pickle is the output (warning file will be overwritten)
       containing nucleotide and quality mismatch frequencies
    """ 
        sys.exit(1)

    for l in srcfile.xreadlines():
        linecount += 1
        if linecount % 100000 == 0:
            print >> sys.stderr, linecount,

        if l.startswith('[map_reads]'):
            break
        if regex.match(l):
            continue #skip over multiple matches, -, N
        process_line(qualities, substitutions, l, trunc)
        processed += 1

    print >> sys.stdout
    print >> sys.stdout, processed, "lines of", linecount, "lines processed"
    print >> sys.stdout
    #
    #print_qualities(qualities)
    #
    #print
    #print
    #
    #print_substitutions(substitutions)
    pickle.dump(((qualities,substitutions), (processed, linecount)), dstfile)
