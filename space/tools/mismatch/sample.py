import pickle
import random
import sys
import pdb

from mismatch import print_substitutions
try:
	qfname=sys.argv[1] #'/fml/ag-raetsch/nobackup/projects/rgasp.tmp/reads/drosophila/CME_W1_CI.8-60-PE/quality_sample.q2'
	ffname=sys.argv[2] #'/fml/ag-raetsch/nobackup/projects/rgasp.tmp/reads/drosophila/CME_W1_CI.8-60-PE/quality_sample.fa'
	rfname=sys.argv[3] #'results.pickle'
	errorinfo=int(sys.argv[4]) # 0 for a normal fastq file, 1 for a fastq file that also contains error information
	generrors=int(sys.argv[5]) # 0 for a normal fastq file, 1 for a fastq file that also contains error information

	print >> sys.stderr, "Producing sample using quality file", qfname, "fasta file", ffname, "and mismatch file", rfname
	((qualities,substitutions), dummy)=pickle.load(file(rfname))

	if not (qfname=='-'):
		qfile=file(qfname)
	ffile=file(ffname)

	fh=ffile.readline()
	fl=ffile.readline()
	if (qfname=='-'):
		ql = 'h' * fl.len()
	else:	
		ql=qfile.readline()
except:
	print >> sys.stderr, \
"""Usage:
	sample.py quality_sample.q reads.fa mismatch.pickle

where
   quality_sample.q are read qualities
   reads.fa are reads
   mismatch.pickle is the output of mismatch.py according to which
   reads are mutated based on nucleotide composition and qualities
""" 
	sys.exit(1)

random.seed(17)

#print_substitutions(substitutions) ;
#c=dict()
#for s1 in ('A', 'C', 'G', 'T', 'N', '-'):
#	c[s1]=0.0 ;
#	for s2 in ('A', 'C', 'G', 'T', 'N', '-'):
#		try:
#			c[s1] += substitutions[(s1,s2)]
#		except:
#			pass
#print c
#p_insert = c['-']/(c['A']+c['C']+c['G']+c['T']+c['N']+c['-'])
p_insert = 0.05 # one percent of all errors are indels
		
while ql and fh and fl:
	if len(ql)!=len(fl):
		print "ql:", ql
		print "fl:", fl
		sys.exit(1)
	print '@'+fh[1:],
	read=''
	qual='' ;
	

	for i in xrange(len(ql)-1):
		
		q=ql[i]
		f=fl[i].upper()
		try:
			(a,b) = qualities[(f,q)]
		except KeyError:
			a=0
			b=1
		if generrors and random.uniform(0,1) < float(a)/float(b) and random.uniform(0,1) < p_insert:
			cum=0
			p=random.uniform(0,1)
			denom=0
			for s in ('A', 'C', 'G', 'T', 'N', '-'):
				try:
					denom+=substitutions[('-',s)]
				except KeyError:
					pass

			for s in ('A', 'C', 'G', 'T', 'N', '-'):
				try:
					cum+=substitutions[('-',s)]
				except KeyError:
					pass
				#print "p", p, "cum", cum/float(denom)
				if p < cum/float(denom):
					if s=='-':
						c=''
					else:
						if errorinfo==1:
							c='[-'+s+']'
						else:
							c=s 
					break
			read+=c
			qual+=ql[i] 

		q=ql[i]
		f=fl[i].upper()
		try:
			(a,b) = qualities[(f,q)]
		except KeyError:
			a=0
			b=1
		if generrors and random.uniform(0,1) < float(a)/float(b):
			cum=0
			p=random.uniform(0,1)
			denom=0
			for s in ('A', 'C', 'G', 'T', 'N', '-'):
				try:
					denom+=substitutions[(f,s)]
				except KeyError:
					pass

			for s in ('A', 'C', 'G', 'T', 'N', '-'):
				try:
					cum+=substitutions[(f,s)]
				except KeyError:
					pass
				#print "p", p, "cum", cum/float(denom)
				if p < cum/float(denom):
					if s=='-':
						if errorinfo==1:
							c='['+f+'-]'
						else:
							c=''
						q='' 
					else:
						if errorinfo==1:
							c='['+f+s+']'
						else:
							c=s
					break
		else:
			c=f
		
		read+=c
		qual+=q ;

	print read
	print '+'+fh[1:],
	print qual

	#ql=qfile.readline()
	fh=ffile.readline()
	fl=ffile.readline()
	if (qfname=='-'):
		ql = 'h' * fl.len()
	else:	
		ql=qfile.readline()

