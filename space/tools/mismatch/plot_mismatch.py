import mismatch
import pickle
import sys
from pylab import plot, legend, title, savefig

rfname=sys.argv[1]
dfname=sys.argv[2]

((qualities,substitutions), (processed, linecount))=pickle.load(file(rfname))

x=mismatch.get_quality_matrix(qualities)
plot(x.T)
#semilogy(x.T)
legend(('A', 'C', 'G', 'T', 'N'))
title(rfname + " proc'd: " + str(processed))
#show()
savefig(dfname)
