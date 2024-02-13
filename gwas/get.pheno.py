# -*- coding:utf-8 -*-
import importlib,sys
import os
import re
#sys.getdefaultencoding()
#reload(sys) 
#sys.setdefaultencoding(')
importlib.reload(sys)
#sys.setdefaultencoding('utf-8')

pheno = sys.argv[3]
pf = open(sys.argv[1],'r')
h = pf.readline().strip().split('\t')
pfdict = {}
for line in pf:
	l = line.strip().split('\t')
	k = l[h.index('SAMPLE_ID')]
	v = l[h.index(pheno)]
	pfdict[k]  = v
pf.close()

fam = open(sys.argv[2],'r')
for line in fam:
	l = line.strip().split('\t')
	if l[0] in pfdict:
		l[-1] = pfdict[l[0]]
	else:
		l[-1] = 'NA'
	r = '\t'.join(l)
	print (r)
fam.close()


