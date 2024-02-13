import sys
import os

vcflist = open(sys.argv[1],'r')
outdir = sys.argv[2]
covar = sys.argv[3]
pheno = sys.argv[4]
phenotable = sys.argv[5]
info = sys.argv[6]
for line in vcflist:
	line = line.strip()
	bn = os.path.basename(line)
	prefix = bn.replace('_withrsid.bychr.vcf.gz','')
	fam = phenotable
	print ('time /software/bin/plink2 --vcf %s dosage=DS --fam %s  --covar %s --covar-variance-standardize --glm --out %s/%s.format.plink --threads 2 --memory 10000 require '%(line,fam,covar,outdir,prefix))
	if sys.argv[7] == 'logistic':
		print('time python ./bin/format3.py %s %s/%s.format.plink.PHENO1.glm.logistic.hybrid %s/%s.format.plink.PHENO1.glm.logistic.format'%(info,outdir,prefix,outdir,prefix))
	elif sys.argv[7] == 'linear':
		print('time python ./bin/format4.py %s %s/%s.format.plink.PHENO1.glm.linear %s/%s.format.plink.PHENO1.glm.linear.format'%(info,outdir,prefix,outdir,prefix))
vcflist.close()
#out1.close()
