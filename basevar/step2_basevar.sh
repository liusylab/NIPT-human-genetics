#!/bin/sh

basevar basetype -R reference.fasta \
    -r chr11:5246595-5248428,chr17:41197764-41276135 \
    -B 200 -t 20 \
    -L bamfile.list \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz && echo "** job done **"
