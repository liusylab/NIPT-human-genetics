#!/bin/sh

basevar basetype -R reference.fasta \
    --regions chr11:5246595-5248428,chr17:41197764-41276135 \
    --batch-count 50 \
    -L bamfile.list \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz \
    --nCPU 4 && echo "** job done **"
