# VCF_annotation

For this challenge, you are asked to prototype a variant annotation tool. We will provide you with a VCF file, and you will create a small software program to annotate each variant in the file.
Each variant must be annotated with the following pieces of information:
1. Type of variation (substitution, insertion, CNV, etc.) and their effect (missense, silent,
intergenic, etc.). If there are multiple effects, annotate with the most deleterious
possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from ExAC API (API documentation is available here:
http://exac.hms.harvard.edu/).
6. Any additional annotations that you feel might be relevant.


# Usage

python anno.py -i <yourvcffile.vcf>

# Output files

vcf_output.csv
1. Chrom: Chromosome
2. Pos: Position
3. Ref: Reference allele
4. Alt: alternate allele
5. Variant_Type:Variant type extract from INFO['TYPE']
6. Variant_Effect: Variant effect from ExAC 'major_consequence'
7. DP: Depth of sequence coverage at the site of variation
8. AO: Depth of sequence coverage at the site of variation
9. Ratio: Percentage of reads supporting the variant versus those supporting reference read, if -1 means RO equals to zero
10. AF_ExAC:Allele frequency of variant from ExAC API
