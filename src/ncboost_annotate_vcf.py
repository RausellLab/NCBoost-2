import polars as pl
import os, sys
import vcfpy
import tabix
import subprocess
from tqdm import tqdm

input_path = sys.argv[1]
db_path = sys.argv[2]

## Stream-load and annotate variants
# 
# This scripts receives a vcf file as input and add the NCBoost score at the end of the INFO field.
# Please only input vcf files with bi-allelic variant representation (one variant per line).
# If you want to annotate a multiallelic vcf file (more than one variant per line), please use bcftools to convert it first.

name, extension = os.path.splitext(input_path)
output_path = f'{name}_scored.vcf'

def get_nvar_in_file(input_path):
    command = f"wc -l {input_path}"
    proc = subprocess.Popen([command], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell = True)
    out, _ = proc.communicate()
    l_var = int(out.decode("utf-8").split(' ')[0])
    return(l_var)

from src.ncboost_functions import ncboost_query_score

q_chr = '0'
reader = vcfpy.Reader.from_path(input_path)
reader.header.add_info_line({'ID' : 'NCBoost', 'Type' : 'Float', 'Description' : 'NCBoost hg38 score', 'Number' : 1})
writer = vcfpy.Writer.from_path(output_path, reader.header)

n_var = get_nvar_in_file(input_path)
n_annotated = 0
for record in tqdm(reader, desc='Variants', total=n_var):
    l_chr = record.CHROM
    l_pos = record.POS
    l_ref = record.REF
    l_alt = str(record.ALT[0]).replace("')", "").split(", value='")[1]

    if q_chr != l_chr:
        q_chr = l_chr
        ncboost_path = f"{db_path}/WG_chr{l_chr}.tsv.gz"
        tb = tabix.open(ncboost_path)
    else:
        tb = tb
    out = []
    ncb_score = ncboost_query_score(l_chr=l_chr, l_pos=l_pos, l_ref=l_ref, l_alt=l_alt, tb=tb)
    if ncb_score != None:
        n_annotated = n_annotated + 1
    record.INFO['NCBoost'] = ncb_score
    writer.write_record(record)

print(f'A total of {n_annotated} among {n_var} variants were succesfully annotated. Results have been saved in {output_path}')

