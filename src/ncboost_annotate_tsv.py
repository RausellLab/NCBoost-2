import polars as pl
import os, sys
from tqdm import tqdm
import subprocess
import tabix
from src.ncboost_functions import get_chr_list, ncboost_query_score

input_path = sys.argv[1]
db_path = sys.argv[2]


name, extension = os.path.splitext(input_path)
output_path = f'{name}_scored.tsv'


## Stream-load and annotate variants
# This scripts receives a tsv file as input and add a new column containing NCBoost score chromosome rank percentile.

ncboost_path = f"{db_path}/ncboost_v2_hg38_20260202_light.tsv.gz"
tb = tabix.open(ncboost_path)

n_annotated = 0
n_var = 0

# Loading variants chromosome per chromosome
for l_chr in tqdm(get_chr_list(), desc='Chromosome', total=len(get_chr_list()), leave=True):
    df = (
        pl.scan_csv(input_path, separator="\t",schema_overrides={'chr':str, 'pos':int})
        .filter(pl.col("chr") == l_chr)
        .collect(streaming=True)
    )
    ncb_score_list = []
    # iterating on variants and annotate with NCBoost chr rank percentile
    for l_var in tqdm(df.iter_rows(named=True), desc='Variants', total=df.shape[0], leave=False):
        l_chr = l_var['chr']
        l_pos = l_var['pos']
        l_ref = l_var['ref']
        l_alt = l_var['alt']
        out = []
        ncb_score = ncboost_query_score(l_chr=l_chr, l_pos=l_pos, l_ref=l_ref, l_alt=l_alt, tb=tb)
        ncb_score_list.append(ncb_score)
        if ncb_score != None:
            n_annotated = n_annotated + 1
    df = df.with_columns(pl.Series('NCBoost_chr_rank_percentile', ncb_score_list))
    n_var = n_var + df.shape[0]
    if l_chr == '1':
        with open(output_path, "w") as f:
            df.write_csv(f,
            separator = "\t", 
            include_header = True
            )
    else:
        with open(output_path, "a") as f:
            df.write_csv(f,
            separator = "\t",
            include_header = False
            )

print(f'A total of {n_annotated} among {n_var} variants were succesfully annotated. Results have been saved in {output_path}')





