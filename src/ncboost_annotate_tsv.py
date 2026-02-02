import polars as pl
import os, sys
from tqdm import tqdm

input_path = sys.argv[1]
db_path = sys.argv[2]


name, extension = os.path.splitext(input_path)
output_path = f'{name}_scored.tsv'


## Stream-load and annotate variants
# This scripts receives a tsv file as input and add a new column containing NCBoost screo rank percentile per chromosome.

def get_nvar_in_file(input_path):
    command = f"wc -l {input_path}"
    proc = subprocess.Popen([command], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell = True)
    out, _ = proc.communicate()
    l_var = int(out.decode("utf-8").split(' ')[0])
    return(l_var)

from src.ncboost_functions import get_chr_list, ncboost_query_score
import tabix

q_chr = 0
n_annotated = 0
n_var = 0
for l_chr in tqdm(get_chr_list(), desc='Chromosome', total=len(get_chr_list()), leave=True):
    if q_chr != l_chr:
        q_chr = l_chr
        ncboost_path = f"{db_path}/WG_chr{l_chr}.tsv.gz"
        tb = tabix.open(ncboost_path)
    else:
        tb = tb

    df = (
        pl.scan_csv(input_path, separator="\t",schema_overrides={'chr':str, 'pos':int})
        .filter(pl.col("chr") == l_chr)
        .collect(streaming=True)
    )

    ncb_score_list = []

    for l_var in tqdm(df.iter_rows(named=True), desc='Variants', total=df.shape[0], leave=False):
        l_chr = l_var['chr']
        l_pos = l_var['pos']
        l_ref = l_var['ref']
        l_alt = l_var['alt']

        if q_chr != l_chr:
            q_chr = l_chr
            ncboost_path = f"{db_path}/WG_chr{l_chr}.tsv.gz"
            tb = tabix.open(ncboost_path)
        else:
            tb = tb
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





