#!/usr/bin/env python3
# make_tx_tpm_matrix.py
import glob, os, re, gzip, pandas as pd, sys, io

gtf_dir   = 'expr_gtf'              # adjust if needed
out_file  = 'matrices/tx_tpm.csv'   # will create matrices/ automatically

os.makedirs(os.path.dirname(out_file), exist_ok=True)
rx_tid = re.compile(r'transcript_id "([^"]+)"')
rx_tpm = re.compile(r'TPM "([^"]+)"')
frames = []

for fp in sorted(glob.glob(f'{gtf_dir}/*.gtf*')):
    sample = os.path.basename(fp).split('.')[0]          # A1_001_sorted → A1_001_sorted
    openf  = gzip.open if fp.endswith('.gz') else open
    rows = []
    with openf(fp, 'rt') as fh:
        for line in fh:
            if '\ttranscript\t' not in line:         # skip exon etc.
                continue
            tid = rx_tid.search(line).group(1)
            tpm = float(rx_tpm.search(line).group(1))
            rows.append((tid, tpm))
    df = pd.DataFrame(rows, columns=['transcript_id', sample]).set_index('transcript_id')
    frames.append(df)

tpm = pd.concat(frames, axis=1).fillna(0.0)
tpm.to_csv(out_file)
print(f'{tpm.shape[0]:,} transcripts × {tpm.shape[1]} samples ➜ {out_file}')

