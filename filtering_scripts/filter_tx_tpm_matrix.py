import pandas as pd
tpm = pd.read_csv('matrices/tx_tpm.csv', index_col=0)
keep = tpm[(tpm >= 1).sum(axis=1) >= 3].index           # boolean mask
keep.to_series().to_csv('pass_tx.txt', index=False)
print(len(keep), 'transcripts kept')
