import pandas as pd 
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="""PhaTYP is a python library for bacteriophages' lifestyles prediction. 
                                 PhaTYP is a BERT-based model and rely on protein-based vocabulary to convert DNA sequences into sentences for prediction.""")
parser.add_argument('--midfolder', help='folder to store the intermediate files', type=str, default='phatyp/')
parser.add_argument('--out', help='folder to store the intermediate files', type=str, default='train.txt')
parser.add_argument('--mode', help='pretrain or finetune', type=str, default='pretrain')
inputs = parser.parse_args()



feat_df = pd.read_csv(f'{inputs.midfolder}/bert_feat.csv')
feat_df = feat_df.drop(columns = ['label'])

if inputs.mode == 'finetune':
    feat_df.to_csv(inputs.out, index=False)
else:
    feat_df.to_csv(inputs.out, index=False, header=False)
