#!/usr/bin/python3

import sys
import numpy as np
import pandas as pd

def txLength(i, df):
    return sum(map(np.int, df.loc[i, 'exonEnds'].split(',')[:-1]))\
    - sum(map(np.int, df.loc[i, 'exonStarts'].split(',')[:-1]))

def cdsLength(i, df):
    return sum(filter(lambda x: np.int(df.loc[i, 'cdsStart']) < x < np.int(df.loc[i, 'cdsEnd']), map(np.int, df.loc[i, 'exonEnds'].split(',')[:-1])))\
    - sum(filter(lambda x: np.int(df.loc[i, 'cdsStart']) < x < np.int(df.loc[i, 'cdsEnd']), map(np.int, df.loc[i, 'exonStarts'].split(',')[:-1])))\
    + np.int(df.loc[i, 'cdsEnd']) - np.int(df.loc[i, 'cdsStart'])

def firstUtrLength(i, df):
    return sum(filter(lambda x: x <= np.int(df.loc[i, 'cdsStart']), map(np.int, df.loc[i, 'exonEnds'].split(',')[:-1])))\
    - sum(filter(lambda x: x <= np.int(df.loc[i, 'cdsStart']), map(np.int, df.loc[i, 'exonStarts'].split(',')[:-1])))\
    + np.int(df.loc[i, 'cdsStart'])

def lastUtrLength(i, df):
    return sum(filter(lambda x: x >= np.int(df.loc[i, 'cdsEnd']), map(np.int, df.loc[i, 'exonEnds'].split(',')[:-1])))\
    - sum(filter(lambda x: x >= np.int(df.loc[i, 'cdsEnd']), map(np.int, df.loc[i, 'exonStarts'].split(',')[:-1])))\
    - np.int(df.loc[i, 'cdsEnd'])

def add(i, df):
    df.loc[i, 'tx_len'] = txLength(i, df)
    df.loc[i, 'CDS_len'] = cdsLength(i, df)

    if df.loc[i, 'cdsStart'] == df.loc[i, 'cdsEnd']:
        if df.loc[i, 'strand'] == '+':
            df.loc[i, 'UTR5_len'] = df.loc[i, 'tx_len']
        else:
            df.loc[i, 'UTR3_len'] = df.loc[i, 'tx_len']
    else:
        if df.loc[i, 'strand'] == '+':
            df.loc[i, 'UTR5_len'] = firstUtrLength(i, df)
            df.loc[i, 'UTR3_len'] = lastUtrLength(i, df)

        else:
            df.loc[i, 'UTR5_len'] = lastUtrLength(i, df)
            df.loc[i, 'UTR3_len'] = firstUtrLength(i, df)

def main(args):
    if (args[1] == '--help') or (args[1] == '-h'):
        print('''
        This tool will be calculate the length of CDS and UTR from genePred table of UCSC table browser.
        This was especially designed for NCBI RefSeq.
        
        Usage : python3 CDSruler.py file(input, not gzipped) file(output, not gzipped)
        
        The total length, 5′ UTR, CDS length and 3′ UTR length of each transcript will be calculated in extended columns:
        tx_len, UTR5_len, CDS_len and UTR3_len.

        Please cite my github adress "https://github.com/hyejun18/CDSruler"
        
        -v, --version\tdisplay version information and exit
        -h, --help\tdisplay this help text and exit
        ''')

    elif (args[1] == '--version') or (args[1] == '-v'):
        print('''
        CDSruler 3.1
        Copyright (C) 2021 Hyejun KIM
        
        This is free software, but please cite my github adress "https://github.com/hyejun18/CDSruler"
        ''')

    elif len(args) != 3:
        print('''
        Usage : python3 CDSruler.py file(input, not gzipped) file(output, not gzipped)
        
        Try 'python3 CDSruler.py --help or -h for more information.
        ''')
        exit()

    else:
        df = pd.read_table(str(args[1]))
        dfOut = args[2]

        df['tx_len'] = np.int(0)
        df['UTR5_len'] = np.int(0)
        df['CDS_len'] = np.int(0)
        df['UTR3_len'] = np.int(0)

        for i in range(len(df)):
            if i % 2000 == 0: 
                print(i, ' of total', len(df) + 1, ' transcripts : complete')
            add(i, df)

        df.to_csv(dfOut, sep='\t', index=None)

if __name__ == '__main__':
    main(sys.argv)
