import sys
import glob
from os import path, chdir
import pandas as pd

def _GetTetradGeometry(quadruplex_layers):
    columns = ['Q3', 'Q2', 'Q1', 'Q4']
    rows    = ['top', 'mid', 'bot']
    df_tetrad_geometry = pd.DataFrame(
        quadruplex_layers,
        columns=columns,
        index=rows,
    )
    return df_tetrad_geometry

'''
type_na: g4dna, dsdna
na     : dsdna, propeller, basket, chair, hybrid-i, hybrid-ii
drug   : poh, ewv, ''
'''

## Main
def load_config(confFile=None, verbose=1):
    conf = {}
    if confFile is None:
        confFile = glob.glob('*.conf')[0]
    if verbose: print(f'Loading Config From "{confFile}" ...')
    f = open(confFile, 'r')
    readlines = f.read()
    readlines = readlines.replace('\\\n', ',').split('\n')
    for readline in readlines:
        readlineSplit = readline.split()
        if readline.startswith('#'):
            if verbose: print(readline)
            continue
        elif '=' in readlineSplit:
            if '$' in readline: ## $USER
                conf[readlineSplit[0]] = ' '.join(readlineSplit[2:])
                if verbose: print(f'{readlineSplit[0]} = {conf[readlineSplit[0]]}')               
            elif '/' in readline: ## ./bigtraj_fluctmatch
                filename = path.abspath(' '.join(readlineSplit[2:]))
                filename = filename.replace('/data/hpcUser/chu_02', '/home/chu_02/user') ## symbolic link in crouns
                conf[readlineSplit[0]] = filename
                if verbose: print(f'{readlineSplit[0]} = "{conf[readlineSplit[0]]}"')
            elif 'sequence' in readline:
                conf[readlineSplit[0]] = readlineSplit[2:]
                conf['strandid2sequence'] = {f'STRAND{i}':list(' '+d) for i, d in enumerate(readlineSplit[2:], start=1)}
                if verbose: print(f'{readlineSplit[0]} = {conf[readlineSplit[0]]}')
            elif 'GBAc' in readline:
                gbac = readlineSplit[2:]
                gbac = list(''.join(gbac))
                gbac = [{'A':'anti', 'S':'syn'}[g] for g in gbac]
                conf[readlineSplit[0]] = {resid:g for resid, g in enumerate(gbac, start=1)}
                if verbose: print(f'{readlineSplit[0]} = {conf[readlineSplit[0]]}')
            elif 'quadruplex_layers' in readline:
                quadruplex_layers = ' '.join(readlineSplit[3:]).split(',')[:-1]
                quadruplex_layers = [q.split() for q in quadruplex_layers]
                quadruplex_layers = [[int(q) for q in qu] for qu in quadruplex_layers]
                conf[readlineSplit[0]] = quadruplex_layers
                conf['df_tetrad_geometry'] = _GetTetradGeometry(quadruplex_layers)
                if verbose: print(f"df_tetrad_geometry = {conf['df_tetrad_geometry'].values.tolist()}")
            else:
                data = ' '.join(readlineSplit[2:])
                conf[readlineSplit[0]] = data.split('#')[0] ## exclude comment
                if verbose: print(f'{readlineSplit[0]} = {conf[readlineSplit[0]]}')
                try: ## 5,0
                    if '.' in conf[readlineSplit[0]]:
                        conf[readlineSplit[0]] = float(conf[readlineSplit[0]])
                    else:
                        conf[readlineSplit[0]] = int(conf[readlineSplit[0]]) 
                except:
                    pass
    f.close()
    if 'drug' in conf and 'na' in conf:
        conf['system'] = conf['na']+'_'+conf['drug'] if conf['drug'] else conf['na']
    return conf


class ConfAgent():
    def __init__(self):
        conf = load_config(verbose=0)
        for var, value in conf.items():
            setattr(self, var, value)