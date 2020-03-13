
# coding: utf-8

# In[1]:

import pandas as pd
import yaml
import numpy as np
from collections import OrderedDict
import re


# In[4]:

def isEnglish(s):
    try:
        s.encode(encoding='utf-8').decode('ascii')
    except UnicodeDecodeError:
        return False
    else:
        return True

def check_dimension(s):
    # if the data contains something like
    """
    if the data contains something like ",nm," ", nm," ",nm ," " nm ", we remove the columns
    """
    s_l = s.lower()
    bad_dimensions_base = ['nm', 'ev', 'angstrom', 'td']  # nanometers, electron volts, angstroms, temperature of diss.
    bad_dimensions = []
    bad_dimensions_end = []
    for bad_dim in bad_dimensions_base:
        bad_dimensions.append(',' + bad_dim + ',')
        bad_dimensions.append(', ' + bad_dim + ',')
        bad_dimensions.append(',' + bad_dim + ' ,')
        bad_dimensions.append(' ' + bad_dim + ' ')
        bad_dimensions_end.append(',' + bad_dim)
        bad_dimensions_end.append(', ' + bad_dim)
        bad_dimensions_end.append(' ' + bad_dim)
    for d in bad_dimensions:
        if d in s_l:
            return False
    for d in bad_dimensions_end:
        if s_l.endswith(d):
            return False
    return True  # everything is OK

def do_not_touch(s):
    do_not_touch_cols_bzowski = ['Bzowski, rho* [6]', 'Bzowski, V_0* [6]']
    if s in do_not_touch_cols_bzowski:
        return (True, s[:-4])
    else:
        return (False, s)
    
def get_names(s):
    ionized = s.count('+')
    if ionized == 1:
        return (s.split('+')[0], s.split('+')[1])
    if ionized == 2:
        if s.endswith('+'):
            return (s.split('+')[0], s.split('+')[1]+'+')
        else:
            return (s.split('+')[0]+'+', s.split('+')[2])
    else:
        return (s.split('+')[0]+'+', s.split('+')[2]+'+')


# In[65]:

drop_cols = ['diss,O2,Ea,mTM, K [7]', 'diss,N2,Ea,mTM, K [7]']

db = pd.read_excel('../../../DB/DBinteraction.xlsx', sheetname=None)
db.pop('Interactions', None)
for k in db:
    db[k].dropna(axis=1, inplace=True)
print(len(db), 'interactions in DB')
for k in db:
    for c in db[k].columns:
        dnt = do_not_touch(c)
        if dnt[0]:
            db[k].rename(columns={c: dnt[1]}, inplace=True)
        else:
            if not isEnglish(c):
                db[k].drop(c, axis=1, inplace=True)
                dnt = True
            else:
                dnt = False
            if not dnt:
                if not check_dimension(c):
                    db[k].drop(c, axis=1, inplace=True)
                    dnt = True
            if not dnt:
                if c in drop_cols:
                    db[k].drop(c, axis=1, inplace=True)
                    
endings_leave = ['Park', 'Scanlon', 'mTM']  # so if a column ends with ',Park [4]', it will be renamed to ',Park'
                    
for k in db:
    cols = []
    for c in db[k].columns:
        cnew = c
        for e_l in endings_leave:
            cnew = re.sub(r',' + e_l + '\s*\[\d\]$', ',' + e_l, cnew)
        cnew = re.sub(r',[^,]*\[\d\]$', '', cnew)
        cnew = re.sub(r'\s*\[\d\]$', '', cnew)
        cnew = cnew.strip()
        db[k].rename(columns={c: cnew}, inplace=True)


# In[66]:

def convert_interactions_db_to_dict():
    interactions = {}
    for interaction in db:
        res = []
        for column in db[interaction]:
            tmp = db[interaction][column].values
            res.append([column, tmp.tolist()])
        interactions[get_names(interaction)[0] + ' + ' + get_names(interaction)[1]] = res
    return interactions


# In[67]:

r = convert_interactions_db_to_dict()


# In[69]:

with open('../../../DB/interactions.yaml', 'w') as f:
    for particle in r:
        f.write(particle + ':\n')
        for el in r[particle]:
            if el[0] != 'Particle':
                if len(el[1]) != 0:
                    f.write('  ' + el[0] + ': ')
                    if len(el[1]) == 1:
                        f.write(str(el[1][0]) + '\n')
                    else:
                        f.write('[')
                        for item in el[1][:-1]:
                            f.write(str(item) + ', ')
                        f.write(str(el[1][-1]))
                        f.write(']\n')


# In[ ]:



