
# coding: utf-8

# In[12]:

import pandas as pd
import yaml
import numpy as np
from collections import OrderedDict


# In[13]:

db = pd.read_excel('../../../DB/DBparticles.xlsx', sheetname=None, skiprows=list(range(9)))


# In[14]:

allowed_columns = ['Particle', 'Mass, kg', 'Diameter, m', 'Parameter ε (Lennard-Jones), J',
                   'Dissociation energy, J', 'Formation energy, J', 'Ionization potential, J',
                   'No. electronic  level', 'Statistical weight', 'Electronic energy, J',
                   'Frequency of vibrations (we), m^-1', 'wexe, m^-1',
                   'weye, m^-1', 'weze, m^-1',
                   'Be, m^-1', 'ae, m^-1', 'Moment of Inertia, J*s^2', 'Factor of symmetry',
                   'Parker (zeta^infty)', 'Charge', 'Reduced oscillator mass',
                   'internuclear distance, r_e , m', 'mA/mAB', 'mB/mAB']


# In[15]:

integer_columns = ['Statistical weight', 'No. electronic  level', 'Factor of symmetry', 'Charge',
                   'vibrational level']


# In[16]:

def convert_particles_db_to_dict(df):
    particles = {}
    for particle in df:
        res = []
        for column in db[particle]:
            if column in allowed_columns:
                tmp = db[particle][column].dropna().values
                if column in integer_columns:
                    tmp = np.asarray(tmp, dtype=int)
                if tmp.shape[0] == 1:
                    res.append([column, tmp.tolist()])
                else:
                    res.append([column, tmp.tolist()])
        particles[particle] = res
    return particles


# In[17]:

r = convert_particles_db_to_dict(db)


# In[18]:

with open('../../../DB/particles.yaml', 'w') as f:
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



