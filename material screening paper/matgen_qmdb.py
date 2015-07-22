
# coding: utf-8

# ## Matgen based qmdb & pymatgen
# 
# The Following Gives the code to assimilate and clean data for the Wolverton Paper. 

# Import All libraries

# In[12]:

import requests
import json
import numpy as np
import itertools
import pandas as pd
from __future__ import division

from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.matproj.rest import MPRester 
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition


API_key = 'NrBIvm9wt7Hq1fSD' # Get API key from https://materialsproject.org/open
mp = MPRester(API_key)


elements = pd.read_table('elements.txt', delimiter="\t") # elements data
elements.index = elements.symbol # Convert the index as the elements symbols 
formation_energies = pd.read_table('formation_energies.txt', delimiter="\t")    # formation energy data 
 
# Create a List of proportions for the binary system

prop = np.arange(0.05,1,0.05);
p = np.around(prop,decimals=2);
prop = prop.tolist()


# In[26]:




# In[13]:

# List of a Binary systems possible in the periodic table

binary_sys = []
except_elements = {'Am','At','Bh','Bk','Cf','Cm','Cn','Db','Ds',
                   'Es','Fm','Fr','Hs','Lr','Md','Mt','No','Po',
                   'Ra','Rf','Rg','Rn','Sg'}
element_list = list(set(elements.symbol) - except_elements)
for subset in itertools.combinations(element_list,2):
    binary_sys.append(subset[0]+'-'+subset[1])
binary_sys = sorted (binary_sys)    
len(binary_sys)  #3916 systems; 89 C 2 combination of 89 known elements


# In[14]:

'''
The Following Function takes in a binary system Eg: 'Li-Fe' in the 
same format as shown and returns a list of dictionaries; each giving 
specific properties of the system and various proportions defined by the 
list 'prop'. The returned variable is a List of Dicts.
'''
def material_load_binary(d, sep='-', p = prop):
    return_data = []
     
    d = d.split(sep)
    
    # Create a phase diagram object for the following system:
    entry = mp.get_entries_in_chemsys([d[0],d[1]]) # gets the entries of the chemical system
    pd = PhaseDiagram(entry) # creates a phasediagram object
    pd_analyse = PDAnalyzer(pd) # creates a phase Diagram analysis object

    # Get the features for various proportions Using the get_hull_energy method;
    # (Need to add documentation)
    for i in range(0,len(p)):
        temp_data = {}
        prop_a = p[i]
        prop_b = p[-(i+1)]
        try :
            temp_data['system'] = d[0]+'-'+d[1]
            temp_data['A'] = d[0]
            temp_data['B'] = d[1]
            temp_data[d[0]+'_prop'] = prop_a
            temp_data[d[1]+'_prop'] = prop_b
            temp_data['formation_energy'] = pd_analyse.get_hull_energy(Composition.from_dict({d[0]: prop_a, d[1] : prop_b}))

            # Element Property extraction

            temp_data['avg_atomic_mass'] = prop_a*elements.loc[d[0]].mass + prop_b*elements.loc[d[1]].mass
            temp_data['avg_row'] = prop_a*elements.loc[d[0]].period + prop_b*elements.loc[d[1]].period
            temp_data['avg_col'] = prop_a*elements.loc[d[0]].group + prop_b*elements.loc[d[1]].group
            temp_data['max_z_diff'] = abs (elements.loc[d[0]].z - elements.loc[d[1]].z) # Max Difference in atomic number
            temp_data['avg_z'] = prop_a*elements.loc[d[0]].z + prop_b*elements.loc[d[1]].z
            temp_data['max_radius_diff'] = abs (elements.loc[d[0]].atomic_radii - elements.loc[d[1]].atomic_radii) # Max Difference in atomic radius
            temp_data['avg_radius'] = prop_a*elements.loc[d[0]].atomic_radii + prop_b*elements.loc[d[1]].atomic_radii
            temp_data['max_en_diff'] = abs (elements.loc[d[0]].electronegativity - elements.loc[d[1]].electronegativity) # Max Difference in electronegativity
            temp_data['avg_en'] = prop_a*elements.loc[d[0]].electronegativity + prop_b*elements.loc[d[1]].electronegativity # Avg Difference in electronegativity
            temp_data['avg_s_elec'] = prop_a*elements.loc[d[0]].s_elec +prop_b* elements.loc[d[1]].s_elec
            temp_data['avg_p_elec'] = prop_a*elements.loc[d[0]].p_elec +prop_b* elements.loc[d[1]].p_elec
            temp_data['avg_d_elec'] = prop_a*elements.loc[d[0]].d_elec +prop_b* elements.loc[d[1]].d_elec
            temp_data['avg_f_elec'] = prop_a*elements.loc[d[0]].f_elec +prop_b* elements.loc[d[1]].f_elec
            
            temp_sum = temp_data['avg_s_elec']+temp_data['avg_p_elec']+temp_data['avg_d_elec']+temp_data['avg_f_elec']
            
            temp_data['prop_s_elec'] = temp_data['avg_s_elec']/temp_sum
            temp_data['prop_p_elec'] = temp_data['avg_p_elec']/temp_sum
            temp_data['prop_d_elec'] = temp_data['avg_d_elec']/temp_sum
            temp_data['prop_f_elec'] = temp_data['avg_f_elec']/temp_sum
            
            
            return_data.append(temp_data)
        except :
            pass
    return return_data,temp_data['system']


# In[ ]:

'''
The Following Create a 'master' dataframe (Called so as later the 
dataframe of ternary systems will be appended to this ).  
'''
master = pd.DataFrame()

for i,item in enumerate(binary_sys, 1):
    # item = item.encode('ascii')
    t = material_load_binary(item)
    temp = pd.DataFrame(t[0])
    master = master.append(temp)
    print t[1]

col_names = list(master.filter(regex='_prop').columns)
tern_data[col_names].fillna(0) # replace all NaN's with 0 in the '_prop' columns


master.head()  


# In[ ]:

master.save()
# Check For ternary systems
# The following check for all Stable (stability <=0) and ternary 
# systems in the formation_energies table

is_ternary = formation_energies['composition_id'].str.split().apply(len) == 3
is_stable = formation_energies['stability'] <= 0

ternary_systems = formation_energies[is_stable & is_ternary]
ternary_systems = ternary_systems[['composition_id','stability','delta_e']]

ternary_systems.head()


# In[ ]:

# the Function splits the above composition_id into 
# a list of 3 tuples.  Eg:
#     split_systems('Li1 N1 Zn1')
# >>> [('Li', '1'), ('N', '1'), ('Zn', '1')]

import re
def split_systems(sys):
    ret_list = []
    s = sys.split() 
    r = re.compile("([a-zA-Z]+)([0-9]+)")
    return [r.match(string).groups() for string in s]


# In[ ]:

'''
The following function takes in a composition_id and returns a
pd.Series of all attributes.
'''

def mutate_ts(index,x):
    d = split_systems(x)
    ts = pd.Series()
    ts['index'] = index
    ts['formation_energy'] = ternary_systems['delta_e'][index]
    ts['system'] = d[0][0]+'-'+d[1][0]+'-'+d[2][0]
    ts['A'] = d[0][0]
    ts['B'] = d[1][0]
    ts['C'] = d[2][0]

    sum_prop = int(d[0][1])+int(d[1][1])+int(d[2][1])
    prop_a = int(d[0][1])/sum_prop
    prop_b = int(d[1][1])/sum_prop
    prop_c = int(d[2][1])/sum_prop

    ts[d[0][0]+'_prop'] = prop_a
    ts[d[1][0]+'_prop'] = prop_b
    ts[d[2][0]+'_prop'] = prop_c

    # Element Property extraction

    ts['avg_atomic_mass'] = prop_a*elements.loc[d[0][0]].mass + prop_b*elements.loc[d[1][0]].mass + prop_c*elements.loc[d[2][0]].mass
    ts['avg_row'] = prop_a*elements.loc[d[0][0]].period + prop_b*elements.loc[d[1][0]].period + prop_c*elements.loc[d[2][0]].period
    ts['avg_col'] = prop_a*elements.loc[d[0][0]].group + prop_b*elements.loc[d[1][0]].group + prop_c*elements.loc[d[2][0]].group
    ts['max_z_diff'] = abs (max(elements.loc[d[0][0]].z,elements.loc[d[1][0]].z,elements.loc[d[2][0]].z)
                                -min(elements.loc[d[0][0]].z,elements.loc[d[1][0]].z,elements.loc[d[2][0]].z)) # Max Difference in atomic number
    ts['avg_z'] = prop_a*elements.loc[d[0][0]].z + prop_b*elements.loc[d[1][0]].z +prop_c*elements.loc[d[2][0]].z
    ts['max_radius_diff'] = abs (max(elements.loc[d[0][0]].atomic_radii,elements.loc[d[1][0]].atomic_radii,elements.loc[d[2][0]].atomic_radii)
                                     -min(elements.loc[d[0][0]].atomic_radii,elements.loc[d[1][0]].atomic_radii,elements.loc[d[2][0]].atomic_radii)) # Max Difference in atomic radius
    ts['avg_radius'] = prop_a*elements.loc[d[0][0]].atomic_radii + prop_b*elements.loc[d[1][0]].atomic_radii + prop_c*elements.loc[d[2][0]].atomic_radii
    ts['max_en_diff'] = abs (max(elements.loc[d[0][0]].electronegativity,elements.loc[d[1][0]].electronegativity,elements.loc[d[2][0]].electronegativity)
                                 -min(elements.loc[d[0][0]].electronegativity,elements.loc[d[1][0]].electronegativity,elements.loc[d[2][0]].electronegativity)) # Max Difference in atomic electronegativity
    ts['avg_en'] = prop_a*elements.loc[d[0][0]].electronegativity + prop_b*elements.loc[d[1][0]].electronegativity + prop_c*elements.loc[d[2][0]].electronegativity # Avg Difference in electronegativity
    ts['avg_s_elec'] = prop_a*elements.loc[d[0][0]].s_elec +prop_b* elements.loc[d[1][0]].s_elec + prop_c* elements.loc[d[2][0]].s_elec
    ts['avg_p_elec'] = prop_a*elements.loc[d[0][0]].p_elec +prop_b* elements.loc[d[1][0]].p_elec + prop_c* elements.loc[d[2][0]].p_elec
    ts['avg_d_elec'] = prop_a*elements.loc[d[0][0]].d_elec +prop_b* elements.loc[d[1][0]].d_elec + prop_c* elements.loc[d[2][0]].d_elec
    ts['avg_f_elec'] = prop_a*elements.loc[d[0][0]].f_elec +prop_b* elements.loc[d[1][0]].f_elec + prop_c* elements.loc[d[2][0]].f_elec

    temp_sum = ts['avg_s_elec']+ts['avg_p_elec']+ts['avg_d_elec']+ts['avg_f_elec']

    ts['prop_s_elec'] = ts['avg_s_elec']/temp_sum
    ts['prop_p_elec'] = ts['avg_p_elec']/temp_sum
    ts['prop_d_elec'] = ts['avg_d_elec']/temp_sum
    ts['prop_f_elec'] = ts['avg_f_elec']/temp_sum

    return ts


# In[ ]:

# the following chuck of code takes all entries from the
# ternary_systems['composition_id] table and gives a 
# new table with all the features

tern_data = pd.DataFrame()
for index,row in ternary_systems.iterrows():
    temp = mutate_ts(index,row['composition_id'])
    tern_data = tern_data.append(temp,ignore_index= True)
    
col_names = list(tern_data.filter(regex='_prop').columns)
tern_data[col_names].fillna(0) # replace all NaN's with 0 in the '_prop' columns
tern_data.head()



# In[ ]:

# data = pd.merge(master,tern_data)
tern_data.save()


# In[ ]:

master.tail()


# In[ ]:




# In[ ]:




# In[ ]:



