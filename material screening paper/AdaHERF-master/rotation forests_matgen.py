
# coding: utf-8

# In[ ]:

import pandas as pd
import numpy as np
import scipy as sp

from sklearn.cross_validation import StratifiedKFold

from AdaHERF import AdaHERF


# In[ ]:




# In[123]:

binary_sys = pd.read_pickle('../binary_sys')
binary_sys['type'] = 'binary'
ternary_sys = pd.read_pickle('../tern_data')
ternary_sys['type'] = 'ternary'

data = pd.concat([binary_sys,ternary_sys], ignore_index= True)


# In[124]:

features = list(data.columns)
features = [f for f in features if (f.endswith('prop') 
                                    or f.startswith('prop') 
                                    or f.startswith(('avg','max')) 
                                    or f in 'ABC' 
                                    or f == 'formation_energy'
                                    or f == 'type')]
data = data[features]
data = data.fillna(0)
Y = 'formation_energy'
X = list(set(features) - set(['formation_energy','type','A','B','C']))


# In[125]:

skf = StratifiedKFold(data.type,n_folds=10)


# In[ ]:

results = []
for train,test in skf:
    temp = {}
    i = 1
    fit = rot_forest.fit(data.loc[train,X].values,data.loc[train,Y].values)
    
    pred_y = rot_forest.predict(data.loc[test,X].values)
    temp['predicted'] = pred_y
    temp['actual'] = data.loc[test,Y].values
    results.append(temp)
    print i
    i+=1


# In[ ]:


    
    


# In[ ]:




# In[ ]:




# In[ ]:



