#!/usr/bin/env python
# coding: utf-8

# In[1]:


## generation of SMILEs using the CAS registry number


# In[2]:


import pandas as pd
import cirpy
import numpy as np


# In[3]:


Name='HFO' # name of the sample [HFO, VRO]


# In[4]:


data=pd.read_excel('InstrumentData.xlsx', sheet_name=Name)
data=np.array(data)
CAS1=data[:,3]
CAS2=data[:,13]
CAS3=data[:,16]

FileName="CAS_to_SMILEs_"+Name+".txt"
file=open(FileName, "w")
file.write("SMILE1" + "\t" + "SMILE2"+ "\t" + "SMILE3" + "\n")
SMILE1=list()
SMILE2=list()
SMILE3=list()
for i in range(len(CAS1)):
    SMILE1.append(cirpy.resolve(CAS1[i], 'smiles'))
    SMILE2.append(cirpy.resolve(CAS2[i], 'smiles'))
    SMILE3.append(cirpy.resolve(CAS3[i], 'smiles'))
    print([i])
    file.write(str(SMILE1[i]) + "\t" + str(SMILE2[i])+ "\t" + str(SMILE3[i]) + "\n")
file.close()

