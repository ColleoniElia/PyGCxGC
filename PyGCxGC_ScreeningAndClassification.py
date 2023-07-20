#!/usr/bin/env python
# coding: utf-8

# In[1]:


# prova di plot per dati GCxGC


# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.optimize import fsolve
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


# In[3]:


#settings for screening and classification
Name='VRO' # name of the sample [HFO, VRO]
data=pd.read_excel('InstrumentData.xlsx', sheet_name=Name)
LightGasDistribution=pd.read_excel('InstrumentData.xlsx',sheet_name='LightGasResponse')

OlefinesPreferentially=1 # between the three hints select preferentialy an olefines as pyrolysis product, boolean [0,1]
WriteResults=1 # write txt file containing the results, boolean [0,1]
WriteChromatogramData=1 # write txt file containing functional groups classification, boolean [0,1]


# In[4]:


data=np.array(data)
volume=data[:,9] #before peaks screening
IntensityTOT=np.sum(volume) #before peaks screening


# In[5]:


#load smiles
FileName="SMILE_"+Name+".txt"
fileRead=pd.read_table(FileName,dtype=str,delimiter='\t',skiprows=0)
SMILE1=fileRead["SMILE1"].to_numpy()
SMILE2=fileRead["SMILE2"].to_numpy()
SMILE3=fileRead["SMILE3"].to_numpy()

# convert None into 0, necessary for the next step
for i in range(len(data)):
    if SMILE1[i]=='None':
        SMILE1[i]='0'
    if SMILE2[i]=='None':
        SMILE2[i]='0'
    if SMILE3[i]=='None':
        SMILE3[i]='0'


# In[6]:


# Screening of non plausible molecules

SMILE=SMILE1
InitialMolecules=len(SMILE)
IndexMolecules=list()

for i in range(len(SMILE)): #conversion of smile to molecular formula
    if SMILE[i]!='0':
        mol1 = Chem.MolFromSmiles(SMILE[i])
        molecule1 = CalcMolFormula(mol1)
    else:
        molecule1='0'

    if SMILE2[i]!='0':
        mol2 = Chem.MolFromSmiles(SMILE2[i])
        molecule2 = CalcMolFormula(mol2)
    else:
        molecule2='0'

    if SMILE3[i]!='0':
        mol3 = Chem.MolFromSmiles(SMILE3[i])
        molecule3 = CalcMolFormula(mol3)
    else:
        molecule3='0'
        
    #definition of constarints for non plausible molecules
    if 'Br' in molecule1 or 'Cl' in molecule1 or 'F' in molecule1 or 'Si' in molecule1 or '0' in molecule1 :
        if 'Br' in molecule2 or 'Cl' in molecule2 or 'F' in molecule2 or 'Si' in molecule2 or '0' in molecule2:
            if 'Br' in molecule3 or 'Cl' in molecule3 or 'F' in molecule3 or 'Si' in molecule3 or '0' in molecule3:
                IndexMolecules.append(i)
            else:
                SMILE[i]=SMILE3[i]
        else:
            SMILE[i]=SMILE2[i]

# delete rows where all three ints are not suitable
print('Initial Molecules:', InitialMolecules)
print('deleted rows:',len(IndexMolecules))
print('Final Molecules:', InitialMolecules-len(IndexMolecules))


SMILE=np.delete(SMILE,IndexMolecules, axis=0) # best hint between the three
SMILE2=np.delete(SMILE2,IndexMolecules, axis=0)
SMILE3=np.delete(SMILE3,IndexMolecules, axis=0)

Data=np.delete(data,IndexMolecules, axis=0)
# for i in range(len(SMILE)):
#     print(SMILE[i]) #to check that no atomns but CHSNO are present


# In[7]:


if OlefinesPreferentially==1: # betyween pure hydrocarbon hints select preferentially olefines
    SMILEHC=['0']*len(SMILE)

    # indentify pure hydrocarbons molecules
    SMILE1HC=['0']*len(SMILE1)
    SMILE2HC=['0']*len(SMILE2)
    SMILE3HC=['0']*len(SMILE3)
     
    for i in range(len(SMILE)): #convert SMILE to molecular formula
        if SMILE[i]!='0':
            mol1 = Chem.MolFromSmiles(SMILE[i])
            molecule1 = CalcMolFormula(mol1)
        else:
            molecule1='0'

        if SMILE2[i]!='0':
            mol2 = Chem.MolFromSmiles(SMILE2[i])
            molecule2 = CalcMolFormula(mol2)
        else:
            molecule2='0'

        if SMILE3[i]!='0':
            mol3 = Chem.MolFromSmiles(SMILE3[i])
            molecule3 = CalcMolFormula(mol3)
        else:
            molecule3='0'
            
        # identify only HC molecules
        if 'O' in molecule1 or "S" in molecule1 or 'N' in molecule1 or 'Br' in molecule1 or 'Cl' in molecule1 or 'F' in molecule1 or 'Si' in molecule1:
            SMILE1HC[i]=SMILE1HC[i]
        else:
            SMILE1HC[i]=SMILE1[i]
        
        if 'O' in molecule2 or "S" in molecule2 or 'N' in molecule2 or 'Br' in molecule2 or 'Cl' in molecule2 or 'F' in molecule2 or 'Si' in molecule2:
            SMILE2HC[i]=SMILE2HC[i]
        else:
            SMILE2HC[i]=SMILE2[i]
        
        if 'O' in molecule3 or "S" in molecule3 or 'N' in molecule3 or 'Br' in molecule3 or 'Cl' in molecule3 or 'F' in molecule3 or 'Si' in molecule3:
            SMILE3HC[i]=SMILE3HC[i]
        else:
            SMILE3HC[i]=SMILE3[i]

    # allocation of paraffins and olefines only        
    SMILE1_Par_and_Olef=['0']*len(SMILE1HC)
    SMILE2_Par_and_Olef=['0']*len(SMILE2HC)
    SMILE3_Par_and_Olef=['0']*len(SMILE3HC)
    
    for i in range(len(Data)):
        if 'c1' in SMILE1HC[i] or 'C1' in SMILE1HC[i] or '#' in SMILE1HC[i]:
            SMILE1_Par_and_Olef[i]='0'

        else:
            SMILE1_Par_and_Olef[i]=SMILE1HC[i]
      
        if 'c1' in SMILE2HC[i] or 'C1' in SMILE2HC[i] or '#' in SMILE2HC[i]:
            SMILE2_Par_and_Olef[i]='0'

        else:
            SMILE2_Par_and_Olef[i]=SMILE2HC[i]
        
        if 'c1' in SMILE3HC[i] or 'C1' in SMILE3HC[i] or '#' in SMILE3HC[i]:
            SMILE3_Par_and_Olef[i]='0'
        else:
            SMILE3_Par_and_Olef[i]=SMILE3HC[i]

    #between paraffins and olefines select preferentially olefines
    for i in range(len(Data)):
        if '=' in SMILE1_Par_and_Olef[i]:
            SMILE[i]=SMILE1_Par_and_Olef[i]
            #formula[i]=formula1[i]

        elif '=' in SMILE2_Par_and_Olef[i]:
            SMILE[i]=SMILE2_Par_and_Olef[i]
            #formula[i]=formula2[i]
        elif '=' in SMILE3_Par_and_Olef[i]:
            SMILE[i]=SMILE3_Par_and_Olef[i]
            #formula[i]=formula3[i]

        else:
            SMILE[i]=SMILE[i]


# In[8]:


# data rearrangment according to ligh molecules detected in the first peak
time1=Data[:,6]
time2=Data[:,7]
volume=Data[:,9]

# wrap-around for retention time in second column > 6sec
for i in range(len(time2)):
    if time2[i]>6:
        time2[i]=time2[i]-6
    else:
        time2[i]=time2[i]
            
GasDistribution=np.array(LightGasDistribution)
sample=GasDistribution[:,0]
print(sample)

for i in range(len(sample)):
    if sample[i]==Name:
        Distribution=GasDistribution[i,:]
        limitTime=GasDistribution[i,6]

IndexLight=next(x[0] for x in enumerate(time1) if x[1] > limitTime) #number of componments that have to be summed in the first peak

# print(Distribution)   
# print('index light',IndexLight)
# print(limitTime)
        
# Data rearrangement to account for light gas
time1Light=np.min(time1[0:IndexLight])
time2Light=np.min(time2[0:IndexLight])

volumeLight=np.sum(volume[0:IndexLight])

time1=np.insert(time1[IndexLight:],0,time1Light,axis=0)
time2=np.insert(time2[IndexLight:],0,time2Light,axis=0)
volume=np.insert(volume[IndexLight:],0,volumeLight,axis=0)

smile=SMILE[IndexLight:] #remove smile strings for light gas

percentage=volume/np.sum(volume)
Intensity=np.sum(volume)


# In[9]:


print('intensity remeaning peaks/intensity total peaks',Intensity/IntensityTOT)


# In[10]:


#initialization of array with atoms
C=np.zeros(len(smile))
H=np.zeros(len(smile))
N=np.zeros(len(smile))
S=np.zeros(len(smile))
O=np.zeros(len(smile))
Formula=[]

# creation of array with atomns count
for i in range(len(C)): #conversion of smile to molecular formula
    mol = Chem.MolFromSmiles(smile[i])
    formula = CalcMolFormula(mol)
    Formula.append(formula)
    formula=re.findall('[A-Z][a-z]?|[0-9]+', formula)
    
    if ('C' in formula):
        index=formula.index('C')
        if index+1>len(formula)-1:
            C[i]=1
        elif formula[index+1].isdigit()==True:
            C[i]=(np.array(formula[index+1]))
        else:
            C[i]=1
    else:
        C[i]=0
        
    
    if ('H' in formula):
        index=formula.index('H')
        if index+1>len(formula)-1:
            H[i]=1
        elif formula[index+1].isdigit()==True:
            H[i]=(np.array(formula[index+1]))
        else:
            H[i]=1
    else:
        H[i]=0
        
    
    if ('N' in formula):
        index=formula.index('N')
        if index+1>len(formula)-1:
            N[i]=1
        elif formula[index+1].isdigit()==True:
            N[i]=(np.array(formula[index+1]))
        else:
            N[i]=1
    else:
        N[i]=0

    
    if ('S' in formula):
        index=formula.index('S')
        if index+1>len(formula)-1:
            S[i]=1
        elif formula[index+1].isdigit()==True:
            S[i]=(np.array(formula[index+1]))
        else:
            S[i]=1
    else:
        S[i]=0

        
    
    if ('O' in formula):
        index=formula.index('O')
        if index+1>len(formula)-1:
            O[i]=1
        elif formula[index+1].isdigit()==True:
            O[i]=(np.array(formula[index+1]))
        else:
            O[i]=1
    else:
        O[i]=0   

# # check if formula, atomic matrix and smile correspond 
# for i in range(len(Formula)):
#     print(Formula[i],C[i],H[i],S[i],N[i],O[i],smile[i])


# In[11]:


#initialization of arrays for molecular classes
Hydrocarbons=np.zeros(len(time1)-1)
Scomponent=np.zeros(len(time1)-1)
Ocomponent=np.zeros(len(time1)-1)
Ncomponent=np.zeros(len(time1)-1)
NSO=np.zeros(len(time1)-1)
NS=np.zeros(len(time1)-1)
NO=np.zeros(len(time1)-1)
SO=np.zeros(len(time1)-1)

HydrocarbonsTime1=np.zeros(len(time1)-1) #to allocate time 1
ScomponentTime1=np.zeros(len(time1)-1)
OcomponentTime1=np.zeros(len(time1)-1)
NcomponentTime1=np.zeros(len(time1)-1)
NSO_Time1=np.zeros(len(time1)-1)
NS_Time1=np.zeros(len(time1)-1)
NO_Time1=np.zeros(len(time1)-1)
SO_Time1=np.zeros(len(time1)-1)

HydrocarbonsTime2=np.zeros(len(time1)-1) #to allocate time 2
ScomponentTime2=np.zeros(len(time1)-1)
OcomponentTime2=np.zeros(len(time1)-1)
NcomponentTime2=np.zeros(len(time1)-1)
NSO_Time2=np.zeros(len(time1)-1)
NS_Time2=np.zeros(len(time1)-1)
NO_Time2=np.zeros(len(time1)-1)
SO_Time2=np.zeros(len(time1)-1)

HC_molecules=list() # to allocate SMILE
S_molecules=list()
O_molecules=list()
N_molecules=list()
NSO_molecules=list()
NS_molecules=list()
NO_molecules=list()
SO_molecules=list()

# classification of molecules according to atomic composition
for i in range(len(Formula)):
    if O[i]!=0 and N[i]!=0 and S[i]!=0:
        NSO[i]=percentage[i+1]*100
        NSO_molecules.append(smile[i])
        NSO_Time1[i]=time1[i+1]
        NSO_Time2[i]=time2[i+1]
        
    elif N[i]!=0 and S[i]!=0:
        NS[i]=percentage[i+1]*100
        NS_molecules.append(smile[i])
        NS_Time1[i]=time1[i+1]
        NS_Time2[i]=time2[i+1]
        
    elif O[i]!=0 and N[i]!=0:
        NO[i]=percentage[i+1]*100
        NO_molecules.append(smile[i])
        NO_Time1[i]=time1[i+1]
        NO_Time2[i]=time2[i+1]
        
    elif O[i]!=0 and S[i]!=0:
        SO[i]=percentage[i+1]*100
        SO_molecules.append(smile[i])
        SO_Time1[i]=time1[i+1]
        SO_Time2[i]=time2[i+1]
        
    elif O[i]!=0:
        Ocomponent[i]=percentage[i+1]*100
        O_molecules.append(smile[i])
        OcomponentTime1[i]=time1[i+1]
        OcomponentTime2[i]=time2[i+1]
        
    elif S[i]!=0:
        Scomponent[i]=percentage[i+1]*100
        S_molecules.append(smile[i])
        ScomponentTime1[i]=time1[i+1]
        ScomponentTime2[i]=time2[i+1]
        
    elif N[i]!=0:
        Ncomponent[i]=percentage[i+1]*100
        N_molecules.append(smile[i])
        NcomponentTime1[i]=time1[i+1]
        NcomponentTime2[i]=time2[i+1]
        
    else:
        Hydrocarbons[i]=percentage[i+1]*100
        HC_molecules.append(smile[i])
        HydrocarbonsTime1[i]=time1[i+1]
        HydrocarbonsTime2[i]=time2[i+1]
    
tot=percentage[0]*100+np.sum(Hydrocarbons)+np.sum(Ocomponent)+np.sum(Scomponent)+np.sum(Ncomponent)+np.sum(NSO)+np.sum(NS)+np.sum(NO)+np.sum(SO)

if tot<100.5 and tot>99.5:
    print()
else:
    print('error in classification, the sum of different classes does not close to 100')
    exit()
    

# plt.figure(figsize=(15, 10))
# plt.bar(['lightgas','HC','O','S','N','NSO','NS','NO','SO'],
#          [percentage[0]*100,np.sum(Hydrocarbons),np.sum(Ocomponent),np.sum(Scomponent),np.sum(Ncomponent),
#          np.sum(NSO),np.sum(NS),np.sum(NO),np.sum(SO)])
# plt.grid()

# print('HC %',np.sum(Hydrocarbons))
# print('total molecules', len(smile))
# print('Hydrocarbon molecules', len(HC_molecules))
# print('O molecules', len(O_molecules))
# print('S molecules', len(S_molecules))
# print('N molecules', len(N_molecules))
# print('NSO molecules', len(NSO_molecules))
# print('NS molecules', len(NS_molecules))
# print('NO molecules', len(NO_molecules))
# print('SO molecules', len(SO_molecules))


# In[12]:


#classification of SMILEs for different molecular classes
smileHC=list() 
smileNS=list()
smileNO=list()
smileSO=list()
smileO=list()
smileS=list()
smileN=list()
smileNSO=list()

percentageHC=[]
percentageNS=[]
percentageNO=[]
percentageSO=[]
percentageO=[]
percentageS=[]
percentageN=[]
percentageNSO=[]

time1HC=list() 
time1NS=list()
time1NO=list()
time1SO=list()
time1O=list()
time1S=list()
time1N=list()
time1NSO=list()

time2HC=list() 
time2NS=list()
time2NO=list()
time2SO=list()
time2O=list()
time2S=list()
time2N=list()
time2NSO=list()

for i in range(len(smile)):
    if Hydrocarbons[i]!=0:
        smileHC.append(smile[i])
        percentageHC.append(percentage[i+1])
        time1HC.append(HydrocarbonsTime1[i])
        time2HC.append(HydrocarbonsTime2[i])
        
    if NSO[i]!=0:
        smileNSO.append(smile[i])
        percentageNSO.append(percentage[i+1])
        time1NSO.append(NSO_Time1[i])
        time2NSO.append(NSO_Time2[i])
        
    if NS[i]!=0:
        smileNS.append(smile[i])
        percentageNS.append(percentage[i+1])
        time1NS.append(NS_Time1[i])
        time2NS.append(NS_Time2[i])
        
    if NO[i]!=0:
        smileNO.append(smile[i])
        percentageNO.append(percentage[i+1])
        time1NO.append(NO_Time1[i])
        time2NO.append(NO_Time2[i])

    if SO[i]!=0:
        smileSO.append(smile[i])
        percentageSO.append(percentage[i+1])
        time1SO.append(SO_Time1[i])
        time2SO.append(SO_Time2[i])
        
    if Ocomponent[i]!=0:
        smileO.append(smile[i])
        percentageO.append(percentage[i+1])
        time1O.append(OcomponentTime1[i])
        time2O.append(OcomponentTime2[i])
        
    if Ncomponent[i]!=0:
        smileN.append(smile[i])
        percentageN.append(percentage[i+1])
        time1N.append(NcomponentTime1[i])
        time2N.append(NcomponentTime2[i])
        
    if Scomponent[i]!=0:
        smileS.append(smile[i])
        percentageS.append(percentage[i+1])
        time1S.append(ScomponentTime1[i])
        time2S.append(ScomponentTime2[i])


# In[13]:


# classification of hydrocarbons according to functional groups
paraffins=[]
aromatics=[]
naphtenes=[]
cycloalkenes=[]
olefines=[]
alkini=[]
other=[]

SMILE_paraffins=list()
SMILE_aromatics=list()
SMILE_naphtenes=list()
SMILE_cycloalkenes=list()
SMILE_olefines=list()
SMILE_alkini=list()

time1_paraffins=[] 
time2_paraffins=[] 
time1_aromatics=[] 
time2_aromatics=[] 
time1_naphtenes=[] 
time2_naphtenes=[] 
time1_cycloalkenes=[] 
time2_cycloalkenes=[] 
time1_olefines=[] 
time2_olefines=[] 
time1_alkini=[] 
time2_alkini=[] 

typeHC=list()

for i in range(len(smileHC)):
    if smileHC[i]!=0:
            if 'c1' in smileHC[i]:
                time1_aromatics.append(time1HC[i])
                time2_aromatics.append(time2HC[i])
                aromatics.append(percentageHC[i]*100)
                SMILE_aromatics.append(smileHC[i])
                typeHC.append('aromatics')
                
            elif '#' in smileHC[i]:
                time1_alkini.append(time1HC[i])
                time2_alkini.append(time2HC[i])
                alkini.append(percentageHC[i]*100)
                SMILE_alkini.append(smileHC[i])
                typeHC.append('alkini')
                
            elif 'C1' in smileHC[i] and '=' in smileHC[i]:
                time1_cycloalkenes.append(time1HC[i])
                time2_cycloalkenes.append(time2HC[i])
                cycloalkenes.append(percentageHC[i]*100) 
                SMILE_cycloalkenes.append(smileHC[i])
                typeHC.append('cycloalkenes')
                
            elif 'C1' in smileHC[i]:
                time1_naphtenes.append(time1HC[i])
                time2_naphtenes.append(time2HC[i])
                naphtenes.append(percentageHC[i]*100)
                SMILE_naphtenes.append(smileHC[i])
                typeHC.append('naphtenes')
                
            elif '=' in smileHC[i]:
                time1_olefines.append(time1HC[i])
                time2_olefines.append(time2HC[i])
                olefines.append(percentageHC[i]*100)
                SMILE_olefines.append(smileHC[i])
                typeHC.append('olefines')
        
            elif 'C' in smileHC[i]:
                time1_paraffins.append(time1HC[i])
                time2_paraffins.append(time2HC[i])
                paraffins.append(percentageHC[i]*100)
                SMILE_paraffins.append(smileHC[i])
                typeHC.append('paraffin')
            
            else:
                other.append(percentageHC[i]*100) #to check if something is not classified

percentage_paraffins=paraffins
percentage_aromatics=aromatics
percentage_naphtenes=naphtenes
percentage_cycloalkenes=cycloalkenes
percentage_olefines=olefines
percentage_alkini=alkini

paraffins=np.sum(paraffins)
aromatics=np.sum(aromatics)
naphtenes=np.sum(naphtenes)
cycloalkenes=np.sum(cycloalkenes)
olefines=np.sum(olefines)
alkini=np.sum(alkini)

# plt.figure(figsize=(15, 10))
# plt.bar(['lightgas','paraffins','aromatics','naphtenes','cycloalkenes','olefines','alkini',
#          'O','S','N','NSO','NS','NO','SO'],
#          [percentage[0]*100,paraffins,aromatics,naphtenes,cycloalkenes,olefines,alkini,
#           np.sum(Ocomponent),np.sum(Scomponent),np.sum(Ncomponent),
#          np.sum(NSO),np.sum(NS),np.sum(NO),np.sum(SO)]) 
# plt.grid()
# plt.xticks(rotation='vertical')

tot=np.absolute(paraffins+aromatics+naphtenes+cycloalkenes+olefines+alkini-np.sum(Hydrocarbons))
if tot<0.5:
    print()
else:
    print('error in classification of HC molecules')
    exit()
    


# In[14]:


# # Hydrocarbons distribution
# print('Hydrocarbon molecules', len(smileHC))
# print('Paraffin molecules', len(SMILE_paraffins))
# print('Olefine molecules', len(SMILE_olefines))
# print('Aromatics molecules', len(SMILE_aromatics))
# print('Alkin molecules', len(SMILE_alkini))
# print('Naphtenes molecules', len(SMILE_naphtenes))
# print('Cycloalkenes molecules', len(SMILE_cycloalkenes))

# HCdistribution="HC_distribution"+title+".txt"
# file=open(HCdistribution, "w")
# file.write("type" + "\t" + "smile"+ "\t" + "percentage" + "\n")
# for i in range(len(smileHC)):
#     print(typeHC[i],smileHC[i],percentageHC[i]*100)
#     file.write(str(typeHC[i]) + "\t" + str(smileHC[i])+ "\t" +str(percentageHC[i]*100) +"\n")
# file.close()


# In[15]:


# classification of oxygen moleculse according to functional groups
Oxydril=[]
Carbonil=[]
Ester=[]
Acid=[]
Oxydril_plus_Carbonil=[]
Furan=[]

SMILE_oxydril=list()
SMILE_carbonyl=list()
SMILE_ester=list()
SMILE_acid=list()
SMILE_oxydril_plus_carbonil=list()
SMILE_furan=list()

typeO=list()

for i in range(len(smileO)):
            if 'CO' in smileO[i] and '=O' in smileO[i]:
                Oxydril_plus_Carbonil.append(percentageO[i]*100)
                SMILE_oxydril_plus_carbonil.append(smileO[i])
                typeO.append('Oxydril plus carbonil')
                
            elif 'o' in smileO[i]:
                Furan.append(percentageO[i]*100)
                SMILE_furan.append(smileO[i])
                typeO.append('Furan')
        
            elif '(O)=O' in smileO[i]:
                Acid.append(percentageO[i]*100)
                SMILE_acid.append(smileO[i])
                typeO.append('Acid')
                
            elif '(=O)O' in smileO[i]:
                Ester.append(percentageO[i]*100)
                SMILE_ester.append(smileO[i])
                typeO.append('Ester')
        
            elif '=O' in smileO[i]:
                Carbonil.append(percentageO[i]*100)
                SMILE_carbonyl.append(smileO[i])
                typeO.append('Carbonyl')
        
            elif 'O' in smileO[i]:
                Oxydril.append(percentageO[i]*100)
                SMILE_oxydril.append(smileO[i])
                typeO.append('Oxydril')
                

      
            


Oxydril_plus_Carbonil=np.sum(Oxydril_plus_Carbonil)
Ester=np.sum(Ester)
Acid=np.sum(Acid)
Carbonil=np.sum(Carbonil)
Oxydril=np.sum(Oxydril)
Furan=np.sum(Furan)

# plt.figure(figsize=(15, 10))
# plt.bar(['lightgas','paraffins','aromatics','naphtenes','cycloalkenes','olefines','alkini',
#          'Oxydril_plus_Carbonil','Ester','Acid','Carbonil','Oxydril','Furan','S','N','NSO','NS','NO','SO'],
#          [percentage[0]*100,paraffins,aromatics,naphtenes,cycloalkenes,olefines,alkini,
#           Oxydril_plus_Carbonil,Ester,Acid,Carbonil,Oxydril,Furan,
#           np.sum(Scomponent),np.sum(Ncomponent),np.sum(NSO),np.sum(NS),np.sum(NO),np.sum(SO)]) 
# plt.xticks(rotation='vertical')
# plt.grid()

tot=np.absolute(Oxydril_plus_Carbonil+Ester+Acid+Carbonil+Oxydril+Furan-np.sum(Ocomponent))
if tot<0.5:
    print()
else:
    print('error in classification of oxygen molecules')
    exit()


# In[16]:


# # Oxygen molecules distribution
# print('O molecules', len(smileO))
# print('Oxydril plus carbonil molecules', len(SMILE_oxydril_plus_carbonil))
# print('Ester molecules', len(SMILE_ester))
# print('Acid molecules', len(SMILE_acid))
# print('Carbonil molecules', len(SMILE_carbonyl))
# print('Oxydril molecules', len(SMILE_oxydril))
# print('Furan molecules', len(SMILE_furan))

# Odistribution="O_distribution"+title+".txt"
# file=open(Odistribution, "w")
# file.write("type" + "\t" + "smile"+ "\t" + "percentage" + "\n")
# for i in range(len(smileO)):
#     print(typeO[i],smileO[i],percentageO[i]*100)
#     file.write(str(typeO[i]) + "\t" + str(smileO[i])+ "\t" +str(percentageO[i]) +"\n")
# file.close()


# In[17]:


# classification of S molecules according to functional groups
Sulfide=[]
Thiophene=[]
Sulfide_plus_Thiophene=[]

SMILE_sulfides=list()
SMILE_thiophene=list()
SMILE_sulfide_plus_thiophene=list()

typeS=list()

for i in range(len(smileS)):
            if 's' in smileS[i] and 'S' in smileS[i]:
                Sulfide_plus_Thiophene.append(percentageS[i]*100)
                SMILE_sulfide_plus_thiophene.append(smileS[i])
                typeS.append('Sulfide and thiophene')
                
            elif 's' in smileS[i]:
                Thiophene.append(percentageS[i]*100)
                SMILE_thiophene.append(smileS[i])
                typeS.append('thiophene')
        
            elif 'S' in smileS[i]:
                Sulfide.append(percentageS[i]*100)
                SMILE_sulfides.append(smileS[i])
                typeS.append('Sulfide')


Sulfide_plus_Thiophene=np.sum(Sulfide_plus_Thiophene)
Sulfide=np.sum(Sulfide)
Thiophene=np.sum(Thiophene)

tot=np.absolute(Sulfide_plus_Thiophene+Sulfide+Thiophene-np.sum(Scomponent))
if tot<0.5:
    print()
else:
    print('error in classification of sulfur molecules')
    exit()


# In[18]:


# # Sulfur molecules distribution
# print('S molecules', len(smileS))
# print('Thiophene plus sulfides molecules', len(SMILE_sulfide_plus_thiophene))
# print('Thiophene molecules', len(SMILE_thiophene))
# print('Sulfide molecules', len(SMILE_sulfides))
    
# Sdistribution="S_distribution"+title+".txt"
# file=open(Sdistribution, "w")
# file.write("type" + "\t" + "smile"+ "\t" + "percentage" + "\n")
# for i in range(len(smileS)):
#     print(typeS[i],smileS[i])
#     file.write(str(typeS[i]) + "\t" + str(smileS[i])+ "\t" +str(percentageS[i]) +"\n")
# file.close()


# In[19]:


# classification of N molecules according to functional groups
Pyridine=[]
Ammine=[]
Pyridine_plus_Ammine=[]

SMILE_pyridine=list()
SMILE_ammine=list()
SMILE_pyridine_plus_ammine=list()

typeN=list()

for i in range(len(smileN)):
            if 'n' in smileN[i] and 'N' in smileN[i]:
                Pyridine_plus_Ammine.append(percentageN[i]*100)
                SMILE_pyridine_plus_ammine.append(smileN[i])
                typeN.append('pyridine and ammine')
                
            elif 'n' in smileN[i]:
                Pyridine.append(percentageN[i]*100)
                SMILE_pyridine.append(smileN[i])
                typeN.append('pyridine')
        
            elif 'N' in smileN[i]:
                Ammine.append(percentageN[i]*100)
                SMILE_ammine.append(smileN[i])
                typeN.append('ammine')
                
Pyridine_plus_Ammine=np.sum(Pyridine_plus_Ammine)
Ammine=np.sum(Ammine)
Pyridine=np.sum(Pyridine)

tot=np.absolute(Pyridine_plus_Ammine+Ammine+Pyridine-np.sum(Ncomponent))

if tot<0.5:
    print()
else:
    print('error in classification of nitrogen molecules')
    exit()


# In[20]:


# # Nitrogen molecules distribution
# print('N molecules', len(smileN))
# print('Pyridine plus ammine molecules', len(SMILE_pyridine_plus_ammine))
# print('Pyridine molecules', len(SMILE_pyridine))
# print('Ammine molecules', len(SMILE_ammine))

# #print(SMILE_ammine) #check
# Ndistribution="NCA distribution/N_distribution"+title+".txt"
# file=open(Ndistribution, "w")
# file.write("type" + "\t" + "smile"+ "\t" + "percentage" + "\n")
# for i in range(len(smileN)):
#     print(typeN[i],smileN[i])
#     file.write(str(typeN[i]) + "\t" + str(smileN[i])+ "\t" +str(percentageN[i]) +"\n")
# file.close()


# In[21]:


#definition of light gas
C2int=Distribution[2]
CH4int=Distribution[3]
H2Sint=Distribution[4]
CO2int=Distribution[5]

ratioC2_CH4=C2int/CH4int
ratioC2_H2S=C2int/H2Sint
ratioC2_CO2=C2int/CO2int

LightGas=percentage[0]*100

def equations(variables,tot,ratioC2_CH4,ratioC2_H2S,ratioC2_CO2):
    CH4,C2,H2S,CO2=variables
    eq1=CH4+C2+H2S+CO2-tot
    eq2=C2-ratioC2_CH4*CH4
    eq3=C2-ratioC2_H2S*H2S
    eq4=C2-ratioC2_CO2*CO2 
    return [eq1,eq2,eq3,eq4]

firstguess=[6,0,0,0]
[CH4,C2,H2S,CO2]=fsolve(equations,firstguess,args=(LightGas,ratioC2_CH4,ratioC2_H2S,ratioC2_CO2,))


tot=np.absolute(CH4+C2+H2S+CO2-LightGas)
if tot<0.5:
    print()
else:
    print('error in classification of light gas')
    exit()
    


# In[22]:


plt.rcParams['font.size'] = 38
plt.figure(figsize=(15, 10))

# plt.bar(['CH4','C2+CO','H2S','CO2'],[CH4,C2,H2S,CO2],color = "red")   
# plt.bar(['paraffins','aromatics','naphtenes','cycloalkenes','olefines','alkini'],
#         [paraffins,aromatics,naphtenes,cycloalkenes,olefines,alkini], color='black')
# plt.bar(['Oxydril_plus_Carbonil','Ester','Acid','Carbonil','Oxydril','Furan'],
#         [Oxydril_plus_Carbonil,Ester,Acid,Carbonil,Oxydril,Furan], color='blue')
# plt.bar(['Thiophene','Sulfide','Thiophene+Sulfide'],
#         [Thiophene,Sulfide,Sulfide_plus_Thiophene], color='orange')
# plt.bar(['N_aromatic','N_chain','N_aromatic+chain'],
#         [Pyridine,Ammine,Pyridine_plus_Ammine], color='green')
# plt.bar(['NSO','NS','NO','SO'],[np.sum(NSO),np.sum(NS),np.sum(NO),np.sum(SO)],color='purple')

plt.bar(['1','2','3','4'],[CH4,C2,H2S,CO2],color = "red")   
plt.bar(['5','6','7','8','9','10'],
        [paraffins,aromatics,naphtenes,cycloalkenes,olefines,alkini], color='black')
plt.bar(['11','12','13','14','15','16'],
        [Oxydril_plus_Carbonil,Ester,Acid,Carbonil,Oxydril,Furan], color='blue')
plt.bar(['17','18','19'],
        [Thiophene,Sulfide,Sulfide_plus_Thiophene], color='orange')
plt.bar(['20','21','22'],
        [Pyridine,Ammine,Pyridine_plus_Ammine], color='green')
plt.bar(['23','24','25','26'],[np.sum(NSO),np.sum(NS),np.sum(NO),np.sum(SO)],color='purple')


plt.xticks(rotation='vertical')
plt.legend(['Light Molecules','Hydrocarbons','Oxygenated molecules','Sulfur Molecules','Nitrogen Molecules','Heteroatoms'],fontsize=25)
plt.ylabel('mass %')
plt.xlabel('Class Index')
#plt.ylim([0,50])
#plt.title(title)
plt.grid()

#plt.savefig('Output/Distribution_'+Name+'.jpg',bbox_inches="tight")

tot=np.absolute(paraffins+aromatics+naphtenes+cycloalkenes+olefines+alkini+
             Oxydril_plus_Carbonil+Ester+Acid+Carbonil+Oxydril+Furan+
             Thiophene+Sulfide+Sulfide_plus_Thiophene+
             Pyridine+Ammine+Pyridine_plus_Ammine+
             np.sum(NSO)+np.sum(NS)+np.sum(NO)+np.sum(SO)+
             CH4+C2+H2S+CO2)
if tot<100.5 and tot>99.5:
    print()
else:
    print('error in classification')
    exit()

print('CH4',CH4)
print('C2+CO',C2)
print('H2S',H2S)
print('CO2',CO2)
print('paraffins',paraffins)
print('aromatics',aromatics)
print('naphtenes',naphtenes)
print('cycloalkenes',cycloalkenes)
print('olefines',olefines)
print('alkini',alkini)
print('Oxydril_plus_Carbonil',Oxydril_plus_Carbonil)
print('Ester',Ester)
print('Acid',Acid)
print('Carbonil',Carbonil)
print('Oxydril',Oxydril)
print('Furan',Furan)
print('Thiophene',Thiophene)
print('Sulfide',Sulfide)
print('Sulfide_plus_Thiophene',Sulfide_plus_Thiophene)
print('Pyridine',Pyridine)
print('Ammine',Ammine)
print('Pyridine_plus_Ammine',Pyridine_plus_Ammine)
print('NSO',np.sum(NSO))
print('NS',np.sum(NS))
print('NO',np.sum(NO))
print('SO',np.sum(SO))


# In[23]:


# write results
if WriteResults==1:
    FileResults="Output/Distribution"+Name+".txt"
    file=open(FileResults, "w")
    file.write('CH4' + "\t" + str(CH4) +"\n")
    file.write('C2_CO' + "\t" + str(C2) +"\n")
    file.write('H2S' + "\t" + str(H2S) +"\n")
    file.write('CO2' + "\t" + str(CO2) +"\n")
    file.write('paraffins' + "\t" + str(paraffins) +"\n")
    file.write('aromatics' + "\t" + str(aromatics) +"\n")
    file.write('naphtenes' + "\t" + str(naphtenes) +"\n")
    file.write('cycloalkenes' + "\t" + str(cycloalkenes)+ "\n")
    file.write('olefines' + "\t" + str(olefines) +"\n")
    file.write('alkini' + "\t" + str(alkini) +"\n")
    file.write('Oxydril_plus_Carbonil' + "\t" + str(Oxydril_plus_Carbonil)+ "\n")
    file.write('Ester' + "\t" + str(Ester) +"\n")
    file.write('Acid' + "\t" + str(Acid) +"\n")
    file.write('Carbonil' + "\t" + str(Carbonil)+ "\n")
    file.write('Oxydril' + "\t" + str(Oxydril) +"\n")
    file.write('Furan' + "\t" + str(Furan) +"\n")
    file.write('Thiophene' + "\t" + str(Thiophene)+ "\n")
    file.write('Sulfide' + "\t" + str(Sulfide) +"\n")
    file.write('Sulfide_plus_Thiophene' + "\t" + str(Sulfide_plus_Thiophene)+ "\n")
    file.write('N_aromatic' + "\t" + str(Pyridine) +"\n")
    file.write('N_chain' + "\t" + str(Ammine) +"\n")
    file.write('N_aromatic+chain' + "\t" + str(Pyridine_plus_Ammine) +"\n")
    file.write('NSO' + "\t" + str(np.sum(NSO))+ "\n")
    file.write('NS' + "\t" + str(np.sum(NS)) +"\n")
    file.write('NO' + "\t" + str(np.sum(NO)) +"\n")
    file.write('SO' + "\t" + str(np.sum(SO)) +"\n")
    file.close()


# In[24]:


if WriteChromatogramData==1:
    
    FileResults="Output/ChromatogramData"+Name+".txt"
    file=open(FileResults, "w")
    file.write('class' + "\t" +'time1' + "\t" +'time2' + "\t" + 'percentage' +"\n")
    
    for i in range(len(smile)):
        file.write('screened GCxGC'+ "\t" +str(time1[i]) + "\t" + str(time2[i])+ "\t" + str(percentage[i])+"\n")
    
    for i in range(len(time1_paraffins)):
        file.write('paraffin'+ "\t" +str(time1_paraffins[i]) + "\t" + str(time2_paraffins[i])+ "\t" + str(percentage_paraffins[i])+"\n")

    for i in range(len(time1_olefines)):
        file.write('olefin'+ "\t" +str(time1_olefines[i]) + "\t" + str(time2_olefines[i])+ "\t" + str(percentage_olefines[i])+"\n")

    for i in range(len(time1_naphtenes)):
        file.write('naphthene'+ "\t" +str(time1_naphtenes[i]) + "\t" + str(time2_naphtenes[i])+ "\t" + str(percentage_naphtenes[i])+"\n")
    
    for i in range(len(time1_aromatics)):
        file.write('aromatic'+ "\t" +str(time1_aromatics[i]) + "\t" + str(time2_aromatics[i])+ "\t" + str(percentage_aromatics[i])+"\n")
    
    for i in range(len(time1_cycloalkenes)):
        file.write('cycloalkene'+ "\t" +str(time1_cycloalkenes[i]) + "\t" + str(time2_cycloalkenes[i])+ "\t" + str(percentage_cycloalkenes[i])+"\n")

    for i in range(len(time1_alkini)):
        file.write('alkini'+ "\t" +str(time1_alkini[i]) + "\t" + str(time2_alkini[i])+ "\t" + str(percentage_alkini[i])+"\n")
    
    for i in range(len(time1S)):
        file.write('sulfur'+ "\t" +str(time1S[i]) + "\t" + str(time2S[i])+ "\t" + str(percentageS[i])+"\n")
     
    for i in range(len(time1O)):
        file.write('oxygen'+ "\t" +str(time1O[i]) + "\t" + str(time2O[i])+ "\t" + str(percentageO[i])+"\n")
     
    for i in range(len(time1N)):
        file.write('nitrogen'+ "\t" +str(time1N[i]) + "\t" + str(time2N[i])+ "\t" + str(percentageN[i])+"\n")
    
    for i in range(len(time1NO)):
        file.write('NO'+ "\t" +str(time1NO[i]) + "\t" + str(time2NO[i])+ "\t" + str(percentageNO[i])+"\n")
     
    for i in range(len(time1NS)):
        file.write('NS'+ "\t" +str(time1NS[i]) + "\t" + str(time2NS[i])+ "\t" + str(percentageNS[i])+"\n")
     
    for i in range(len(time1SO)):
        file.write('SO'+ "\t" +str(time1SO[i]) + "\t" + str(time2SO[i])+ "\t" + str(percentageSO[i])+"\n")
     
    for i in range(len(time1NSO)):
        file.write('NSO'+ "\t" +str(time1NSO[i]) + "\t" + str(time2NSO[i])+ "\t" + str(percentageNSO[i])+"\n")
     
    file.close()


# In[ ]:




