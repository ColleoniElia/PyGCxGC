#!/usr/bin/env python
# coding: utf-8

# In[1]:


# GCxGC Plotter


# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# In[3]:


# settings
Name='HFO' # name of the sample [HFO, VRO]
PlotChromatogramRawData=1 # boolean [0,1]
PlotChromatogramScreenedData=1 # boolean [0,1]
PlotChromatogramClass=1 ## boolean [0,1]
ClassName=['paraffin','aromatic'] # [paraffin, olefin, naphthene, aromatic, cycloalkene,alkini, sulfur, oxygen, nitrogen, NO,NS,SO,NSO]


# In[4]:


# load data
string='Output/ChromatogramData'+Name+'.txt'
ChromatogramData=pd.read_csv(string,delimiter='\t')
ChromatogramData=np.array(ChromatogramData)
Class=ChromatogramData[:,0]
time1=ChromatogramData[:,1]
time2=ChromatogramData[:,2]
percentage=ChromatogramData[:,3]


# In[5]:


if PlotChromatogramRawData==1:
    #load data
    data=pd.read_excel('InstrumentData.xlsx', sheet_name=Name)
    data=np.array(data)
    Time1=data[:,6]
    Time2=data[:,7]
    volume=data[:,9]
    Percentage=volume/np.sum(volume)

    plt.figure(figsize=(15, 10))
    size= 40
    sns.scatterplot(x=Time1, y=Time2, hue=Percentage/np.max(Percentage), size=Percentage/np.max(Percentage),sizes=(50, 1000),palette='inferno')
    plt.grid()
    plt.title('Raw data 2D Chromatogram')
    plt.xlabel("first column [min]", fontsize=size)
    plt.ylabel("second column [sec]", fontsize=size)
    plt.ylim([0,6])
    plt.legend(title="normalized intensity",framealpha=0.2, loc='best', fontsize=15)
    plt.tick_params(axis='x', labelsize=size)
    plt.tick_params(axis='y', labelsize=size)
    plt.show()


# In[6]:


if PlotChromatogramScreenedData==1:
    Time1=[]
    Time2=[]
    Percentage=[]
    for i in range(len(Class)):
        if Class[i]=='screened GCxGC':
            Time1.append(time1[i])
            Time2.append(time2[i])
            Percentage.append(percentage[i])
            
    plt.figure(figsize=(15, 10))
    size= 40
    sns.scatterplot(x=Time1, y=Time2, hue=Percentage/np.max(Percentage), size=Percentage/np.max(Percentage),sizes=(50, 1000),palette='inferno')
    plt.grid()
    plt.title('Screened data 2D Chromatogram')
    plt.xlabel("first column [min]", fontsize=size)
    plt.ylabel("second column [sec]", fontsize=size)
    plt.ylim([0,6])
    plt.legend(title="normalized intensity",framealpha=0.2, loc='best', fontsize=15)
    plt.tick_params(axis='x', labelsize=size)
    plt.tick_params(axis='y', labelsize=size)
    plt.show()


# In[7]:


if PlotChromatogramClass==1:
    for j in range(len(ClassName)):
        Time1=[]
        Time2=[]
        Percentage=[]

        for i in range(len(Class)):
            if Class[i]==ClassName[j]:
                Time1.append(time1[i])
                Time2.append(time2[i])
                Percentage.append(percentage[i])

        plt.figure(figsize=(15, 10))
        size= 40
        sns.scatterplot(x=Time1, y=Time2, hue=Percentage/np.max(Percentage), size=Percentage/np.max(Percentage),sizes=(50, 1000),palette='inferno')
        plt.grid()
        plt.title(ClassName[j]+' 2D Chromatogram')
        plt.xlabel("first column [min]", fontsize=size)
        plt.ylabel("second column [sec]", fontsize=size)
        plt.ylim([0,6])
        plt.legend(title="normalized intensity",framealpha=0.2, loc='best', fontsize=15)
        plt.tick_params(axis='x', labelsize=size)
        plt.tick_params(axis='y', labelsize=size)
        plt.show()
    

