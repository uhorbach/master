
# coding: utf-8

# In[2]:


import pandas as pd
import math as mt
import numpy as np
import sys
import matplotlib, matplotlib.pyplot as plt


# In[ ]:


if len(sys.argv) < 4:
    print("You need to pass four files as arguments! Three inputfiles and one outputfile")
    exit()


# In[ ]:


inputfiles = sys.argv[1:3]
outputfile = sys.argv[4]
print(inputfiles)

# In[]:
data=list()
i=1
while i<=3: #for dose1 dose2 dose3

	print("Processing files: ", sys.argv[i], "into", outputfile)
	data1=np.loadtxt(sys.argv[i],delimiter=',',skiprows=9, usecols=3, max_rows=1) 
	data.append(data1)
	i+=1
print(data)
mean=np.mean(data)
sta=np.std(data)/mt.sqrt(3)
err_rel=sta/mean
print(mean,sta)
mean=np.array2string(mean) #array2string, to write it in a file
sta=np.array2string(sta) #array2string, to write it in a file
err_rel=np.array2string(err_rel) #array2string, to write it in a file
op=open(outputfile,'a') #open the file you want it in 
op.write(mean) #write the string in the new line, first the average
op.write(',') # then a comma as a delimiter
op.write(sta) # the the absolur error
op.write(',') # then a comma as a delimiter
op.write(err_rel) # then the relative error
op.write('\n') # get a new line for the next argument
op.close() 

