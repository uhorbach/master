
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit




data_C40=np.loadtxt('C40.txt', delimiter=',')
data_C60=np.loadtxt('C60.txt', delimiter=',')
data_C80=np.loadtxt('C80.txt', delimiter=',')
data_C100=np.loadtxt('C100.txt', delimiter=',')
data_C150=np.loadtxt('C150.txt', delimiter=',')
data_C200=np.loadtxt('C200.txt', delimiter=',')
data_C250=np.loadtxt('C250.txt', delimiter=',')
data_662=np.loadtxt('662.txt', delimiter=',')



particle=np.linspace(5000000,25000000,5) #number of particles

Luftkerma=np.array([particle,
                data_C40[:,0],
                data_C60[:,0],
                data_C80[:,0],
                data_C100[:,0],
                data_C150[:,0],
                data_C200[:,0],
                data_C250[:,0],
                data_662[:,0]])

label=np.array(['C40','C60','C80','C100','C150','C200','C250','662'])


# In[4]:


#Farben
r=np.linspace(1,0,12)
g=np.linspace(0,1,12)
b=np.linspace(1,0,12)
rgb=np.array([r,g,b])



# In[8]:


m_ges = np.empty(shape=[8])
m_err = np.empty(shape=[8])
m_rel = np.empty(shape=[8])

i=0
while i<=7:
    x=particle
    y=Luftkerma[i+1,:]
    def f(x, m):
        y = m * x
        return y

    params, covariance_matrix = curve_fit(f, x, y)
    m = params
    m_ges[i]=m
    m_err[i]=np.sqrt(covariance_matrix[0][0])
    m_rel[i]=m_err[i]/m_ges[i]

    plt.plot(x, y,"x", color=rgb[:,i])
    plt.plot(x, m*x, label=label[i], color=rgb[:,i])
    plt.legend(loc="upper left")
    plt.xlabel('AnzahlTeilchen')
    plt.ylabel('Dosis in Gy')
    plt.savefig('Luftkerma.pdf')
    i+=1


# In[26]:


datages=pd.DataFrame({'energy': label,
                     'm': m_ges,
                     'm_err': m_err,
                     'm_rel': m_rel})
print(datages)
datages.to_excel('Luftkerma.xlsx')


