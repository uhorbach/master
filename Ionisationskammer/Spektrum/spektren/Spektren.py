import numpy as np
import sys
import matplotlib, matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker


files = sys.argv[1:]

for file in files:
  print("Processing file: ", file)
  data=np.loadtxt(file,skiprows=0)

  plt.bar(data[:,0], data[:,1])
  plt.xlabel('Energy in keV')
  plt.ylabel('Counts')
  plt.savefig(file.replace('.txt', '.png'))
  plt.close()

