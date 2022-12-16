#!/usr/bin/env python
# coding: utf-8

# In[93]:


import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager

dataframe = pd.read_csv("WR_plot_res.csv",sep=',')
#print(dataframe.head(10))

#dataframe.columns  # show DataFrame with dataframe.columns

x = dataframe.dt
y = dataframe.WR 

plt.xlim(0,80)     # set x range
plt.ylim(-15,15)     # set y range

# Set the linewidth 
plt.plot(x, y, linewidth=2)

# Add x and y lables, and set their font size
plt.xlabel("dt (s)", fontsize=10)
plt.ylabel("WR", fontsize=10)
#plt.title("Jacobi_1d_WR_dt")


#top and right show ticks or values
ax = plt.gca()  # get current axes, --> move axes
ax.tick_params(top=True,labeltop=False,right=True,labelright = True)

plt.title("Jacobi_1d_WR_dt")
#plt.scatter(x, y, )
plt.show() 


# In[ ]:




