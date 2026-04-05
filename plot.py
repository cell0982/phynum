import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np

df = pd.read_csv('mirage_data.csv')

fig = plt.figure()
ax = fig.add_subplot(projection='3d')


ax.plot(df['x'], df['y'], df['z'], c='blue')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
plt.figure()
#plt.plot(df['x'], df['y'])
plt.plot(df['z'], df['x'])
