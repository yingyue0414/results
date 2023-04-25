#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Raw surface point data & parameters
# =====

# In[2]:


timeStep = 0.1 #us
diffConst = 0.001 #um^2/us
numStep = 10000
samplingRate = 100000/timeStep
membraneBendingConst = 4.17

df = pd.read_csv("surfacepoint.csv", header = None)
df = df.iloc[:, :-1]
df


# Extract $h(x,y)$ 
# ====

# In[3]:


import matplotlib.pyplot as plt
df_z = df.iloc[:,2::3]
df_z.columns = range(df_z.columns.size)
df_z


# h(x,y)~t
# =

# In[4]:


#plot h(x,y)~t
#pick (x,y,z)*
index = 245

timeScale = np.linspace(0, numStep*timeStep, numStep+1)

plt.plot(timeScale, df.iloc[:,index])
plt.xlabel("Time(us)")
plt.ylabel("Displacement on z-axis(um)")
plt.show()


# In[5]:


#plot h(x,y)~t
#pick (x,y,z)*
index = 2

timeScale = np.linspace(0, numStep*timeStep, numStep+1)

plt.plot(timeScale, df.iloc[:,index])
plt.xlabel("Time(us)")
plt.ylabel("Displacement on z-axis(um)")
plt.show()


# The periodic BC only allows for drifting of the membrane but not rotating. This means that the displacement can be largely correlated to membrane drifting over long period of simulation. To eliminate the effect of overall drifitng, side length is used for FFT.

# l(x1,y1,x2,y2)~t
# =

# In[6]:


#plot l(x1,y1,x2,y2)~t
#first 500 us
index1 = 241
index2 = 242

timeScale = np.linspace(0, numStep*timeStep/100, int(numStep/100)+1)

def getPointTrimmed(df, index):
    df_pt = df.iloc[0:int(numStep/100)+1,index*3:index*3+3]
    df_pt.columns = range(df_pt.columns.size)
    return df_pt

df_disp = (np.power(getPointTrimmed(df, index1) - getPointTrimmed(df, index2),2)).sum(axis=1)-1

plt.plot(timeScale, df_disp)
plt.xlabel("Time(us)")
plt.ylabel("Length of mesh side(nm)")
plt.show()


# In[7]:


#fft of h~t
from scipy.fft import fft

def getPoint(df, index):
    df_pt = df.iloc[0:int(numStep)+1,index*3:index*3+3]
    df_pt.columns = range(df_pt.columns.size)
    return df_pt

#calculate time vs vibration
timeScale = np.linspace(0, timeStep * numStep, numStep + 1) #0.5 sec, 500 samples
vibrationOutput = (np.power(getPoint(df, index1) - getPoint(df, index2),2)).sum(axis=1)-1

vibrationOutput -= vibrationOutput.mean()
fftResult = fft(np.array(vibrationOutput))
fftResult


# FFT~t
# =

# In[8]:


shiftedResult = np.abs(np.fft.fftshift(fftResult))[:numStep // 2] * 1 / numStep #second half is only imaging of first half
#fft result
plt.plot((samplingRate * np.linspace(0,4999,5000)/numStep), shiftedResult)
plt.xlabel("Freq (Hz)")
plt.ylabel("Arbitrary Amplitude")
plt.show()


# In[9]:


#plot z over given y
df_slice = df.iloc[5000:7000:500,102*3:152*3]
df_slice.columns = range(df_slice.columns.size)
df_slice


# FFT(h(x) given t) ~ y
# =

# In[10]:


#average of FFT results
fftStep = 10

df_slice = df.iloc[0:10000:fftStep, 102*3:152*3]
df_slice.columns = range(df_slice.columns.size)
#fft of z~y
df_slice_z = df_slice.iloc[:,2::3]
df_slice_z_shifted = df_slice_z.sub(df_slice_z.mean(axis=1), axis=0)
df_slice_z_shifted.columns = range(df_slice_z_shifted.columns.size)
df_slice_z_shifted


# In[11]:


fftResult = fft(df_slice_z_shifted)
shiftedResult = np.abs(fftResult) * 1 / numStep #second half is only imaging of first half
avgfft = np.zeros(len(fftResult[0]) // 2)
for i in range(0,1000):
#fft result
    avgfft += shiftedResult[i][:len(fftResult[0]) // 2]
avgfft /= 1000
avgfft


# In[12]:


plt.plot((np.linspace(0,24,25)), avgfft)
plt.xlabel("Average y-Step vibration period")
plt.ylabel("Arbitrary Amplitude")
plt.show()


# Resolving Translation
# ====
# use the difference in z between a given point and the average of mesh points to get rid of the effect of translation of the whole membrane
# 

# In[13]:


z_mean = df_z.mean(axis=1)
z_mean


# FFT2
# ===
# Plot the relation ship between h(q) and |q|

# $$L = \text{side length of membrane}$$
# $$(x,y) = \text{position of surface points}$$
# 
# $$\vec{q}={L}\begin{pmatrix}{n_{x}}\\{n_{y}}\\\end{pmatrix}$$
# 
# Therefore,
# 
# $$\vec{q}^2=\frac{1}{L^2}({n_{x}^2}+{n_{y}^2})$$

# In[32]:


#h(q) over |q| (using FFT2)
x_scale = 41 # n_x
y_scale = 47 # n_y
time = 5000
L = 40

# extract x y info to np.array at time t
def extract_xy_at_t(df, time, x_scale, y_scale, zmean = 0.0):
    arr_xy = np.zeros([y_scale,x_scale])
    for y in range(0, y_scale):
        for x in range(0, x_scale):
            arr_xy[y][x] = df.iloc[time,(x + y * x_scale)*3+2] - zmean
    return arr_xy

# vector_q = (1/L * row, 1/L * col)
# q2 = ((1/L * row) ** 2 + (1/L * col) ** 2) ** 4 * pi ^ 2
# 4pi^2  = 39.4784176044
def extract_q2(nrow, ncol): #l in um
    arr_q2 = np.zeros([nrow,ncol])
    for row in range(0, nrow):
        for col in range(0, ncol):
            arr_q2[row][col] = ((1/L * row) ** 2 + (1/L * col) ** 2)# 4pi^2
    return arr_q2

arr_xy = extract_xy_at_t(df, time, x_scale, y_scale) #z matrices at x,y
#arr_xy
arr_fft2 = np.abs(np.fft.fft2(arr_xy))
arr_fft2_trimmed = arr_fft2[0:len(arr_fft2)//2][0:len(arr_fft2[0])//2]# rest 3/4 of the matrices are just mirror image
arr_fft2_trimmed 


# In[33]:


arr_sqr_fft2 = (np.power(arr_fft2_trimmed, 2)).ravel() #flatten for plotting
arr_sqr_fft2


# In[34]:


qspace_len = len(arr_sqr_fft2) # get number of points in |q|
arr_sqr_fft2_sum = np.zeros(qspace_len)
numitr = 3000
for time in range(3000,3000+numitr): # average h(q)h*(q) over given time scale: <h(q)h*(q)>
    arr_xy = extract_xy_at_t(df, time, x_scale, y_scale, zmean = z_mean[time]) # resolve translation
    arr_fft2 = np.abs(np.fft.fft2(arr_xy))
    arr_fft2_trimmed = arr_fft2[0:len(arr_fft2)//2][0:len(arr_fft2[0])//2]
    arr_sqr_fft2 = (np.power(arr_fft2_trimmed, 2)).ravel() # For complex number h(q): h(q)h*(q) = |h(q)|^2
    arr_sqr_fft2_sum += arr_sqr_fft2
    if (time % 200 == 0):
        print("Complete:" + str(time))
arr_sqr_fft2_average = arr_sqr_fft2_sum/numitr
arr_sqr_fft2_average


# In[35]:


#flatten q2 for plotting
q2 = extract_q2(len(arr_fft2_trimmed),len(arr_fft2_trimmed[0]))
arr_q2_flatten = q2.ravel()
arr_q2_flatten


# In[36]:


#plot in q space
plt.scatter(np.power(arr_q2_flatten, 0.5), arr_sqr_fft2_average) # take sqrt of q^2 to get |q|
plt.ylabel('$<h(q)h^*(q)> (nm^3)$')
plt.xlabel('$|q| (nm^{-1})$')
plt.xlim([0,2])


# FF2 - Trim out the periodic region
# ===

# In[51]:


#h(q) over |q| (using FFT2)
x_scale = 41
y_scale = 47
time = 5000
L = (x_scale + y_scale * 0.86602540378)/2 -5

# extract x y info to np.array at time t
def extract_xy_at_t(df, time, x_scale, y_scale, zmean = 0.0):
    arr_xy = np.zeros([y_scale,x_scale])
    for y in range(0, y_scale):
        for x in range(0, x_scale):
            arr_xy[y][x] = df.iloc[time,(x + y * x_scale)*3+2] - zmean
    return arr_xy

# vector_q = (1/L * row, 1/L * col)
# q2 = (1/L * row) ** 2 + (1/L * col) ** 2
def extract_q2(nrow, ncol): #l in um
    arr_q2 = np.zeros([nrow,ncol])
    for row in range(0, nrow):
        for col in range(0, ncol):
            arr_q2[row][col] = ((1/L * row) ** 2 + (1/L * col) ** 2) * 39.4784176044
    return arr_q2


# In[52]:


arr_periodic_trimmed = []
for y in range(0,y_scale-4):
    for x in range(0,x_scale-4):
        for i in range(0,3):
            arr_periodic_trimmed.append((x+y*x_scale)*3+i)
arr_periodic_trimmed


# In[53]:


df_periodic_trimmed = df.iloc[:,arr_periodic_trimmed]
df_periodic_trimmed.columns = range(0,len(df_periodic_trimmed.columns))
df_periodic_trimmed


# In[54]:


qspace_len = 666 # get number of points in |q|
arr_sqr_fft2_sum = np.zeros(qspace_len)
numitr = 2000
for time in range(5000,5000+numitr): # average h(q)h*(q) over given time scale: <h(q)h*(q)>
    arr_xy = extract_xy_at_t(df_periodic_trimmed, time, x_scale-4, y_scale-4, zmean = z_mean[time]) # resolve translation
    arr_fft2 = np.abs(np.fft.fft2(arr_xy))
    arr_fft2_trimmed = arr_fft2[0:len(arr_fft2)//2][0:len(arr_fft2[0])//2]
    arr_sqr_fft2 = (np.power(arr_fft2_trimmed, 2)).ravel() # For complex number h(q): h(q)h*(q) = |h(q)|^2
    arr_sqr_fft2_sum += arr_sqr_fft2
    if (time % 200 == 0):
        print("Complete:" + str(time))
arr_sqr_fft2_average = arr_sqr_fft2_sum/numitr
arr_sqr_fft2_average


# In[55]:


#flatten q2 for plotting
q2 = extract_q2(len(arr_fft2_trimmed),len(arr_fft2_trimmed[0]))
arr_q2_flatten = q2.ravel()
arr_q2_flatten


# In[56]:


len(arr_q2_flatten)


# In[57]:


#plot in q space
plt.scatter(np.power(arr_q2_flatten, 0.5), arr_sqr_fft2_average) # take sqrt of q^2 to get |q|
plt.ylabel('$<h(q)h^*(q)> (nm^3)$')
plt.xlabel('$|q| (nm^{-1})$')


# In[59]:


#plot in q space
x = np.power(arr_q2_flatten, -1)[1:]
y = arr_sqr_fft2_average[1:]
m, b = np.polyfit(x, y, 1)
plt.plot(x, m*x+b, "r")
plt.scatter(x, y) # take sqrt of q^2 to get |q|
plt.ylabel('$<h(q)h^*(q)> (nm^3)$')
plt.xlabel('$|q|^{-4} (nm^{4})$')
print(m)


# In[60]:


#binned statistic
from scipy.stats import binned_statistic

x = np.power(arr_q2_flatten, 0.5)[1:]
y = arr_sqr_fft2_average[1:]
mean_stat = binned_statistic(x, y, 
                             statistic='mean', 
                             bins=30)

print(mean_stat.statistic)
# array([0.198,   nan, 0.28 , 0.355, 0.265])
print(mean_stat.bin_edges)
# array([0. , 0.5, 1. , 1.5, 2. , 2.5])
#mean_stat.binnumber
# array([1, 1, 1, ..., 4, 5, 5])


# In[61]:


def getBinMidpoint(bin_edges):
    step = bin_edges[1] - bin_edges[0]
    print(step)
    return (bin_edges + step)[0:len(bin_edges)-1]


# In[62]:


#plot in q space
plt.plot(getBinMidpoint(mean_stat.bin_edges), mean_stat.statistic) # take sqrt of q^2 to get |q|
plt.ylabel('$<h(q)h^*(q)> (nm^3)$')
plt.xlabel('$|q| (nm^{-1})$')


# In[48]:


df_s = pd.DataFrame({"q":getBinMidpoint(mean_stat.bin_edges),"hq":mean_stat.statistic})


# In[49]:


df_s


# In[50]:


df_s.to_pickle("4040.pkl")


# In[ ]:




