#!/usr/bin/env python
# coding: utf-8

# In[1]:


#things to do 
#1 comment and review

#imports
import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
import tqdm as tq
from scipy.interpolate import interp1d


# In[ ]:





# In[2]:


#downloading the lightcurve file for our example star KIC 10685175
lcfc = lk.search_lightcurvefile("KIC 10685175",mission="Kepler").download_all()
lc = lcfc.PDCSAP_FLUX.stitch().remove_nans()


# In[3]:


#noise threshhold function, determines at what noise levels the frequency crosses the .99 percent correct recovery line
def noise_threshold(time,noise ,frequency , max_frequency, min_frequency = 0, fap = .01, max_runs= 1000):
    
    #creating sinusoidal light curve 
    flux = np.sin(2*np.pi*frequency*time)
    lc = lk.LightCurve(time,flux)
    
    #lightcurve to frequency periodogram
    per = lc.to_periodogram(nyquist_factor = 0.01)
    nyquist = per.nyquist.value
    
    frequency_candidates = []
    
    min_factor = int(np.floor(min_frequency/nyquist))
    max_factor = int(np.ceil(max_frequency/nyquist))

    #picking out which peaks to sample in each nyquist 'zone'
    for i in range(min_factor, max_factor):
        per = lc.to_periodogram(minimum_frequency=i*nyquist,maximum_frequency=(i+1)*nyquist)
        frequency_candidates.append(per.frequency_at_max_power.value)
    
    frequency_candidates = np.array(frequency_candidates)
    frequency_candidates = frequency_candidates[(frequency_candidates > min_frequency)&(frequency_candidates < max_frequency)]
    
    #sampling only near the peaks, increasing effeciency
    frequency_resolution = 1/(time[-1]-time[0])
    n_samples = 41
    local_sampling = np.linspace(-.1*frequency_resolution,.1*frequency_resolution,n_samples)
    
    frequency_sample = np.zeros(n_samples*len(frequency_candidates))

    for i in range(len(frequency_candidates)):
        frequency_sample[i*n_samples:(i+1)*n_samples] = local_sampling + frequency_candidates[i]
        
    results = np.zeros(max_runs)
    max_incorrect_runs = (fap*max_runs)
    incorrect_results = 0
    percentage = 0     
    for i in tq.tqdm(range(max_runs)):
        flux = np.sin(2*np.pi*(frequency*time + np.random.rand())) + np.random.randn(len(time))*noise

        lc = lk.LightCurve(time,flux)

        per = lc.to_periodogram(frequency=frequency_sample,ls_method="slow")

        frequency_dif = np.abs(per.frequency_at_max_power.value - frequency)

        results[i] = (frequency_dif < frequency_resolution)
        
        if(frequency_dif > frequency_resolution):
                incorrect_results += 1
        if(incorrect_results > max_incorrect_runs):
            break
    
        percentage = np.sum(results)/(i+1)
        
    return percentage


# In[4]:


lc = lcfc.PDCSAP_FLUX.stitch().remove_nans()
time = lc.time
print(time)

print(noise_threshold(time,60.3, 289.39094, 300, min_frequency = 50, max_runs=100))


# In[ ]:




