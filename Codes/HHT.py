'''
Author: LetMeFly
Date: 2022-08-25 20:21:32
LastEditors: LetMeFly
LastEditTime: 2022-08-25 21:40:43
'''
# %matplotlib inline
from scipy.signal import hilbert

def HHT(EEG):
    instantaneousFrequency = hilbert(EEG)
    print(instantaneousFrequency)
    print(instantaneousFrequency.shape)

    

    return instantaneousFrequency