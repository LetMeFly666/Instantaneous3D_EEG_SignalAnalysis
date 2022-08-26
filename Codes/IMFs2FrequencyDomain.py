'''
Author: LetMeFly
Date: 2022-08-25 16:25:58
LastEditors: LetMeFly
LastEditTime: 2022-08-26 09:34:50
'''
import numpy as np
import matplotlib.pyplot as plt
from Visualize import showIMFs
def IMFs2FrequencyDomain(IMFs):
    for i in range(len(IMFs)):
        IMFs[i] = np.fft.fft(IMFs[i])
    
    showIMFs(IMFs, 0, 40, "Change IMFs into frequency domain")
    return IMFs
