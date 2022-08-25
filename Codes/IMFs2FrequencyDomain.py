'''
Author: LetMeFly
Date: 2022-08-25 16:25:58
LastEditors: LetMeFly
LastEditTime: 2022-08-25 19:25:13
'''
import numpy as np
import matplotlib.pyplot as plt
def IMFs2FrequencyDomain(IMFs):
    for i in range(len(IMFs)):
        IMFs[i] = np.fft.fft(IMFs[i])
    fig, axes = plt.subplots(IMFs.shape[0], 1)
    if IMFs.shape[0] == 1:
        axes = list(axes)
    axes[0].set_title("Change IMFs into frequency domain")
    for num, IMF in enumerate(IMFs):
        axes[num].plot(np.arange(0, 40, 0.01), IMF)
        axes[num].set_ylabel("IMF " + str(num + 1))
    return IMFs
