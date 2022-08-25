'''
Author: LetMeFly
Date: 2022-08-25 10:50:54
LastEditors: LetMeFly
LastEditTime: 2022-08-25 20:13:51
'''
import numpy as np
from matplotlib import pyplot as plt
def cutoffNoise(IMFs):
    return IMFs  # Try result
    IMFs[IMFs < 0.5] = 0
    IMFs[IMFs > 32] = 0
    fig, axes = plt.subplots(IMFs.shape[0], 1)
    if IMFs.shape[0] == 1:
        axes = list(axes)
    axes[0].set_title("Cut off noise")
    for num, IMF in enumerate(IMFs):
        axes[num].plot(np.arange(0, 40, 0.01), IMF)
        axes[num].set_ylabel("IMF " + str(num + 1))
    plt.show()
    return IMFs
