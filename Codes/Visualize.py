'''
Author: LetMeFly
Date: 2022-08-25 20:02:22
LastEditors: LetMeFly
LastEditTime: 2022-08-25 20:09:53
'''
import numpy as np
from matplotlib import pyplot as plt

def showIMFs(IMFs, startTime, endTime, title="IMFs"):
    fig, axes = plt.subplots(IMFs.shape[0], 1)
    if IMFs.shape[0] == 1:
        axes = list(axes)
    axes[0].set_title(title)
    for num, IMF in enumerate(IMFs):
        ax = axes[num]        
        ax.plot(np.arange(0, endTime - startTime, 0.01), IMF)
        ax.set_ylabel("IMF " + str(num + 1))
    plt.show()
