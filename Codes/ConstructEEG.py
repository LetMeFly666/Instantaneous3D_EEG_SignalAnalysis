'''
Author: LetMeFly
Date: 2022-08-25 13:34:22
LastEditors: LetMeFly
LastEditTime: 2022-08-25 14:42:34
'''
import numpy as np
import matplotlib.pyplot as plt
def ConstructEEG(IMFs):
    EEG = IMFs.sum(axis=0)
    EEG = EEG[1000:]
    print(EEG)
    print(EEG.shape)
    ax = plt.subplot()
    ax.set_title('Construct EEG')
    ax.plot(np.arange(0, 30, 0.01), EEG)
    plt.show()
    return EEG
    