'''
Author: LetMeFly
Date: 2022-08-25 17:01:48
LastEditors: LetMeFly
LastEditTime: 2022-08-25 19:47:57
'''
import numpy as np
from matplotlib import pyplot as plt

def showOriginalEEGdata(EEG) -> None:
    ax = plt.subplot()
    ax.set_title('Original EEG')
    ax.plot(np.arange(0, 40, 0.01), EEG)
    plt.show()
