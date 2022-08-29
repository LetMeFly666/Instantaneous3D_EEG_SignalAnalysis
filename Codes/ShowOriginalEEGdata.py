'''
Author: LetMeFly
Date: 2022-08-25 17:01:48
LastEditors: LetMeFly
LastEditTime: 2022-08-29 10:48:20
'''
import numpy as np
from matplotlib import pyplot as plt
from BaseClass import Data

def showOriginalEEGdata(EEG: Data) -> None:
    ax = plt.subplot()
    ax.set_title('Original EEG')
    ax.plot(EEG.getTimeRangeArray(), EEG.getData())
    plt.show()
