'''
Author: LetMeFly
Date: 2022-08-25 10:50:54
LastEditors: LetMeFly
LastEditTime: 2022-08-29 20:38:17
'''
import numpy as np
from matplotlib import pyplot as plt
from Visualize import showIMFs
from BaseClass import Data

def cutoffNoise(IMFs: Data) -> Data:
    # return IMFs  # Try result
    originalData = IMFs.data.copy()
    for th, IMF in enumerate(IMFs.data):
        maxVal = IMF.max()
        IMFs.data[th] = IMF * 32 / maxVal
    IMFs.data[IMFs.data < 0.5] = 0
    IMFs.data[IMFs.data > 32] = 0
    showIMFs(IMFs, "Cut off noise")
    IMFs.data = originalData
    return IMFs
