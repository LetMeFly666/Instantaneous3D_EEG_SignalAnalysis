'''
Author: LetMeFly
Date: 2022-08-25 10:50:54
LastEditors: LetMeFly
LastEditTime: 2022-08-29 20:48:14
'''
import numpy as np
from matplotlib import pyplot as plt
from Visualize import showIMFs
from BaseClass import Data

def cutoffNoise(IMFs: Data) -> Data:
    # return IMFs  # Try result
    # originalData = IMFs.data.copy()
    originalMaxVal = []
    resizeRate = 35
    for th, IMF in enumerate(IMFs.data):
        maxVal = IMF.max()
        originalMaxVal.append(maxVal)
        IMFs.data[th] = IMF * resizeRate / maxVal
    IMFs.data[IMFs.data < 0.5] = 0
    IMFs.data[IMFs.data > 32] = 0
    showIMFs(IMFs, "Cut off noise")
    for th, IMF in enumerate(IMFs.data):
        IMFs.data[th] = IMF * originalMaxVal[th] / resizeRate
    # IMFs.data = originalData
    return IMFs
