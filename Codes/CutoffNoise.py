'''
Author: LetMeFly
Date: 2022-08-25 10:50:54
LastEditors: LetMeFly
LastEditTime: 2022-08-29 10:52:58
'''
import numpy as np
from matplotlib import pyplot as plt
from Visualize import showIMFs
from BaseClass import Data

def cutoffNoise(IMFs: Data) -> Data:
    return IMFs  # Try result
    IMFs.data[IMFs < 0.5] = 0
    IMFs.data[IMFs > 32] = 0
    showIMFs(IMFs, "Cut off noise")
    return IMFs
