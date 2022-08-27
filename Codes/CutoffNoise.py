'''
Author: LetMeFly
Date: 2022-08-25 10:50:54
LastEditors: LetMeFly
LastEditTime: 2022-08-27 21:50:36
'''
import numpy as np
from matplotlib import pyplot as plt
from Visualize import showIMFs

def cutoffNoise(IMFs):
    # return IMFs  # Try result
    IMFs[IMFs < 0.5] = 0
    IMFs[IMFs > 32] = 0
    showIMFs(IMFs, 0, 40, "Cut off noise")
    return IMFs
