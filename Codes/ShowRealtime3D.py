'''
Author: LetMeFly
Date: 2022-08-26 21:46:23
LastEditors: LetMeFly
LastEditTime: 2022-08-26 22:00:23
'''
from Visualize import show3d
import numpy as np

def showRealtime3D(EEG, instantaneousFrequency):
    t = np.arange(2, 38, 0.01)
    show3d(EEG, instantaneousFrequency, t)
