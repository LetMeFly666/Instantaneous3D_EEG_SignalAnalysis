'''
Author: LetMeFly
Date: 2022-08-26 21:46:23
LastEditors: LetMeFly
LastEditTime: 2022-08-29 21:05:31
'''
from Visualize import show3d
import numpy as np
from BaseClass import Data

def showRealtime3D(EEG: Data, instantaneousFrequency: Data) -> None:
    show3d(instantaneousFrequency, EEG)
    show3d(instantaneousFrequency, EEG, ColorfulAndResize=True)
