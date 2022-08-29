'''
Author: LetMeFly
Date: 2022-08-25 10:28:32
LastEditors: LetMeFly
LastEditTime: 2022-08-29 10:38:54
'''
import numpy as np
from BaseClass import Data

def txt2numpy(file="../Data/case7.txt") -> Data:
# def txt2numpy(file="../Data/CA001F9S_1-1+.txt"):
    arr = np.loadtxt(file)
    return Data(arr, 0, 100)
