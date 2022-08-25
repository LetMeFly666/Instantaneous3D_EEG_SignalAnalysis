'''
Author: LetMeFly
Date: 2022-08-25 10:28:32
LastEditors: LetMeFly
LastEditTime: 2022-08-25 10:33:27
'''
import numpy as np
def txt2numpy(file="../Data/CA001F9S_1-1+.txt"):
    arr = np.loadtxt(file)
    return arr
