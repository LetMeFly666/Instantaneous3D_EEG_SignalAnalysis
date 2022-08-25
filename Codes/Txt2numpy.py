'''
Author: LetMeFly
Date: 2022-08-25 10:28:32
LastEditors: LetMeFly
LastEditTime: 2022-08-25 16:06:22
'''
import numpy as np
def txt2numpy(file="../Data/case7.txt"):
    arr = np.loadtxt(file)
    return arr
