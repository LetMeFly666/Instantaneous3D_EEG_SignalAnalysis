'''
Author: LetMeFly
Date: 2022-08-25 10:28:32
LastEditors: LetMeFly
LastEditTime: 2022-08-25 18:59:18
'''
import numpy as np
def txt2numpy(file="../Data/case7.txt"):
# def txt2numpy(file="../Data/CA001F9S_1-1+.txt"):
    arr = np.loadtxt(file)
    return arr
