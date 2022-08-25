'''
Author: LetMeFly
Date: 2022-08-25 10:34:12
LastEditors: LetMeFly
LastEditTime: 2022-08-25 10:36:31
'''
from Txt2numpy import txt2numpy
from Data2IMFs import data2IMFs
data = txt2numpy()
data2IMFs(data)
