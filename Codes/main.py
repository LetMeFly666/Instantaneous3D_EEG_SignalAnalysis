'''
Author: LetMeFly
Date: 2022-08-25 10:34:12
LastEditors: LetMeFly
LastEditTime: 2022-08-25 14:39:28
'''
from Txt2numpy import txt2numpy
from Data2IMFs import data2IMFs
from CutoffNoice import cutoffNoice
from IMF2TimeDomain import IMF2TimeDomain
from ConstructEEG import ConstructEEG

data = txt2numpy()
IMFs = data2IMFs(data)
IMFs = cutoffNoice(IMFs=IMFs)
IMFs = IMF2TimeDomain(IMFs=IMFs)
EEG = ConstructEEG(IMFs=IMFs)

