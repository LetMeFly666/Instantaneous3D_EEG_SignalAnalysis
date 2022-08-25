'''
Author: LetMeFly
Date: 2022-08-25 10:34:12
LastEditors: LetMeFly
LastEditTime: 2022-08-25 16:40:26
'''
import matplotlib.pyplot as plt
from Txt2numpy import txt2numpy
from Data2IMFs import data2IMFs
from IMFs2FrequencyDomain import IMFs2FrequencyDomain
from CutoffNoice import cutoffNoice
from IMF2TimeDomain import IMF2TimeDomain
from ConstructEEG import ConstructEEG

import os
# 取消QT警告
os.environ["QT_DEVICE_PIXEL_RATIO"] = "0"
os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
os.environ["QT_SCREEN_SCALE_FACTORS"] = "1"
os.environ["QT_SCALE_FACTOR"] = "1"

data = txt2numpy()  # 从文本文件读入数据
# TODO: show original EEG figure
IMFs = data2IMFs(data)  # 将EEG数据转换为IMF
IMFs = IMFs2FrequencyDomain(IMFs=IMFs)  # 将IMF转换到频域
IMFs = cutoffNoice(IMFs=IMFs)
IMFs = IMF2TimeDomain(IMFs=IMFs)
EEG = ConstructEEG(IMFs=IMFs)

plt.show()

