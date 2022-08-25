'''
Author: LetMeFly
Date: 2022-08-25 10:34:12
LastEditors: LetMeFly
LastEditTime: 2022-08-25 19:24:19
'''
import matplotlib.pyplot as plt
from Txt2numpy import txt2numpy
from ShowOriginalEEGdata import showOriginalEEGdata
from Data2IMFs import data2IMFs
from IMFs2FrequencyDomain import IMFs2FrequencyDomain
from CutoffNoise import cutoffNoise
from IMF2TimeDomain import IMF2TimeDomain
from ConstructEEG import ConstructEEG

import os
# 取消QT警告
os.environ["QT_DEVICE_PIXEL_RATIO"] = "0"
os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
os.environ["QT_SCREEN_SCALE_FACTORS"] = "1"
os.environ["QT_SCALE_FACTOR"] = "1"

EEG = txt2numpy()  # 从文本文件读入数据
showOriginalEEGdata(EEG)  # 显示原本EEG文件
IMFs = data2IMFs(EEG)  # step1.将EEG数据转换为IMF
IMFs = IMFs2FrequencyDomain(IMFs=IMFs)  # step2.将IMF转换到频域
IMFs = cutoffNoise(IMFs=IMFs)  # step3.去除每个IMF中的噪声
IMFs = IMF2TimeDomain(IMFs=IMFs)
EEG = ConstructEEG(IMFs=IMFs)

plt.show()

