'''
Author: LetMeFly
Date: 2022-08-25 10:34:12
LastEditors: LetMeFly
LastEditTime: 2022-08-26 21:59:42
'''
from Txt2numpy import txt2numpy
from ShowOriginalEEGdata import showOriginalEEGdata
from Data2IMFs import data2IMFs
from IMFs2FrequencyDomain import IMFs2FrequencyDomain
from CutoffNoise import cutoffNoise
from IMF2TimeDomain import IMF2TimeDomain
from ConstructEEG import ConstructEEG
from HHT import HHT
from ShowRealtime3D import showRealtime3D

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
IMFs = IMF2TimeDomain(IMFs=IMFs)  # step4.转回到时域
EEG = ConstructEEG(IMFs=IMFs)  # step5.选取有效范围的IMF并构建EEG
instantaneousFrequency = HHT(EEG)  # step6.通过HHT获得实时频率
showRealtime3D(EEG, instantaneousFrequency)  # step7.显示为实时三维图

