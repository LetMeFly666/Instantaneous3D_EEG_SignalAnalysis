'''
Author: LetMeFly
Date: 2022-08-25 13:34:22
LastEditors: LetMeFly
LastEditTime: 2022-08-29 16:33:32
'''
import numpy as np
import matplotlib.pyplot as plt
from Visualize import showIMFs, showEEG
from BaseClass import Data

def ConstructEEG(IMFs: Data) -> Data:
    # 丢掉国小的IMF
    # 后面当我没说 # 这里比论文中小改进的一点是：此方法不是绝对定义的切除几秒，而是自动计算出应该切掉的边缘部分
    # start = 0
    # end = 40
    # for i in range(len(IMFs)):
    #     avg = IMFs[i].mean()  # 均值
    #     var = IMFs[i].var()  # 方差
    #     print(avg, var)

    # start = IMFs.getStartTime() + 2
    # end = IMFs.getEndTime() - 2
    # IMFs.startTime = start
    # IMFs.endTime = end
    # IMFs.data = IMFs.data[:, 2 * IMFs.getFPS() : IMFs.getDataLength() - 2 * IMFs.getFPS()]
    # IMFs.dataLength = IMFs.getDataLength()
    
    IMFs.setSubDataByPoints(200, 200)

    showIMFs(IMFs, "Remove edge effects data")

    # 丢弃过小的IMF
    IMFs.data = IMFs.data[:6]
    showIMFs(IMFs, "Desert small IMFs")



    EEG = Data(IMFs.data.sum(axis=0), IMFs.getStartTime(), IMFs.getFPS())
    showEEG(EEG, "Construct EEG")
    # EEG = EEG[1000:]
    # ax = plt.subplot()
    # ax.set_title('Construct EEG')
    # ax.plot(np.arange(0, end - start, 0.01), EEG)
    # plt.show()
    return EEG
