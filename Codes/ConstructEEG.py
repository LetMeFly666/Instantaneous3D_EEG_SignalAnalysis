'''
Author: LetMeFly
Date: 2022-08-25 13:34:22
LastEditors: LetMeFly
LastEditTime: 2022-08-25 20:12:04
'''
import numpy as np
import matplotlib.pyplot as plt
from ShowIMFs import showIMFs

def ConstructEEG(IMFs):
    # 丢掉国小的IMF
    # 后面当我没说 # 这里比论文中小改进的一点是：此方法不是绝对定义的切除几秒，而是自动计算出应该切掉的边缘部分
    # start = 0
    # end = 40
    # for i in range(len(IMFs)):
    #     avg = IMFs[i].mean()  # 均值
    #     var = IMFs[i].var()  # 方差
    #     print(avg, var)
    start = 2
    end = 40 - 2
    IMFs = IMFs[:, start * 100 : end * 100]
    showIMFs(IMFs, start, end, "Remove edge effects data")

    # 丢弃过小的IMF
    IMFs = IMFs[:6]
    showIMFs(IMFs, start, end, "Desert small IMFs")



    EEG = IMFs.sum(axis=0)
    # EEG = EEG[1000:]
    print(EEG)
    print(EEG.shape)
    ax = plt.subplot()
    ax.set_title('Construct EEG')
    ax.plot(np.arange(0, end - start, 0.01), EEG)
    plt.show()
    return EEG
    