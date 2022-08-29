'''
Author: LetMeFly
Date: 2022-08-25 20:02:22
LastEditors: LetMeFly
LastEditTime: 2022-08-29 11:12:59
'''
import numpy as np
from matplotlib import pyplot as plt
from BaseClass import Data

def showIMFs(IMFs: Data, title="IMFs") -> None:
    fig, axes = plt.subplots(IMFs.data.shape[0], 1)
    if IMFs.data.shape[0] == 1:
        axes = list(axes)
    axes[0].set_title(title)
    for num, IMF in enumerate(IMFs.data):
        ax = axes[num]        
        ax.plot(IMFs.getTimeRangeArray(), IMF)
        ax.set_ylabel("IMF " + str(num + 1))
    plt.show()


""" show EEG or something that has only 1 graph to show """
def showEEG(EEG: Data, title="EEG") -> None:
    ax = plt.subplot()
    ax.set_title(title)
    ax.plot(EEG.getTimeRangeArray(), EEG.data)
    plt.show()


def show3d(x: Data, y: Data):
    x_A = x
    y_time = x.getTimeRangeArray()
    z_F = y

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    im = ax.scatter(x_A.data, y_time, z_F.data, s=1)  # s: 点的大小
    ax.set_xlabel("x_A")
    ax.set_ylabel("y_time")
    ax.set_zlabel("z_F")
    plt.show()
