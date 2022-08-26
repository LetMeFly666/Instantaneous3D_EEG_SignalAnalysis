'''
Author: LetMeFly
Date: 2022-08-25 20:02:22
LastEditors: LetMeFly
LastEditTime: 2022-08-26 22:00:34
'''
import numpy as np
from matplotlib import pyplot as plt

def showIMFs(IMFs, startTime, endTime, title="IMFs") -> None:
    fig, axes = plt.subplots(IMFs.shape[0], 1)
    if IMFs.shape[0] == 1:
        axes = list(axes)
    axes[0].set_title(title)
    for num, IMF in enumerate(IMFs):
        ax = axes[num]        
        ax.plot(np.arange(0, endTime - startTime, 0.01), IMF)
        ax.set_ylabel("IMF " + str(num + 1))
    plt.show()


""" show EEG or something that has only 1 graph to show """
def showEEG(EEG, startTime, endTime, title="EEG") -> None:
    ax = plt.subplot()
    ax.set_title(title)
    ax.plot(np.arange(0, endTime - startTime, 0.01), EEG)
    plt.show()


def show3d(x, y, z):
    x_time = np.arange(2, 38, 0.01)
    y_A = x
    z_F = y

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    im = ax.scatter(x_time, y_A, z_F, s=1)  # s: 点的大小
    ax.set_xlabel("x_time")
    ax.set_ylabel("y_A")
    ax.set_zlabel("z_F")
    plt.show()
