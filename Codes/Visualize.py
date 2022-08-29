'''
Author: LetMeFly
Date: 2022-08-25 20:02:22
LastEditors: LetMeFly
LastEditTime: 2022-08-29 21:05:00
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


# https://www.jb51.net/article/258747.htm
def show3d(F: Data, A: Data, ColorfulAndResize=False):
    x_F = F
    y_time = A.getTimeRangeArray()
    z_A = A

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if ColorfulAndResize:
        colors = []
        for f in x_F.data:
            if f < 0.5:
                colors.append("red")
            elif f < 4:
                colors.append("blue")
            elif f < 8:
                colors.append("yellow")
            elif f < 16:
                colors.append("green")
            elif f < 32:
                colors.append("pink")
            else:
                colors.append("red")
        im = ax.scatter(np.abs(x_F.data), y_time, np.abs(z_A.data), c=colors, s=1)  # s: 点的大小
        ax.set_zlim(0, 1000)
        ax.set_title("Colorful & Resize")
    else:
        im = ax.scatter(np.abs(x_F.data), y_time, np.abs(z_A.data), s=1)  # s: 点的大小
        

    ax.set_xlabel("x_F")
    ax.set_ylabel("y_time")
    ax.set_zlabel("z_A")
    plt.show()
