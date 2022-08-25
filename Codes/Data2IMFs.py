'''
Author: LetMeFly
Date: 2022-08-24 20:49:20
LastEditors: LetMeFly
LastEditTime: 2022-08-25 10:36:04
'''
import numpy as np
from PyEMD import EMD, Visualisation
import os


def data2IMFs(data):
    os.environ["QT_DEVICE_PIXEL_RATIO"] = "0"
    os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
    os.environ["QT_SCREEN_SCALE_FACTORS"] = "1"
    os.environ["QT_SCALE_FACTOR"] = "1"

    t = np.arange(0, 40, 0.01)
    print(t)
    print(t.shape)
    # S = np.sin(13*t + 0.2*t**1.4) - np.cos(3*t)
    print(data.shape)

    # Extract imfs and residue
    # In case of EMD
    emd = EMD()
    emd.emd(data)
    imfs, res = emd.get_imfs_and_residue()

    # In general:
    #components = EEMD()(S)
    #imfs, res = components[:-1], components[-1]

    vis = Visualisation()
    vis.plot_imfs(imfs=imfs, residue=res, t=t, include_residue=True)
    vis.plot_instant_freq(t, imfs=imfs)
    vis.show()

