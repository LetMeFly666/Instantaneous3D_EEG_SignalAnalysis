'''
Author: LetMeFly
Date: 2022-08-24 20:49:20
LastEditors: LetMeFly
LastEditTime: 2022-08-25 16:52:05
'''
import numpy as np
from PyEMD import EMD  # , Visualisation
from matplotlib import pyplot as plt


def data2IMFs(data):

    t = np.arange(0, 40, 0.01)
    # print(t)
    # print(t.shape)
    # print(data.shape)

    # Extract imfs and residue
    # In case of EMD
    emd = EMD()
    emd.emd(data)
    imfs, res = emd.get_imfs_and_residue()

    # In general:
    #components = EEMD()(S)
    #imfs, res = components[:-1], components[-1]

    fig, axes = plt.subplots(imfs.shape[0], 1)
    if imfs.shape[0] == 1:
        axes = list(axes)
    axes[0].set_title("The IMFs.")
    for num, IMF in enumerate(imfs):
        axes[num].plot(np.arange(0, 40, 0.01), IMF)
        axes[num].set_ylabel("IMF " + str(num + 1))


    # vis = Visualisation()
    # vis.plot_imfs(imfs=imfs, residue=res, t=t, include_residue=True)
    # vis.plot_instant_freq(t, imfs=imfs)
    # vis.show()

    return imfs

