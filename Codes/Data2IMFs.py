'''
Author: LetMeFly
Date: 2022-08-24 20:49:20
LastEditors: LetMeFly
LastEditTime: 2022-08-29 21:19:10
'''
import numpy as np
from LetEMD import EMD  # , Visualisation
from matplotlib import pyplot as plt
from BaseClass import Data


def data2IMFs(data: Data):

    # t = np.arange(0, 40, 0.01)
    # print(t)
    # print(t.shape)
    # print(data.shape)

    # Extract imfs and residue
    # In case of EMD
    emd = EMD()
    emd.emd(data.getData())
    imfs, res = emd.getIMFsAndResidue()

    # In general:
    #components = EEMD()(S)
    #imfs, res = components[:-1], components[-1]

    fig, axes = plt.subplots(imfs.shape[0] + 1, 1)
    if imfs.shape[0] == 1:
        axes = list(axes)
    axes[0].set_title("The IMFs")
    for num, IMF in enumerate(imfs):
        axes[num].plot(data.getTimeRangeArray(), IMF)
        axes[num].set_ylabel("IMF " + str(num + 1))
    axes[imfs.shape[0]].plot(data.getTimeRangeArray(), res)
    axes[imfs.shape[0]].set_ylabel("Residue")

    plt.show()

    # vis = Visualisation()
    # vis.plot_imfs(imfs=imfs, residue=res, t=t, include_residue=True)
    # vis.plot_instant_freq(t, imfs=imfs)
    # vis.show()

    return Data(imfs, data.getStartTime(), data.getFPS())

