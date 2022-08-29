'''
Author: LetMeFly
Date: 2022-08-25 13:13:28
LastEditors: LetMeFly
LastEditTime: 2022-08-29 10:55:03
'''
import numpy as np
# import matplotlib.pyplot as plt
from BaseClass import Data
from Visualize import showIMFs

def IMF2TimeDomain(IMFs: Data) -> Data:
    for i in range(len(IMFs.data)):
        IMFs.data[i] = np.fft.ifft(IMFs.data[i])
    
    showIMFs(IMFs, "Back to time domain")
    # ax = plt.subplot()
    # ax.set_title('Back to time domain')
    # for thisIMF in IMFs:
    #     ax.plot(np.arange(0, 40, 0.01), thisIMF)
    # plt.show()
    
    # fig, axes = plt.subplots(IMFs.shape[0], 1)
    # if IMFs.shape[0] == 1:
    #     axes = list(axes)
    # axes[0].set_title("Back to time domain")
    # for num, IMF in enumerate(IMFs):
    #     ax = axes[num]        
    #     ax.plot(np.arange(0, 40, 0.01), IMF)
    #     ax.set_ylabel("IMF " + str(num + 1))
    #
    # plt.show()

    return IMFs
