'''
Author: LetMeFly
Date: 2022-08-25 13:13:28
LastEditors: LetMeFly
LastEditTime: 2022-08-25 19:49:01
'''
import numpy as np
import matplotlib.pyplot as plt
def IMF2TimeDomain(IMFs):
    for i in range(len(IMFs)):
        IMFs[i] = np.fft.ifft(IMFs[i])
    # ax = plt.subplot()
    # ax.set_title('Back to time domain')
    # for thisIMF in IMFs:
    #     ax.plot(np.arange(0, 40, 0.01), thisIMF)
    # plt.show()
    
    fig, axes = plt.subplots(IMFs.shape[0], 1)
    if IMFs.shape[0] == 1:
        axes = list(axes)
    axes[0].set_title("Back to time domain")
    for num, IMF in enumerate(IMFs):
        ax = axes[num]        
        ax.plot(np.arange(0, 40, 0.01), IMF)
        ax.set_ylabel("IMF " + str(num + 1))
    
    plt.show()
    return IMFs
