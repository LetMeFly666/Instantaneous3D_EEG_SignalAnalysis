'''
Author: LetMeFly
Date: 2022-08-25 16:25:58
LastEditors: LetMeFly
LastEditTime: 2022-08-29 10:52:29
'''
from scipy.fftpack import fft
import numpy as np
from Visualize import showIMFs
from BaseClass import Data

def IMFs2FrequencyDomain(IMFs: Data):
    for num, IMF in enumerate(IMFs.getData()):
        frequency = fft(IMF)
        magnitude = np.abs(frequency)
        phase = np.angle(frequency)
        IMFs.data[num] = magnitude

    showIMFs(IMFs, "Change IMFs into frequency domain")    
    return IMFs

# import numpy as np
# from Visualize import showIMFs
# def IMFs2FrequencyDomain(IMFs):
#     for i in range(len(IMFs)):
#         IMFs[i] = np.fft.fft(IMFs[i])
    
#     showIMFs(IMFs, 0, 40, "Change IMFs into frequency domain")
#     return IMFs


# from scipy.fftpack import fft
# import numpy as np
# from Visualize import showIMFs

# def IMFs2FrequencyDomain(IMFs):
#     # for num, IMF in enumerate(IMFs):
#     #     IMFs[num] = fft(IMF)
#     # data_freq = IMFs
#     data_freq = fft(IMFs)
#     showIMFs(IMFs, 0, 40, "Try FFT")
#     mdata = np.abs(data_freq)  # magnitude
#     pdata = np.angle(data_freq)  # phase
#     showIMFs(mdata, 0, 40, "magnitude")
#     showIMFs(pdata, 0, 40, "phase")
#     return data_freq

