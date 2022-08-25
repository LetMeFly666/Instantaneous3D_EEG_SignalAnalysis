'''
Author: LetMeFly
Date: 2022-08-25 20:21:32
LastEditors: LetMeFly
LastEditTime: 2022-08-25 21:58:18
'''
# %matplotlib inline
from scipy.signal import hilbert
import numpy as np
from Visualize import showEEG

fs = 38 - 2

def HHT(EEG):
    analytic_signal = hilbert(EEG)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    instantaneousFrequency = (np.diff(instantaneous_phase) / (2.0 * np.pi) * fs)
    # 上面instantaneousFrequency.shape = (3599,)
    instantaneousFrequency = np.append(instantaneousFrequency, 0)
    print(instantaneousFrequency)
    print(instantaneousFrequency.shape)

    showEEG(instantaneousFrequency, 2, 38, "HHT")

    return instantaneousFrequency