'''
Author: LetMeFly
Date: 2022-08-25 20:21:32
LastEditors: LetMeFly
LastEditTime: 2022-08-29 11:03:40
'''
# %matplotlib inline
from scipy.signal import hilbert
import numpy as np
from Visualize import showEEG
from BaseClass import Data

def HHT(EEG: Data) -> Data:
    analytic_signal = hilbert(EEG.data)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    instantaneousFrequency = (np.diff(instantaneous_phase) / (2.0 * np.pi) * (EEG.getEndTime() - EEG.getStartTime()))  # TODO: 这里是否发生了错误
    # 上面instantaneousFrequency.shape = (3599,)
    instantaneousFrequency = np.append(instantaneousFrequency, 0)
    print(instantaneousFrequency)
    print(instantaneousFrequency.shape)

    instantaneousFrequency = Data(instantaneousFrequency, EEG.getStartTime(), EEG.getFPS())
    showEEG(instantaneousFrequency, "HHT")

    return instantaneousFrequency