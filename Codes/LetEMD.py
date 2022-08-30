import numpy as np
# https://blog.csdn.net/The_Time_Runner/article/details/100882454
# https://blog.csdn.net/weixin_39781323/article/details/115808594
from scipy.interpolate import interp1d
import BaseFunction


class EMD():
    def __init__(self):
        self.nbsym = 2
        self.dtype = np.float64
        self.IMFs = None
        self.residue = None
    
    def __call__(self, S):
        return self.emd(S, T=None, max_imf=-1)

    def emd(self, S, T=None, max_imf: int=-1):
        T = np.arange(0, len(S), dtype=S.dtype)
        T = BaseFunction.formatTimeArray(T)
        S, T = BaseFunction.ChangeToSameType(S, T)
        oldIMF = np.nan
        IMFNum = 0
        EXTNum = -1
        self.dtype = S.dtype
        residue = S.astype(self.dtype)
        IMF = np.zeros(len(S), dtype=self.dtype)
        IMF2 = np.empty((IMFNum, len(S)))
        
        done = False
        while not done:
            residue[:] = S - np.sum(IMF2[:IMFNum], axis=0)
            IMF = residue.copy()
            mean = np.zeros(len(S), dtype=self.dtype)
            n = 0
            while True:
                n += 1
                if n >= 1000:  # 设置为最多执行1000次
                    print("数据过大或发生了死循环！")
                    break
                residueEXT = BaseFunction.peekDetection(T, IMF)
                maxLoc, minLoc, indexer = residueEXT[0], residueEXT[2], residueEXT[4]
                EXTNum = len(minLoc) + len(maxLoc)
                nzm = len(indexer)
                if EXTNum > 2:
                    envMax, envMin, eMax, eMin = self.extract_max_min_spline(T, IMF)
                    mean[:] = 0.5 * (envMax + envMin)
                    oldIMF = IMF.copy()
                    IMF[:] = IMF - mean
                    # -------------
                    residueEXT = BaseFunction.peekDetection(T, IMF)
                    maxLoc, temp, minLoc, temp, indexer_ = residueEXT
                    EXTNum = len(maxLoc) + len(minLoc)
                    nzm = len(indexer_)
                    if oldIMF is np.nan:
                        continue
                    if self.check_imf(IMF, oldIMF, eMax, eMin) and abs(EXTNum - nzm) < 2:
                        break
                else:
                    done = True
                    break
            IMF2 = np.vstack((IMF2, IMF.copy()))
            IMFNum += 1
            if self.timeToStop(S, IMF2) or IMFNum == max_imf:
                done = True
                break
        if EXTNum <= 2:
            IMF2 = IMF2[:-1]
        self.IMFs = IMF2.copy()
        self.residue = S - np.sum(self.IMFs, axis=0)
        if not np.allclose(self.residue, 0):
            IMF2 = np.vstack((IMF2, self.residue))
        return IMF2

    def extract_max_min_spline(self, T: np.ndarray, S: np.ndarray):
        
        ext_res = BaseFunction.peekDetection(T, S)
        max_pos, max_val = ext_res[0], ext_res[1]
        min_pos, min_val = ext_res[2], ext_res[3]

        if len(max_pos) + len(min_pos) < 3:
            return [-1] * 4
        max_extrema, min_extrema = self._prepare_points_parabol(T, S, max_pos, max_val, min_pos, min_val)

        _, max_spline = self.cubicSpline(T, max_extrema)
        _, min_spline = self.cubicSpline(T, min_extrema)

        return max_spline, min_spline, max_extrema, min_extrema

    def _prepare_points_parabol(self, T, S, max_pos, max_val, min_pos, min_val):
        """
        Performs mirroring on signal which extrema do not necessarily
        belong on the position array.

        See :meth:`EMD.prepare_points`.
        """

        # Need at least two extrema to perform mirroring
        max_extrema = np.zeros((2, len(max_pos)), dtype=self.dtype)
        min_extrema = np.zeros((2, len(min_pos)), dtype=self.dtype)

        max_extrema[0], min_extrema[0] = max_pos, min_pos
        max_extrema[1], min_extrema[1] = max_val, min_val

        # Local variables
        nbsym = self.nbsym
        end_min, end_max = len(min_pos), len(max_pos)

        # Left bound
        d_pos = max_pos[0] - min_pos[0]
        left_ext_max_type = d_pos < 0  # True -> max, else min

        # Left extremum is maximum
        if left_ext_max_type:
            if (S[0] > min_val[0]) and (np.abs(d_pos) > (max_pos[0] - T[0])):
                # mirror signal to first extrema
                expand_left_max_pos = 2 * max_pos[0] - max_pos[1 : nbsym + 1]
                expand_left_min_pos = 2 * max_pos[0] - min_pos[0:nbsym]
                expand_left_max_val = max_val[1 : nbsym + 1]
                expand_left_min_val = min_val[0:nbsym]
            else:
                # mirror signal to beginning
                expand_left_max_pos = 2 * T[0] - max_pos[0:nbsym]
                expand_left_min_pos = 2 * T[0] - np.append(T[0], min_pos[0 : nbsym - 1])
                expand_left_max_val = max_val[0:nbsym]
                expand_left_min_val = np.append(S[0], min_val[0 : nbsym - 1])

        # Left extremum is minimum
        else:
            if (S[0] < max_val[0]) and (np.abs(d_pos) > (min_pos[0] - T[0])):
                # mirror signal to first extrema
                expand_left_max_pos = 2 * min_pos[0] - max_pos[0:nbsym]
                expand_left_min_pos = 2 * min_pos[0] - min_pos[1 : nbsym + 1]
                expand_left_max_val = max_val[0:nbsym]
                expand_left_min_val = min_val[1 : nbsym + 1]
            else:
                # mirror signal to beginning
                expand_left_max_pos = 2 * T[0] - np.append(T[0], max_pos[0 : nbsym - 1])
                expand_left_min_pos = 2 * T[0] - min_pos[0:nbsym]
                expand_left_max_val = np.append(S[0], max_val[0 : nbsym - 1])
                expand_left_min_val = min_val[0:nbsym]

        if not expand_left_min_pos.shape:
            expand_left_min_pos, expand_left_min_val = min_pos, min_val
        if not expand_left_max_pos.shape:
            expand_left_max_pos, expand_left_max_val = max_pos, max_val

        expand_left_min = np.vstack((expand_left_min_pos[::-1], expand_left_min_val[::-1]))
        expand_left_max = np.vstack((expand_left_max_pos[::-1], expand_left_max_val[::-1]))

        # Right bound
        d_pos = max_pos[-1] - min_pos[-1]
        right_ext_max_type = d_pos > 0

        # Right extremum is maximum
        if not right_ext_max_type:
            if (S[-1] < max_val[-1]) and (np.abs(d_pos) > (T[-1] - min_pos[-1])):
                # mirror signal to last extrema
                idx_max = max(0, end_max - nbsym)
                idx_min = max(0, end_min - nbsym - 1)
                expand_right_max_pos = 2 * min_pos[-1] - max_pos[idx_max:]
                expand_right_min_pos = 2 * min_pos[-1] - min_pos[idx_min:-1]
                expand_right_max_val = max_val[idx_max:]
                expand_right_min_val = min_val[idx_min:-1]
            else:
                # mirror signal to end
                idx_max = max(0, end_max - nbsym + 1)
                idx_min = max(0, end_min - nbsym)
                expand_right_max_pos = 2 * T[-1] - np.append(max_pos[idx_max:], T[-1])
                expand_right_min_pos = 2 * T[-1] - min_pos[idx_min:]
                expand_right_max_val = np.append(max_val[idx_max:], S[-1])
                expand_right_min_val = min_val[idx_min:]

        # Right extremum is minimum
        else:
            if (S[-1] > min_val[-1]) and len(max_pos) > 1 and (np.abs(d_pos) > (T[-1] - max_pos[-1])):
                # mirror signal to last extremum
                idx_max = max(0, end_max - nbsym - 1)
                idx_min = max(0, end_min - nbsym)
                expand_right_max_pos = 2 * max_pos[-1] - max_pos[idx_max:-1]
                expand_right_min_pos = 2 * max_pos[-1] - min_pos[idx_min:]
                expand_right_max_val = max_val[idx_max:-1]
                expand_right_min_val = min_val[idx_min:]
            else:
                # mirror signal to end
                idx_max = max(0, end_max - nbsym)
                idx_min = max(0, end_min - nbsym + 1)
                expand_right_max_pos = 2 * T[-1] - max_pos[idx_max:]
                expand_right_min_pos = 2 * T[-1] - np.append(min_pos[idx_min:], T[-1])
                expand_right_max_val = max_val[idx_max:]
                expand_right_min_val = np.append(min_val[idx_min:], S[-1])

        if not expand_right_min_pos.shape:
            expand_right_min_pos, expand_right_min_val = min_pos, min_val
        if not expand_right_max_pos.shape:
            expand_right_max_pos, expand_right_max_val = max_pos, max_val

        expand_right_min = np.vstack((expand_right_min_pos[::-1], expand_right_min_val[::-1]))
        expand_right_max = np.vstack((expand_right_max_pos[::-1], expand_right_max_val[::-1]))

        max_extrema = np.hstack((expand_left_max, max_extrema, expand_right_max))
        min_extrema = np.hstack((expand_left_min, min_extrema, expand_right_min))

        return max_extrema, min_extrema

    """
    三次样条
    """
    def cubicSpline(self, T: np.ndarray, extrema: np.ndarray):
        t = T[np.r_[T >= extrema[0, 0]] & np.r_[T <= extrema[0, -1]]]
        if extrema.shape[1] > 3:  # 大于三个点，使用内置库
            return t, interp1d(extrema[0], extrema[1], kind="cubic")(t)
        else:  # 否则还得手动实现一波
            return BaseFunction.MyCubicSpline(extrema[0], extrema[1], t)

    """
    是否停止分解
    对于剩下的数据，最大值和最小值之差小于0.001 或 绝对值之和小于0.005时，停止分解
    """
    def timeToStop(self, S: np.ndarray, IMF: np.ndarray) -> bool:
        remain = S - np.sum(IMF, axis=0)
        if np.max(remain) - np.min(remain) < 0.001:  # 这里设置为0.001 #--------------------
            return True
        if np.sum(np.abs(remain)) < 0.005:  # 这里设置为0.005 #--------------------
            return True
        return False

    """
    如果连续筛选不影响信号，则信号为IMF
    """
    def check_imf(self, newIMF: np.ndarray, oldIMF: np.ndarray, eMax: np.ndarray, eMin: np.ndarray) -> bool:
        if np.any(eMax[1] < 0) or np.any(eMin[1] > 0):
            return False

        if np.sum(newIMF ** 2) < 1e-10:
            return False

        IMFDiff = newIMF - oldIMF
        imf_diff_sqrd_sum = np.sum(IMFDiff ** 2)

        svar = imf_diff_sqrd_sum / (max(oldIMF) - min(oldIMF))
        if svar < 0.001:  # svar_thr，定义为0.001
            return True

        # Standard deviation test
        std = np.sum((IMFDiff / newIMF) ** 2)
        if std < 0.2:  # std_thr，定义为0.2
            return True

        energy_ratio = imf_diff_sqrd_sum / np.sum(oldIMF * oldIMF)
        if energy_ratio < 0.2:  # energy_ratio_thr，定义为0.2
            return True

        return False

    def getIMFsAndResidue(self):
        return self.IMFs, self.residue
