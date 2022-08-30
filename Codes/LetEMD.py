import numpy as np
# https://blog.csdn.net/The_Time_Runner/article/details/100882454
from scipy.interpolate import interp1d
import BaseFunction


class EMD():
    def __init__(self):
        self.energy_ratio_thr = 0.2
        self.std_thr = 0.2
        self.svar_thr = 0.001
        self.nbsym = 2
        self.scale_factor = 1.0
        self.DTYPE = np.float64
        self.FIXE = 0
        self.FIXE_H = 0
        self.MAX_ITERATION = 1000
        self.imfs = None
        self.residue = None

    def emd(self, S: np.ndarray, T = None, max_imf: int = -1) -> np.ndarray:
        T = np.arange(0, len(S), dtype=S.dtype)
        T = BaseFunction.formatTimeArray(T)
        S, T = BaseFunction.ChangeToSameType(S, T)
        self.DTYPE = S.dtype
        N = len(S)

        residue = S.astype(self.DTYPE)
        imf = np.zeros(len(S), dtype=self.DTYPE)
        imf_old = np.nan

        imfNo = 0
        extNo = -1
        IMF = np.empty((imfNo, N))
        finished = False

        while not finished:

            residue[:] = S - np.sum(IMF[:imfNo], axis=0)
            imf = residue.copy()
            mean = np.zeros(len(S), dtype=self.DTYPE)

            # Counters
            n = 0  # All iterations for current imf.
            n_h = 0  # counts when |#zero - #ext| <=1

            while True:
                n += 1
                if n >= self.MAX_ITERATION:
                    break

                ext_res = BaseFunction.peekDetection(T, imf)
                max_pos, min_pos, indzer = ext_res[0], ext_res[2], ext_res[4]
                extNo = len(min_pos) + len(max_pos)
                nzm = len(indzer)

                if extNo > 2:

                    max_env, min_env, eMax, eMin = self.extract_max_min_spline(T, imf)
                    mean[:] = 0.5 * (max_env + min_env)

                    imf_old = imf.copy()
                    imf[:] = imf - mean

                    # Fix number of iterations
                    if self.FIXE:
                        if n >= self.FIXE:
                            break

                    # Fix number of iterations after number of zero-crossings
                    # and extrema differ at most by one.
                    elif self.FIXE_H:
                        tmp_residue = BaseFunction.peekDetection(T, imf)
                        max_pos, min_pos, ind_zer = (tmp_residue[0], tmp_residue[2], tmp_residue[4])
                        extNo = len(max_pos) + len(min_pos)
                        nzm = len(ind_zer)

                        if n == 1:
                            continue

                        n_h = n_h + 1 if abs(extNo - nzm) < 2 else 0

                        if n_h >= self.FIXE_H:
                            break

                    else:
                        ext_res = BaseFunction.peekDetection(T, imf)
                        max_pos, _, min_pos, _, ind_zer = ext_res
                        extNo = len(max_pos) + len(min_pos)
                        nzm = len(ind_zer)

                        if imf_old is np.nan:
                            continue

                        f1 = self.check_imf(imf, imf_old, eMax, eMin)
                        f2 = abs(extNo - nzm) < 2

                        if f1 and f2:
                            break

                else:
                    finished = True
                    break

            IMF = np.vstack((IMF, imf.copy()))
            imfNo += 1

            if self.end_condition(S, IMF) or imfNo == max_imf:
                finished = True
                break

        if extNo <= 2:
            IMF = IMF[:-1]

        self.imfs = IMF.copy()
        self.residue = S - np.sum(self.imfs, axis=0)

        if not np.allclose(self.residue, 0):
            IMF = np.vstack((IMF, self.residue))

        return IMF


    def __call__(self, S):
        return self.emd(S, T=None, max_imf=-1)

    def extract_max_min_spline(self, T: np.ndarray, S: np.ndarray):
        
        ext_res = BaseFunction.peekDetection(T, S)
        max_pos, max_val = ext_res[0], ext_res[1]
        min_pos, min_val = ext_res[2], ext_res[3]

        if len(max_pos) + len(min_pos) < 3:
            return [-1] * 4
        max_extrema, min_extrema = self._prepare_points_parabol(T, S, max_pos, max_val, min_pos, min_val)

        _, max_spline = self.spline_points(T, max_extrema)
        _, min_spline = self.spline_points(T, min_extrema)

        return max_spline, min_spline, max_extrema, min_extrema

    def _prepare_points_parabol(self, T, S, max_pos, max_val, min_pos, min_val):
        """
        Performs mirroring on signal which extrema do not necessarily
        belong on the position array.

        See :meth:`EMD.prepare_points`.
        """

        # Need at least two extrema to perform mirroring
        max_extrema = np.zeros((2, len(max_pos)), dtype=self.DTYPE)
        min_extrema = np.zeros((2, len(min_pos)), dtype=self.DTYPE)

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

    def spline_points(self, T: np.ndarray, extrema: np.ndarray):
        t = T[np.r_[T >= extrema[0, 0]] & np.r_[T <= extrema[0, -1]]]
        if extrema.shape[1] > 3:  # 大于三个点，使用内置库
            return t, interp1d(extrema[0], extrema[1], kind="cubic")(t)
        else:  # 否则还得手动实现一波
            return BaseFunction.cubic_spline_3pts(extrema[0], extrema[1], t)

    def end_condition(self, S: np.ndarray, IMF: np.ndarray) -> bool:
        # 对于剩下的数据，最大值和最小值之差小于0.001 或 绝对值之和小于0.005时，停止分解
        remain = S - np.sum(IMF, axis=0)

        if np.max(remain) - np.min(remain) < 0.001:  # 这里设置为0.001 #--------------------
            return True

        if np.sum(np.abs(remain)) < 0.005:  # 这里设置为0.005 #--------------------
            return True

        return False

    def check_imf(self, imf_new: np.ndarray, imf_old: np.ndarray, eMax: np.ndarray, eMin: np.ndarray) -> bool:
        """
        Huang criteria for **IMF** (similar to Cauchy convergence test).
        Signal is an IMF if consecutive siftings do not affect signal
        in a significant manner.
        """
        if np.any(eMax[1] < 0) or np.any(eMin[1] > 0):
            return False

        if np.sum(imf_new ** 2) < 1e-10:
            return False

        # Precompute values
        imf_diff = imf_new - imf_old
        imf_diff_sqrd_sum = np.sum(imf_diff * imf_diff)

        # Scaled variance test
        svar = imf_diff_sqrd_sum / (max(imf_old) - min(imf_old))
        if svar < self.svar_thr:
            return True

        # Standard deviation test
        std = np.sum((imf_diff / imf_new) ** 2)
        if std < self.std_thr:
            return True

        energy_ratio = imf_diff_sqrd_sum / np.sum(imf_old * imf_old)
        if energy_ratio < self.energy_ratio_thr:
            return True

        return False

    def getIMFsAndResidue(self):
        return self.imfs, self.residue
