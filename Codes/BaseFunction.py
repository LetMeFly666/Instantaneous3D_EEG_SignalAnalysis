'''
Author: LetMeFly
Date: 2022-08-30 09:44:58
LastEditors: LetMeFly
LastEditTime: 2022-08-30 20:20:52
'''
import numpy as np
from scipy.interpolate import Akima1DInterpolator

def ChangeToSameType(x: np.ndarray, y: np.ndarray):
        dtype = np.find_common_type([x.dtype, y.dtype], [])
        if x.dtype != dtype:
            x = x.astype(dtype)
        if y.dtype != dtype:
            y = y.astype(dtype)
        return x, y

def formatTimeArray(t):
        """
        规范化时间数组，使其不会在微小值上爆炸。返回的数组从0开始，最小增量为1。
        """
        d = np.diff(t)
        return (t - t[0]) / np.min(d)


# https://blog.csdn.net/bodybo/article/details/77335129
def cubic_spline_3pts(x, y, T):
    """scipy.interpolate.interp1d不支持少于4个点的三次样条曲线"""
    x0, x1, x2 = x
    y0, y1, y2 = y

    x1x0, x2x1 = x1 - x0, x2 - x1
    y1y0, y2y1 = y1 - y0, y2 - y1
    _x1x0, _x2x1 = 1.0 / x1x0, 1.0 / x2x1

    m11, m12, m13 = 2 * _x1x0, _x1x0, 0
    m21, m22, m23 = _x1x0, 2.0 * (_x1x0 + _x2x1), _x2x1
    m31, m32, m33 = 0, _x2x1, 2.0 * _x2x1

    v1 = 3 * y1y0 * _x1x0 * _x1x0
    v3 = 3 * y2y1 * _x2x1 * _x2x1
    v2 = v1 + v3

    M = np.array([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])
    v = np.array([v1, v2, v3]).T
    k = np.array(np.linalg.inv(M).dot(v))

    a1 = k[0] * x1x0 - y1y0
    b1 = -k[1] * x1x0 + y1y0
    a2 = k[1] * x2x1 - y2y1
    b2 = -k[2] * x2x1 + y2y1

    t = T[np.r_[T >= x0] & np.r_[T <= x2]]
    t1 = (T[np.r_[T >= x0] & np.r_[T < x1]] - x0) / x1x0
    t2 = (T[np.r_[T >= x1] & np.r_[T <= x2]] - x1) / x2x1
    t11, t22 = 1.0 - t1, 1.0 - t2

    q1 = t11 * y0 + t1 * y1 + t1 * t11 * (a1 * t11 + b1 * t1)
    q2 = t22 * y1 + t2 * y2 + t2 * t22 * (a2 * t22 + b2 * t2)
    q = np.append(q1, q2)

    return t, q


def peekDetection(T: np.ndarray, S: np.ndarray):
    """极值检测"""

    # Finds indexes of zero-crossings
    S1, S2 = S[:-1], S[1:]
    indzer = np.nonzero(S1 * S2 < 0)[0]
    if np.any(S == 0):
        indz = np.nonzero(S == 0)[0]
        if np.any(np.diff(indz) == 1):
            zer = S == 0
            dz = np.diff(np.append(np.append(0, zer), 0))
            debz = np.nonzero(dz == 1)[0]
            finz = np.nonzero(dz == -1)[0] - 1
            indz = np.round((debz + finz) / 2.0)

        indzer = np.sort(np.append(indzer, indz))

    # Finds local extrema
    d = np.diff(S)
    d1, d2 = d[:-1], d[1:]
    indmin = np.nonzero(np.r_[d1 * d2 < 0] & np.r_[d1 < 0])[0] + 1
    indmax = np.nonzero(np.r_[d1 * d2 < 0] & np.r_[d1 > 0])[0] + 1

    # When two or more points have the same value
    if np.any(d == 0):

        imax, imin = [], []

        bad = d == 0
        dd = np.diff(np.append(np.append(0, bad), 0))
        debs = np.nonzero(dd == 1)[0]
        fins = np.nonzero(dd == -1)[0]
        if debs[0] == 1:
            if len(debs) > 1:
                debs, fins = debs[1:], fins[1:]
            else:
                debs, fins = [], []

        if len(debs) > 0:
            if fins[-1] == len(S) - 1:
                if len(debs) > 1:
                    debs, fins = debs[:-1], fins[:-1]
                else:
                    debs, fins = [], []

        lc = len(debs)
        if lc > 0:
            for k in range(lc):
                if d[debs[k] - 1] > 0:
                    if d[fins[k]] < 0:
                        imax.append(np.round((fins[k] + debs[k]) / 2.0))
                else:
                    if d[fins[k]] > 0:
                        imin.append(np.round((fins[k] + debs[k]) / 2.0))

        if len(imax) > 0:
            indmax = indmax.tolist()
            for x in imax:
                indmax.append(int(x))
            indmax.sort()

        if len(imin) > 0:
            indmin = indmin.tolist()
            for x in imin:
                indmin.append(int(x))
            indmin.sort()

    local_max_pos = T[indmax]
    local_max_val = S[indmax]
    local_min_pos = T[indmin]
    local_min_val = S[indmin]

    return local_max_pos, local_max_val, local_min_pos, local_min_val, indzer
