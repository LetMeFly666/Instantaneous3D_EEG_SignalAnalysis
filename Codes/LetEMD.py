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
    
    def getIMFsAndResidue(self):
        return self.IMFs, self.residue

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
                    envMax, envMin, eMax, eMin = self.extractMaxAndMinSpline(T, IMF)
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
                    if self.ifIsIMF(IMF, oldIMF, eMax, eMin) and abs(EXTNum - nzm) < 2:
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

    def extractMaxAndMinSpline(self, T: np.ndarray, S: np.ndarray):
        
        extremaLoc = BaseFunction.peekDetection(T, S)
        maxLoc, maxVal = extremaLoc[0], extremaLoc[1]
        minLoc, minVal = extremaLoc[2], extremaLoc[3]
        if len(maxLoc) + len(minLoc) < 3:
            return [-1] * 4
        maxExtrema, minExtrema = self.preparePoints(T, S, maxLoc, maxVal, minLoc, minVal)

        temp, maxSpline = self.cubicSpline(T, maxExtrema)
        temp, minSpline = self.cubicSpline(T, minExtrema)
        return maxSpline, minSpline, maxExtrema, minExtrema

    """
    镜像操作
    """
    def preparePoints(self, T, S, maxLoc, maxVal, minLoc, minVal):
        # 至少需要两个极值
        maxExtrema = np.zeros((2, len(maxLoc)), dtype=self.dtype)
        minExtrema = np.zeros((2, len(minLoc)), dtype=self.dtype)
        maxExtrema[0], minExtrema[0] = maxLoc, minLoc
        maxExtrema[1], minExtrema[1] = maxVal, minVal
        nbsym = self.nbsym
        minEnd, maxEnd = len(minLoc), len(maxLoc)
        # 左边界
        locD = maxLoc[0] - minLoc[0]
        ifLeftExtremaIsMax = locD < 0  # 如果locD < 0，那么leftExtremaMaxType就为True，就说明为左极值最大值；否则为最小值

        # 左极值为最大值
        if ifLeftExtremaIsMax:
            if (S[0] > minVal[0]) and (np.abs(locD) > (maxLoc[0] - T[0])):
                # 第一个极值
                expandLeftMaxPos = 2 * maxLoc[0] - maxLoc[1 : nbsym + 1]
                expandLeftMinPos = 2 * maxLoc[0] - minLoc[0:nbsym]
                expandLeftMaxVal = maxVal[1 : nbsym + 1]
                expandLeftMinVal = minVal[0:nbsym]
            else:
                # 起始位置
                expandLeftMaxPos = 2 * T[0] - maxLoc[0:nbsym]
                expandLeftMinPos = 2 * T[0] - np.append(T[0], minLoc[0 : nbsym - 1])
                expandLeftMaxVal = maxVal[0:nbsym]
                expandLeftMinVal = np.append(S[0], minVal[0 : nbsym - 1])
        # 左极值为最小值
        else:
            if (S[0] < maxVal[0]) and (np.abs(locD) > (minLoc[0] - T[0])):
                # 第一个极值
                expandLeftMaxPos = 2 * minLoc[0] - maxLoc[0:nbsym]
                expandLeftMinPos = 2 * minLoc[0] - minLoc[1 : nbsym + 1]
                expandLeftMaxVal = maxVal[0:nbsym]
                expandLeftMinVal = minVal[1 : nbsym + 1]
            else:
                # 起始位置
                expandLeftMaxPos = 2 * T[0] - np.append(T[0], maxLoc[0 : nbsym - 1])
                expandLeftMinPos = 2 * T[0] - minLoc[0:nbsym]
                expandLeftMaxVal = np.append(S[0], maxVal[0 : nbsym - 1])
                expandLeftMinVal = minVal[0:nbsym]

        if not expandLeftMinPos.shape:
            expandLeftMinPos, expandLeftMinVal = minLoc, minVal
        if not expandLeftMaxPos.shape:
            expandLeftMaxPos, expandLeftMaxVal = maxLoc, maxVal

        expandLeftMin = np.vstack((expandLeftMinPos[::-1], expandLeftMinVal[::-1]))
        expandLeftMax = np.vstack((expandLeftMaxPos[::-1], expandLeftMaxVal[::-1]))

        # 右边界
        locD = maxLoc[-1] - minLoc[-1]
        rightExtremaMaxYype = locD > 0
        # 右极值是最大值
        if not rightExtremaMaxYype:
            if (S[-1] < maxVal[-1]) and (np.abs(locD) > (T[-1] - minLoc[-1])):
                # 最后一个极值
                maxIndex = max(0, maxEnd - nbsym)
                minIndex = max(0, minEnd - nbsym - 1)
                expandRightMaxLoc = 2 * minLoc[-1] - maxLoc[maxIndex:]
                expandRightMinPos = 2 * minLoc[-1] - minLoc[minIndex:-1]
                expandRightMaxVal = maxVal[maxIndex:]
                expandRightMinVal = minVal[minIndex:-1]
            else:
                # 终止位置
                maxIndex = max(0, maxEnd - nbsym + 1)
                minIndex = max(0, minEnd - nbsym)
                expandRightMaxLoc = 2 * T[-1] - np.append(maxLoc[maxIndex:], T[-1])
                expandRightMinPos = 2 * T[-1] - minLoc[minIndex:]
                expandRightMaxVal = np.append(maxVal[maxIndex:], S[-1])
                expandRightMinVal = minVal[minIndex:]
        # 右极值是最小值
        else:
            if (S[-1] > minVal[-1]) and len(maxLoc) > 1 and (np.abs(locD) > (T[-1] - maxLoc[-1])):
                # 最后一个极值
                maxIndex = max(0, maxEnd - nbsym - 1)
                minIndex = max(0, minEnd - nbsym)
                expandRightMaxLoc = 2 * maxLoc[-1] - maxLoc[maxIndex:-1]
                expandRightMinPos = 2 * maxLoc[-1] - minLoc[minIndex:]
                expandRightMaxVal = maxVal[maxIndex:-1]
                expandRightMinVal = minVal[minIndex:]
            else:
                # 终止位置
                maxIndex = max(0, maxEnd - nbsym)
                minIndex = max(0, minEnd - nbsym + 1)
                expandRightMaxLoc = 2 * T[-1] - maxLoc[maxIndex:]
                expandRightMinPos = 2 * T[-1] - np.append(minLoc[minIndex:], T[-1])
                expandRightMaxVal = maxVal[maxIndex:]
                expandRightMinVal = np.append(minVal[minIndex:], S[-1])

        if not expandRightMinPos.shape:
            expandRightMinPos, expandRightMinVal = minLoc, minVal
        if not expandRightMaxLoc.shape:
            expandRightMaxLoc, expandRightMaxVal = maxLoc, maxVal

        expand_right_min = np.vstack((expandRightMinPos[::-1], expandRightMinVal[::-1]))
        expand_right_max = np.vstack((expandRightMaxLoc[::-1], expandRightMaxVal[::-1]))

        maxExtrema = np.hstack((expandLeftMax, maxExtrema, expand_right_max))
        minExtrema = np.hstack((expandLeftMin, minExtrema, expand_right_min))

        return maxExtrema, minExtrema

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
    def ifIsIMF(self, newIMF: np.ndarray, oldIMF: np.ndarray, eMax: np.ndarray, eMin: np.ndarray) -> bool:
        if np.any(eMax[1] < 0) or np.any(eMin[1] > 0):
            return False
        if np.sum(newIMF ** 2) < 1e-10:
            return False

        IMFDiff = newIMF - oldIMF
        PingFang = np.sum(IMFDiff ** 2)
        svar = PingFang / (max(oldIMF) - min(oldIMF))
        if svar < 0.001:  # svar_thr，定义为0.001
            return True

        # 标准差检验
        std = np.sum((IMFDiff / newIMF) ** 2)
        if std < 0.2:  # std_thr，定义为0.2
            return True
        energyRatio = PingFang / np.sum(oldIMF * oldIMF)
        if energyRatio < 0.2:  # energy_ratio_thr，定义为0.2
            return True
        return False
