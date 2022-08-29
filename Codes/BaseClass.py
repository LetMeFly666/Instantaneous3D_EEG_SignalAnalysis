'''
Author: LetMeFly
Date: 2022-08-29 09:39:46
LastEditors: LetMeFly
LastEditTime: 2022-08-29 10:46:30
'''
import numpy as np


class Data():
    """
    数据类
    ====
        支持一维数据和二维数据
        ----
            一维数据：[1, 2, 1, 2, 3, ...]    Example: data.shape = (500, )
            二维数据：[[1, 2, 1, ..], [2, 6, 3, ..], ..]  Example: data.shape = (5, 500)
    
    成员
    ====
        变量：
        ----
            data: np.ndarray，主要数据
            dataLength: 数据的长度。这里时指单条数据共采样了多少次
            startTime: 数据的开始采样时间     --+
            endTime: 数据的结束采样时间       --.\__--> [startTime, endTime)
            FPS: 每秒采样多少次
        
        方法：
        ----
            getData() -> np.ndarray:  获取采样数据
            getDataLength(data: np.ndarray) -> int:  获取数据的长度（采样次数）
            getStartTime() -> float: 获取数据的开始采样时间
            getEndTime() -> float: 获取数据的结束采样时间
            getFPS() -> int: 获取数据的采样频率
            getTimeRangeArray() -> np.ndarray:  获取时间的数组（数组中为所有的采样时刻）
    
    构造函数的参数
    ====
        startTime
        endTime
        FPS
    """
    
    def __init__(self, data: np.ndarray, startTime: float, FPS: int) -> None:
        self.data = data
        self.dataLength = self.getDataLength()
        self.startTime = startTime
        self.endTime = self.startTime + self.dataLength / FPS
        self.FPS = FPS
    
    def getDataLength(self) -> int:
        data = self.data
        if len(data.shape) == 1:
            return data.shape[0]
        else:
            return data.shape[1]
    
    def getData(self) -> np.ndarray:
        return self.data
    
    def getStartTime(self) -> float:
        return self.startTime
    
    def getEndTime(self) -> float:
        return self.endTime
    
    def getFPS(self) -> int:
        return self.FPS
    
    def getTimeRangeArray(self) -> int:
        # TODO: 是否会有由于精度问题导致的BUG
        return np.arange(self.getStartTime(), self.getEndTime(), 1 / self.getFPS())


if __name__ == "__main__":
    data1 = Data(np.array([1, 2, 3, 5, 7, 8, 9, 9, 0, 2]), 1, 5)
    print("data1.data:", data1.getData())
    print("data1.dataLength:", data1.getDataLength())
    print(f"data1.timeRange: [{data1.getStartTime()}, {data1.getEndTime()})", )
    print("data1.FPS:", data1.getFPS())
    print("data1.timeRangeArray:", data1.getTimeRangeArray())

    print("-" * 50)

    data2 = Data(np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [10, 9, 8, 7, 6, 5, 4, 3, 2, 1], [5, 6, 7, 8, 9, 0, 1, 2, 3, 4]]), 1.5, 10)
    print("data2.data:", data2.getData())
    print("data2.dataLength:", data2.getDataLength())
    print(f"data2.timeRange: [{data2.getStartTime()}, {data2.getEndTime()})", )
    print("data2.FPS:", data2.getFPS())
    print("data2.timeRangeArray:", data2.getTimeRangeArray())
