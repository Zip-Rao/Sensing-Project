'''
创建信号类，定义传感中的虚拟信号，可用于构建脉冲序列，或模拟外部磁场信号
支持：零信号，恒定信号，正弦信号，高斯脉冲信号，非对称脉冲信号等，返回列表形式的信号序列
'''
import numpy as np
from qutip import *

class Signal:
    def __init__(self, type = 0, t_list=None, **kwargs):
        '''
        初始化Signal对象
        :param type: 信号类型，0-零信号，1-恒定信号，2-正弦信号，3-高斯脉冲信号 ，4-非对称脉冲信号 , 5-slepian脉冲信号 
        :param t_list: 时间列表（ns），用于生成信号序列
        :param kwargs: 信号参数，可以包括幅度，频率，相位，中心位置，宽度，直流偏置等
        '''
        self.type = type
        self.t_list = t_list
        self.params = kwargs

        self.init()

        self.signal = self.generate()
    def init(self):
        '''
        根据type初始化信号参数
        '''
        default_params = {
            'amplitude': 1.0,
            'frequency': 0.1,   #(GHZ)
            'phase': 0.0,       #(rad)
            'center': 0.0,      #(ns)
            'width': 1.0,       #(ns)
            'offset': 0.0,      #(DC offset)
            't_rabi':np.linspace(0,40,100), #(ns)
            'tau':20.0,           #(ns)
            'rise':5.0,         #(ns)
            'fall':2.0,         #(ns)
        }
        # 补全参数，用默认值代替
        for key, value in default_params.items():
            if key not in self.params:
                self.params[key] = value
    
    def generate(self):
        '''
        在给定时序上生成信号序列
        '''
        t_array = np.array(self.t_list)
        
        signal = self.back_signal()

        match self.type:
            case 0: # 零信号
                return signal
            case 1: # 恒定信号
                signal += self.params['amplitude'] * np.ones_like(t_array)
                return signal+self.params['offset']
            case 2: # 正弦信号
                omega = 2 * np.pi * self.params['frequency']
                signal += self.params['amplitude'] * np.sin(omega * t_array + self.params['phase'])
                return signal+self.params['offset']
            case 3: # 高斯脉冲信号
                sigma = self.params['width'] / 4
                signal += self.params['amplitude'] * np.exp(-((t_array - self.params['center']) ** 2) / (2 * sigma ** 2))
                return signal+self.params['offset']
            case 4: # 冲击信号
                rise = self.params['rise']
                fall = self.params['fall']
                t0 = self.params['center']
                signal += self.params['amplitude'] * np.where(t_array < t0, 
                                                              0, 
                                                              (1 - np.exp(-(t_array - t0)/rise)) * np.exp(-(t_array - t0)/fall))
                return signal+self.params['offset']
            case 5: # slepian脉冲信号
                pass



    def back_signal(self, t_list=None):
        '''
        背景信号（可加噪声）
        '''
        if t_list is None:
            t_list = self.t_list
        return np.zeros_like(t_list)
    
    def value_at(self, t):
        '''
        获取信号在时间t处的值，如果t不在t_list外，则返回0
        :param t: 时间点（ns）
        '''
        signal = self.signal
        t_array = np.array(self.t_list)
        idx = (np.abs(t_array - t)).argmin()
        if t_array[0] <= t <= t_array[-1]:  # 间隔小于0.1ns则认为在范围内，视具体的t_list间隔修改
            return signal[idx]
        elif t < t_array[0] or t > t_array[-1]:
            return 0.0
    
    def plot(self):
        '''
        绘制信号波形
        '''
        import matplotlib.pyplot as plt

        signal = self.signal
        plt.figure(figsize=(10,4))
        plt.plot(self.t_list, signal)
        plt.xlabel('Time (ns)')
        plt.ylabel('Signal Amplitude')
        plt.title('Signal Waveform')
        plt.grid(True)
        plt.show()
        

