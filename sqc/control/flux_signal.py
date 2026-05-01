'''
创建信号类，定义传感中的虚拟信号，可用于构建脉冲序列，或模拟外部磁场信号
支持：零信号，恒定信号，正弦信号，高斯脉冲信号，非对称脉冲信号等，返回列表形式的信号序列
'''
import numpy as np
from qutip import *
from sqc.control.waveform import Waveform

class FluxSignal(Waveform):
    def __init__(self, type = 0, t_list=None, **kwargs):
        '''
        初始化Signal对象
        :param type: 信号类型，0-零信号，1-恒定信号，2-正弦信号，3-高斯脉冲信号 ，4-非对称脉冲信号 , 5-双峰信号 ， 6-基展开信号 7-复杂信号 8-用户自定义信号
        :param t_list: 时间列表（ns），用于生成信号序列
        :param kwargs: 信号参数，可以包括幅度，频率，相位，中心位置，宽度，直流偏置等
        '''
        self.type = type
        self.t_list = t_list
        self.params = kwargs

        self.init()
        self.generate_basis_signal()

        self.signal = self.generate()

    @property
    def samples(self):
        """Alias used by the new Waveform-oriented API."""
        return self.signal

    @samples.setter
    def samples(self, value):
        self.signal = np.asarray(value, dtype=float)

    def copy(self):
        '''
        创建当前Signal对象的副本
        '''
        new_signal = Signal(type=self.type, t_list=self.t_list.copy(), **self.params)

        return new_signal
    
    def init(self):
        '''
        根据type初始化信号参数
        '''
        default_params = {
            'amplitude': 1.0,
            'frequency': 0.01,   #(GHZ)
            'phase': 0.0,       #(rad)
            'center': 0.0,      #(ns)
            'width': 30.0,       #(ns)
            'offset': 0.0,      #(DC offset)
            't_rabi':np.linspace(0,40,100), #(ns)
            'tau':20.0,           #(ns)
            'rise':5.0,         #(ns)
            'fall':2.0,         #(ns)
            'noise_level': 0.0,   #(信号叠加的高斯噪声标准差)
            'seed': 42,
        }
        # 补全参数，用默认值代替
        for key, value in default_params.items():
            if key not in self.params:
                self.params[key] = value
        if self.type == 6:
            if "n_basis" not in self.params:
                self.params["n_basis"] = 10  # 默认使用10个基函数
            if "b" not in self.params:
                self.params["b"] = np.ones(self.params["n_basis"])  # 默认权重为1
            if "basis_type" not in self.params:
                self.params["basis_type"] = "bspline"  # 默认使用B样条基函数
    def generate(self):
        '''
        在给定时序上生成信号序列
        '''
        t_array = np.array(self.t_list)
        
        signal = self.back_signal(noise_level=self.params['noise_level'], seed=self.params['seed'], t_list=t_array)

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
                # 光滑的非对称脉冲信号，使用双指数函数
                signal += self.params['amplitude'] * (np.exp(-(t_array - t0)/fall - np.exp(-(t_array - t0)/rise))) 
                # 有突变的非对称脉冲信号
                #signal += self.params['amplitude'] * np.where(t_array < t0, 
                #                                              0, 
                #                                              (1 - np.exp(-(t_array - t0)/rise)) * np.exp(-(t_array - t0)/fall))
                return signal+self.params['offset']
            case 5: # 双峰信号
                sigma = self.params['width'] / 4
                center1 = self.params['center'] - self.params['width'] / 2
                center2 = self.params['center'] + self.params['width'] / 2
                signal += self.params['amplitude'] * (np.exp(-((t_array - center1) ** 2) / (2 * sigma ** 2)) + 
                                                      np.exp(-((t_array - center2) ** 2) / (2 * sigma ** 2)))
                return signal+self.params['offset']
            case 6: # 基展开信号
                b = self.params['b']
                for i, phi_i in enumerate(self.basis_functions):
                    signal += b[i] * phi_i(t_array)
                return signal+self.params['offset']
            case 7: # 复杂信号
                def wave_packet(t, A, t0, w, f, phi):
                    """
                    t: 时间数组
                    A: 振幅
                    t0: 波包中心时间
                    w: 高斯包络的宽度 (标准差)
                    f: 频率
                    phi: 相位
                    """
                    return 0.0001 *A * np.exp(-(t - t0)**2 / (2 * w**2)) * np.cos(2 * np.pi * f * (t - t0) + phi)
                tlist = self.t_list
               # p1: 左侧小波包
                p1 = wave_packet(tlist, A=8, t0=50, w=12, f=0.06, phi=0)

                # p2: 中间主波包 (振幅最大)
                p2 = wave_packet(tlist, A=25, t0=90, w=18, f=0.03, phi=np.pi/4)

                # p3: 右侧次级波包
                p3 = wave_packet(tlist, A=15, t0=120, w=20, f=0.02, phi=np.pi/2)
                signal += p1  + p2  + p3
                return signal+self.params['offset']
            case 8: # 用户自定义信号
                if "signal" not in self.params:
                    raise ValueError("Custom signal not provided.")
                return self.params["signal"] + self.params['offset']
    def truncate(self, t_start, t_end):
        '''
        截取信号在[t_start, t_end]时间范围内的部分，其他部分暂置零，未来可以考虑振铃等边界效应
        '''
        t_array = np.array(self.t_list)
        signal = self.signal.copy()
        signal[(t_array < t_start) | (t_array > t_end)] = 0.0
        self.signal = signal


    
    def generate_basis_signal(self):
        '''
        产生基函数
        '''
        if not self.type == 6:
            return
        basis_type = self.params['basis_type']
        n_basis = self.params['n_basis']
        t_list = self.t_list
        if basis_type == "bspline":
            from scipy.interpolate import BSpline
            degree = 3  # B样条的阶数
            # 生成均匀节点
            knots = np.linspace(t_list[0], t_list[-1], n_basis - degree + 1)
            # 添加边界条件
            knots = np.r_[[t_list[0]] * degree, knots, [t_list[-1]] * degree]
            self.basis_functions = []
            for i in range(n_basis):
                coeffs = np.zeros(n_basis)
                coeffs[i] = 1.0
                spline = BSpline(knots, coeffs, degree)
                t_lo, t_hi = t_list[0], t_list[-1]
                def make_safe_spline(spl, lo, hi):
                    def f(t):
                        t = np.asarray(t, dtype=float)
                        result = spl(t)
                        result[(t < lo) | (t > hi)] = 0.0
                        return result
                    return f
                self.basis_functions.append(make_safe_spline(spline, t_lo, t_hi))
        elif basis_type == "fourier":
            T = t_list[-1] - t_list[0]
            self.basis_functions = []
            self.basis_functions.append(lambda t: np.ones_like(t) / np.sqrt(T))  # 添加常数项
            n_max = (n_basis - 1) // 2
            for n in range(1, n_max + 1):
                self.basis_functions.append(lambda t, n=n: np.sqrt(2/T) * np.sin(2 * np.pi * n * (t - t_list[0]) / T))
                self.basis_functions.append(lambda t, n=n: np.sqrt(2/T) * np.cos(2 * np.pi * n * (t - t_list[0]) / T))
        elif basis_type == "legendre":
            from scipy.special import legendre
            self.basis_functions = []
            for n in range(n_basis):
                Pn = legendre(n)
                Pn_normalized = lambda t, Pn=Pn: Pn((2 * (t - t_list[0]) / (t_list[-1] - t_list[0]) - 1))  # 归一化到[-1, 1]
                self.basis_functions.append(Pn_normalized)
        else:
            raise ValueError("Unsupported basis type.")
    
    def update_signal(self, **kwargs):
        '''
        更新信号参数并重新生成信号
        '''
        for key, value in kwargs.items():
            self.params[key] = value
        self.signal = self.generate()
        
    def back_signal(self, noise_level, seed, t_list=None):
        '''
        背景信号（可加噪声）
        '''
        if t_list is None:
            t_list = self.t_list
        # 添加高斯噪声
        np.random.seed(seed)
        noise = np.random.normal(0, noise_level, size=len(t_list))
        return noise
    
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
        

Signal = FluxSignal


class CompositeSignal(Signal):
    def __init__(self, signals):
        '''
        初始化CompositeSignal对象
        :param signals: Signal对象列表，表示要叠加的多个信号
        '''
        self.signals = signals
        self.t_list = self.get_t_list()
        self.signal = self.get_signal()

    def get_t_list(self):
        '''
        获取所有信号的时间列表的并集，并排序
        '''
        t_list = []
        curr = 0.0
        for signal in self.signals:
            pulse_list = [t + curr for t in signal.t_list]
            t_list.extend(pulse_list)
            if pulse_list:
                # 保证t_list没有重复的时间点
                curr = pulse_list[-1] + 1e-9  # 在最后一个时间点基础上加一个小的时间间隔，避免重复，同时两个脉冲之间有一个小的间隔，从而避开coeff边界的处理
        return np.array(t_list)
    
    def get_signal(self):
        '''
        在CompositeSignal的时间列表上叠加所有子信号的值
        '''
        signal = np.concatenate([signal.signal for signal in self.signals])
        return np.array(signal)
    
    def plot(self):
        '''
        绘制复合信号波形
        '''
        import matplotlib.pyplot as plt

        plt.figure(figsize=(10,4))
        plt.plot(self.t_list, self.signal)
        plt.xlabel('Time (ns)')
        plt.ylabel('Composite Signal Amplitude')
        plt.title('Composite Signal Waveform')
        plt.grid(True)
        plt.show()
