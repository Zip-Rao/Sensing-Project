# 仿真产出与结论
按照脉冲扫描的方式，可以将时变磁场测量分为三种协议：
- 基于脉冲延迟时间扫描的瞬态磁场协议
- 基于演化时间扫描的Ramsey协议和差分回波协议
- 基于频谱扫描的CPMG协议

为了达到最高的磁场灵敏度，Transmon qubit的偏置工作磁场设置在斜率最大处。为了能测量正负磁场，这种偏置是必要的
## 延迟时间扫描
通过等效时间采样得到$p_{meas}$，随后采用了两种方法进行磁场重建
- Wiener反卷积
- [Levenberg-Marquardt优化算法](./numerical_inversion_scheme.md)


### 两种重建方法对比：输入信号
![alt text](../result/different_signals/comparison.png)
![alt text](image-1.png)
测试了四种磁场信号：
- 正弦信号：0.01 Phi_0, 0.01 GHz,0.0005
- 脉冲信号：用alpha函数表示，噪声幅度0.0001
- 双峰信号：两个重叠的高斯函数，噪声幅度0.0001
- 复杂信号：三个重叠的高斯波包，噪声幅度0.00005

其中，LM优化用的是fourier基，最大迭代次数不超过20次。

可以看到，对于窄带宽信号（sine，step），LM优化表现更好。宽带宽信号（double peak，complex）中，Wiener反卷积表现更好。

理论上，LM优化最终可以将残差降低到任意小，但是考虑优化时间，迭代次数以及基函数个数都有一定限制，因此对于宽带宽信号，LM优化的表现可能不如Wiener反卷积。

两种重建方法对于噪声都有一定的鲁棒性，因为选用了较大的正则化参数，另外由于在采样过程得到的是激发态概率，其是信号的积分，因此对于高频噪声有一定的滤波作用。

### 两种重建方法对比：幅值
由于Transmon qubit的响应函数为非线性函数，因此当信号幅值较大时，Wiener反卷积的线性假设不在成立，而LM优化算法则不受此限制，因此在信号幅值较大时，LM优化算法的表现更好。

可以通过hammerstein-Wiener反卷积来改进Wiener反卷积算法，使其适用于非线性系统。
![alt text](image-3.png)
从图中可知，当幅度超过0.03时，Wiener反卷积的重建误差迅速增加，而LM优化算法的重建误差相对较小。
![alt text](../result/signal_amp/amplitude_scan/check_amplitude_0.08.png)
![alt text](../result/signal_amp/amplitude_scan/check_amplitude_0.04.png)
![alt text](../result/signal_amp/amplitude_scan/check_amplitude_0.02.png)
![alt text](../result/signal_amp/amplitude_scan/check_amplitude_0.01.png)
### Wiener反卷积：正则化参数lambda
lambda为10时，重建结果较好，lambda过大导致全局频率被压低，信号的高频部分也被抑制，导致重建信号展宽变宽，幅度变小，lambda过小导致高频噪声被放大，从而出现振荡
![alt text](../result/lambda_scan/wiener_lambda_scan.png)

### Levenberg-Marquardt优化：正则化参数lambda
LM算法对于lambda的选择相对不敏感，因为这里的lambda主要用于高频信号的惩罚，而LM算法通过自适应调节阻尼因子，来补偿lamda选择的不当。不过，具体来看，可以看到，较小的lambda由于对高频信号惩罚不足，导致优化算法过拟合了高频噪声，从而在低频区域出现了一定的高频振荡，而较大的lambda则过度惩罚了高频信号，导致重建信号的高频部分被抑制，从而出现了展宽和幅度降低的现象。
![alt text](../result/lambda_scan/lm_lambda_scan.png)

### Levenberg-Marquardt优化：迭代次数与收敛性
可以看到，残差在前几次迭代中迅速下降，之后缓慢下降，在最后收敛到0.0146。这是LM算法的典型行为，优化前期是梯度下降主导的快速下降，后期为Gauss-Newton主导的精细调整。实际应用中，发现在残差小于0.1时重建的波形已经非常接近输入信号，因此LM算法可以在预期的时间内收敛到满意的结果。
![alt text](../result/convergence/lm_convergence.png)
## 演化时间扫描
验证结论
### Ramsey协议
改进后的Ramsey协议可以实现对任意时变磁场的重建，但是由于重建过程依赖数值微分，因此对于噪声的耐受性较差。

重建信号的整体位移和边界效应来源于脉冲信号的宽度，当时间尺度在ns级别时，脉冲信号的宽度不能忽略。上述问题可以通过后续的数据处理来解决。
![alt text](../result/ramsey/field_reconstruction(0.001Hz).png)
磁场信号：幅度：0.001，频率：0.001，噪声幅度：0.0001
![alt text](../result/ramsey/field_reconstruction(0.004Hz).png)
磁场信号：幅度：0.001，频率：0.004，噪声幅度：0.0001

### 差分回波协议
TODO...
## 频谱扫描

### CPMG协议作为带通滤波器
TODO...
## 其他实现

[Transmon qubit能级的数值求解，单比特门，两比特门](../note/Numerical%20calculation.ipynb)


## TODO
- 扩展Ramsey协议，支持相位unwrap(IQ调制)，
- 基函数比较：Fourier，B-spline，Legendre
- [差分回波协议](./_sensing%20theory.md)
- [CPMG协议作为带通滤波器](./_sensing%20theory.md)
- tomo
- 加噪声，加耗散
- 并行优化，代码结构优化：
    - 统一时间轴t_global，统一时间间隔，优化函数接口与扩展性
    - 全局统一t_rabi
    - 超导量子比特的实际架构：参考Practical Guide for Building Superconducting Quantum Devices
- demo
- 当脉冲时间不可忽略时，完善以上协议


- 研究量子传感协议，是为了探索以下几种应用：
    - 高灵敏高分辨磁场传感
    - calibration：qubit频率
    - 波形矫正：predistortion

实际的calibration包含从qubit频率标定，单比特门标定，读出标定，多比特门标定等单向依赖的多个步骤，可包含在一个自适应的DAG图中。本项目主要关注单比特门标定前的标定工作，包括qubit频率标定和波形矫正。

qubit频率标定有两种，一种是标定$f(\Phi)$曲线，一种是在确定工作点后标定频率$f_{01}$

瞬态磁场协议的应用：
- qubit频率标定（不太适用，推导？）
- 预失真标定（适用，推导？）


