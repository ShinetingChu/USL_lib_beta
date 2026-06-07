#!/usr/bin/env python
# coding: utf-8

"""温度传递函数模块 - 计算温度波动对频率的影响。"""

from typing import List, Tuple
import numpy as np
from ._psd import calc_psd_single


def transfer_temp(
    t_a: List[float],
    T_a: List[float],
    t_0: float,
    tau: float,
    D_t: float,
    fs: float,
) -> Tuple[float, float, np.ndarray, np.ndarray]:
    """计算温度到频率的传递函数和Allan偏差。

    通过温度变化PSD和热时间常数计算腔体频率波动Allan偏差。

    Args:
        t_a: 时间序列
        T_a: 温度数据
        t_0: 过零点温度 (K)
        tau: 热时间常数 (s)
        D_t: 温度偏置
        fs: 采样率 (Hz)

    Returns:
        (D_t, alpha, t_x, sigma_y) 温度偏置, CTE, 平均时间数组, Allan偏差
    """
    num_t, dt, nu, S_out = calc_psd_single(t_a, T_a, fs)
    # 热膨胀系数 (CTE)
    beta = 1.8E-8  # /K^2
    alpha = beta * D_t  # /K, coefficient of thermal expansion(CTE)

    # 计算腔体PSD
    S_cav = np.array(S_out[1:]) / (1 + (np.array(nu[1:]) * tau) ** 2)

    # 由PSD计算频率Allan偏差
    num_dot = int(np.log2(num_t))
    t_x = 2 ** np.linspace(0, num_dot - 1, num_dot) * dt
    df = nu[1] - nu[0]
    a_t_x = np.array([t_x]).T
    a_nu = np.array([nu[1:]])
    a_S_cav = np.array([S_cav])
    sigma_y = 2 * df * np.sum(
        (np.sin(np.pi * a_t_x * a_nu) ** 4 / (np.pi * a_t_x * a_nu) ** 2) * a_S_cav,
        axis=1,
    )

    return D_t, alpha, t_x, sigma_y
