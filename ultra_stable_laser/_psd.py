#!/usr/bin/env python
# coding: utf-8

"""PSD/CSD功率谱密度计算模块。"""

from typing import List, Tuple
import numpy as np
from scipy import signal
from scipy.fftpack import fft


def psd_welch(
    t_x: List[float],
    f_x: List[float],
    fs: float,
    nfft_0: int = 1024 * 16 * 8,
) -> Tuple[np.ndarray, np.ndarray]:
    """使用Welch方法计算功率谱密度。

    Args:
        t_x: 时间序列（仅用于长度参考）
        f_x: 数据序列
        fs: 采样率 (Hz)
        nfft_0: FFT点数

    Returns:
        (Pxx, nu) 功率谱密度和对应频率
    """
    nu, Pxx = signal.welch(f_x, fs, window='hann', scaling='density',
                           nperseg=0.5 * nfft_0, nfft=nfft_0)
    return Pxx, nu  # unit: (f_x's unit)^2/Hz, Hz


def psd_int_allan(
    nu: List[float],
    Pxx: List[float],
) -> Tuple[np.ndarray, List[float]]:
    """从PSD积分计算Allan方差秒稳的影响。

    从频域起点开始积分计算对Allan方差秒稳的影响。

    Args:
        nu: 频率数组 (Hz)
        Pxx: 功率谱密度 ((f_x's unit)^2/Hz)

    Returns:
        (nu[2:], sigma_f) 截断频率和对应的积分Allan偏差
    """
    df = nu[1] - nu[0]
    t_x = 1  # @1s
    sigma_f: List[float] = []
    m=10
    for i in range(2, int(len(nu)/m)):
        a_nu = np.array([nu[1:i*m]])  # unit: Hz
        a_S = np.array([(Pxx[1:i*m])])  # unit: (Pxx's unit)^2/Hz
        sigma_f.append(2 * df * np.sum(
            np.sin(np.pi * t_x * a_nu) ** 4 / (np.pi * t_x * a_nu) ** 2 * a_S
        ))
    return nu[2:len(nu)], sigma_f


def plot_csd(
    f_13: List[float],
    f_23: List[float],
    fs: float,
    f0: float,
    nfft: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """计算交叉功率谱密度 (CSD)。

    Args:
        f_13: 第一通道数据
        f_23: 第二通道数据
        fs: 采样率 (Hz)
        f0: 参考频率
        nfft: FFT点数

    Returns:
        (f, csd_12) 频率和归一化CSD
    """
    f, Pxy1 = signal.csd(f_13, f_23, fs, window='hann', nperseg=nfft)
    csd_12 = np.array(np.abs(Pxy1)) / f0 ** 2
    return f, csd_12


def calc_psd_single(
    t_a: List[float],
    T_a: List[float],
    fs: float,
) -> Tuple[int, float, np.ndarray, np.ndarray]:
    """计算单边功率谱密度（FFT方法）。

    原PSD函数，使用FFT直接计算单边功率谱密度。
    被transfer_temp调用。

    Args:
        t_a: 时间序列
        T_a: 数据序列（如温度、频率等）
        fs: 采样率 (Hz)

    Returns:
        (num_t, dt, nu, S_out) 点数, 时间间隔, 频率数组, 单边PSD
    """
    tot = t_a[-1] - t_a[0]
    num_t = len(t_a)
    dt = tot / num_t
    nu = np.linspace(0, fs / 2, int(len(t_a) / 2) + 1)
    yy = fft(T_a)
    S_out = np.power(abs(yy[0:int(num_t / 2) + 1] * dt), 2) / tot
    S_out[1:] = 2 * S_out[1:]
    return num_t, dt, nu, S_out
