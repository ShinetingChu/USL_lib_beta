#!/usr/bin/env python
# coding: utf-8

"""画图函数模块 - 温度、频率稳定性、PSD、Allan偏差等数据可视化。"""

from datetime import datetime
from typing import List, Tuple, Optional, Dict, Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import signal, optimize
from scipy.fftpack import fft

import allantools as alt

from ._io import (
    KK_data_read_single,
    keysight_data_read,
    sim_keysight_data_read,
    labview_data_read,
    _datetime_to_epoch,
)
from ._drift import move_long_drift
from ._psd import psd_welch
from ._allan import allan_adev, allan_mdev
from ._transfer import transfer_temp
from ._utils import sin_func

# ────────────────────────────────────────────────────────────
# 温度数据绘图：合并 plot_data_labview 和 temp_read_psd_allan_2
# ────────────────────────────────────────────────────────────


def plot_temp_stability(
    path: str,
    fs: float = 2,
    label: str = '123',
    i: int = 1,
    nu_0: float = 429.228E12 / 2,
    start: int = 0,
    end: int = -1,
    label2: str = 'inloop',
    nfft_n: int = 1024,
    switch: int = 0,
    plot_transfer: bool = False,
    scale_by_nu0: bool = True,
) -> None:
    """绘制腔体温度稳定性（合并自 plot_data_labview 和 temp_read_psd_allan_2）。

    读取Labview温度数据，绘制时域曲线、功率谱密度和Allan偏差。
    可选：长漂补偿、温度传递函数分析。

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)
        label: 图例标签
        i: 数据列索引
        nu_0: 光频率 (Hz)
        start: 起始数据点
        end: 结束数据点
        label2: 子图标签
        nfft_n: Welch FFT点数
        switch: 长漂补偿开关 (0=关闭, 1=打开)
        plot_transfer: 是否绘制温度传递函数图（图4）
        scale_by_nu0: 是否除以 nu_0^2 转换为分数频率
    """
    # 读取数据
    t_a: List[float] = []
    t_2: List[datetime] = []
    f_1_d_raw: List[float] = []
    with open(path, 'r') as file:
        next(file)
        userlines = file.readlines()
    for line in userlines:
        datetime_obj = datetime.strptime(
            line.split('\t')[0], "%Y-%m-%d %H:%M:%S.%f"
        )
        t_a.append(
            time.mktime(datetime_obj.timetuple()) + datetime_obj.microsecond / 1E6
        )
        t_2.append(datetime_obj)
        f_1_d_raw.append(float(line.split('\t')[i]))

    t_a = t_a[start:end]
    t_2 = t_2[start:end]
    f_1_d_raw = f_1_d_raw[start:end]

    # 长漂补偿
    if switch >= 1:
        t_a, f_1_d, d_13 = move_long_drift(t_a, f_1_d_raw, switch)
        print(d_13)
    else:
        f_1_d = f_1_d_raw

    # 图1：时域
    plt.figure(1)
    plt.title('Cavity drift (' + label2 + ')')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.plot(t_2, np.array(f_1_d) - np.mean(f_1_d), label=label)
    plt.xticks(rotation=30)
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    # 根据scale_by_nu0选择是否除以光频率
    if scale_by_nu0:
        adev_data = np.array(f_1_d) / nu_0
        psd_scale = nu_0 ** 2
    else:
        adev_data = np.array(f_1_d)
        psd_scale = 1.0

    # Allan偏差
    taus_1, adevs_1, error_1 = allan_adev(adev_data, fs)
    f_psd, Pxx_1 = signal.welch(
        np.array(f_1_d) - np.mean(f_1_d), fs,
        window='hann', nperseg=nfft_n, nfft=nfft_n
    )

    # 图2：PSD
    plt.figure(2)
    plt.plot(f_psd, np.array(Pxx_1) / psd_scale, label=label)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Temperature stability (' + label2 + ')')
    plt.xlabel("Frequency(Hz)")
    plt.ylabel("Power Spectral Density($K^2/Hz$)")
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    # 图3：Allan偏差
    plt.figure(3)
    plt.errorbar(
        taus_1, np.array(adevs_1), np.array(error_1),
        fmt='o--', ecolor='r', elinewidth=2, capsize=4, label=label
    )
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Averaging Time(s)')
    plt.ylabel('Allan Deviation $\\sigma_y$')
    plt.title('Temperature stability (' + label2 + ')')
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    # 图4：温度传递函数（可选）
    if plot_transfer:
        D_t, alpha, t_x, sigma_y = transfer_temp(t_a, f_1_d, 123.85, 5.4E5, 0.1, fs)
        plt.figure(4)
        plt.plot(
            t_x, np.sqrt(sigma_y) * alpha, '.-',
            label=label + "_$\\Delta$T={:.2f}".format(D_t)
        )
        plt.xlim([1, 1000])
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Time interval(s)')
        plt.ylabel('Allan Deviation')
        plt.title('Allan Deviation of Frequency instability by Temperature')
        plt.grid(which='both', linestyle='dashed')
        plt.legend()


def plot_data_labview(*args, **kwargs):
    """(保持向后兼容) 调用合并后的 plot_temp_stability。"""
    plot_temp_stability(*args, **kwargs)


def temp_read_psd_allan_2(
    path: str,
    fs: float = 2,
    label: str = '123',
    i: int = 1,
    start: int = 0,
    end: int = -1,
    label2: str = 'inloop',
    nfft_n: int = 1024,
) -> None:
    """(保持向后兼容) 调用 plot_temp_stability (无漂移补偿, 带传递函数图, 不缩放)。"""
    plot_temp_stability(
        path, fs=fs, label=label, i=i, nu_0=429.228E12 / 2,
        start=start, end=end, label2=label2, nfft_n=nfft_n,
        switch=0, plot_transfer=True, scale_by_nu0=False,
    )


# ────────────────────────────────────────────────────────────
# 通用绘图辅助：稳定度三图（时域、PSD、ADEV）
# ────────────────────────────────────────────────────────────


def _plot_stability_triple(
    t_a: List[float],
    f_1_d: List[float],
    fs: float,
    k_p: float = 1,
    title: str = '',
    label: str = '',
    nfft_n: int = 1024,
    xlim_psd: Optional[Tuple[float, float]] = None,
    ylabel_time: str = '$\\Delta f/f$',
    ylabel_psd: str = 'Power Spectral Density($Hz^2/Hz$)',
    divide_nu0: bool = False,
    nu_0: float = 1.0,
    psd_scale: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """绘制稳定性三图：时域、PSD、Allan偏差。

    Args:
        t_a: 时间序列
        f_1_d: 数据序列
        fs: 采样率 (Hz)
        k_p: 比例因子
        title: 图表标题
        label: 图例标签
        nfft_n: FFT点数
        xlim_psd: PSD图的x轴范围
        ylabel_time: 时域图y轴标签
        ylabel_psd: PSD图y轴标签
        divide_nu0: 是否除以光频率
        nu_0: 光频率
        psd_scale: PSD缩放因子

    Returns:
        (taus, adevs, f_psd, Pxx)
    """
    # 图1：时域
    plt.figure(1)
    plt.title('Stability (' + title + ')')
    plt.xlabel('Time (s)')
    plt.ylabel(ylabel_time)
    plt.plot(
        np.array(t_a) - t_a[0],
        np.array(f_1_d) - f_1_d[0],
        label=label
    )
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    # 计算Allan偏差和PSD
    taus, adevs, error = allan_adev(f_1_d, fs)
    f_psd, Pxx = signal.welch(
        np.array(f_1_d) - np.mean(f_1_d), fs,
        window='hann', nperseg=nfft_n, nfft=nfft_n
    )

    # 图2：PSD
    plt.figure(2)
    plt.plot(f_psd, np.array(Pxx) * psd_scale, label=label)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Stability (' + title + ')')
    plt.xlabel("Frequency(Hz)")
    plt.ylabel(ylabel_psd)
    if xlim_psd:
        plt.xlim(xlim_psd)
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    # 图3：Allan偏差
    div = nu_0 if divide_nu0 else 1
    plt.figure(3)
    plt.errorbar(
        taus, np.array(adevs) * k_p / div,
        np.array(error) * k_p / div,
        fmt='o--', ecolor='r', elinewidth=2, capsize=4, label=label
    )
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Averaging Time(s)')
    plt.ylabel('Allan Deviation $\\sigma_y$')
    plt.title('Stability (' + title + ')')
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    return taus, adevs, f_psd, Pxx


# ────────────────────────────────────────────────────────────
# 各USB/示波器绘图函数 (简化为调用辅助函数)
# ────────────────────────────────────────────────────────────


def plot_pico_USB_err(
    path: str,
    switch: int,
    label: str,
    fs: float,
    k_p: float,
    title: str,
    start: int,
    end: int,
    nfft_n: int = 1024,
) -> Tuple[List[float], List[float]]:
    """绘制Pico示波器USB采集的光强稳定性数据。

    Args:
        path: 数据文件路径
        switch: 长漂补偿开关
        label: 图例标签
        fs: 采样率 (Hz)
        k_p: 比例因子
        title: 图表标题
        start: 起始数据点
        end: 结束数据点
        nfft_n: FFT点数

    Returns:
        (t_a, f_1_d) 时间和补偿后的数据
    """
    t: List[float] = []
    f: List[float] = []
    with open(path, 'r') as file:
        userlines = file.readlines()
    for line in userlines[3:-1]:
        t.append(float(line.split(',')[0]))
        f.append(float(line.split(',')[1]))

    t_a, f_1_d, d_13 = move_long_drift(t[start:end], f[start:end], switch)
    print(d_13)

    _plot_stability_triple(
        t_a, f_1_d, fs, k_p=k_p,
        title=title, label=label, nfft_n=nfft_n,
        xlim_psd=(0.001, fs / 2),
        ylabel_time='$\\Delta f/f$',
        ylabel_psd='Power Spectral Density($V^2/Hz$)',
    )

    return t_a, f_1_d


def plot_keysight_USB_power(
    path: str,
    switch: int,
    label: str,
    fs: float,
    k_p: float,
    title: str,
    start: int,
    end: int,
    factor: float = 1,
    nfft_n: int = 1024,
    switch2: int = 1,
    width: float = 100,
) -> Tuple[List[float], np.ndarray, np.ndarray, np.ndarray]:
    """绘制Keysight频率计数器USB采集的拍频稳定性数据。

    Args:
        path: 数据文件路径
        switch: 长漂补偿开关
        label: 图例标签
        fs: 采样率 (Hz)
        k_p: 比例因子
        title: 图表标题
        start: 起始数据点
        end: 结束数据点
        factor: 缩放因子
        nfft_n: FFT点数
        switch2: 离群值过滤开关
        width: 离群过滤带宽

    Returns:
        (t_a, f_detrend, f_psd, Pxx_scaled)
    """
    t, f = keysight_data_read(path, fs)
    t_a_m, f_1_d_m, d_13 = move_long_drift(
        t[start:end], f[start:end], switch
    )

    t_a, f_1_d = _filter_outliers(t_a_m, f_1_d_m, switch2, width)

    plt.figure(1)
    plt.title('Frequency changes with time (' + title + ')')
    plt.xlabel('Time (s)')
    plt.ylabel('Beat frequency (Hz)')
    plt.plot(np.array(t_a), np.array(f_1_d) - f_1_d[0], label=label)
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    taus_1, adevs_1, error_1 = allan_adev(f_1_d, fs)
    f_1_arr, Pxx_1 = signal.welch(
        np.array(f_1_d) - np.mean(f_1_d), fs,
        window='hann', nperseg=nfft_n, nfft=nfft_n
    )

    plt.figure(2)
    plt.plot(f_1_arr, np.array(Pxx_1) * k_p ** 2, label=label)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Stability (' + title + ')')
    plt.xlabel("Frequency(Hz)")
    plt.ylabel("Power Spectral Density($Hz^2/Hz$)")
    plt.xlim([0.001, fs / 2])
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    plt.figure(3)
    plt.errorbar(
        taus_1, np.array(adevs_1) * k_p, np.array(error_1) * k_p,
        fmt='o--', ecolor='r', elinewidth=2, capsize=4, label=label
    )
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Averaging Time(s)')
    plt.ylabel('Allan Deviation $\\sigma_y$')
    plt.title('Stability (' + title + ')')
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    return t_a, np.array(f_1_d) - f_1_d[0], f_1_arr, np.array(Pxx_1) * k_p ** 2


def _filter_outliers(
    t_a_m: List[float],
    f_1_d_m: List[float],
    switch2: int,
    width: float,
) -> Tuple[List[float], List[float]]:
    """根据中值附近的带宽过滤离群值。"""
    t_a: List[float] = []
    f_1_d: List[float] = []
    if switch2 == 1:
        f_mean = np.mean(f_1_d_m[0:300])
        for i in range(0, len(f_1_d_m)):
            if np.abs(f_mean - f_1_d_m[i]) <= width:
                t_a.append(t_a_m[i])
                f_1_d.append(f_1_d_m[i])
    else:
        t_a = t_a_m
        f_1_d = f_1_d_m
    return t_a, f_1_d


def plot_keysight_six_half_USB_power(
    path: str,
    switch: int,
    label: str,
    fs: float,
    k_p: float,
    title: str,
    start: int,
    end: int,
    factor: float = 1,
    nfft_n: int = 1024,
    switch2: int = 1,
    width: float = 100,
) -> Tuple[np.ndarray, np.ndarray]:
    """绘制Keysight六位半万用表USB采集的电压稳定性数据。

    Args:
        path: 数据文件路径
        switch: 长漂补偿开关
        label: 图例标签
        fs: 采样率 (Hz)
        k_p: 比例因子
        title: 图表标题
        start: 起始数据点
        end: 结束数据点
        factor: 缩放因子
        nfft_n: FFT点数
        switch2: 离群值过滤开关
        width: 离群过滤带宽

    Returns:
        (f_psd, Pxx_scaled)
    """
    t, f = keysight_data_read(path, fs)
    t_a_m, f_1_d_m, d_13 = move_long_drift(
        t[start:end], f[start:end], switch
    )

    t_a, f_1_d = _filter_outliers(t_a_m, f_1_d_m, switch2, width)

    plt.figure(1)
    plt.title('Frequency changes with time (' + title + ')')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.plot(
        np.array(t_a) - t_a[0],
        (np.array(f_1_d) - f_1_d[0]) * factor,
        label=label
    )
    plt.grid(which='both', linestyle='dashed')

    taus_1, adevs_1, error_1 = allan_adev(t_a, f_1_d)
    f_1_arr, Pxx_1 = signal.welch(
        np.array(f_1_d) - np.mean(f_1_d), fs,
        window='hann', nperseg=nfft_n, nfft=nfft_n
    )

    plt.figure(2)
    plt.plot(f_1_arr, np.array(Pxx_1) * k_p ** 2, label=label)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Stability (' + title + ')')
    plt.xlabel("Voltage (V)")
    plt.ylabel("Power Spectral Density($V^2/Hz$)")
    plt.xlim([0.001, fs / 2])
    plt.grid(which='both', linestyle='dashed')

    plt.figure(3)
    plt.errorbar(
        taus_1, np.array(adevs_1) * k_p, np.array(error_1) * k_p,
        fmt='o--', ecolor='r', elinewidth=2, capsize=4, label=label
    )
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Averaging Time(s)')
    plt.ylabel('Allan Deviation $\\sigma_y$')
    plt.title('Stability (' + title + ')')
    plt.grid(which='both', linestyle='dashed')
    print(adevs_1[3] * k_p)

    return f_1_arr, np.array(Pxx_1) * k_p ** 2


def plot_sim_keysight_USB_power(
    path: str,
    switch: int,
    label1: str = '123',
    fs: float = 1,
    k_p: float = 1,
    title: str = '123',
    start: int = 0,
    end: int = -1,
    nfft_n: int = 1024,
) -> np.ndarray:
    """绘制模拟Keysight USB采集的数据稳定性。

    Args:
        path: 数据文件路径
        switch: 长漂补偿开关
        label1: 图例标签
        fs: 采样率 (Hz)
        k_p: 比例因子
        title: 图表标题
        start: 起始数据点
        end: 结束数据点
        nfft_n: FFT点数

    Returns:
        adevs_a allantools计算的Allan偏差
    """
    t, f = sim_keysight_data_read(path, fs)
    t_a, f_1_d, d_13 = move_long_drift(
        t[start:end], f[start:end], switch
    )

    plt.figure(0)
    plt.title('Light intensity stability (' + title + ')')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.plot(
        np.array(t_a) - t_a[0],
        np.array(f_1_d),
        label=label1
    )
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    plt.figure(1)
    plt.title('Light intensity stability (' + title + ')')
    plt.xlabel('Time (s)')
    plt.ylabel('$\\Delta f/f$')
    plt.plot(
        np.array(t_a) - t_a[0],
        np.array(f_1_d) * k_p,
        label=label1
    )
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    taus_1, adevs_1, error_1 = allan_adev(t_a, f_1_d)
    f_1_arr, Pxx_1 = signal.welch(
        np.array(f_1_d) - np.mean(f_1_d), fs,
        window='hann', nperseg=nfft_n, nfft=nfft_n
    )
    (taus_a, adevs_a, errors_1, ns_1) = alt.adev(
        f_1_d, rate=fs, data_type="freq", taus=1
    )

    plt.figure(2)
    plt.plot(f_1_arr, np.array(Pxx_1) * k_p ** 2, label=label1)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Stability (' + title + ')')
    plt.xlabel("Frequency(Hz)")
    plt.ylabel("Power Spectral Density($V^2/Hz$)")
    plt.xlim([0.01, 50])
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    plt.figure(3)
    plt.errorbar(
        taus_1, np.array(adevs_1) * k_p, np.array(error_1) * k_p,
        fmt='o--', ecolor='r', elinewidth=2, capsize=4, label=label1
    )
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Averaging Time(s)')
    plt.ylabel('Allan Deviation $\\sigma_y$')
    plt.title('stability (' + title + ')')
    plt.xlim([1E-2, 1E2])
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=(1, 0))

    return adevs_a


# ────────────────────────────────────────────────────────────
# K+K频率计数器绘图 (单段 / 多段 / 两文件合并)
# ────────────────────────────────────────────────────────────

# 通道配置：差异化标题、x轴范围、拟合线系数、图例位置
_CHANNEL_CONFIGS: Dict[str, Dict[str, Any]] = {
    'CH1': {
        'title1': 'Beat Frequency of U3 and Si1',
        'title2': 'Allan Deviation of U3 and Si1',
        'title3': 'Modified Allan Deviation of U3 and Si1',
        'xlim_psd': [1E-3, None],
        'fit_line_coeff': None,
        'legend_loc_psd': (1, 0),
        'legend_loc_ad': (1, 0),
        'legend_loc_mad': (1, 0),
        'xticks_rotation': False,
        'plot_time': True,
    },
    'CH1-2': {
        'title1': 'Beat Frequency of U3 and Si1',
        'title2': 'Allan Deviation of U3 and Si1',
        'title3': 'Modified Allan Deviation of U3 and Si1',
        'xlim_psd': [1E-3, None],
        'fit_line_coeff': None,
        'legend_loc_psd': (1, 0),
        'legend_loc_ad': (1, 0),
        'legend_loc_mad': (1, 0),
        'xticks_rotation': False,
        'plot_time': True,
    },
    'CH2': {
        'title1': 'Beat Frequency of U1 and U3',
        'title2': 'Allan Deviation of U1 and U3',
        'title3': 'Modified Allan Deviation of U1 and U3',
        'xlim_psd': [1E-2, None],
        'fit_line_coeff': 1.6E-31,
        'legend_loc_psd': 'best',
        'legend_loc_ad': 'lower left',
        'legend_loc_mad': 'best',
        'xticks_rotation': True,
        'plot_time': True,
    },
    'CH3': {
        'title1': 'Beat Frequency of U1 and Si1',
        'title2': 'Allan Deviation of U1 and Si1',
        'title3': 'Modified Allan Deviation of U1 and Si1',
        'xlim_psd': [1E-2, None],
        'fit_line_coeff': 1.6E-31,  # K_K_plot uses 1.6E-31; path1_path2 overrides to 3E-33
        'legend_loc_psd': (1, 0),
        'legend_loc_ad': (1, 0),
        'legend_loc_mad': (1, 0),
        'xticks_rotation': True,
        'plot_time': True,
    },
}


def _plot_kk_core(
    t_a: List[float],
    f_1_d: List[float],
    fs: float,
    ch: str,
    label1: str,
    nu_0: float,
    nfft_0: int,
    plot_mdev: bool = True,
    fig_offset: int = 1,
    legend_label_override: Optional[str] = None,
    fit_line_coeff: Optional[float] = None,
) -> Dict[str, Any]:
    """K+K数据核心绘图（时域 + PSD + ADEV + MDEV）。

    Args:
        t_a: 时间序列
        f_1_d: 频率数据
        fs: 采样率
        ch: 通道名
        label1: 图例标签
        nu_0: 光频率
        nfft_0: FFT点数
        plot_mdev: 是否绘制MDEV
        fig_offset: 图号偏移
        legend_label_override: 图例标签覆盖

    Returns:
        计算结果字典
    """
    cfg = _CHANNEL_CONFIGS.get(ch, _CHANNEL_CONFIGS['CH1'])
    leg = legend_label_override or label1

    # 时域图
    if cfg['plot_time']:
        plt.figure(fig_offset)
        plt.title(cfg['title1'] + ' (' + leg + ')')
        plt.xlabel('Time (s)')
        plt.ylabel('Frequency(Hz)')
        plt.plot(
            np.array(t_a) - t_a[0],
            np.array(f_1_d) - f_1_d[0],
            label=leg
        )
        if cfg['xticks_rotation']:
            plt.xticks(rotation=30)
        plt.grid(which='both', linestyle='dashed')
        plt.legend(loc=(1, 0))

    # Allan偏差
    taus_a3, adevs_a3, error_a3 = allan_adev(f_1_d, fs)
    taus_m3, adevs_m3, error_m3 = allan_mdev(f_1_d, fs)
    (taus_1, adevs_1, errors_1, ns_1) = alt.adev(
        f_1_d, rate=fs, data_type="freq", taus=1
    )

    P_m_13, freq_m_13 = psd_welch(t_a, f_1_d, fs, nfft_0)

    # PSD图
    plt.figure(fig_offset + 1)
    plt.plot(freq_m_13, P_m_13 / nu_0 ** 2, label=leg)
    plt.xscale('log')
    plt.yscale('log')
    plt.title(cfg['title1'])
    plt.xlabel("Frequency(Hz)")
    plt.ylabel("Power Spectral Density($\\sigma_y/\\sqrt{Hz}$)")
    xlim_low = cfg['xlim_psd'][0]
    xlim_high = cfg['xlim_psd'][1] if cfg['xlim_psd'][1] is not None else fs
    plt.xlim([xlim_low, xlim_high])
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=cfg['legend_loc_psd'])

    # 拟合线（优先使用传入系数，否则使用通道配置）
    fc = fit_line_coeff if fit_line_coeff is not None else cfg['fit_line_coeff']
    if fc is not None:
        fit_line = [fc / i for i in freq_m_13]
        plt.plot(freq_m_13, fit_line, 'r', linestyle='--')

    # ADEV图
    plt.figure(fig_offset + 2)
    plt.errorbar(
        taus_a3, np.array(adevs_a3) / nu_0, np.array(error_a3) / nu_0,
        fmt='o--', ecolor='r', elinewidth=2, capsize=4, label=leg
    )
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Averaging Time(s)')
    plt.ylabel('Allan Deviation $\\sigma_y$')
    plt.title(cfg['title2'])
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc=cfg['legend_loc_ad'])

    # MDEV图
    if plot_mdev:
        plt.figure(fig_offset + 3)
        plt.errorbar(
            taus_m3, np.array(adevs_m3) / nu_0, np.array(error_m3) / nu_0,
            fmt='o--', ecolor='r', elinewidth=2, capsize=4, label=leg
        )
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('Averaging Time(s)')
        plt.ylabel('Modified Allan Deviation $\\sigma_y$')
        plt.title(cfg['title3'])
        plt.grid(which='both', linestyle='dashed')
        plt.legend(loc=cfg['legend_loc_mad'])

    result = {
        'freq_psd': freq_m_13,
        'psd': P_m_13 / nu_0 ** 2,
        'taus_m': taus_m3,
        'adevs_m': adevs_m3 / nu_0,
        'error_m': error_m3 / nu_0,
        'taus_a': taus_a3,
        'adevs_a': adevs_a3 / nu_0,
        'adevs_1': adevs_1 / nu_0,
        'mdevs_1': None,
    }
    if plot_mdev:
        (taus_1_m, mdevs_1_m, errors_1_m, ns_1_m) = alt.mdev(
            f_1_d, rate=fs, data_type="freq", taus=1
        )
        result['mdevs_1'] = mdevs_1_m / nu_0
        result['taus_1'] = taus_1_m

    return result


def K_K_plot(
    path: str,
    fs: float,
    CH: str = 'CH1',
    label1: str = '123',
    start: int = 1,
    end: int = 10,
    nu_0: float = 429.228E12,
    nfft_0: int = 1024 * 4,
    switch1: int = 0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """K+K频率计数器数据绘图（单段数据）。

    读取数据，绘制时域 + PSD + ADEV + MDEV四图，按通道配置差异化标题/范围/拟合线。

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)
        CH: 通道名 ('CH1', 'CH1-2', 'CH2', 'CH3')
        label1: 图例标签
        start: 起始数据点
        end: 结束数据点
        nu_0: 光频率 (Hz)
        nfft_0: FFT点数
        switch1: 长漂补偿开关

    Returns:
        (freq_psd, psd, taus_m, adevs_m, error_m)
    """
    t, t_data, f_1 = KK_data_read_single(path, fs, start, end, channel=CH)
    t_a, f_1_d, d_13 = move_long_drift(t, f_1, switch1)
    print("Drift is: {} Hz/s.".format(d_13))

    result = _plot_kk_core(t_a, f_1_d, fs, CH, label1, nu_0, nfft_0)

    # 按原代码的差异处理打印信息和返回值
    if CH == 'CH1':
        if result.get('adevs_1') is not None and result.get('mdevs_1') is not None:
            print(result['adevs_1'], result['mdevs_1'])
    elif CH == 'CH1-2':
        n = 8
        if n < len(result.get('taus_a', [])):
            idx = min(n, len(result['taus_a']) - 1)
            print(result['adevs_a'][idx], result.get('adevs_m', [])[idx]
                  if len(result.get('adevs_m', [])) > idx else 0)
    elif CH == 'CH2':
        n = 8
        if n < len(result.get('taus_a', [])):
            idx = min(n, len(result['taus_a']) - 1)
            print(result['adevs_a'][idx], result.get('adevs_m', [])[idx]
                  if len(result.get('adevs_m', [])) > idx else 0)
    elif CH == 'CH3':
        if result.get('adevs_1') is not None:
            print(result['adevs_1'])

    return (result['freq_psd'], result['psd'],
            result['taus_m'], result['adevs_m'], result['error_m'])


def K_K_plot_path1_path2(
    path1: str,
    path2: str,
    fs: float,
    CH: str = 'CH1',
    label1: str = '123',
    start: int = 1,
    end: int = 10,
    nu_0: float = 429.228E12,
    nfft_0: int = 1024 * 4,
    switch1: int = 0,
    buttom: int = 0,
    width: float = 30,
) -> np.ndarray:
    """K+K频率计数器双文件合并绘图。

    读取两个数据文件并拼接，绘制时域 + PSD + ADEV + MDEV四图。

    Args:
        path1: 第一数据文件路径
        path2: 第二数据文件路径
        fs: 采样率 (Hz)
        CH: 通道名 ('CH1', 'CH2', 'CH3')
        label1: 图例标签
        start: 起始数据点
        end: 结束数据点
        nu_0: 光频率 (Hz)
        nfft_0: FFT点数
        switch1: 长漂补偿开关
        buttom: 离群值过滤开关
        width: 离群过滤带宽

    Returns:
        adevs_1 / nu_0 ADEV@1s值
    """
    t_a, t_data_a, f_1_a = KK_data_read_single(path1, fs, 0, -1, channel=CH)
    t_b, t_data_b, f_1_b = KK_data_read_single(path2, fs, 0, -1, channel=CH)

    t = t_a + t_b
    t = t[start:end]
    f_1 = f_1_a + f_1_b
    f_1 = f_1[start:end]

    t_a_m, f_1_d_m, d_13 = move_long_drift(t, f_1, switch1)
    print('Drift is {:.3f}mHz/s'.format(-d_13 * 1000))

    t_a, f_1_d = _filter_outliers(t_a_m, f_1_d_m, buttom, width)

    # path1_path2的CH3使用3E-33系数（与单段K_K_plot的1.6E-31不同）
    _path12_fit_coeff: Optional[float] = 3E-33 if CH == 'CH3' else None
    result = _plot_kk_core(
        t_a, f_1_d, fs, CH, label1, nu_0, nfft_0,
        legend_label_override=_get_path12_legend(ch=CH, label=label1),
        fit_line_coeff=_path12_fit_coeff,
    )

    if result.get('adevs_1') is not None:
        print(result['adevs_1'])
    return result.get('adevs_1', np.array([]))


def _get_path12_legend(ch: str, label: str) -> str:
    """获取双文件合并的图例标签。"""
    if ch == 'CH1':
        return label
    elif ch == 'CH2':
        return 'U1&U3' + label
    elif ch == 'CH3':
        return 'Si1&U1' + label
    return label


def K_K_single_plot(
    path: str,
    fs: float,
    CH: str = 'CH1',
    label1: str = '123',
    start: int = 1,
    end: int = 10,
    nu_0: float = 429.228E12,
    nfft_0: int = 1024 * 4,
    switch1: int = 0,
) -> None:
    """K+K单段数据简易绘图（替代旧版K_K_single_plot，保持向后兼容）。

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)
        CH: 通道名 ('CH1', 'CH2', 'CH3', 'CH1-2')
        label1: 图例标签
        start: 起始数据点
        end: 结束数据点
        nu_0: 光频率 (Hz)
        nfft_0: FFT点数
        switch1: 长漂补偿开关
    """
    t, t_data, f_1 = KK_data_read_single(path, fs, start, end, channel=CH)
    t_a, f_1_d, d_13 = move_long_drift(t, f_1, switch1)

    result = _plot_kk_core(t_a, f_1_d, fs, CH, label1, nu_0, nfft_0, plot_mdev=True)
    if result.get('adevs_1') is not None and result.get('mdevs_1') is not None:
        print(result['adevs_1'], result['mdevs_1'])


# ────────────────────────────────────────────────────────────
# 鉴频斜率 / SR780 / 相噪 绘图
# ────────────────────────────────────────────────────────────


def freq_disc_slope(
    path: str,
    fs: float,
    V: float,
    f_mod: float,
    channel: str = 'CH1',
) -> float:
    """鉴频斜率拟合计算。

    读取拍频数据，通过正弦拟合获取调制幅度，计算鉴频斜率。

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)
        V: 调制电压 (V)
        f_mod: 调制频率 (Hz)
        channel: K+K通道名

    Returns:
        slope 鉴频斜率 (V/Hz)
    """
    t_0, t_data_0, f_3_0 = KK_data_read_single(
        path, fs, 0 * fs, 10 * fs, channel=channel
    )
    t_a, f_3_d, d_23 = move_long_drift(t_0, f_3_0, 0)
    t_b = [i / fs for i in range(len(f_3_d))]

    A_0, omega_0, phi_0, C_0 = optimize.curve_fit(
        sin_func, t_b, np.array(f_3_d) - f_3_d[0],
        bounds=([0, 2 * np.pi * (f_mod - 0.01), 0, -600],
                [8000, 2 * np.pi * (f_mod + 0.01), 2 * np.pi, 600])
    )[0]

    t_c = np.linspace(0, t_b[-1], num=fs * 100)
    plt.scatter(t_b, np.array(f_3_d) - f_3_d[0])
    plt.plot(np.array(t_c), sin_func(np.array(t_c), A_0, omega_0, phi_0, C_0), 'r')

    print(A_0)
    return V / (2 * A_0)


def SR780_data_concatenate(
    nu_a: List[float] = None,
    nu_b: List[float] = None,
    nu_c: List[float] = None,
    nu_d: List[float] = None,
    psd_a: List[float] = None,
    psd_b: List[float] = None,
    psd_c: List[float] = None,
    psd_d: List[float] = None,
    label1: str = 'up_light_far_off_resonant',
) -> None:
    """SR780频谱分析仪数据拼接绘图。

    将多个频段的PSD数据拼接，计算Allan偏差。

    Args:
        nu_a: 频段A频率
        nu_b: 频段B频率
        nu_c: 频段C频率
        nu_d: 频段D频率
        psd_a: 频段A PSD
        psd_b: 频段B PSD
        psd_c: 频段C PSD
        psd_d: 频段D PSD
        label1: 图例标签
    """
    if nu_a is None:
        nu_a = [0, 1]
    if nu_b is None:
        nu_b = [0, 1]
    if nu_c is None:
        nu_c = [0, 1]
    if nu_d is None:
        nu_d = [0, 1]
    if psd_a is None:
        psd_a = [0, 1]
    if psd_b is None:
        psd_b = [0, 1]
    if psd_c is None:
        psd_c = [0, 1]
    if psd_d is None:
        psd_d = [0, 1]

    nu_all = nu_a + nu_b + nu_c + nu_d
    psd_all = np.concatenate((psd_a, psd_b, psd_c, psd_d))

    t_x = np.linspace(1 / (2 * nu_all[-1]), 1 / (2 * nu_all[0]), len(nu_all))

    df_1 = nu_a[1] - nu_a[0]
    df_2 = nu_b[1] - nu_b[0]
    df_3 = nu_c[1] - nu_c[0]
    df_4 = nu_d[1] - nu_d[0]

    sigma_f_2 = []
    for t_x_i in t_x:
        s = (2 * df_1 * np.dot(
            np.sin(np.pi * t_x_i * np.array(nu_a)) ** 4 / (np.pi * t_x_i * np.array(nu_a)) ** 2,
            np.array(psd_a)
        ) + 2 * df_2 * np.dot(
            np.sin(np.pi * t_x_i * np.array(nu_b)) ** 4 / (np.pi * t_x_i * np.array(nu_b)) ** 2,
            np.array(psd_b)
        ) + 2 * df_3 * np.dot(
            np.sin(np.pi * t_x_i * np.array(nu_c)) ** 4 / (np.pi * t_x_i * np.array(nu_c)) ** 2,
            np.array(psd_c)
        ) + 2 * df_4 * np.dot(
            np.sin(np.pi * t_x_i * np.array(nu_d)) ** 4 / (np.pi * t_x_i * np.array(nu_d)) ** 2,
            np.array(psd_d)
        ))
        sigma_f_2.append(np.sqrt(s))

    plt.figure(5)
    plt.plot(nu_all, psd_all, label=label1)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude ($Hz^{-1}$)')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc='best')

    plt.figure(7)
    plt.plot(t_x, sigma_f_2, label=label1)
    plt.xlabel('Average time (s)')
    plt.ylabel('Fractional frequency')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc='best')
