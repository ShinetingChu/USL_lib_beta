#!/usr/bin/env python
# coding: utf-8

"""数据读取模块 - 各种频率计数器和仪器的数据读取函数。"""

from datetime import datetime
from typing import List, Tuple, Optional, Union
import numpy as np
import pandas as pd


def KK_data_read(
    path: str,
    fs: float,
    begin: int,
    end: int,
    CH1: int = 3,
    CH2: int = 4,
    CH3: int = 5,
    CH4: int = 6,
    CH5: int = 8,
) -> Tuple[List[float], List[datetime], List[float], List[float], List[float], List[float], List[float]]:
    """K+K频率计数器多通道数据读取。

    使用pandas读取，自动解析时间戳和多列频率数据。

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)
        begin: 开始行
        end: 结束行
        CH1-CH5: 各通道在文件中的列索引

    Returns:
        (t, t_data, f_1, f_2, f_3, f_4, f_5)
    """
    nrows = end - (begin + 1)
    if nrows <= 0:
        return [], [], [], [], [], [], []

    df = pd.read_csv(path, sep='\s+', header=None, skiprows=begin + 1, nrows=nrows)

    # 从两列拼接解析日期时间："20" + column0 + column1 → "%Y%m%d%H%M%S.%f"
    datetime_str = '20' + df.iloc[:, 0].astype(str) + df.iloc[:, 1].astype(str)
    t_data: List[datetime] = pd.to_datetime(datetime_str, format="%Y%m%d%H%M%S.%f").to_list()

    t: List[float] = [td.timestamp() for td in t_data]

    f_1: List[float] = df.iloc[:, CH1].astype(float).tolist()
    f_2: List[float] = df.iloc[:, CH2].astype(float).tolist()
    f_3: List[float] = df.iloc[:, CH3].astype(float).tolist()
    f_4: List[float] = df.iloc[:, CH4].astype(float).tolist()
    f_5: List[float] = df.iloc[:, CH5].astype(float).tolist()

    return t, t_data, f_1, f_2, f_3, f_4, f_5


def KK_data_read_single(
    path: str,
    fs: float,
    begin: int,
    end: int,
    channel: str = 'CH1',
) -> Tuple[np.ndarray, Union[List[datetime], np.ndarray], np.ndarray]:
    """K+K频率计数器单通道数据读取。

    使用pandas读取，支持低采样率（fs<=50）时解析时间戳，
    高采样率时直接生成等间隔时间序列。

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)
        begin: 开始行
        end: 结束行
        channel: 通道名 ('CH1', 'CH2', 'CH3', 'CH1-2')

    Returns:
        (t, t_data, f_1)
    """
    channel_map = {'CH1': 3, 'CH2': 4, 'CH3': 5, 'CH1-2': 6}
    col_idx = channel_map[channel]

    nrows = end - (begin + 1)
    if nrows <= 0:
        return np.array([]), [], np.array([])

    df = pd.read_csv(path, sep='\s+', header=None, skiprows=begin + 1, nrows=nrows)

    if fs <= 50:
        # 低采样率：解析每一行的时间戳
        datetime_str = '20' + df.iloc[:, 0].astype(str) + df.iloc[:, 1].astype(str)
        t_data: List[datetime] = pd.to_datetime(datetime_str, format="%Y%m%d%H%M%S.%f").to_list()
        t: List[float] = [td.timestamp() for td in t_data]
        f_1: List[float] = df.iloc[:, col_idx].astype(float).tolist()
    else:
        # 高采样率：等间隔时间序列，不解析时间戳
        f_1: List[float] = df.iloc[:, col_idx].astype(float).tolist()
        t: List[float] = list(np.linspace(0, len(f_1) / fs, len(f_1)))
        t_data = t

    return np.array(t) - t[0], t_data, np.array(f_1) - f_1[0]


def keysight_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """Keysight频率计数器数据读取。

    使用pandas读取，每行一个频率值。

    Args:
        path: 数据文件路径（每行一个频率值）
        fs: 采样率 (Hz)

    Returns:
        (t, f) 时间序列和频率数据
    """
    df = pd.read_csv(path, header=None)
    f: List[float] = df.iloc[:, 0].astype(float).tolist()
    t: List[float] = [i / fs for i in range(len(f))]
    return t, f


def keysight_six_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """Keysight六位半万用表数据读取（制表符分隔）。

    使用pandas读取，两列数据：时间、电压。

    Args:
        path: 数据文件路径（两列：时间、电压）
        fs: 采样率 (Hz)

    Returns:
        (t, f) 时间和电压数据
    """
    df = pd.read_csv(path, sep='\t', header=None)
    t: List[float] = df.iloc[:, 0].astype(float).tolist()
    f: List[float] = df.iloc[:, 1].astype(float).tolist()
    return t, f


def sim_keysight_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """模拟Keysight数据读取（逗号分隔，第二列数据）。

    使用pandas读取，跳过第一行（表头），读取第二列数据。

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)

    Returns:
        (t, f) 时间和数据序列
    """
    df = pd.read_csv(path, skiprows=1, header=None)
    f: List[float] = df.iloc[:, 1].astype(float).tolist()
    t: List[float] = [(i + 1) / fs for i in range(len(f))]
    return t, f


def labview_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """LabView自编程序数据读取。

    使用pandas读取，第一列时间戳，第二列温度。

    Args:
        path: 数据文件路径（第一列时间戳，第二列温度）
        fs: 采样率 (Hz)

    Returns:
        (t, T) 时间和温度数据
    """
    df = pd.read_csv(path, sep='\t', skiprows=1, header=None)
    t_data = pd.to_datetime(df.iloc[:, 0], format="%Y-%m-%d %H:%M:%S.%f")
    t: List[float] = t_data.astype('int64') / 1e9  # nanoseconds → seconds
    T: List[float] = df.iloc[:, 1].astype(float).tolist()
    return t, T


def data_read_SR780_dbm(
    path: str,
    label1: str = '123',
    d_k: float = 1E-3,
) -> Tuple[List[float], List[float]]:
    """SR780频谱分析仪数据读取，纵坐标为dBm。

    使用pandas读取，两列：频率、dBm。

    Args:
        path: 数据文件路径（两列：频率、dBm）
        label1: 图例标签
        d_k: 系数

    Returns:
        (nu_0, p_v2_Hz) 频率和功率谱密度
    """
    df = pd.read_csv(path, sep='\t', skiprows=3, header=None)
    nu_0: List[float] = df.iloc[:, 0].astype(float).tolist()
    p_dbm: List[float] = df.iloc[:, 1].astype(float).tolist()

    df_val = nu_0[1] - nu_0[0]
    p_v2_Hz = np.array(p_dbm) / df_val / (d_k * 429.228E12) ** 2

    return nu_0, p_v2_Hz


def oscilloscope_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """示波器频率数据读取。

    使用pandas读取，逗号分隔，第二列为频率。

    Args:
        path: 数据文件路径（逗号分隔，第二列为频率）
        fs: 采样率 (Hz)

    Returns:
        (t, f) 时间和频率数据
    """
    df = pd.read_csv(path, skiprows=1, encoding='gbk', header=None)
    # 限制到前500万行（与原代码 userlines[0:50e5] 一致）
    df = df.iloc[:int(50e5)]
    f: List[float] = df.iloc[:, 1].astype(float).tolist()
    t: List[float] = [i / fs for i in range(len(f))]
    return t, f


def pn_plot_to_psd(path: str) -> Tuple[List[float], List[float], List[float]]:
    """相噪仪相噪读取及PSD转换。

    使用pandas读取，逗号分隔：频率, 相噪dBc/Hz。

    Args:
        path: 数据文件路径（逗号分隔：频率, 相噪dBc/Hz）

    Returns:
        (nu, pn, psd) 频率, 相噪(dBc/Hz), 功率谱密度
    """
    df = pd.read_csv(path, encoding='gbk', header=None)
    nu: List[float] = df.iloc[:, 0].astype(float).tolist()
    pn: List[float] = df.iloc[:, 1].astype(float).tolist()

    pn_rad = 2 * 10 ** (np.array(pn) / 10)
    psd = (pn_rad * np.array(nu) ** 2).tolist()
    return nu, pn, psd
