#!/usr/bin/env python
# coding: utf-8

"""数据读取模块 - 各种频率计数器和仪器的数据读取函数（基于pandas优化）。"""

from datetime import datetime
import time
from typing import List, Tuple, Optional, Union
import numpy as np
import pandas as pd
from pathlib import Path


def _parse_kk_datetime(date_str: str, time_str: str) -> datetime:
    """解析K+K频率计数器的日期时间格式"yyMMdd HHmmss.ffffff"。

    原格式无世纪前缀，统一添加'20'前缀补全为四位年份。

    Args:
        date_str: 日期字符串，如 "230101"
        time_str: 时间字符串，如 "120000.123456"

    Returns:
        解析后的datetime对象
    """
    return datetime.strptime('20' + date_str + time_str, "%Y%m%d%H%M%S.%f")


def datetime_to_epoch(dt: datetime) -> float:
    """将datetime对象转换为Unix时间戳（浮点数，含微秒）。

    Args:
        dt: datetime对象

    Returns:
        Unix时间戳（秒，含微秒小数部分）
    """
    return time.mktime(dt.timetuple()) + dt.microsecond / 1E6



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

    使用pandas快速读取，按行号范围[begin+1, end)选取数据。

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)
        begin: 开始行（0-indexed）
        end: 结束行
        CH1-CH5: 各通道在文件中的列索引（0-indexed）

    Returns:
        (t, t_data, f_1, f_2, f_3, f_4, f_5)
    """
    # 读取所有列（空格分隔），跳过表头行之前的数据
    df = pd.read_csv(
        path,
        sep=r'\s+',
        header=None,
        skiprows=begin + 1,
        nrows=end - begin - 1 if end > begin + 1 else None,
        engine='python',
    )

    # 解析日期时间（第0列=日期，第1列=时间）
    t_data = [_parse_kk_datetime(str(df.iloc[i, 0]), str(df.iloc[i, 1]))
              for i in range(len(df))]
    t = [datetime_to_epoch(dt) for dt in t_data]

    # 提取各通道频率
    f_1 = df.iloc[:, CH1].astype(float).tolist()
    f_2 = df.iloc[:, CH2].astype(float).tolist()
    f_3 = df.iloc[:, CH3].astype(float).tolist()
    f_4 = df.iloc[:, CH4].astype(float).tolist()
    f_5 = df.iloc[:, CH5].astype(float).tolist()

    return t, t_data, f_1, f_2, f_3, f_4, f_5


def KK_data_read_single(
    path: str,
    fs: float,
    begin: int,
    end: int,
    channel: str = 'CH1',
) -> Tuple[np.ndarray, Union[List[datetime], np.ndarray], np.ndarray]:
    """K+K频率计数器单通道数据读取。

    使用pandas快速读取，支持低速(fs≤50)和高速(fs>50)两种模式：
    - 低速模式：解析日期时间戳
    - 高速模式：仅读取频率值，时间等间隔生成

    Args:
        path: 数据文件路径
        fs: 采样率 (Hz)
        begin: 开始行（0-indexed）
        end: 结束行
        channel: 通道名 ('CH1'→col3, 'CH2'→col4, 'CH3'→col5, 'CH1-2'→col6)

    Returns:
        (t, t_data, f_1)
    """
    # 列索引映射
    channel_col = {'CH1': 3, 'CH2': 4, 'CH3': 5, 'CH1-2': 6}.get(channel, 3)

    # 读取数据（空格分隔，无表头），日期/时间列保持为字符串以保留前导零
    df = pd.read_csv(
        path,
        sep=r'\s+',
        header=None,
        skiprows=begin + 1,
        nrows=end - begin - 1 if end > begin + 1 else None,
        dtype={0: str, 1: str},
        engine='python',
    )

    # 提取频率值
    f_1 = df.iloc[:, channel_col].astype(float).to_numpy()

    if fs <= 50:
        # 低速模式：解析日期时间
        t_data = [_parse_kk_datetime(str(df.iloc[i, 0]), str(df.iloc[i, 1]))
                  for i in range(len(df))]
        t = np.array([datetime_to_epoch(dt) for dt in t_data])
        t = t - t[0]
    else:
        # 高速模式：等间隔时间
        t = np.linspace(0, len(f_1) / fs, len(f_1))
        t_data = t

    f_1 = f_1 - f_1[0]

    return t, t_data, f_1


def keysight_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """Keysight频率计数器数据读取。

    每行一个频率值，使用pandas快速读取。

    Args:
        path: 数据文件路径（每行一个频率值）
        fs: 采样率 (Hz)

    Returns:
        (t, f) 时间序列和频率数据
    """
    f = pd.read_csv(path, header=None, dtype=float).iloc[:, 0].tolist()
    t = [i / fs for i in range(len(f))]
    return t, f


def keysight_six_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """Keysight六位半万用表数据读取（制表符分隔）。

    Args:
        path: 数据文件路径（两列：时间、电压，制表符分隔）
        fs: 采样率 (Hz)

    Returns:
        (t, f) 时间和电压数据
    """
    df = pd.read_csv(path, sep='\t', header=None, dtype=float)
    t = df.iloc[:, 0].tolist()
    f = df.iloc[:, 1].tolist()
    return t, f


def sim_keysight_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """模拟Keysight数据读取（逗号分隔，第二列数据，跳过表头）。

    Args:
        path: 数据文件路径（逗号分隔，第二列为数据）
        fs: 采样率 (Hz)

    Returns:
        (t, f) 时间和数据序列
    """
    df = pd.read_csv(path, header=0, dtype=float)
    f = df.iloc[:, 1].tolist()
    t = [(i + 1) / fs for i in range(len(f))]  # skip header → from index 1
    return t, f


def labview_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """LabView自编程序数据读取。

    制表符分隔，第一列时间戳(格式: "%Y-%m-%d %H:%M:%S.%f")，第二列温度。
    第一行为表头。

    Args:
        path: 数据文件路径（第一列时间戳，第二列温度）
        fs: 采样率 (Hz)

    Returns:
        (t, T) 时间和温度数据
    """
    df = pd.read_csv(path, sep='\t', header=0)

    # 解析时间戳
    time_col = df.columns[0]
    t_data = pd.to_datetime(df.iloc[:, 0], format="%Y-%m-%d %H:%M:%S.%f")
    t = [datetime_to_epoch(dt.to_pydatetime()) for dt in t_data]

    T = df.iloc[:, 1].astype(float).tolist()

    return t, T

def DAQ970_data_read(path: str, fs: float,channel:int) -> Tuple[List[float], List[float], List[int]]:
    """LabView自编程序数据读取。

    制表符分隔，第一列时间戳(格式: "%Y-%m-%d %H:%M:%S.%f")，第二列温度。
    第一行为表头。

    Args:
        path: 数据文件路径（第一列时间戳，第二列温度）
        fs: 采样率 (Hz)

    Returns:
        (t, T) 时间和温度数据
    """
    df = pd.read_csv(path, sep=',', header=0)

    # 解析时间戳
    time_col = df.columns[0]
    t_data = pd.to_datetime(df.iloc[:, 0], format="%Y-%m-%d %H:%M:%S.%f")
    t = [datetime_to_epoch(dt.to_pydatetime()) for dt in t_data]

    T = df.iloc[:, channel].astype(float).tolist()

    return t, T


def data_read_SR780_dbm(
    path: str,
    label1: str = '123',
    d_k: float = 1E-3,
) -> Tuple[List[float], List[float]]:
    """SR780频谱分析仪数据读取，纵坐标为dBm。

    制表符分隔，前3行为元数据，从第4行开始为数据（两列：频率、dBm）。

    Args:
        path: 数据文件路径（两列：频率、dBm）
        label1: 图例标签
        d_k: 系数

    Returns:
        (nu_0, p_v2_Hz) 频率和功率谱密度
    """
    df = pd.read_csv(path, sep='\t', header=None, skiprows=3, dtype=float)
    nu_0 = df.iloc[:, 0].tolist()
    p_dbm = df.iloc[:, 1].tolist()

    df_ = nu_0[1] - nu_0[0]
    p_v2_Hz = np.array(p_dbm) / df_ / (d_k * 429.228E12) ** 2

    return nu_0, p_v2_Hz


def oscilloscope_data_read(path: str, fs: float) -> Tuple[List[float], List[float]]:
    """示波器频率数据读取。

    逗号分隔，第一行为表头，第一列为时间索引，第二列为频率。
    最大读取5E6行。

    Args:
        path: 数据文件路径（逗号分隔，第二列为频率，gbk编码）
        fs: 采样率 (Hz)

    Returns:
        (t, f) 时间和频率数据
    """
    df = pd.read_csv(path, encoding='gbk', header=0, dtype=float, nrows=int(50E5))
    f = df.iloc[:, 1].tolist()
    t = [i / fs for i in range(len(f))]
    return t, f


def pn_plot_to_psd(path: str) -> Tuple[List[float], List[float], List[float]]:
    """相噪仪相噪读取及PSD转换。

    逗号分隔，无表头，两列：频率(Hz)、相噪(dBc/Hz)，gbk编码。

    Args:
        path: 数据文件路径（逗号分隔：频率, 相噪dBc/Hz, gbk编码）

    Returns:
        (nu, pn, psd) 频率(Hz), 相噪(dBc/Hz), 功率谱密度
    """
    df = pd.read_csv(path, encoding='gbk', header=None, dtype=float)
    nu = df.iloc[:, 0].tolist()
    pn = df.iloc[:, 1].tolist()

    pn_rad = 2 * 10 ** (np.array(pn) / 10)
    psd = (pn_rad * np.array(nu) ** 2).tolist()

    return nu, pn, psd
