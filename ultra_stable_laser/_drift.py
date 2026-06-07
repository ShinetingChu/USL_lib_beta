#!/usr/bin/env python
# coding: utf-8

"""漂移补偿模块 - 长漂去除和线性补偿。"""

from typing import List, Tuple
import numpy as np


def move_long_drift(
    t_a: List[float],
    f_1: List[float],
    switch: int,
    length: int = 30,
) -> Tuple[List[float], List[float], float]:
    """去除长漂（线性漂移补偿）。

    switch=1时启用长漂补偿，根据首尾段的平均值斜率进行线性修正。

    Args:
        t_a: 时间序列
        f_1: 原始数据序列
        switch: 漂移补偿开关 (0=关闭, 1=打开)
        length: 首尾用于计算漂移的数据点数

    Returns:
        (t_a, f_1_d, d_13) 时间序列, 补偿后数据, 漂移率 (Hz/s)
    """
    f_1_d: List[float] = []
    d_13 = -(np.mean(f_1[-1 * length:-1]) - np.mean(f_1[0:length])) / (
        np.mean(t_a[-1 * length:-1]) - np.mean(t_a[0:length])
    )
    for l in range(0, len(t_a)):
        if switch >= 1:
            f_1_d.append(f_1[l] + d_13 * (t_a[l] - t_a[0]))
        else:
            f_1_d.append(f_1[l])
    return t_a, f_1_d, d_13
