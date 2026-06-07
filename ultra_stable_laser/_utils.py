#!/usr/bin/env python
# coding: utf-8

"""工具函数模块 - 正弦拟合、线性拟合、斜率计算等。"""

from typing import List, Tuple, Any
import numpy as np
from scipy import optimize


def sin_func(
    t: np.ndarray,
    A: float,
    omega: float,
    phi: float,
    C: float,
) -> np.ndarray:
    """正弦函数（用于curve_fit拟合）。

    Args:
        t: 自变量
        A: 振幅
        omega: 角频率
        phi: 相位
        C: 常数偏置

    Returns:
        正弦函数值
    """
    return A * np.sin(omega * t + phi) + C


def line_fit(
    x: np.ndarray,
    a: float,
    b: float,
) -> np.ndarray:
    """一元一次线性函数（用于curve_fit拟合）。

    Args:
        x: 自变量
        a: 斜率
        b: 截距

    Returns:
        线性函数值
    """
    return a * x + b


def calculate_slope(
    v_in: List[float],
    k_p_1: List[float],
) -> Tuple[List[float], np.ndarray, float, float]:
    """线性拟合并计算斜率。

    Args:
        v_in: 自变量
        k_p_1: 因变量

    Returns:
        (v_in, k_p_1_fit, A1, B1) 自变量, 拟合值, 斜率, 截距
    """
    A1, B1 = optimize.curve_fit(line_fit, v_in, k_p_1)[0]
    v_in_0 = np.arange(
        np.min(v_in) - (np.max(v_in) - np.min(v_in)) / 100,
        np.max(v_in) + (np.max(v_in) - np.min(v_in)) / 100,
        (np.max(v_in) - np.min(v_in)) / 100,
    )
    k_p_1_fit = line_fit(v_in_0, A1, B1)
    return v_in, k_p_1_fit, A1, B1
