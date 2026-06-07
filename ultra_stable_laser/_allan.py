#!/usr/bin/env python
# coding: utf-8

"""Allan方差计算模块 - 提供多种Allan偏差计算和三角帽方法。"""

from typing import List, Tuple
import numpy as np
import allantools as alt
import sympy as sp


def allan_adev(
    f_x: List[float],
    fs: float,
) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray]]:
    """计算Allan偏差 (ADEV)，带置信区间。

    Args:
        f_x: 频率数据序列
        fs: 采样率 (Hz)

    Returns:
        (taus, adevs, error) 平均时间, Allan偏差, 置信区间误差
    """
    dt = 1 / fs
    t_n = np.array([1, 2, 5, 8, 10, 20, 50, 80, 100, 200, 500, 800, 1000,
                    2000, 5000, 8000, 10000, 20000, 50000, 80000]) * dt
    (taus, adevs, errors, ns) = alt.adev(f_x, rate=fs, data_type="freq", taus=t_n)
    cis = []
    for (t, dev) in zip(taus, adevs):
        edf = alt.edf_greenhall(alpha=0, d=2, m=round(t / taus[0]),
                                N=len(f_x), overlapping=False, modified=False)
        (lo, hi) = alt.confidence_interval(dev=dev, edf=edf)
        cis.append((lo, hi))
    err_lo = np.array([d - ci[0] for (d, ci) in zip(adevs, cis)])
    err_hi = np.array([ci[1] - d for (d, ci) in zip(adevs, cis)])
    error = [err_lo, err_hi]
    return taus, adevs, error


def allan_oadev(
    f_x: List[float],
    fs: float,
) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray]]:
    """计算重叠Allan偏差 (OADEV)，带置信区间。

    Args:
        f_x: 频率数据序列
        fs: 采样率 (Hz)

    Returns:
        (taus, adevs, error) 平均时间, 重叠Allan偏差, 置信区间误差
    """
    dt = 1 / fs
    t_n = np.array([1, 2, 5, 8, 10, 20, 50, 80, 100, 200, 500, 800, 1000,
                    2000, 5000, 8000, 10000, 20000, 50000, 80000]) * dt
    (taus, adevs, errors, ns) = alt.oadev(f_x, rate=fs, data_type="freq", taus=t_n)
    cis = []
    for (t, dev) in zip(taus, adevs):
        edf = alt.edf_greenhall(alpha=0, d=2, m=round(t / taus[0]),
                                N=len(f_x), overlapping=True, modified=False)
        (lo, hi) = alt.confidence_interval(dev=dev, edf=edf)
        cis.append((lo, hi))
    err_lo = np.array([d - ci[0] for (d, ci) in zip(adevs, cis)])
    err_hi = np.array([ci[1] - d for (d, ci) in zip(adevs, cis)])
    error = [err_lo, err_hi]
    return taus, adevs, error


def allan_mdev(
    f_x: List[float],
    fs: float,
) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray]]:
    """计算修正Allan偏差 (MDEV)，带置信区间。

    Args:
        f_x: 频率数据序列
        fs: 采样率 (Hz)

    Returns:
        (taus, adevs, error) 平均时间, 修正Allan偏差, 置信区间误差
    """
    dt = 1 / fs
    t_n = np.array([1, 2, 5, 8, 10, 20, 50, 80, 100, 200, 500, 800, 1000,
                    2000, 5000, 8000, 10000, 20000, 50000, 80000]) * dt
    (taus, adevs, errors, ns) = alt.mdev(f_x, rate=fs, data_type="freq", taus=t_n)
    cis = []
    for (t, dev) in zip(taus, adevs):
        edf = alt.edf_greenhall(alpha=0, d=2, m=round(t / taus[0]),
                                N=len(f_x), overlapping=False, modified=True)
        (lo, hi) = alt.confidence_interval(dev=dev, edf=edf)
        cis.append((lo, hi))
    err_lo = np.array([d - ci[0] for (d, ci) in zip(adevs, cis)])
    err_hi = np.array([ci[1] - d for (d, ci) in zip(adevs, cis)])
    error = [err_lo, err_hi]
    return taus, adevs, error


def allan_hdev(
    f_x: List[float],
    fs: float,
) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray]]:
    """计算Hadamard偏差 (HDEV)，带置信区间。

    Args:
        f_x: 频率数据序列
        fs: 采样率 (Hz)

    Returns:
        (taus, adevs, error) 平均时间, Hadamard偏差, 置信区间误差
    """
    dt = 1 / fs
    t_n = np.array([1, 2, 5, 8, 10, 20, 50, 80, 100, 200, 500, 800, 1000,
                    2000, 5000, 8000, 10000, 20000, 50000, 80000]) * dt
    (taus, adevs, errors, ns) = alt.hdev(f_x, rate=fs, data_type="freq", taus=t_n)
    cis = []
    for (t, dev) in zip(taus, adevs):
        edf = alt.edf_greenhall(alpha=0, d=2, m=round(t / taus[0]),
                                N=len(f_x), overlapping=False, modified=True)
        (lo, hi) = alt.confidence_interval(dev=dev, edf=edf)
        cis.append((lo, hi))
    err_lo = np.array([d - ci[0] for (d, ci) in zip(adevs, cis)])
    err_hi = np.array([ci[1] - d for (d, ci) in zip(adevs, cis)])
    error = [err_lo, err_hi]
    return taus, adevs, error


def allan_psd(
    t_x: List[float],
    f_x: List[float],
    label: str,
    nfft: int = 1024,
) -> Tuple[np.ndarray, np.ndarray]:
    """从时域数据计算PSD（使用plt.psd），用于Allan方差分析。

    Args:
        t_x: 时间序列
        f_x: 频率数据
        label: 图例标签
        nfft: FFT点数

    Returns:
        (Pxx, nu) 功率谱密度和对应频率
    """
    import matplotlib.pyplot as plt
    tot = t_x[-1] - t_x[0]
    num_t = len(t_x)
    dt = tot / num_t
    fs = 1 / dt
    Pxx, nu = plt.psd(f_x, NFFT=nfft, Fs=fs, detrend='mean',
                       window=np.hanning(nfft), noverlap=int(nfft * 3 / 4),
                       sides='onesided', label=label)
    return Pxx, nu


def way_one_plus_and_minus(
    taus: List[float],
    U12: List[float],
    U23: List[float],
    U13: List[float],
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """通过加/减法分离三个噪声源的Allan偏差。

    Args:
        taus: 平均时间序列
        U12: 通道1-2的Allan偏差
        U23: 通道2-3的Allan偏差
        U13: 通道1-3的Allan偏差

    Returns:
        (taus, U1, U2, U3) 各通道分离后的Allan偏差
    """
    U1: List[float] = []
    U2: List[float] = []
    U3: List[float] = []
    for i in range(0, len(U23)):
        U1.append(np.sqrt((U12[i] ** 2 + U13[i] ** 2 - U23[i] ** 2) / 2))
        U2.append(np.sqrt((U12[i] ** 2 + U23[i] ** 2 - U13[i] ** 2) / 2))
        U3.append(np.sqrt((U13[i] ** 2 + U23[i] ** 2 - U12[i] ** 2) / 2))
    return taus, U1, U2, U3


def solve_equ(
    c1: float, c2: float, c3: float, c4: float, c5: float, c6: float,
) -> List:
    """解一元六次方程。

    Args:
        c1-c6: 方程系数

    Returns:
        方程的解
    """
    x = sp.Symbol('x')
    f = c1 * x + c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6
    x = sp.solve(f)
    return x


def three_cornered_hat(
    f_1_d: List[float],
    f_2_d: List[float],
    fs: float,
) -> Tuple[List[float], List[float], List[float]]:
    """三角帽方法计算各通道独立噪声。

    以U3为参考，通过协方差分析分离三个通道的本底噪声。

    Args:
        f_1_d: 通道1频率数据
        f_2_d: 通道2频率数据
        fs: 采样率 (Hz)

    Returns:
        (u1, u2, u3) 各通道本底噪声的Allan偏差
    """
    dt = 1 / fs
    tau_s = np.array([1, 2, 5, 8, 10, 20, 50, 80, 100, 200, 500, 800, 1000,
                      2000, 5000, 8000, 10000, 20000, 50000, 80000]) * dt
    s_11: List[float] = []
    s_22: List[float] = []
    s_12: List[float] = []
    S: List[float] = []
    for i in range(0, len(tau_s)):
        d_mean_U13: List[float] = []
        d_mean_U23: List[float] = []
        s_s_1: List[float] = []
        s_s_2: List[float] = []
        s_s_12: List[float] = []
        for j in range(1, int(len(f_1_d) / int(tau_s[i] * fs))):
            d_mean_U13.append(
                np.mean(f_1_d[int(tau_s[i] * fs) * j:(j + 1) * int(tau_s[i] * fs)]) -
                np.mean(f_1_d[(j - 1) * int(tau_s[i] * fs):j * int(tau_s[i] * fs)])
            )
            d_mean_U23.append(
                np.mean(f_2_d[int(tau_s[i] * fs) * j:(j + 1) * int(tau_s[i] * fs)]) -
                np.mean(f_2_d[(j - 1) * int(tau_s[i] * fs):j * int(tau_s[i] * fs)])
            )
        for k in range(0, len(d_mean_U13)):
            s_s_1.append(d_mean_U13[k] ** 2)
            s_s_2.append(d_mean_U23[k] ** 2)
            s_s_12.append(d_mean_U13[k] * d_mean_U23[k])
        s_11.append(np.mean(s_s_1) / 2)
        s_22.append(np.mean(s_s_2) / 2)
        s_12.append(np.mean(s_s_12) / 2)
        S.append(s_11[i] * s_22[i] - s_12[i] ** 2)

    # 计算c
    c_1: List[float] = []
    c_2: List[float] = []
    c_3: List[float] = []
    c_4: List[float] = []
    c_5: List[float] = []
    c_6: List[float] = []
    for i in range(0, len(s_11)):
        c_1.append(3 * ((S[i]) ** 0.5) * s_12[i] * (s_11[i] - s_12[i]) * (s_22[i] - s_12[i]))
        c_2.append(2.25 * (S[i] ** 2) + 2 * (s_11[i] + s_22[i] + s_12[i]) * c_1[i] / (3 * (S[i] ** 0.5)))
        c_3.append(3 * (S[i] ** 1.5) * (s_11[i] + s_22[i]) + c_1[i] / 3)
        c_4.append(S[i] * (1.5 * S[i] + (s_11[i] + s_22[i] - s_12[i]) * (s_11[i] + s_22[i] + s_12[i])))
        c_5.append((S[i] ** 1.5) * (s_11[i] + s_22[i]))
        c_6.append((S[i] ** 2) / 4)

    # 解方程
    f_solve: List = []
    for i in range(0, len(c_1)):
        f_solve.append(solve_equ(c_1[i], c_2[i], c_3[i], c_4[i], c_5[i], c_6[i]))

    f_solve_deal = np.zeros((len(f_solve), len(f_solve[0])))
    for i in range(0, len(f_solve)):
        for j in range(0, len(f_solve[i])):
            if f_solve[i][j] > 0:
                f_solve_deal[i][j] = f_solve[i][j]
            else:
                f_solve_deal[i][j] = 0
    # 求最小正根
    f_min_plus: List[float] = []
    for line in f_solve_deal:
        if np.max(line) > 0:
            temp: List[float] = []
            for i in line:
                if i > 0:
                    temp.append(i)
            f_min_plus.append(np.min(temp))
        else:
            f_min_plus.append(0)

    # 计算b
    b_0: List[float] = []
    b_1: List[float] = []
    b_2: List[float] = []
    for i in range(0, len(f_min_plus)):
        b_0.append((S[i] ** 0.5) * (s_12[i] ** 2) +
                   (s_12[i] ** 2) * (s_11[i] + s_22[i]) * f_min_plus[i] +
                   (S[i] ** 0.5) * (s_12[i] ** 2) * (f_min_plus[i] ** 2))
        b_1.append(-(S[i] ** 0.5) * s_12[i] -
                   (2 * s_12[i] ** 2 + 1.5 * S[i]) * f_min_plus[i] -
                   (S[i] ** 0.5) * (s_11[i] + s_22[i]) * (f_min_plus[i] ** 2) -
                   S[i] * (f_min_plus[i] ** 3) / 2)
        b_2.append((S[i] ** 0.5) +
                   2 * (s_11[i] + s_22[i] - s_12[i]) * f_min_plus[i] +
                   3 * (S[i] ** 0.5) * f_min_plus[i] ** 2)

    # 计算a, r
    a_20: List[float] = []
    a_02: List[float] = []
    a_11: List[float] = []
    a_10: List[float] = []
    r_33: List[float] = []
    r_13: List[float] = []
    r_23: List[float] = []
    for i in range(0, len(b_0)):
        a_20.append(2 * S[i] ** 0.5 + f_min_plus[i] * s_22[i])
        a_02.append(2 * S[i] ** 0.5 + f_min_plus[i] * s_11[i])
        a_11.append(S[i] ** 0.5 - f_min_plus[i] * s_12[i])
        r_33.append(-b_1[i] / b_2[i])
        a_10.append((S[i] ** 0.5) * (2 * r_33[i] + s_12[i]))
        r_13.append(r_33[i] - a_10[i] * (a_02[i] - a_11[i]) / (a_20[i] * a_02[i] - a_11[i] ** 2))
        r_23.append(r_33[i] - a_10[i] * (a_20[i] - a_11[i]) / (a_20[i] * a_02[i] - a_11[i] ** 2))

    r_11: List[float] = []
    r_12: List[float] = []
    r_22: List[float] = []
    for i in range(0, len(r_13)):
        r_11.append(s_11[i] - r_33[i] + 2 * r_13[i])
        r_12.append(s_12[i] - r_33[i] + r_13[i] + r_23[i])
        r_22.append(s_22[i] - r_33[i] + 2 * r_23[i])

    # 计算u1, u2, u3
    u1 = [(i) ** 0.5 for i in r_11]
    u2 = [(i) ** 0.5 for i in r_22]
    u3 = [(i) ** 0.5 for i in r_33]

    return u1, u2, u3
