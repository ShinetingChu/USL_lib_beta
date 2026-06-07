#!/usr/bin/env python
# coding: utf-8

"""
超稳激光系统主要函数集合
=========================

整合了Allan方差、PSD、线漂补偿、三角帽计算、温度传递函数等。
重构为结构化包，提供向后兼容的公共API。

更新日志:
  - 20241024: 添加控温层温度读取函数，绘制控温层温度随时间变化，PSD，adev曲线
  - 20241022: 修改adev，mdev，hdev，oadev函数输入参数，修改三角帽计算adev函数，添加fs=50Hz和100Hz在1s处取点
  - 20240901: 修改画图比例，图例，字体大小等
  - 20240101: plot_keysight_USB_power, plot_pico_USB_err
  - 20231023: ad function pn_plot_to_psd, oscilloscope_data_read
  - 20231009: ad function K_K_plot
  - 20230607: 重构为结构化包（__init__.py 仅作导出，功能拆分到子模块）

子模块:
  _io        — 数据读取
  _allan     — Allan偏差计算
  _psd       — PSD/CSD
  _drift     — 漂移补偿
  _utils     — 工具函数
  _transfer  — 温度传递函数
  _plot      — 画图函数
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


# ────────────────────────────────────────────────────────────
# 全局配置
# ────────────────────────────────────────────────────────────

def configure_plotting() -> None:
    """配置matplotlib绘图样式和全局参数。

    应用 science 风格（若可用），设置字体、大小、dpi等参数。
    """
    try:
        plt.style.use(['science'])
    except (OSError, IOError):
        # science style not available - use default matplotlib settings
        pass

    plt.rcParams.update({
        "text.usetex": False,
        "font.family": "serif",
        "font.serif": "arial",
        "mathtext.fontset": "stix",
        "font.size": 12,
        "savefig.bbox": "standard",
    })
    plt.rcParams['figure.figsize'] = (6.0, 4.0)
    plt.rcParams['image.interpolation'] = 'nearest'
    plt.rcParams['image.cmap'] = 'gray'
    plt.rcParams['savefig.dpi'] = 100
    plt.rcParams['figure.dpi'] = 200


# 默认执行一次配置
configure_plotting()

# 自定义色图
upper = mpl.cm.Blues(np.arange(256))
lower = np.ones((int(256 / 4), 4))
for i in range(3):
    lower[:, i] = np.linspace(1, upper[0, i], lower.shape[0])
cmap0 = np.vstack((lower, upper))
cmap0 = mpl.colors.ListedColormap(cmap0, name='myColorMap0', N=cmap0.shape[0])

# 全局绘图尺寸和字体
col_width = 3.375  # inch (半个A4宽度)
fs = np.array([10, 9, 6.7]) * 2


# ────────────────────────────────────────────────────────────
# 导入子模块，导出所有公共函数
# ────────────────────────────────────────────────────────────

from ._io import (
    KK_data_read,
    KK_data_read_single,
    keysight_data_read,
    keysight_six_data_read,
    sim_keysight_data_read,
    labview_data_read,
    data_read_SR780_dbm,
    oscilloscope_data_read,
    pn_plot_to_psd,
)

from ._allan import (
    allan_adev,
    allan_oadev,
    allan_mdev,
    allan_hdev,
    allan_psd,
    three_cornered_hat,
    way_one_plus_and_minus,
    solve_equ,
)

from ._psd import (
    psd_welch,
    psd_int_allan,
    plot_csd,
    calc_psd_single,
)
from ._drift import (
    move_long_drift,
)
from ._utils import (
    sin_func,
    line_fit,
    calculate_slope,
)
from ._transfer import (
    transfer_temp,
)
from ._plot import (
    plot_data_labview,
    temp_read_psd_allan_2,
    plot_temp_stability,
    plot_pico_USB_err,
    plot_keysight_USB_power,
    plot_keysight_six_half_USB_power,
    plot_sim_keysight_USB_power,
    K_K_single_plot,
    K_K_plot,
    K_K_plot_path1_path2,
    freq_disc_slope,
    SR780_data_concatenate,
)

# 向后兼容：保留 PSD 函数别名
from ._psd import calc_psd_single as PSD

# __all__ 导出列表
__all__ = [
    # 配置
    'configure_plotting', 'cmap0', 'col_width', 'fs',
    # IO
    'KK_data_read', 'KK_data_read_single',
    'keysight_data_read', 'keysight_six_data_read',
    'sim_keysight_data_read', 'labview_data_read',
    'data_read_SR780_dbm', 'oscilloscope_data_read',
    'pn_plot_to_psd',
    # Allan
    'allan_adev', 'allan_oadev', 'allan_mdev', 'allan_hdev',
    'allan_psd', 'three_cornered_hat', 'way_one_plus_and_minus',
    'solve_equ',
    # PSD
    'psd_welch', 'psd_int_allan', 'plot_csd', 'PSD', 'calc_psd_single',
    # Drift
    'move_long_drift',
    # Utils
    'sin_func', 'line_fit', 'calculate_slope',
    # Transfer
    'transfer_temp',
    # Plot
    'plot_data_labview', 'temp_read_psd_allan_2', 'plot_temp_stability',
    'plot_pico_USB_err', 'plot_keysight_USB_power',
    'plot_keysight_six_half_USB_power', 'plot_sim_keysight_USB_power',
    'K_K_single_plot', 'K_K_plot', 'K_K_plot_path1_path2',
    'freq_disc_slope', 'SR780_data_concatenate',
]
