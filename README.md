# ultra-stable-laser (USL_lib_beta)

**超稳激光系统数据分析工具箱**

Data analysis toolbox for **ultra-stable laser (USL) systems**. Provides functions for reading various frequency counter / instrument data formats, computing Allan deviations, power spectral density (PSD), drift compensation, temperature transfer functions, and publication-ready plotting.

---

## 功能概述 / Overview

- **数据读取** — 支持 K+K 频率计数器、Keysight 频率计数器、Keysight 万用表、LabView、SR780 频谱分析仪、示波器、相噪仪等多种数据格式
- **Allan 方差** — 计算 ADEV, OADEV, MDEV, HDEV 及置信区间；三角帽法分离多通道本底噪声
- **PSD 分析** — Welch 法 PSD 估计、FFT 法单边 PSD、交叉功率谱密度 (CSD)
- **漂移补偿** — 长漂移线性去除 (linear drift removal)
- **温度传递函数** — 温度波动到腔体频率稳定性的传递函数分析
- **可视化** — 时域、PSD、Allan 偏差三图联动，支持多种图例配置和拟合线

## 模块说明 / Module Structure

| 模块 | 说明 | 主要函数 |
|------|------|----------|
| `ultra_stable_laser._io` | 数据读取 | `KK_data_read`, `KK_data_read_single`, `keysight_data_read`, `keysight_six_data_read`, `sim_keysight_data_read`, `labview_data_read`, `data_read_SR780_dbm`, `oscilloscope_data_read`, `pn_plot_to_psd` |
| `ultra_stable_laser._allan` | Allan偏差计算 | `allan_adev`, `allan_oadev`, `allan_mdev`, `allan_hdev`, `allan_psd`, `three_cornered_hat`, `way_one_plus_and_minus`, `solve_equ` |
| `ultra_stable_laser._psd` | 功率谱密度 | `psd_welch`, `psd_int_allan`, `plot_csd`, `calc_psd_single` |
| `ultra_stable_laser._drift` | 漂移补偿 | `move_long_drift` |
| `ultra_stable_laser._utils` | 工具函数 | `sin_func`, `line_fit`, `calculate_slope` |
| `ultra_stable_laser._transfer` | 温度传递函数 | `transfer_temp` |
| `ultra_stable_laser._plot` | 可视化 | `plot_temp_stability`, `K_K_plot`, `K_K_single_plot`, `K_K_plot_path1_path2`, `freq_disc_slope`, `SR780_data_concatenate`, `plot_pico_USB_err`, `plot_keysight_USB_power`, `plot_keysight_six_half_USB_power`, `plot_sim_keysight_USB_power` |

---

## 主要函数列表 / Key Functions

### 数据读取 (IO)

| 函数 | 用途 |
|------|------|
| `KK_data_read(path, fs, begin, end, CH1..CH5)` | K+K频率计数器多通道数据读取（5通道） |
| `KK_data_read_single(path, fs, begin, end, channel)` | K+K频率计数器单通道数据读取 |
| `keysight_data_read(path, fs)` | Keysight频率计数器单列频率数据读取 |
| `keysight_six_data_read(path, fs)` | Keysight六位半万用表两列（时间、电压）数据读取 |
| `sim_keysight_data_read(path, fs)` | 模拟Keysight数据读取（逗号分隔，第二列为数据） |
| `labview_data_read(path, fs)` | LabView温度数据读取（时间戳、温度） |
| `data_read_SR780_dbm(path, label1, d_k)` | SR780频谱分析仪dBm数据读取及PSD转换 |
| `oscilloscope_data_read(path, fs)` | 示波器频率数据读取（逗号分隔，第二列为频率） |
| `pn_plot_to_psd(path)` | 相噪仪数据读取及PSD转换 |

### Allan 偏差 (Allan)

| 函数 | 用途 |
|------|------|
| `allan_adev(f_x, fs)` | 计算Allan偏差 (ADEV)，带置信区间 |
| `allan_oadev(f_x, fs)` | 计算重叠Allan偏差 (OADEV)，带置信区间 |
| `allan_mdev(f_x, fs)` | 计算修正Allan偏差 (MDEV)，带置信区间 |
| `allan_hdev(f_x, fs)` | 计算Hadamard偏差 (HDEV)，带置信区间 |
| `allan_psd(t_x, f_x, label, nfft)` | 从时域数据计算PSD（用于Allan分析） |
| `three_cornered_hat(f_1_d, f_2_d, fs)` | 三角帽法分离三通道本底噪声 |
| `way_one_plus_and_minus(taus, U12, U23, U13)` | 加减法分离噪声源的Allan偏差 |

### PSD 分析 (PSD)

| 函数 | 用途 |
|------|------|
| `psd_welch(t_x, f_x, fs, nfft_0)` | Welch法功率谱密度估计 |
| `psd_int_allan(nu, Pxx)` | 从PSD积分计算Allan方差秒稳影响 |
| `plot_csd(f_13, f_23, fs, f0, nfft)` | 交叉功率谱密度 (CSD) 计算 |
| `calc_psd_single(t_a, T_a, fs)` | FFT法单边功率谱密度计算 |

### 漂移补偿 (Drift)

| 函数 | 用途 |
|------|------|
| `move_long_drift(t_a, f_1, switch, length)` | 线性漂移补偿 / 去除长漂 |

### 工具函数 (Utils)

| 函数 | 用途 |
|------|------|
| `sin_func(t, A, omega, phi, C)` | 正弦函数（用于curve_fit拟合） |
| `line_fit(x, a, b)` | 一元一次线性函数（用于curve_fit拟合） |
| `calculate_slope(v_in, k_p_1)` | 线性拟合并计算斜率 |

### 温度传递函数 (Transfer)

| 函数 | 用途 |
|------|------|
| `transfer_temp(t_a, T_a, t_0, tau, D_t, fs)` | 温度到频率的传递函数和Allan偏差计算 |

### 可视化 (Plot)

| 函数 | 用途 |
|------|------|
| `plot_temp_stability(...)` | 腔体温度稳定性三图（时域、PSD、Allan偏差） |
| `K_K_plot(...)` | K+K频率计数器四图（时域、PSD、ADEV、MDEV） |
| `K_K_single_plot(...)` | K+K单段数据简易绘图 |
| `K_K_plot_path1_path2(...)` | K+K双文件合并绘图 |
| `plot_pico_USB_err(...)` | Pico示波器光强稳定性绘图 |
| `plot_keysight_USB_power(...)` | Keysight频率计数器拍频稳定性绘图 |
| `plot_keysight_six_half_USB_power(...)` | Keysight六位半万用表电压稳定性绘图 |
| `plot_sim_keysight_USB_power(...)` | 模拟Keysight数据稳定性绘图 |
| `freq_disc_slope(...)` | 鉴频斜率拟合计算 |
| `SR780_data_concatenate(...)` | SR780频谱分析仪多段数据拼接绘图 |

---

## 安装方式 / Installation

### 从 PyPI 安装

```bash
pip install ultra-stable-laser
```

### 从源码安装

```bash
git clone git@github.com:ShinetingChu/USL_lib_beta.git
cd USL_lib_beta
pip install -e .
```

### 导入使用

```python
import ultra_stable_laser as usl
```

---

## 快速使用示例 / Quick Start

以下展示 K+K 频率计数器数据的完整处理流程：读取数据 → 去除漂移 → 计算Allan偏差 → 绘制稳定度三图。

```python
import ultra_stable_laser as usl

# 1. 读取 K+K 频率计数器单通道数据
fs = 100  # 采样率 100 Hz
path = "path/to/your/data.txt"
t, t_data, f = usl.KK_data_read_single(path, fs, begin=1, end=100000, channel='CH1')

# 2. 去除线性漂移
t_d, f_d, drift_rate = usl.move_long_drift(t.tolist(), f.tolist(), switch=1)
print(f"漂移率: {drift_rate:.2e} Hz/s")

# 3. 计算Allan偏差 (ADEV)
taus, adevs, errors = usl.allan_adev(f_d.tolist(), fs)

# 4. 计算功率谱密度 (PSD)
f_psd, Pxx = usl.psd_welch(t_d.tolist(), f_d.tolist(), fs)

# 5. 一键绘图（时域、PSD、ADEV、MDEV）
result = usl.K_K_plot(path, fs, CH='CH1', label1='Sample', start=1, end=100000)
```

更多示例请参考各函数文档字符串。

---

## 依赖 / Dependencies

- `numpy` — 数组计算
- `scipy` — FFT、信号处理、曲线拟合
- `matplotlib` — 数据可视化
- `pandas` — 数据读取 (IO模块)
- `allantools>=2024.4` — Allan方差计算内核
- `sympy` — 符号计算（三角帽法解方程）

可选依赖：
- `scienceplots` — 科学论文风格绘图 (已内建 fallback)

---

## 许可证 / License

[MIT License](LICENSE)
