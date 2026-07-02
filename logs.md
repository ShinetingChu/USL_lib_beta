# 更新日志 / Changelog

## [0.1.1] - 2026-07-02

### Bugfix

1. **修复 `three_cornered_hat()` 调用 `solve_equ` 传参错误**: 修复了三棱锥法中方程求解的传参顺序问题。
2. **移除 `__pycache__/` 目录**: 删除版本控制中的编译缓存文件，添加 `.gitignore`。
3. **包文件微调**: 重构部分代码结构，优化导入方式。

## [0.1.0] - 2025-06-07

### 重构为结构化包 (Major Refactor)

将原有的单文件 `USL_lib_beta/__init__.py`（~58KB 单体模块）重构为按功能拆分的子模块包，存入 `ultra_stable_laser/` 目录下。

#### 子模块划分

| 模块 | 功能 | 包含函数 |
|------|------|----------|
| `_io.py` | 数据读取 | `KK_data_read`, `KK_data_read_single`, `keysight_data_read`, `keysight_six_data_read`, `sim_keysight_data_read`, `labview_data_read`, `data_read_SR780_dbm`, `oscilloscope_data_read`, `pn_plot_to_psd` |
| `_allan.py` | Allan方差计算 | `allan_adev`, `allan_oadev`, `allan_mdev`, `allan_hdev`, `allan_psd`, `three_cornered_hat`, `way_one_plus_and_minus`, `solve_equ` |
| `_psd.py` | PSD/CSD功率谱密度 | `psd_welch`, `psd_int_allan`, `plot_csd`, `calc_psd_single` |
| `_drift.py` | 漂移补偿 | `move_long_drift` |
| `_utils.py` | 工具函数 | `sin_func`, `line_fit`, `calculate_slope` |
| `_transfer.py` | 温度传递函数 | `transfer_temp` |
| `_plot.py` | 画图函数 | `plot_data_labview`, `temp_read_psd_allan_2`, `plot_temp_stability`, `plot_pico_USB_err`, `plot_keysight_USB_power`, `plot_keysight_six_half_USB_power`, `plot_sim_keysight_USB_power`, `K_K_single_plot`, `K_K_plot`, `K_K_plot_path1_path2`, `freq_disc_slope`, `SR780_data_concatenate` |

#### 主要改动

1. **包名变更**: `USL_lib_beta` → `ultra-stable-laser` (PyPI: `ultra_stable_laser`)
2. **代码重构**: 所有功能从单体 `__init__.py` 拆分到独立子模块
3. **类型注解**: 所有函数添加完整类型注解 (type hints)
4. **文档字符串**: 所有函数添加详细 docstring（中英文参数说明）
5. **性能优化**: `_io.py` 中全部数据读取函数从逐行 `for` 循环改写为 `pandas.read_csv()` 向量化读取
6. **重复代码合并**:
   - `plot_data_labview` + `temp_read_psd_allan_2` → `plot_temp_stability()`，保留向后兼容包装
   - `K_K_plot` / `K_K_plot_path1_path2` 通道分支参数化 (`_CHANNEL_CONFIGS` + `_plot_kk_core`)，大幅精简重复代码
7. **死代码清理**: 删除被注释的 `plt.savefig`、`plt.ylim`、`plt.xlim` 等代码块
8. **无用导入清理**: 删除 `from cProfile import label`、重复的 `numpy` / `pyplot` / `time` / `datetime` 导入
9. **全局配置优化**: rcParams 封装到 `configure_plotting()` 函数，支持 `scienceplots` 样式 fallback
10. **向后兼容**: `__init__.py` 导入所有子模块函数，保持原有 `import ultra_stable_laser as usl` 用法不变

#### 新增公共函数

- `configure_plotting()` — 配置 matplotlib 绘图参数
- `datetime_to_epoch()` — 将 datetime 对象转换为 Unix 时间戳
- `plot_temp_stability()` — 合并后的温度稳定性绘图函数
- `calc_psd_single()` — FFT 法单边 PSD 计算（原 `PSD` 函数，保留别名）

#### 数据读取性能提升 (pandas 改造)

| 函数 | 原实现 | 新实现 | 预期性能 |
|------|--------|--------|----------|
| `KK_data_read` | 逐行 `readlines` + `strptime` | `pandas.read_csv` + 列索引 | 10x-50x (大文件) |
| `KK_data_read_single` | 逐行 + 条件分支 | `pandas.read_csv` + 通道映射 | 10x-50x |
| `keysight_data_read` | 逐行 + `for i in range` | `pandas.read_csv` | 5x-20x |
| `keysight_six_data_read` | 逐行 + `split('\t')` | `pandas.read_csv(sep='\t')` | 5x-20x |
| `sim_keysight_data_read` | 逐行 + `split(',')` | `pandas.read_csv` | 5x-20x |
| `labview_data_read` | 逐行 `strptime` | `pandas.to_datetime` 批量 | 10x-50x |
| `data_read_SR780_dbm` | 逐行从第4行起 | `pandas.read_csv(skiprows=3)` | 5x-20x |
| `oscilloscope_data_read` | 逐行 + range 限制 | `pandas.read_csv(nrows=5E6)` | 5x-20x |
| `pn_plot_to_psd` | 逐行2列 + 手动转换 | `pandas.read_csv` + 向量化 | 5x-20x |

#### 绘图函数 pandas 改造

- `plot_temp_stability`: 内联文件读取改用 `pd.read_csv` + `pd.to_datetime`
- `plot_pico_USB_err`: 内联文件读取改用 `pd.read_csv(skiprows=3)`
