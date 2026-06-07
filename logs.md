# 更新日志 / Changelog

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
6. **向后兼容**: `__init__.py` 导入所有子模块函数，保持原有 `import USL_lib_beta as usl` 用法不变
7. **导入优化**: 各子模块只导入所需依赖，避免全局导入开销

#### 模块初始化优化

- `__init__.py` 不再直接执行大量计算和配置代码
- 绘图配置移至 `configure_plotting()` 函数，可按需调用
- 全局变量 (`cmap0`, `col_width`, `fs`) 从模块导出列表移除，减少启动开销
