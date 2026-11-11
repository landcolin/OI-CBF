# OI-CBF 源码测试报告

测试日期：2026-11-08。结论：两个示例的默认仿真可以完成，默认轨迹未检测到安全越界、输入超限或非有限数值；扩展测试发现输入限幅处理、参数一致性和数值稳定性等问题，不能据此认为源码全面正确。

## 环境与范围

- Windows，Python 3.12.14，独立环境 `.test-environment/`。
- 已安装 `requirements.txt` 的全部锁定依赖，`pip check` 无依赖冲突。
- 原有 Python、MATLAB、README 和依赖文件未修改；新增诊断脚本和结果文件。
- 直接调用原始 `Simulation.sim()`，分别运行两个模型的 bCBF、blending 和 OI，保持原始初始条件、步长、备份时域和积分容差。
- 数值检查包括输入界限、采样安全函数、非有限数值、QP 回退、有限差分梯度/雅可比，以及闭式 OI 与独立 QP 求解器的对照。
- 绘图检查使用已保存的完整状态、控制及插值系数；未重建备份轨迹装饰线。没有运行 MATLAB 脚本。

## 默认仿真结果

| 模型 | 控制器 | 实际推进步数 | 仿真时长 | 最小采样安全函数 h | 输入超限 | QP 失败回退 |
|---|---|---:|---:|---:|---|---:|
| 双积分器 | bCBF | 351 | 7.02 s | 0.001800 | 无 | 0 |
| 双积分器 | blending | 351 | 7.02 s | 8.56789e-7 | 无 | 不适用 |
| 双积分器 | OI | 351 | 7.02 s | 0.001800 | 无 | 不适用 |
| 飞机 | bCBF | 734 | 110.10 s | 50.3381 m | 无 | 0 |
| 飞机 | blending | 734 | 110.10 s | 12.5525 m | 无 | 不适用 |
| 飞机 | OI | 734 | 110.10 s | 47.9702 m | 无 | 不适用 |

六次默认仿真全部完成，无运行期警告，状态与控制量均为有限数。数组控制列 0 是源码初始化占位列，输入检查排除了该列。

其他验证：

- 双积分器 OI 与 bCBF 的最大状态差：`3.33e-16`，符合该示例中两者等价的预期。
- 两个模型各抽取 12 个默认 OI 轨迹状态，将源码生成的全部约束交给 quadprog 独立求解。闭式解与 QP 的最大插值系数差分别为 `1.11e-16`、`5.55e-17`；均可行。
- 每个模型在初始状态附近固定随机种子抽取 25 个状态，与中心有限差分对比。双积分器雅可比/梯度误差为 0；飞机最大雅可比绝对误差 `7.59e-9`、安全函数梯度误差 `2.26e-7`。
- 双积分器 OI 将步长改为 `0.01 s` 和 `0.04 s` 后分别运行 7 s，输入均未超限。按恒定控制区间的解析抛物线检查连续时间位置，最小安全裕量分别为 `0.00095`、`0.0032`。

## 已复现的问题

### 1. 高优先级：OI 在不可行约束下仍输出控制量，可能违反输入界限

位置：`double_integrator/safety.py:335`，飞机中相同求解模式位于 `aircraft/safety.py:457`。

双积分器默认参数下，直接调用：

```python
env = Simulation(verbose=False, bCBF_QP=False, OI_QP_flag=True)
u, _, _ = env.OI_QP(np.array([-0.128, 0.5]), 1.0)
```

实际得到 `mu=1.020454520509091`、`u=-1.0409090410181818`，超出了源码定义的 `[-1,1]` 输入界限。

该状态并非已经无法制动的状态：施加备份控制 `u=-1` 时，位置为 `p(t)=-0.128+0.5t-t²/2`，在 `t=0.5` 达到最大值 `-0.003`，1 s 时速度为 `-0.5`，因此整段备份轨迹安全且终点满足备份约束。最小安全裕量 `0.003` 还大于源码收紧量 `0.002`。

问题在于：给定的增益和收紧参数并不保证此处的全部 QP 约束可行，而实现只计算正系数约束给出的最大下界，没有检查 `mu<=1`、零/负系数约束或最终约束残差。论文闭式解的可行性前提不能代替运行期检查。默认轨迹上闭式公式与 QP 一致，并不能排除此问题。

建议：增加可行性/有限性检查和明确的失败处理；重新验证参数的适用状态域。单纯将 `mu` 裁剪到 `[0,1]` 会掩盖原约束不可行，不能当作安全性修复。飞机同样缺少这些检查，但本次未复现飞机默认状态的输入超限。

### 2. 中优先级：飞机 OI 的终端约束增益与 bCBF 不一致

位置：`aircraft/safety.py:448` 使用 `self.alpha(hb_phi)`，但 bCBF 终端约束在第 285 行使用 `self.alpha_b(hb_phi)`；两者分别为 `0.072*x` 与 `0.5*x`。

例如默认初始状态的备份终点 `hb≈0.00991204`，源码 OI 终端右端为 `-0.000713667`，改用 `alpha_b` 则为 `-0.00495602`。这会改变控制器可行域及与 bCBF 的比较条件。

这是一处需要确认作者意图的参数不一致。两个函数都属于可用的增益形式，不能仅凭名称差异认定它破坏了安全定理。在本次抽取的 3 个状态中，终端约束没有决定最优值，改用 `alpha_b` 未改变当步 `mu`。

### 3. 中优先级：飞机非默认工况下的指数运算会溢出

位置：`aircraft/control.py:104` 的 softplus 饱和，以及 `aircraft/safety.py:117`、第 131–134 行的 soft-min 和梯度。

- 从默认状态仅将滚转角速度 `P` 改为 `-30`，备份控制第一个分量为 `nan`；这是极端压力测试，不属于默认飞行轨迹。
- 从默认状态仅将高度 `H` 改为 `12000`，soft-min 为 `-inf`，梯度含非有限数，触发 `overflow encountered in exp`。该状态在地理围栏内，但不在备份集内；此测试针对函数自身的数值稳定性，并不声称此处具有论文的闭环安全保证。

建议使用稳定 softplus（例如 `logaddexp`）及 log-sum-exp/归一化权重，并检查非有限数传播。默认 6 次仿真没有出现此问题。

### 4. 低优先级：控制器配置缺少校验

位置：两个 `main_sim.py` 的控制器分支。

`safety_flag=True` 且 `bCBF_QP=False, regFun=False, OI_QP_flag=False` 时，第一次控制迭代抛出 `UnboundLocalError`，因为 `sdt` 未赋值。多个标志同时为真时也只会隐式选择第一个分支。

建议在初始化时要求启用安全过滤时恰好选择一种模式，并给出明确的配置异常。

### 5. 低优先级：耗时统计没有计算或返回

位置：`double_integrator/main_sim.py:68`、`aircraft/main_sim.py:97`。

虽然收集了求解耗时，两种模型返回的 `avg_solver_t` 和 `max_solver_t` 在六次仿真中始终为 `[]`，无法用于复现平均/最大求解耗时比较。诊断脚本的整段墙钟耗时也不能替代论文所比较的求解器耗时。

## 绘图与环境限制

- 两个对比入口都显式设置 `latex_plots=True`。当前系统缺少 LaTeX 可执行程序，实际绘图测试均报：`Failed to process string with tex because latex could not be found`。
- 使用同一份完整仿真数组，将参数设为 `latex_plots=False`，两个模型的默认图类型均成功生成 PDF。这是环境依赖问题，不是仿真控制逻辑错误。
- Python 3.12 对原有若干绘图标签中的反斜杠字符串给出 `SyntaxWarning`；目前不阻止解析及关闭 LaTeX 后的绘图。

## 重现

在项目根目录执行 PowerShell 命令：

```powershell
.test-environment/Scripts/python.exe tests/validate_source.py double_integrator
.test-environment/Scripts/python.exe tests/validate_source.py aircraft
.test-environment/Scripts/python.exe tests/validate_edges.py double_integrator
.test-environment/Scripts/python.exe tests/validate_edges.py aircraft
```

`validate_edges.py` 需要先运行对应的 `validate_source.py`。这些是诊断脚本：记录已知失败和通过结果，不把“发现缺陷”转换为进程非零退出码；应读取 JSON 和报告判断结果。

- `double_integrator.json`、`aircraft.json`：默认仿真及数值检查。
- `double_integrator_edges.json`、`aircraft_edges.json`：边界、参数差异和绘图测试。
- `*_bCBF.npz`、`*_blending.npz`、`*_OI.npz`：完整状态、控制、插值系数数组。
- `double_integrator_plots/plots/`、`aircraft_plots/plots/`：关闭 LaTeX 后生成的 PDF。

本次未修改控制算法或自动修复问题，也未穷举状态域、证明所有条件下的安全性，或逐项复现论文的全部实验。
