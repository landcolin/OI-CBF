# OI-CBF MATLAB 实现与数学审查

原项目：David E. J. van Wijk 等人的 OI-CBF Python 开源代码。对应论文为
*Safety-Critical Control with Bounded Inputs: A Closed-Form Solution for Backup Control Barrier Functions*，arXiv:2510.05436v1。

本目录实现双积分器和固定翼飞机的动力学、标称/备份控制、bCBF、函数混合、OI、仿真和绘图。
原有 Python 源文件与 `aircraft/matlab/` 中的符号脚本保持原样。

## 运行

已实际验证环境：Windows + MATLAB R2026b。仿真只使用基础 MATLAB，无需 Python、Optimization Toolbox、Symbolic Math Toolbox 或外部 LaTeX。
源码使用 `tiledlayout`、`exportgraphics`、`jsonencode` 等现代 MATLAB 功能；未实测其他发行版或 GNU Octave。

将 MATLAB 当前目录切到项目根目录，然后执行：

```matlab
addpath('matlab');
resultsDI = run_comparison('double_integrator', 'source');
resultsAircraft = run_comparison('aircraft', 'source');
```

每个入口运行 bCBF、blending、OI 三组仿真，默认保存到 `matlab/results/<model>/<profile>/`：

- `comparison.mat`：状态、控制、混合系数、备份轨迹、耗时与配置。
- `controls_and_safety.png`：安全函数、混合系数、全部输入及其界限。
- `trajectories.png`：相平面或水平飞行轨迹。
- `state_histories.png`：全部状态的时间历程。

只运行一种控制器：

```matlab
p = oicbf.config('aircraft', 'source');
r = oicbf.simulate(p, 'OI');
plot(r.t, r.h); yline(0);
```

### 使用 Optimization Toolbox

已确认本机 `quadprog` 可用，并已完成独立交叉验证。默认继续使用无需工具箱的 `projection`；可显式切换 bCBF 求解器：

```matlab
p = oicbf.config('aircraft', 'source');
p.qpSolver = 'quadprog';
r = oicbf.simulate(p, 'bCBF');
```

一次运行三种控制器并使用 `quadprog` 求解 bCBF：

```matlab
run_comparison('double_integrator', 'source', [], 'quadprog');
run_comparison('aircraft', 'source', [], 'quadprog');
```

OI 在两种设置下都使用论文的闭式解；`quadprog` 用作其独立验证参照。`solveQP` 将本项目的 `C*u>=d` 转成 MATLAB 的 `-C*u<=-d`，目标取 `H=I, f=-up`，采用 active-set 算法并检查退出状态和残差。

`validate_toolbox` 检查 80 组独立 QP、OI 闭式解和两个完整 bCBF 仿真；`run_validation` 检测到工具箱时自动包含这项验证。结果见 `test_results/matlab/toolbox_validation.json`。

### 使用 Symbolic Math Toolbox

```matlab
symbolicReport = validate_symbolic();
```

该检查独立符号求导飞机闭环动力学、两种配置的六项备份约束，并核对数值雅可比/梯度；另外精确验证 soft-min 梯度恒等式与双积分器雅可比。结果为 `test_results/matlab/symbolic_validation.json`，同时生成可检查的 `symbolic_aircraft_jacobian.m`。

符号推导在验证阶段使用，运行仿真仍用现有稳定数值函数。符号核验不等于已经证明备份集不变性或修复了论文的半径问题。

### 使用 Parallel Computing Toolbox

```matlab
parallelReport = validate_parallel();
batch = run_batch_tests(6, 'parallel', 2);
```

`run_batch_tests` 每个模型生成指定数量的初始状态，每个状态运行三种控制器。例如上面的 6 表示每模型 6 个状态、合计 36 个场景。每场景推进 40 个控制步，是短时批量探测，不代替完整时域的 `validate_simulations`。状态中显式包含已知边界反例，随机部分固定种子，以便重复比较。

`validate_parallel` 对 18 个场景先串行后用两个进程 worker 并行，比较完成/不可行/错误分类及数值结果。报告的 `passed` 指串行、并行的一致性，不表示所有场景安全。本次结果为 14 个完成、4 个不可行；完成场景中有 1 个飞机 blending 已知边界反例发生围栏越界，已单独记录 `safetyViolations=1`。

结果见 `batch_serial.json`、`batch_parallel.json`、`parallel_validation.json`。首次创建并行池有启动开销；本次小批量含启动后的总耗时高于串行。批量规模增大后可再测实际收益。函数会关闭自己创建的池，保留用户已经存在的池。

`run_validation` 检测到这两个工具箱时自动运行对应验证。

## 两种公式配置

`source` 是默认值，保留 Python 原始参数与数学选择，以便跨语言复现。它包含原代码的航向误差绝对值、高度约束归一化、OI 终端增益选择、离散收紧、跳过预测时刻零等。

`paper` 用于核对论文的字面公式：使用有符号航向误差、未归一化的高度约束、终端 `alphaB`、包含预测时刻零、不加入论文未写出的源码收紧项。论文没有给出所有仿真增益、阈值与 softplus 的精确表达式；这些仍取原项目的数值实现。因此该选项不是论文所有数值细节的独立重建，也不意味着安全定理已被证明。

```matlab
p = oicbf.config('aircraft', 'paper');
[C, d] = oicbf.constraintSystem(p, p.x0, oicbf.control(p,p.x0), 'OI');
```

**两种配置都保留论文式 (33) 的原始半径公式**。审查发现：给定负转弯角会产生负 `rho`，使论文写出的 `C_B` 不一定包含于安全集；这不是换语言能够解决的问题。详细反例见完整报告。

显式尝试几何公式修订时：

```matlab
p.radiusConvention = 'signed-corrected';
```

它将 `h6` 改为 `h + rho*(n_g'*n(psi)) - abs(rho) - c6`。这是审查提出的修订，未冒充原论文公式，也未证明有限邻域备份集或整个采样控制器的安全性。

## 数值与异常处理

- 飞机使用 `ode45`；双积分器使用同一 ODE 的解析解。
- MATLAB 按列展开 STM；Python 的 STM 按行展开。两者通过各自一致的重塑方式表达相同的数学矩阵，测试已验证矩阵及方向导数。
- bCBF 的一维/二维欧氏投影与原二次规划等价。二维实现对输入矩形依次做半平面裁剪，再求最近边界点；已用 80 组 Python quadprog 参考问题对照。
- softplus 与 soft-min 使用代数等价的稳定表达式，避免原代码的指数溢出；保留备份饱和函数中的 `1e-6` 正则项。
- 默认 `p.failurePolicy='error'`：约束不可行时明确抛出 `oicbf:Infeasible`。这与原 Python OI 不检查可行性、可能输出超限控制的行为不同。
- `p.failurePolicy='backup'` 返回备份输入并标记失败；此回退不等于恢复论文的安全保证。
- 仅为复现原始缺陷，可显式设置 `p.failurePolicy='legacy'`。务必检查 `info.feasible` 和残差。
- 耗时统计在 MATLAB 中实际计算；没有沿用 Python 的空统计数组。全程序墙钟耗时不能直接当作论文算法性能指标。

## 验证

项目当前已提供本次生成的参考 MAT 文件。直接运行：

```matlab
addpath('matlab');
report = run_validation();
```

包含公式/梯度/雅可比检查、原 Python 单点和备份流对照、独立 QP 对照、论文配置的独立约束重建、六组完整仿真与轨迹比较。结果保存在 `test_results/matlab/`。

需要重新从 Python 生成参考文件时，在项目根目录按顺序运行：

```powershell
.test-environment/Scripts/python.exe tests/validate_source.py double_integrator
.test-environment/Scripts/python.exe tests/validate_source.py aircraft
.test-environment/Scripts/python.exe tests/audit_paper_math.py double_integrator
.test-environment/Scripts/python.exe tests/audit_paper_math.py aircraft
.test-environment/Scripts/python.exe tests/export_qp_fixtures.py
```

`.test-environment/` 是当前机器的独立 Python 测试环境。换机器可自行创建环境并安装根目录 `requirements.txt` 后使用相应的 Python 路径。

完整逐式核查、偏差解释、反例与实测数值见 [数学一致性与 MATLAB 报告](../test_results/MATH_AND_MATLAB_REPORT.md)。

## 文件对应

| MATLAB | 原 Python / 论文 |
|---|---|
| `+oicbf/config.m` | setupDynamics、setupControl、setupASIF 的参数 |
| `+oicbf/dynamics.m` | dynamics.py 的 f_x/g_x；式 (16)、(30) |
| `+oicbf/control.m`、`softSaturate.m` | control.py；式 (32) |
| `+oicbf/jacobian.m` | computeJacobianSTM/CL；式 (12)、(29) |
| `+oicbf/barriers.m` | safety.py 的 h、hb 和梯度；式 (31)、(33) 后的软最小值 |
| `+oicbf/backupFlow.m` | 三种备份流积分；式 (8)、(12)、(29) |
| `+oicbf/constraintSystem.m` | bCBF/OI 约束；式 (10)、(17)、(19)、(20) |
| `+oicbf/filter.m` | bCBF、blending、OI；式 (13)–(15)、(18)、(21) |
| `+oicbf/projectPolytope.m` | quadprog 的等价一维/二维求解 |
| `+oicbf/solveQP.m` | 投影求解器或 MATLAB quadprog 的统一接口 |
| `+oicbf/simulate.m` | main_sim.py；每步零阶保持控制 |
| `run_comparison.m`、`plotComparison.m` | comparison.py 和原绘图脚本的原生图形替代 |

图形用于展示相同仿真量，不逐像素复制论文排版；未重新生成原有符号雅可比脚本，也不需要它参与运行。
