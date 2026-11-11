# Python—论文—MATLAB 数学一致性审查报告

日期：2026-11-08。论文依据：用户提供的 8 页 PDF，arXiv:2510.05436v1（2025-10-06）。本报告针对该版本与当前本地源码，不推断后续论文版本是否已修正。

## 三项任务的结论

1. **Python 与论文的核心动力学、混合控制结构、敏感度方程和 OI 闭式公式对应，但不是逐项完全一致。** 飞机标称航向控制、高度备份约束、OI 终端增益以及离散化处理有明确差异；闭式解的可行性前提在源码中未检查。
2. **已完成可独立运行的 MATLAB 实现**，覆盖两个模型及 bCBF、blending、OI 三种控制器，包含数值积分、全部状态/控制记录、耗时统计和原生 MATLAB 绘图。原有 Python 源文件保持不变。
3. **MATLAB 的公式及跨语言数值检查通过，但不能据此宣称论文的全部安全性结论成立。** 本次明确复现了论文式 (33) 负转弯半径的反例，并通过独立数值构造核对 MATLAB 论文公式配置。实现正确、复现正确与安全定理前提成立，是不同的结论。

## 验证方法

- 阅读论文相关完整页面；人工核对了第 7 页式 (31)–(34)、邻域备份集与标称控制说明的实际 PDF 排版。
- 动力学使用独立抄写的式 (16)、(30) 作为参照。
- 梯度与闭环雅可比使用中心有限差分校验，避免只比较两个相同的解析实现。
- 分别积分完整 STM 与方向导数，检验式 (12)、(29) 的等价关系。
- MATLAB `paper` 配置的 OI 约束由另一套“字面安全函数 + 有限差分梯度 + 完整 STM”的代码重建，独立于生产代码的方向导数实现。
- 用原 Python 的 quadprog 解作为独立参考，检验 MATLAB 的一维/二维投影求解器。
- 完整运行两种语言的两个模型、三种控制器，逐个状态分量比较轨迹与控制输入。
- Python 源文件 SHA-256 记录在各模型的 `*_paper_audit.json`，验证期间未修改它们。

## 逐式对应表

“一致”表示核对的代数/数值实现相符，不自动包含全状态域安全证明。“近似”表示有明确数值离散化。“偏差”表示数学表达或参数选择不同。

| 论文公式/内容 | Python 对应 | 结论 | MATLAB 对应 |
|---|---|---|---|
| (1)、(6)：仿射动力学、备份闭环 | 两个 dynamics.py 的 f_x/g_x 与备份传播 | 一致 | dynamics、control、backupFlow |
| (2)–(4)：传统 CBF 与无界输入闭式解 | 本项目示例未单独实现传统 CBF 控制器 | 背景公式，非缺失的示例功能 | 未新增无关控制器 |
| (5)、(7)：备份集、隐式安全集 | 安全函数沿有限备份轨迹取样，终端 hb 检查 | 数值近似；集合前提需要另行核实 | barriers、backupFlow、constraintSystem |
| (8)：备份流初值问题 | integrateStateBackup | 一致，solve_ivp 数值积分 | backupFlow：飞机 ode45、双积分器解析解 |
| (9)、Lemma 3：存在合适增益以满足可行性 | 固定 alpha、alpha_b 参数 | 不能从“存在”直接推出给定增益全域可行；已发现反例 | 默认显式检查不可行性 |
| (10)、(bCBF-QP)：安全/终端约束和最小干预目标 | vanilla_bCBF | 链式法则和 QP 目标一致；源码额外收紧并跳过预测零时刻 | constraintSystem + projectPolytope |
| (11)：总导数的链式法则 | gradh @ S @ (f+g*u) | 一致 | constraintSystem |
| (12)：STM 方程和单位初值 | STM 传播、computeJacobianSTM/CL | 数值核验一致 | jacobian、backupFlow('stm') |
| (13)–(15)：混合控制、指数混合函数 | vanilla_blending、h_I | 混合与指数公式一致；飞机 h_I 额外收紧 | filter('blending') |
| (16)：双积分器 | A、B、primary/backupControl | 一致；终端速度条件结合路径安全条件表达所需交集 | dynamics、control |
| (17)、(18)：导数插值与最优混合控制 | OI_QP | 一致 | constraintSystem、filter |
| (19c)、(OI-QP) 终端增益 alpha_b | 飞机 OI_QP 使用 alpha，bCBF 使用 alpha_b | 参数不一致；两者分别为 0.072 与 0.5 | source 保留，paper 用 alphaB |
| (19a,b,d)：OI 系数 a、b、b_b | hdot_up、hdot_ub 及梯度 | 一致，包含源码特有收紧项 | constraintSystem |
| (20)：mu ∈ [0,1] 和线性约束 | G、h 包含上下界 | 约束构建有上界，但求解后未验证满足它 | filter 验证全部残差、范围与有限性 |
| (21)–(26)：闭式解和 KKT | 只取正系数约束的最大下界 | 在 QP 可行时一致；不可行时不能继续将该值当作最优可行解 | 相同闭式解 + 可行性检查；另有 legacy 模式 |
| (27)：连续备份时域极限 | 有限离散时刻 | 未实现无限精细连续极限 | 有限网格，明确保留数值近似 |
| (28)、(29)：方向导数与其初值 | q_up、q_ub 及传播 | 数值核验一致 | backupFlow('directional') |
| (30)：八状态固定翼飞机模型和输入界限 | aircraft/dynamics.py 与 control.py | 一致 | dynamics、config |
| (31)：线性地理围栏 | h_x_geo、grad_h_geo | 一致 | barriers |
| (32)：备份转弯控制结构 | backupControl | 结构一致；具体 softplus 和 1e-6 项是源码选定细节 | control、softSaturate |
| 第 V-B 节标称控制 | primaryControl 中 abs(epsi) | 与论文有符号航向误差不同；atan2/角度环绕是合理的圆周解释 | source 保留 abs，paper 去掉 abs |
| (33) 后邻域约束 h3 | 1-((Hstar-H)/cH)^2 | 与 c3²-(Hstar-H)² 不同；单约束零集相同，soft-min 集合不同 | source 归一化，paper 字面形式 |
| (33)、h6、Proposition 2 | hb6_x、rho 计算 | 源码照抄论文，但负 rho 使论文的集合包含结论存在反例 | 两种配置均显式保留原式，另提供有名称的修订选项 |
| (34)：理想协调转弯的圆轨迹 | 源码数值积分，未直接用圆轨迹解析式 | 理想常值转弯情形可对应；softplus 精确平衡还需核实 | 同样数值积分，不把近似平衡当作恒等式 |

方向导数系统维数为 `3n`，STM 系统维数为 `n+n²`：飞机分别为 24 与 72；双积分器两者都为 6，双积分器并不减少 ODE 数量，但闭式解仍省去在线 QP 求解。

## 需要重视的数学差异和反例

### A. 飞机标称控制增加了绝对值

`aircraft/control.py:68` 用 `Kpsi*abs(epsi)`，论文标称转弯角说明使用有符号的 `Kpsi*(psi_star-psi)`。在航向误差为 `-0.5 rad` 时，同样用源码的 tanh 饱和：

- Python 期望滚转角约 `+0.785384887 rad`。
- 有符号版本约 `-0.785384887 rad`。

这属于控制器行为变化，不能当作实现细节忽略。源码注释说明绝对值是为了持续转向围栏。`atan2` 和角度环绕则是对圆周角度的合理处理；论文的 `arctan` 简写没有给出所有象限处理细节。

### B. 高度约束的归一化会改变软最小值

`aircraft/safety.py:61` 使用 `h3=1-(DeltaH/400)²`。论文写的是 `h3=c3²-DeltaH²`。当 `DeltaH=150`，取 `c3=400`：前者为 `0.859375`，后者为 `137500`。

两者各自的非负集合均为 `|DeltaH|<=400`，但只有第三项被缩放，其他五项和 kappa 不变。因此放入六项 soft-min 后，函数、梯度权重及最终备份集都发生改变，不能说完全等价。这里沿用 400 作为对照参数；论文并未公布所有阈值的具体数值。

### C. 飞机 OI 的终端增益与 bCBF 不同

论文 (19c) 使用 `alpha_b`。`aircraft/safety.py:448` 的 OI 分支使用 `alpha(hb_phi)`，第 285 行的 bCBF 分支使用 `alpha_b(hb_phi)`。本地参数为 `alpha(s)=0.072s`、`alpha_b(s)=0.5s`。

两种增益都可以属于 class-K 函数，因此这本身不等于定理必然失效；但它使源码比较的两个控制器具有不同终端条件。MATLAB `source` 保留此差异；`paper` 使用 `alphaB`。

### D. OI 的可行性前提没有由代码保证

双积分器状态 `x=[-0.128,0.5]` 的备份轨迹为 `p(t)=-0.128+0.5t-t²/2`，在 `t=0.5` 时位置最大为 `-0.003`，1 s 后速度为 `-0.5`。它可被备份控制安全制动，且路径裕量超过源码收紧量 `0.002`。

原 Python 却返回：

```text
mu = 1.020454520509091
u  = -1.0409090410181818
允许的输入范围 = [-1,1]
全部线性约束的最小残差 = -0.020454520509090912
```

在 `mu=1` 处约束残差也有负值，说明给定增益、收紧与状态组合下，论文可行性引理所需前提没有落实。公式 (21) 的正系数最大下界不能单独解决不可行 QP。

MATLAB 默认抛出 `oicbf:Infeasible`，不会把超限值当作正常输出；`failurePolicy='legacy'` 可准确复现上述原始结果。简单裁剪 mu 无法修复原约束不可行。

### E. 离散化收紧不是论文所列原式，并且参数缺少上界证据

Python 在 bCBF 与 OI 的路径约束中额外使用 `alpha(h-delta)`，其中 `delta=Delta/2*Lh*sup_fcl`；飞机 blending 的 h_I 也有该项，而双积分器 blending 没有。论文给出的相应公式未包含此项。

如果 `sup_fcl` 的含义是备份向量场范数上界，源码数值在默认初始状态就不满足该含义：

| 模型 | 配置 sup_fcl | 实际 f_cl(x0) 的范数 | 实际 abs(dh/dt) |
|---|---:|---:|---:|
| 双积分器 | 0.2 | 1.11803399 | 0.5 |
| 飞机 | 150 | 200.006237 | 200 |

因此，即使把它解释为安全函数变化率上界，也存在默认状态反例。源码没有给出其他可验证的适用域说明，不能用这些常数作为严格的连续时间收紧证书。

两种语言还都逐步保持输入，而论文主要给出连续状态反馈结果；有限网格、跳过预测时刻零以及零阶保持都需要各自的误差条件。默认仿真安全不能代替这些条件的证明。

### F. 论文式 (33) 的负转弯半径存在明确反例

第 7 页写 `rho=VT²/(g*tan(phi_star))`，并在仿真中取 `phi_star=-pi/4`，所以 `rho=-4077.471967 m`。同页使用：

```text
h6 = h + rho*(n_g'*n(psi)-1) - c6
```

其证明用到了 `rho*(n_g'*n(psi)-1)<=0`，但该不等式对负 rho 不成立。取：

```text
phi=-pi/4, theta=0, psi=0,
pN=26000, pE=0, H=10000, P=0, Nz=sqrt(2)
```

该状态满足 (33) 的五个等式，且 `h6=3072.471967>0`，但 `h=-1000<0`。按给定字面公式，它属于论文备份集却不属于安全集，直接反驳此参数下的 `C_B subset C_S` 断言。Python 和 MATLAB 均复现此代数反例。

更强的动态检验：只将上述位置改为 `pN=24990`，初始 `h=10>0`、`h6=4082.471967>0`，源码 soft-min `hb≈0.00991431>0`。沿原 Python 备份控制积分 60 s，采样得到最小 `h≈-4067.474267 m`。因此“备份集合内且初始在围栏内”的该状态也会沿备份控制越界。这个反例不依赖 OI 求解器。

由 (34) 写出轨迹圆心 `c=p+rho*n(psi)`，对单位围栏法向量，整个圆上的最小围栏裕量应为：

```text
min h = h(x) + rho*n_g'*n(psi) - abs(rho)
```

可考虑将 h6 修订为 `h + rho*n_g'*n(psi) - abs(rho) - c6`。MATLAB 仅在显式选择 `radiusConvention='signed-corrected'` 时使用该式。此修订的几何依据不等于对近似备份邻域和整个采样控制系统的完整安全证明；建议与作者确认转弯半径的符号约定。

### G. 精确平衡与源码 softplus 的微小偏差

在理想转弯平衡状态，原 Python 备份输入的第二分量比 `Nz_star` 小约 `1.414212e-6`。所以由源码的具体 softplus 与 `1e-6` 正则项，不能将论文证明中的“这些状态严格恒定”直接当作精确代码恒等式。该量很小，默认仿真仍正常，但与严格不变性论证有关。

## MATLAB 实现及可追踪的变更

目录：`matlab/`；使用说明：`matlab/README.md`。

- `source` 配置保留原始数学选择，以验证转换是否忠实。
- `paper` 配置显式调整上述标称控制、高度约束、终端增益、预测零时刻和未公布收紧项，以进行字面公式核验；未捏造论文未给出的具体增益或精确 sat 函数。
- softplus 和 soft-min 改为数值稳定的代数等价形式；不保留浮点溢出这一实现缺陷。
- 某些原 Python 当前状态的 hb/梯度会溢出，但仿真控制器评估的是备份终点，二者不矛盾。跨语言软最小值检查对非有限原值使用导出的六项安全函数独立稳定重算，没有把 NaN 的相互比较当作通过。
- OI 默认检查可行性；可显式请求 legacy 行为，便于复现反例。
- 使用解析闭环雅可比，已经有限差分与原始生成雅可比双重对照。
- bCBF 使用与原目标等价的一维/二维欧氏投影，不依赖优化工具箱。
- 备份 STM 采用 MATLAB 的按列展开；Python 按行展开。状态重塑始终内部一致，跨语言比较的是数学矩阵，不是未经转换的内存布局。
- 实际计算耗时，修复原始返回空统计数组的接口不足；这些计时不作为论文速度比的复现结论。
- 双积分器源代码实际推进 351 步（7.02 s），飞机 734 步（110.10 s）；MATLAB 刻意保留，而非无说明地裁剪到标题中的整数秒。

## 实际验证结果

实际使用 MATLAB R2026b 运行了全部测试，无需优化工具箱。验证入口为 `run_validation`，失败断言会使 MATLAB 批处理以错误结束。

### 数学与求解器

- `validate_math`：57 项检查全部通过，包括动力学、控制、梯度、雅可比、STM、方向导数、闭式解与 QP、异常处理及已知反例复现。
- `validate_projection`：80 组固定随机种子参考问题，72 组可行、8 组不可行；与 Python quadprog 的可行性判断全部一致，最大输入差 `1.05e-15`，最大约束违反量约 `4.44e-16`。
- `validate_paper_profile`：两个模型各 3 个状态，用独立的字面公式和完整 STM 重建 OI 约束。双积分器最大系数差 `1.49e-11`，飞机最大系数差 `3.86e-8`；均通过预设 `1e-5` 系数误差阈值。
- 飞机闭环雅可比与原 Python 生成雅可比的最大差 `3.55e-15`；与中心有限差分的最大差 `3.70e-9`。
- 在收紧积分容差到 `1e-10` 后，飞机方向导数与完整 STM 乘初始方向的最大差 `2.71e-8`。原始 `1e-6` 容差下两个独立积分的误差更大，不能将其混作代数公式错误。
- 最终新增 MATLAB 文件经 Code Analyzer 检查，诊断条目为 0。

“反例复现检查通过”表示测试成功检测到了被审查对象的问题，不表示对应论文性质成立。JSON 中明确保留 `paperSafetyClaimsValidated=false`。

### 六组完整 source 配置仿真

| 模型 | 控制器 | MATLAB 最小采样 h | 最大状态差（混合单位，详见 JSON 的逐分量数据） | 最大输入差 | 输入/安全检查 |
|---|---|---:|---:|---:|---|
| 双积分器 | bCBF | 0.0018 | 6.66e-16 | 5.55e-15 | 通过 |
| 双积分器 | blending | 8.56789e-7 | 2.11e-15 | 2.04e-14 | 通过 |
| 双积分器 | OI | 0.0018 | 6.66e-16 | 7.11e-15 | 通过 |
| 飞机 | bCBF | 50.338061 m | 2.96e-4 | 4.94e-6 | 通过 |
| 飞机 | blending | 12.552515 m | 1.91e-5 | 1.01e-8 | 通过 |
| 飞机 | OI | 47.970165 m | 2.07e-4 | 4.76e-6 | 通过 |

飞机上表最大状态差均发生在以米计的位置/高度分量；全六组中位置分量最大差约为 `0.000296 m`。MATLAB 与 SciPy 的自适应积分和 QP 算法不同，未要求逐位相同。

轨迹比较的阈值在查看结果前写入测试：双积分器各分量 `1e-8`；飞机位置/高度 `2 m`，角度、滚转角速度及载荷因子各 `0.002`（按各自单位）。实际差异远小于这些阈值。所有默认运行均产生有限状态和输入，未发生 QP 不可行，采样安全函数未变负，输入违反界限最多为约 `2.22e-16` 的舍入量。

六组完整仿真使用 `source` 配置。`paper` 配置通过的是独立逐式数值核验，未声称它已经复现论文全部未公开数值细节或验证了所有安全定理。

## 复现与交付位置

在 MATLAB 中将当前目录设为项目根目录：

```matlab
addpath('matlab');
run_comparison('double_integrator', 'source');
run_comparison('aircraft', 'source');
report = run_validation();
```

MATLAB 结果：

- `test_results/matlab/math_validation.json`
- `test_results/matlab/projection_validation.json`
- `test_results/matlab/paper_profile_validation.json`
- `test_results/matlab/simulation_validation.json`
- `test_results/matlab/code_analyzer.json`
- `test_results/matlab/*_bCBF.mat`、`*_blending.mat`、`*_OI.mat`
- `test_results/matlab/double_integrator/`、`aircraft/` 的全部状态、控制与轨迹图。

Python 审查及跨语言参考：

- `tests/audit_paper_math.py`、`tests/export_qp_fixtures.py`
- `test_results/aircraft_paper_audit.json`、`double_integrator_paper_audit.json`
- `test_results/aircraft_reference.mat`、`double_integrator_reference.mat`、`projection_reference.mat`

参考文件生成顺序与机器环境说明见 `matlab/README.md`。这些数值结果证明本次检查范围内的实现对应与复现质量，不构成对任意初始条件、任意离散步长、任意增益或全部备份集的安全认证。

## 补充：安装 Optimization Toolbox 后的独立验证

用户安装工具箱后，已确认 MATLAB R2026b 可以实际调用 `quadprog`。新增 `solveQP` 接口及 `validate_toolbox`，保留默认投影求解器，并支持显式指定 `p.qpSolver='quadprog'`。OI 生产控制仍采用闭式解。

- MATLAB quadprog 与 Python quadprog 的 80 组问题：可行性判断全部相同，最大控制差 `8.05e-16`。
- MATLAB quadprog 与本项目投影求解器：最大控制差 `5.00e-16`。
- OI 闭式解与 MATLAB quadprog：两个模型合计 6 个抽样状态的最大 mu 差 `1.11e-16`。
- 两个模型均用 MATLAB quadprog 完整重跑 bCBF：所有 QP 可行，输入未超限，最小采样安全函数分别为 `0.0018`、`50.338061 m`。
- 与投影算法的完整轨迹最大差：双积分器 `1.11e-16`，飞机 `1.64e-11`；飞机最大控制差 `7.54e-12`。
- 新增 MATLAB 文件的 Code Analyzer 检查无诊断。

详细数值见 `test_results/matlab/toolbox_validation.json`。该验证进一步支持两种 QP 实现等价；不会消除前文已复现的论文公式和可行性前提问题。

## 补充：符号数学与并行计算验证

用户安装 Symbolic Math Toolbox、Parallel Computing Toolbox 后，已在 MATLAB R2026b Update 1 确认实际可调用。

`validate_symbolic` 直接从动力学与备份控制表达式符号求导，使用精确十进制常数而非调用已有数值雅可比。飞机 12 个参考状态的最大雅可比误差为 `2.84e-14`，两种配置的分量安全函数梯度最大误差为 `4.55e-13`。soft-min 梯度公式和双积分器雅可比的符号恒等式均通过。生成的飞机雅可比保存在 `test_results/matlab/symbolic_aircraft_jacobian.m`，数值报告为 `symbolic_validation.json`。

`run_batch_tests` 提供可重复的短时初始条件批量测试；`validate_parallel` 使用两进程池与串行执行相同的 18 个场景，每场景推进 40 个控制步。完成/不可行/其他错误分类完全相同，最大数值差为 0；本次创建的并行池已自动关闭。

场景结果为 14 个完成、4 个检测到约束不可行、0 个其他错误。14 个完成场景中，显式加入的飞机 blending 边界反例发生围栏越界，最小 h 为约 `-1172.7524 m`；该场景初始状态为前文已讨论的 `phi=phi_star, pN=24990, Nz=Nz_star` 等理想转弯状态。它再次暴露原备份集假设问题，不是并行计算造成的数值差异。报告中的 `passed` 表示并行可重复性通过，不表示这些压力测试场景全部安全。

本次串行计算约 `3.10 s`，并行计算约 `2.36 s`，但池启动约 `14.35 s`，因此小批量总耗时更长；不据此声称现有单次仿真已经加速。

新入口、参数及输出位置见 `matlab/README.md`，结果为 `test_results/matlab/parallel_validation.json` 和两份 `batch_*.json`。代码分析器对新增代码未报告诊断。`run_validation` 已接入这两项可选工具箱验证。
