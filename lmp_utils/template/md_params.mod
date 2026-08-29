# md_params.mod — 升温 MD 的全部可调参数（相列表 + 温度阶梯 + 采样节奏）
#
# 运行期由 mymetal.ml.n2p2.workflow_md.PeiN2p2MD.submit_md 在每个作业目录重新写入；
# 此处这份是**默认值 + 文档**，单独手跑 md_heating_template.in 时直接可用。
#
# 相列表用 index 风格变量：LAMMPS 的 clear 不清除 input script variables
# （见 doc/src/clear.rst），所以 md_heating_template.in 才能用一个 in 文件、
# 靠 next + jump 跨 clear 依次跑完 FCC 与 HCP 两个体相。
# 六个 index 变量必须**等长且一一对应**，next 一次性推进它们。

# 相名（只用于输出文件名）/ stretch_model.mod 的 lat 编号（1=hcp 2=fcc 3=bcc）
variable phase index fcc  hcp
variable latid index 2    1

# 每相的初始晶格常数猜测（Å）：升温前先做 box/relax 弛豫到该势函数自己的平衡值，
# 故这里只需大致正确。Au: FCC a0≈4.08、HCP a0≈2.88（与 pei_lmp_run_properties 的
# ini_aa 同源）。HCP 的 c/a 由 stretch_model.mod 按理想值 1.633 给出。
variable aa0   index 4.08 2.88

# 超胞复制倍数：FCC 常规胞 4 原子 -> 2x2x2 = 32 原子；HCP 六方胞 2 原子 -> 3x3x2 = 36 原子。
# 这些帧最终要送 DFT 静态标记，原子数直接决定 DFT 代价，也要与训练集口径相称
# （stage1/0 把训练结构规整到 [20,30] 逼近 24）。
variable nx    index 2    3
variable ny    index 2    3
variable nz    index 2    2

# 温度阶梯：T_k = t_start + (k-1)*t_step，共 n_temp 档，逐档升温（不重新初始化速度）。
# 默认沿用 chap_4_TL.pdf Stage I 的 25 -> 450 K / 25 K 一档 / 共 18 档。
# ⚠️ 该设置来自 Mg（Tm≈923 K）；Au 的 Tm≈1337 K，同比例应到 ~650 K。
#    要按 Au 重标定就改这三个数，模板其余部分无需动。
variable t_start equal 25.0
variable t_step  equal 25.0
variable n_temp  equal 18

# 采样节奏（步）：每档先 n_equil 步平衡（不采样），再 n_prod 步采样，每 n_dump 步落一帧。
# 默认每档 20 帧 -> 单相 18*20 = 360 帧 -> 单个作业（FCC+HCP）720 帧。
variable n_equil equal 2000
variable n_prod  equal 10000
variable n_dump  equal 500

# 时间步长（metal 单位为 ps）与 NPT 阻尼时间；阻尼取 100*dt / 1000*dt 的常规值。
variable dt      equal 0.001
variable t_damp  equal 0.1
variable p_damp  equal 1.0
