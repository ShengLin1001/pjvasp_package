"""
planar submodule

平面缺陷物性的**唯一配置来源**：算哪些面的表面能、算哪些层错、每个层错的初猜矢量。

三个下游各自读这里，避免同一张表在三处漂移：
    - ``lmp_utils/lmp_universal/pei_lmp_run_properties``（bash）通过
      :func:`format_planar_table` 取纯文本表后按行 ``read`` 进 case 循环；
    - ``lmp_utils/post/{ssf,surface,seqsf}.py`` 直接 import 字典遍历；
    - ``mymetal/ml/n2p2/workflow.py`` 用 :func:`get_planar_columns` 生成 epoch 扫描列名。

Functions:
    - format_planar_table: 把某一类的配置渲染成每行一条记录的纯文本（供 bash 读）。
    - get_planar_columns: 某一类在 epoch 扫描表里的列名顺序。
"""

# 稳定层错（tilted-cell 单点全弛豫）：(相, 滑移面模型, 层错名, t1, t2)。
#
# t1/t2 是层错矢量在面内基矢 a1/a2 上的分数，只是**初猜**：ssf.mod 放开 xz/yz，
# 弛豫会把它推到真实位置（post 写出的 t1/t2 与 dt1/dt2 就是弛豫后的结果）。
# 取值来源：Yin et al., Acta Mater. 2017 表 3 与图 2 的对称性分析（表 3 列的是
# Mg/Ti/Zr 等的数值，对 Au 只作初猜）。对称性精确定死的位置（Basal I2、Prism SF1/SF2）
# 直接写死；Pyr 上随 c/a 漂移的位置取 Mg 值起步。
#
# FCC_111 / HCP_basal 的 (t1, t2) 与 gsfe.mod 的 (bp1, bp2) 保持同号，
# 这样 ssf 的结果可以直接和 γ-line 端点对拍 —— 对称性禁止面内弛豫的体系（如 FCC ISF）
# 两者必须重合，这是最好的回归基准。
# t1/t2 是 Cartesian Δxz/lx 与 Δyz/ly，不是斜 a2 基矢上的分数。
# Pyr. I 窄面 e2 与 Yin 图 2 同向，表 3 的 SF1 第二分量为负（原初猜 -0.11）。
# 2026-10-11 PJ：DFT（A51-2，从 t2=±0.11 两个起点弛豫）都落到 t2≈-0.17，初猜改为 -0.17；
# 此前（stage2、stage3 round 0）的 Pyr. I-SF1 物性结果是 -0.11 起点，跨轮比较时注明。
# t1/t2 写成 LAMMPS equal-style 可求值的表达式字符串，交给 LAMMPS 求值，避免这里先转浮点掉精度。
LIST_SSF_FAULT = [
    ('fcc', 'FCC_111',     'ISF', '-0.5',     '1.0/3.0'),
    ('fcc', 'FCC_100',     'SF',  '0.5',      '0.5'),
    ('hcp', 'HCP_basal',   'I2',  '-0.5',     '1.0/3.0'),
    ('hcp', 'HCP_prism1w', 'SF1', '0.5',      '0.0'),
    ('hcp', 'HCP_prism1w', 'SF2', '0.5',      '0.5'),
    ('hcp', 'HCP_pyr1n',   'SF1', '0.5',      '-0.17'),
    ('hcp', 'HCP_pyr1w',   'SF2', '0.0',      '0.43'),
    ('hcp', 'HCP_pyr1w',   'SF3', '0.5',      '0.49'),
    ('hcp', 'HCP_pyr2',    'SF1', '0.0',      '0.46'),
]

# 表面能：(相, 滑移面模型)。真空插在超胞的 z 边界上，切开的就是该模型 z 偏移所选的那个面，
# 所以 Prism / Pyr 一律用**宽面**（wide）变体 —— 层间距大的那个面才是物理上的解理面。
LIST_SURF_TYPE = [
    ('fcc', 'FCC_111'),
    ('fcc', 'FCC_100'),
    ('hcp', 'HCP_basal'),
    ('hcp', 'HCP_prism1w'),
    ('hcp', 'HCP_pyr1w'),
    ('hcp', 'HCP_pyr2'),
]

# 堆垛序列层错（单剪切生成不了，必须改层序）：(相, 层错名, n_repeat)。
# 序列本身定义在 mymetal.build.bulk.stacking.DICT_SEQ_FAULT。
LIST_SEQSF_FAULT = [
    ('hcp', 'HCP_basal_I1', 10),
    ('hcp', 'HCP_basal_I3', 10),
    ('fcc', 'FCC_111_ESF', 7),
    ('fcc', 'FCC_111_TB', 7),
]

# 真空厚度 (Å)：必须 > 2 * rcut，否则两个自由表面隔着周期边界仍能相互作用。
SURF_VACUUM = 20.0

DICT_PLANAR_TABLE = {'ssf': LIST_SSF_FAULT,
                     'surface': LIST_SURF_TYPE,
                     'seqsf': LIST_SEQSF_FAULT}


def fail(msg):
    """打印错误并终止（与包内其它模块的失败口径一致）。"""
    print('❌ ERROR: %s' % msg)
    raise SystemExit(1)


def format_planar_table(group: str = None) -> str:
    """把某一类配置渲染成每行一条记录、空格分隔的纯文本。

    供 bash runner ``while read ... ; do`` 消费，这样 bash 侧不用再抄一份表。

    Args:
        group: ``'ssf'`` / ``'surface'`` / ``'seqsf'``。

    Returns:
        str: 每行一条记录。列依次为
        ssf: ``lat gsfe_type fault t1 t2``；
        surface: ``lat gsfe_type``；
        seqsf: ``lat fault_type n_repeat``。

    Raises:
        SystemExit: ``group`` 不在表内。
    """
    if group not in DICT_PLANAR_TABLE:
        fail('unknown planar group %r, expect one of %s'
             % (group, sorted(DICT_PLANAR_TABLE)))
    return '\n'.join(' '.join(str(field) for field in lrow)
                     for lrow in DICT_PLANAR_TABLE[group])


def get_ssf_name(gsfe_type: str = None, fault: str = None) -> str:
    """SSF 的工作目录名 / 列名后缀，如 ``('HCP_prism1w', 'SF2') -> 'HCP_prism1w_SF2'``。"""
    return '%s_%s' % (gsfe_type, fault)


def get_planar_columns(group: str = None) -> list:
    """某一类在 epoch 扫描表里的列名（不含 ``epoch``），顺序与 :data:`DICT_PLANAR_TABLE` 一致。

    Args:
        group: ``'ssf'`` / ``'surface'`` / ``'seqsf'``。

    Returns:
        list: 列名列表。ssf/seqsf 每项一列 ``ssf_<name>``；
        surface 每个面两列 ``surf_unrel_<type>`` 与 ``surf_rel_<type>``。

    Raises:
        SystemExit: ``group`` 不在表内。
    """
    if group == 'ssf':
        return ['ssf_' + get_ssf_name(gsfe_type, fault)
                for lat, gsfe_type, fault, t1, t2 in LIST_SSF_FAULT]
    if group == 'seqsf':
        return ['ssf_' + fault_type for lat, fault_type, n_repeat in LIST_SEQSF_FAULT]
    if group == 'surface':
        lcol = []
        for lat, gsfe_type in LIST_SURF_TYPE:
            lcol += ['surf_unrel_' + gsfe_type, 'surf_rel_' + gsfe_type]
        return lcol
    fail('unknown planar group %r, expect one of %s'
         % (group, sorted(DICT_PLANAR_TABLE)))
