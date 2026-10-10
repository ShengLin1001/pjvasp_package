"""
stacking submodule

按显式 A/B/C 堆垛序列构建密排面超胞，用于**单剪切无法生成**的层错：
HCP basal 的 I1 与 I3(extrinsic E)、FCC {111} 的 ESF 与孪晶界。

单剪切能生成的层错（Basal I2 / FCC ISF / Prism / Pyr 各 SF）走 tilted-cell 法，
见 :mod:`mymetal.build.bulk.gsfe` 与 lmp_utils 的 ``ssf.mod``；本模块只负责
"必须改层序" 的那一类。

参考：Yin et al., Acta Mater. 2017（references/Yin-2017-Comprehensive first-principles study.pdf）
§3：I1 用 ``(AB)10-(AC)10``、extrinsic E 用 ``(AB)10-C-(AB)9-C``，周期边界下**每胞含 2 个
等价层错**，故层错能除以 ``2A``。

Functions:
    - create_stacking_model: 按 A/B/C 序列串构建密排超胞。
    - get_stacking_sequence: 按层错名查 (perfect, faulted) 序列对。
    - check_stacking_sequence: 校验序列是否为合法密排堆垛（含周期回绕）。
"""

from ase import Atoms
import numpy as np

# 密排面内三个套构位置在 (a1, a2) 基下的分数坐标。
# 由 create_fcc_111() 产出的 data.ini 反解得到：a1=(ax,0,0)、a2=(-ax/2, ax*sqrt(3)/2, 0) 时，
# 相邻 (111) 层依次落在 A -> B -> C，与 FCC 的 ABCABC 顺序一致。改这里会同时改动 fcc/hcp 两支。
DICT_STACKING_OFFSET = {'A': (0.0, 0.0),
                        'B': (2.0 / 3.0, 1.0 / 3.0),
                        'C': (1.0 / 3.0, 2.0 / 3.0)}

# 与 create.py 同源：避免原子正好压在 xlo/ylo 面上，被 LAMMPS minimize 的浮点噪声挤出盒子
# （现象是 write/read_restart 报 "Did not assign all restart atoms correctly"）。
INPLANE_SHIFT = 1e-5

# 层错名 -> (相, 生成 perfect / faulted 序列的函数, 默认重复数 n)。
# 序列写法直接照搬 Yin 2017 §3；FCC 两项是同一构造思路在 ABC 堆垛上的对应物。
# 每个序列在周期回绕后都恰好含 2 个等价层错，post 侧统一按 2A 归一。
DICT_SEQ_FAULT = {
    # HCP：Yin 原文的 (AB)n-(AC)n 与 (AB)n-C-(AB)(n-1)-C，共 4n 层
    'HCP_basal_I1': ('hcp', lambda n: ('AB' * (2 * n),
                                       'AB' * n + 'AC' * n), 10),
    'HCP_basal_I3': ('hcp', lambda n: ('AB' * (2 * n),
                                       'AB' * n + 'C' + 'AB' * (n - 1) + 'C'), 10),
    # FCC：外禀层错 = 插入一层；孪晶 = 上半胞堆垛顺序反向
    'FCC_111_ESF': ('fcc', lambda n: ('ABC' * (2 * n),
                                      'ABC' * n + 'B' + 'ABC' * n + 'B'), 7),
    'FCC_111_TB':  ('fcc', lambda n: ('ABC' * (2 * n),
                                      'ABC' * n + 'BAC' * n), 7),
}


def check_stacking_sequence(sequence: str = None) -> str:
    """校验 A/B/C 序列是否为合法密排堆垛。

    合法的判据只有一条：相邻两层不能同字母（含最后一层与第一层的周期回绕）。
    同字母意味着两层原子面内完全重合、层间距只有 d，是非物理构型。

    Args:
        sequence: 由 'A'/'B'/'C' 组成的堆垛序列串。

    Returns:
        原样返回 ``sequence``，便于 ``seq = check_stacking_sequence(seq)`` 串联。

    Raises:
        SystemExit: 序列为空、含非 ABC 字符、层数少于 2，或存在相邻同层。
    """
    if not sequence:
        fail('empty stacking sequence')

    lbad = sorted(set(sequence) - set(DICT_STACKING_OFFSET))
    if lbad:
        fail('illegal layer symbol(s) %s in sequence %r' % (lbad, sequence))

    nlayer = len(sequence)
    if nlayer < 2:
        fail('stacking sequence needs >= 2 layers, got %d' % nlayer)

    for i in range(nlayer):
        # i = nlayer-1 时比的是周期回绕的那一对，PBC 下同样必须合法
        if sequence[i] == sequence[(i + 1) % nlayer]:
            fail('adjacent identical layers %r at index %d (wrap included) in %r'
                 % (sequence[i], i, sequence))
    return sequence


def get_stacking_sequence(fault_type: str = None, n_repeat: int = None) -> tuple:
    """按层错名取 (perfect, faulted) 序列对。

    Args:
        fault_type: :data:`DICT_SEQ_FAULT` 的键，如 ``'HCP_basal_I1'``。
        n_repeat: 重复单元数 n。None 时取该层错的默认值（HCP 为 10，对应 Yin 的 40 层）。

    Returns:
        tuple: ``(lat, seq_perfect, seq_fault)``，``lat`` 为 ``'fcc'`` / ``'hcp'``。

    Raises:
        SystemExit: 层错名不在表内，或 ``n_repeat`` 不是正整数。
    """
    if fault_type not in DICT_SEQ_FAULT:
        fail('unknown fault_type %r, expect one of %s'
             % (fault_type, sorted(DICT_SEQ_FAULT)))

    lat, func_seq, n_default = DICT_SEQ_FAULT[fault_type]
    n_repeat = n_default if n_repeat is None else int(n_repeat)
    if n_repeat < 2:
        fail('n_repeat must be >= 2 (extrinsic sequences need n-1 >= 1), got %d' % n_repeat)

    seq_perfect, seq_fault = func_seq(n_repeat)
    seq_perfect = check_stacking_sequence(seq_perfect)
    seq_fault = check_stacking_sequence(seq_fault)
    return lat, seq_perfect, seq_fault


def create_stacking_model(sequence: str = None, a: float = None, c: float = None,
                          lat: str = 'fcc', symbol: str = 'Au') -> Atoms:
    """按 A/B/C 序列构建密排面超胞（每层 1 个原子，面外为第 3 轴）。

    胞的约定与 :func:`mymetal.build.bulk.create.create_fcc_111` /
    :func:`~mymetal.build.bulk.create.create_hcp_basal` 保持一致：
    ``a1 = (ax, 0, 0)``、``a2 = (-ax/2, ax*sqrt(3)/2, 0)``、``a3 = (0, 0, nlayer*d)``，
    这样面积 ``A = |a1 x a2|`` 与 y_gsfe 的 ``Asf`` 口径相同，可直接互相比对。

    Args:
        sequence: A/B/C 堆垛序列串，自底向上。
        a: 晶格常数 a（Å）。
        c: HCP 晶格常数 c（Å）；``lat='fcc'`` 时忽略。
        lat: ``'fcc'``（面内 ``ax=a/sqrt(2)``、层间距 ``a/sqrt(3)``）或
            ``'hcp'``（面内 ``ax=a``、层间距 ``c/2``）。
        symbol: 元素符号。

    Returns:
        ase.Atoms: 原子数 = ``len(sequence)`` 的密排超胞，三向周期。

    Raises:
        SystemExit: ``lat`` 非法、``a`` 缺失、``lat='hcp'`` 但 ``c`` 缺失，或序列非法。
    """
    sequence = check_stacking_sequence(sequence)
    if a is None:
        fail('lattice constant a is required')

    if lat == 'fcc':
        ax = a / np.sqrt(2.0)          # (111) 面内最近邻距离
        d = a / np.sqrt(3.0)           # (111) 层间距
    elif lat == 'hcp':
        if c is None:
            fail("lattice constant c is required for lat='hcp'")
        ax = a
        d = c / 2.0                    # (0001) 层间距
    else:
        fail("lat must be 'fcc' or 'hcp', got %r" % lat)

    nlayer = len(sequence)
    cell = np.array([[ax, 0.0, 0.0],
                     [-ax / 2.0, ax * np.sqrt(3.0) / 2.0, 0.0],
                     [0.0, 0.0, nlayer * d]])

    # 面内按分数坐标定位再乘 cell，面外用绝对层号；整体抬 0.25*d 让最底层不贴 zlo 面。
    lfrac = [DICT_STACKING_OFFSET[letter] for letter in sequence]
    frac = np.zeros([nlayer, 3])
    frac[:, 0] = [f[0] for f in lfrac]
    frac[:, 1] = [f[1] for f in lfrac]
    frac[:, 2] = (np.arange(nlayer) + 0.25) / nlayer
    frac[:, 0] += INPLANE_SHIFT
    frac[:, 1] += INPLANE_SHIFT

    atoms = Atoms(symbols=[symbol] * nlayer,
                  scaled_positions=frac,
                  cell=cell,
                  pbc=(True, True, True))
    atoms.wrap()
    return atoms


def fail(msg):
    """打印错误并终止（与包内其它模块的失败口径一致）。"""
    print('❌ ERROR: %s' % msg)
    raise SystemExit(1)
