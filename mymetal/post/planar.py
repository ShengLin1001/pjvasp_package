"""
planar post-processing submodule

平面缺陷的**单点**后处理：稳定层错能（SSF / 堆垛序列层错）与表面能。

与 :mod:`mymetal.post.gsfe` 的分工：
    - ``gsfe`` 处理整条 γ-line（21 步扫描，原子只放开面外、胞只放开 zz）；
    - 本模块处理单点构型 —— 层错矢量已知，直接倾斜 / 直接按层序建模，然后
      **原子 x,y,z 全放开、胞弛豫到 σ3j = 0**（Yin 2017 的 tilted-cell 全弛豫法），
      以及在 slab 里插真空得到的未弛豫 / 弛豫表面能。

LAMMPS 侧的约定（lmp_utils/template 的 ``ssf.mod`` / ``surface.mod`` / ``seqsf.mod``）：
    - ``dump/y_post_data.txt`` 按 jobn 递增记录每个状态的能量与应力，第 0 行必为参考态；
    - ``dump/movie.lammpstrj`` 按同样顺序追加每个状态**弛豫后**的构型；
    - ``dump/ideal.lammpstrj`` 单帧，记录层错**弛豫前**的理想构型，用于算 Yin 式 (1) 的
      原子位移 ``u = X - x``。

Functions:
    - get_layer_gap_change: 弛豫前后层间距变化中位数，识别整胞相变。
    - get_nn_spectrum: 逐原子前 n 个近邻距离（升序），即局域环境的指纹。
    - get_bulk_deviation: 层错胞里最像体相的原子偏离参考胞的程度，识别面内 shuffle 型整胞相变与薄胞。
    - post_ssf: 稳定层错能后处理（tilted-cell 与堆垛序列共用）。
    - post_surface: 表面能后处理（未弛豫 + 弛豫）。
    - read_ssf: 读回 post_ssf 写出的 txt。
    - read_surface: 读回 post_surface 写出的 txt。
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import Any, Dict, Union

from ase.io import read
from ase.neighborlist import neighbor_list
from myvasp import vasp_func as vf


def fail(msg):
    """打印错误并终止（与包内其它 post 模块的失败口径一致）。"""
    print('❌ ERROR: %s' % msg)
    raise SystemExit(1)


def warn(msg):
    """打印警告但继续。"""
    print('⚠️  %s' % msg)


### check —— 读入与一致性校验 ###

def read_planar_frames(post_data_path: Union[str, Path] = None,
                       latoms_path: Union[str, Path] = None,
                       latoms_format: str = 'lammps-dump-text',
                       nframe_expect: int = None) -> tuple:
    """读 y_post_data.txt + movie.lammpstrj，并核对两者帧数一致。

    Args:
        post_data_path: LAMMPS 侧 ``dump/y_post_data.txt``。
        latoms_path: LAMMPS 侧 ``dump/movie.lammpstrj``。
        latoms_format: ASE 格式名。
        nframe_expect: 期望帧数；None 时不校验。

    Returns:
        tuple: ``(Etot, latoms)``，``Etot`` 为 np.ndarray（eV），``latoms`` 为 Atoms 列表。

    Raises:
        SystemExit: 帧数不足 / 与期望不符 / 两个文件帧数对不上 / 参考态 jobn 不是 0。
    """
    jobn, Etot, Eent, pres = vf.vasp_read_post_data(str(post_data_path))
    latoms = read(str(latoms_path), format=latoms_format, index=':')

    if len(jobn) < 2:
        fail('need >= 2 states in %s, got %d' % (post_data_path, len(jobn)))
    if abs(float(jobn[0])) > 1e-5:
        fail('first state in %s is not the reference (jobn=%s)' % (post_data_path, jobn[0]))
    if len(latoms) != len(jobn):
        fail('frame mismatch: %d rows in %s vs %d frames in %s'
             % (len(jobn), post_data_path, len(latoms), latoms_path))
    if nframe_expect is not None and len(jobn) != nframe_expect:
        fail('expect %d states, got %d in %s' % (nframe_expect, len(jobn), post_data_path))

    return np.array(Etot, dtype=float), latoms


def get_inplane_area(atoms=None) -> float:
    """面内面积 ``|a1 x a2|``（Å²），与 y_gsfe 的 ``Asf`` 同一口径。"""
    return float(np.linalg.norm(np.cross(atoms.cell[0, :], atoms.cell[1, :])))


def check_inplane_fixed(latoms: list = None, tol: float = 1e-8) -> None:
    """核对所有帧的面内基矢 (a1, a2) 未被弛豫。

    Yin 2017 §3 强调：面内 ``N1a1, N2a2`` 必须锁死在体相值，否则层错的本征应力会被
    胞内其余部分的反向应力抵消，弹性能随 N3 漂移，层错能失去意义。这是 SSF 计算
    最容易被 ``fix box/relax`` 写错的地方，所以在 post 侧硬校验。

    Args:
        latoms: 各状态的 Atoms 列表。
        tol: 基矢差的范数容差（Å）。

    Raises:
        SystemExit: 任一帧的面内基矢相对第 0 帧有变化。
    """
    for i, atoms in enumerate(latoms):
        dlatt = atoms.cell[:] - latoms[0].cell[:]
        temp = float(np.linalg.norm(dlatt[0:2, :]))
        if temp > tol:
            print('dlatt:', dlatt)
            fail('in-plane lattice relaxed at frame %d (|d(a1,a2)| = %.3e Ang)' % (i, temp))


def get_displacement(atoms_ideal=None, atoms_relaxed=None) -> np.ndarray:
    """Yin 2017 式 (1) 的原子位移 ``u = X - x``（弛豫后 − 理想层错构型）。

    先在分数坐标下取最小镜像，再回到笛卡尔，最后逐分量扣掉平均值 —— 整体平移是
    弛豫的零模（原子全放开时必然出现），不扣掉会把 u 整体抬起来。

    Args:
        atoms_ideal: 弛豫前的理想层错构型。
        atoms_relaxed: 弛豫后的构型。

    Returns:
        np.ndarray: 形状 (natoms, 3) 的位移（Å）。

    Raises:
        SystemExit: 两个构型原子数不一致。
    """
    if len(atoms_ideal) != len(atoms_relaxed):
        fail('ideal/relaxed atom count mismatch: %d vs %d'
             % (len(atoms_ideal), len(atoms_relaxed)))

    cell = np.array(atoms_relaxed.cell[:], dtype=float)
    dpos = atoms_relaxed.positions - atoms_ideal.positions
    dposD = dpos @ np.linalg.inv(cell)
    dposD = dposD - np.around(dposD)          # 最小镜像：跨周期边界的原子不算大位移
    dpos = dposD @ cell
    return dpos - dpos.mean(axis=0)           # 扣掉刚体平移零模


def get_layer_gap_change(atoms_ideal=None, atoms_relaxed=None) -> float:
    """弛豫前后逐个层间距变化量的中位数（Å），用来识别整胞相变。

    层序按弛豫前构型排，原子按编号配对，而不是两边各自排序：末态 z 折回会改变窄/宽层间距
    的排列相位，各自排序会把相位差误当成层间距变化。a1/a2 锁在 xy 面内，所以层间距直接用
    a3 方向的分数坐标乘 c_z（tilted cell 的 a3 带面内分量也成立）。

    正常层错只在层错面附近弛豫，绝大多数层间距不动，中位数近 0；整胞相变（如 HCP prism1w SF2
    的 20 个波纹层全部压平成 10 个平层）每层都变，中位数就是变化幅度本身。

    Args:
        atoms_ideal: 弛豫前的理想层错构型。
        atoms_relaxed: 弛豫后的构型（原子编号与 atoms_ideal 一一对应）。

    Returns:
        float: 层间距变化量绝对值的中位数（Å）。

    Raises:
        SystemExit: 两个构型原子数不一致。
    """
    if len(atoms_ideal) != len(atoms_relaxed):
        fail('ideal/relaxed atom count mismatch: %d vs %d'
             % (len(atoms_ideal), len(atoms_relaxed)))

    s0 = atoms_ideal.get_scaled_positions(wrap=True)[:, 2]
    s1 = atoms_relaxed.get_scaled_positions(wrap=False)[:, 2]
    order = np.argsort(s0)
    gap0 = np.diff(np.append(s0[order], s0[order[0]] + 1.0))
    gap1 = np.diff(np.append(s1[order], s1[order[0]] + 1.0))
    gap1 -= np.round(gap1 - gap0)              # 末态个别原子跨周期边界折回，差整数个周期
    return float(np.median(np.abs(gap1 * atoms_relaxed.cell[2, 2] - gap0 * atoms_ideal.cell[2, 2])))


def get_nn_spectrum(atoms=None, nneigh: int = 12, cutoff: float = None) -> np.ndarray:
    """逐原子前 ``nneigh`` 个近邻的距离（升序），即局域环境的指纹。

    截断内不足 ``nneigh`` 个近邻时用 ``cutoff`` 补齐 —— 近邻被推出截断本身就是“不像体相”，
    补成截断值让它在 :func:`get_bulk_deviation` 里显出大偏差，而不是被悄悄忽略。
    不区分元素（Au 单元素；多元素体系须另按元素配对）。

    Args:
        atoms: 结构。
        nneigh: 取前几个近邻；密排金属第一壳层为 12。
        cutoff: 近邻搜索截断（Å）。

    Returns:
        np.ndarray: 形状 (natoms, nneigh)。
    """
    li, ld = neighbor_list('id', atoms, cutoff)
    spectrum = np.full((len(atoms), nneigh), float(cutoff))
    for i in range(len(atoms)):
        dist = np.sort(ld[li == i])[:nneigh]
        spectrum[i, :len(dist)] = dist
    return spectrum


def get_bulk_deviation(atoms_ref=None, atoms_relaxed=None, nneigh: int = 12) -> float:
    """层错胞里**最像体相**的那个原子，其近邻谱偏离参考胞的程度（Å）。

    逐原子取近邻谱对参考胞各原子近邻谱的最大分量差，再对参考原子取最小（参考胞不必每个原子
    都等价），最后对层错胞所有原子取最小。孤立层错离层错面远的原子回到体相环境，值近 0；
    整个胞没有一个原子像体相时值很大 —— 两种成因：

    - 面内 shuffle 型整胞相变：层间距不变，:func:`get_layer_gap_change` 抓不到。
      DFT 的 HCP pyr2 SF1（A21-3，20 原子）每个原子都变成 10 个近邻 2.79–2.86 Å + 2 个 3.43 Å，
      本值 0.516 Å，而层间距中位数只变 0.002 Å（2026-10-13 审阅）；
    - 胞太薄，层错的弛豫场贯穿整个胞，不存在孤立层错（LAMMPS pyr2 [1,1,5] 为 0.09–0.16 Å）。

    Args:
        atoms_ref: 同尺寸的无层错参考胞（已弛豫）。
        atoms_relaxed: 弛豫后的层错胞（原子数可与参考胞不同，如 FCC ESF）。
        nneigh: 近邻数。

    Returns:
        float: 偏差（Å）。
    """
    # 截断取参考胞第 nneigh 近邻的 1.5 倍：盖住第一壳层，又够大，被推远的近邻落在截断内
    li, ld = neighbor_list('id', atoms_ref, 6.0)
    d_last = max(np.sort(ld[li == i])[nneigh - 1] for i in range(len(atoms_ref)))
    cutoff = 1.5 * float(d_last)
    spectrum_ref = get_nn_spectrum(atoms_ref, nneigh, cutoff)
    spectrum = get_nn_spectrum(atoms_relaxed, nneigh, cutoff)
    ldev = np.abs(spectrum[:, None, :] - spectrum_ref[None, :, :]).max(axis=2).min(axis=1)
    return float(ldev.min())


### main —— 两个对外后处理入口 ###

def post_ssf(post_data_path: Union[str, Path] = None,
             latoms_path: Union[str, Path] = None,
             ideal_path: Union[str, Path] = None,
             latoms_format: str = 'lammps-dump-text',
             nfault: int = 1,
             label: str = 'SSF',
             dt_tol: float = 0.25,
             gamma_tol: float = 1.0,
             gap_tol: float = 0.2,
             bulk_tol: float = 0.05,
             save_txt_path: Union[str, Path] = './y_post_ssf.txt') -> Dict[str, Any]:
    """稳定层错能后处理（tilted-cell 单点 与 堆垛序列层错 共用）。

    两个状态：第 0 帧是同尺寸的**无层错**参考胞（已弛豫），第 1 帧是层错胞（已弛豫）。

    ``gamma = (E_fault - natoms_fault * E0_bulk) / (nfault * Asf)``

    用**每原子**参考能而非整胞能量差，是为了兼容层错胞与参考胞原子数不同的情形
    （FCC 外禀层错 ``ABC*n + B`` 比 ``ABC*n`` 多一层）。两者原子数相同时，该式
    与直接作差完全等价。

    Args:
        post_data_path: ``dump/y_post_data.txt``。
        latoms_path: ``dump/movie.lammpstrj``（2 帧：参考、层错）。
        ideal_path: ``dump/ideal.lammpstrj``（1 帧：层错弛豫前）。None 时跳过位移统计。
        latoms_format: ASE 格式名。
        nfault: 每胞的等价层错数。tilted-cell 为 1；周期堆垛序列胞为 2（Yin 2017 §3）。
        label: 写进 txt 表头的方法名，如 ``'SSF (tilted-cell)'``。
        dt_tol: 层错矢量相对初猜的面内漂移容差（以 a1/a2 为单位）。超过就判定构型
            已经滑出目标层错的势阱（通常是滑回完整晶格），见 ``stable`` 字段。
        gamma_tol: 判定「其实没有层错」的能量阈值（mJ/m²）。
        gap_tol: 整胞相变判据（Å）：弛豫前后层间距变化量的中位数超过它，就判定整个胞变成了
            另一种结构，gamma 随原子数增长、不是层错能。0.2 取自 DFT（A51-2，2026-10-11）：
            正常层错 ≤ 0.063，prism1w SF2 为 0.80；stage2 六势 LAMMPS 正常层错 ≤ 0.03。
            没有 ideal 构型时不判。
        bulk_tol: 体相保留判据（Å）：层错胞里最像体相的原子，其 12 近邻谱偏离参考胞超过它，
            就判定胞里已没有体相原子（面内 shuffle 型整胞相变，或胞太薄容不下孤立层错），
            归入 ``transformed``，见 :func:`get_bulk_deviation`。0.05 取自 2026-10-13 标定：
            正常层错 DFT（A21-3/A51/A52）≤ 0.015、M1 六势 LAMMPS ≤ 0.013；DFT pyr2 SF1 0.516、
            prism1w SF2 ≥ 0.068，LAMMPS pyr2 [1,1,5] 0.09–0.16。
        save_txt_path: 输出 txt 路径。

    Returns:
        Dict[str, Any]: ``gamma``、``Asf``、``E0bulk``、``t1``/``t2``、``dt1``/``dt2``、
        ``da33``、``umax``、``uabs``（逐原子位移模长）、``gap_change``（Å，无 ideal 时为 nan）、
        ``bulk_dev``（Å）、``transformed``（整胞相变或无体相原子，bool）、
        ``stable``（bool，整胞相变也算不稳定）。

    Raises:
        SystemExit: 帧数不对、面内基矢被弛豫、或原子数/元素不一致。
    """
    qe = vf.phy_const('qe')

    Etot, latoms = read_planar_frames(post_data_path, latoms_path,
                                      latoms_format, nframe_expect=2)
    check_inplane_fixed(latoms)

    atoms_ref, atoms_sf = latoms[0], latoms[1]
    Asf = get_inplane_area(atoms_ref)
    a11 = float(atoms_ref.cell[0, 0])
    a22 = float(atoms_ref.cell[1, 1])
    E0bulk = float(Etot[0]) / len(atoms_ref)

    dE = float(Etot[1]) - len(atoms_sf) * E0bulk
    gamma = dE / (nfault * Asf) * qe * 1e23

    # 层错矢量：a3 的面内分量相对参考胞的变化，归一到 a11 / a22（与 y_gsfe 的 da31/a11 同口径）
    da3 = np.array(atoms_sf.cell[2, :], dtype=float) - np.array(atoms_ref.cell[2, :], dtype=float)
    t1, t2, da33 = da3[0] / a11, da3[1] / a22, da3[2]

    # 相对理想倾斜位置的面内弛豫量 —— 这正是「标准法」与「全弛豫法」的差别所在
    dt1 = dt2 = 0.0
    uabs = np.zeros(len(atoms_sf))
    gap_change = float('nan')
    if ideal_path is not None and Path(ideal_path).is_file():
        atoms_ideal = read(str(ideal_path), format=latoms_format, index='-1')
        da3_ideal = (np.array(atoms_ideal.cell[2, :], dtype=float)
                     - np.array(atoms_ref.cell[2, :], dtype=float))
        dt1 = t1 - da3_ideal[0] / a11
        dt2 = t2 - da3_ideal[1] / a22
        uabs = np.linalg.norm(get_displacement(atoms_ideal, atoms_sf), axis=1)
        gap_change = get_layer_gap_change(atoms_ideal, atoms_sf)
    else:
        warn('no ideal config at %s, in-plane relaxation dt1/dt2 and u are reported as 0, '
             'whole-cell transformation not checked' % ideal_path)

    bulk_dev = get_bulk_deviation(atoms_ref, atoms_sf)

    # 这个层错到底存不存在。四种「不存在」成因不同（前两种给出 gamma ~ 0），必须报出来
    # 而不是当成一个层错能：
    #   - 面内滑走（|dt| > dt_tol）：初猜位置不是势阱，弛豫把它推回完整晶格（HCP 窄面 Pyr. I 常见）；
    #   - 停在原地但 gamma ~ 0：该位置本身就等价于一个完整的点阵平移（如 FCC {100} 的 (1/2,1/2)）。
    #   - 整胞相变（gap_change > gap_tol）：dt 可以很小、gamma 也不小，但每一层都变了，
    #     gamma 是两种结构的能量差、随原子数增长（DFT 与 LAMMPS 的 prism1w SF2 都如此）。
    #   - 没有体相原子（bulk_dev > bulk_tol）：层间距没变、面内 shuffle 了整个胞（DFT pyr2 SF1），
    #     或胞太薄、弛豫场贯穿全胞，同样不是孤立层错的能量。
    if_gap = bool(gap_change > gap_tol)               # nan（无 ideal）比较为 False，即不判
    if_bulk_lost = bool(bulk_dev > bulk_tol)
    transformed = if_gap or if_bulk_lost
    stable = (abs(dt1) <= dt_tol) and (abs(dt2) <= dt_tol) and (abs(gamma) > gamma_tol) \
        and not transformed
    if if_gap:
        warn('%s: whole cell transformed (median layer-gap change %.4f Ang > %.2f) — '
             'gamma = %.4f mJ/m^2 is an energy difference between two structures, '
             'not a fault energy' % (label, gap_change, gap_tol, gamma))
    elif if_bulk_lost:
        warn('%s: no bulk-like atom left (NN-spectrum deviation of the most bulk-like atom '
             '%.4f Ang > %.2f) — in-plane whole-cell transformation or cell too thin; '
             'gamma = %.4f mJ/m^2 is not an isolated fault energy'
             % (label, bulk_dev, bulk_tol, gamma))
    elif not stable:
        warn('%s: no stable stacking fault here '
             '(gamma = %.4f mJ/m^2, dt = %.4f, %.4f) — the relaxed cell is the perfect '
             'lattice, do not read gamma as a fault energy' % (label, gamma, dt1, dt2))

    dict_out = {'gamma': gamma, 'Asf': Asf, 'a11': a11, 'a22': a22,
                'E0bulk': E0bulk, 'dE': dE, 'nfault': nfault,
                't1': t1, 't2': t2, 'dt1': dt1, 'dt2': dt2,
                'da33': da33, 'umax': float(uabs.max()), 'uabs': uabs,
                'gap_change': gap_change, 'bulk_dev': bulk_dev,
                'transformed': transformed, 'stable': stable}
    write_ssf_output(dict_out, label=label, save_txt_path=save_txt_path)
    return dict_out


def post_surface(post_data_path: Union[str, Path] = None,
                 latoms_path: Union[str, Path] = None,
                 latoms_format: str = 'lammps-dump-text',
                 label: str = 'surface',
                 save_txt_path: Union[str, Path] = './y_post_surface.txt') -> Dict[str, Any]:
    """表面能后处理：未弛豫 + 弛豫。

    三个状态：第 0 帧是全周期体相超胞（已弛豫），第 1 帧是插真空后的**单点**（未弛豫表面），
    第 2 帧是同一真空胞弛豫后的（弛豫表面）。真空把周期方向切开，一个胞暴露 **2 个**表面：

    ``gamma = (E_slab - E_bulk) / (2A)``

    未弛豫表面能与 chap_4_TL 图 1.9 的 decohesion 曲线收敛平台值等价，代价小得多。

    Args:
        post_data_path: ``dump/y_post_data.txt``。
        latoms_path: ``dump/movie.lammpstrj``（3 帧）。
        latoms_format: ASE 格式名。
        label: 写进 txt 表头的面名，如 ``'FCC_111'``。
        save_txt_path: 输出 txt 路径。

    Returns:
        Dict[str, Any]: ``gamma_unrelaxed``、``gamma_relaxed``、``Asf``、``E0bulk``、
        ``u3max``（弛豫引起的最大面外位移）、``u3``（逐原子面外位移）。

    Raises:
        SystemExit: 帧数不对，或面内基矢被弛豫。
    """
    qe = vf.phy_const('qe')

    Etot, latoms = read_planar_frames(post_data_path, latoms_path,
                                      latoms_format, nframe_expect=3)
    check_inplane_fixed(latoms)

    atoms_bulk, atoms_unrel, atoms_rel = latoms[0], latoms[1], latoms[2]
    Asf = get_inplane_area(atoms_bulk)
    natoms = len(atoms_bulk)
    E0bulk = float(Etot[0]) / natoms

    # 真空帧与体相帧原子数相同（只改胞长、不加减原子），故直接作差即可
    gamma_unrel = (float(Etot[1]) - float(Etot[0])) / (2.0 * Asf) * qe * 1e23
    gamma_rel = (float(Etot[2]) - float(Etot[0])) / (2.0 * Asf) * qe * 1e23
    if gamma_rel > gamma_unrel + 1e-6:
        warn('relaxed surface energy (%.4f) exceeds unrelaxed (%.4f) — minimization suspect'
             % (gamma_rel, gamma_unrel))

    u = get_displacement(atoms_unrel, atoms_rel)
    u3 = u[:, 2]
    z = atoms_rel.positions[:, 2]

    dict_out = {'gamma_unrelaxed': gamma_unrel, 'gamma_relaxed': gamma_rel,
                'Asf': Asf, 'E0bulk': E0bulk, 'natoms': natoms,
                'u3max': float(np.abs(u3).max()), 'u3': u3, 'z': z}
    write_surface_output(dict_out, label=label, save_txt_path=save_txt_path)
    return dict_out


### output —— 写 txt 与读回 ###

def write_ssf_output(dict_out: dict = None, label: str = 'SSF',
                     save_txt_path: Union[str, Path] = './y_post_ssf.txt') -> None:
    """把 :func:`post_ssf` 的结果写成定宽文本表。

    Args:
        dict_out: :func:`post_ssf` 的返回字典。
        label: 方法名，写进表头。
        save_txt_path: 输出路径。
    """
    with open(str(save_txt_path), 'w+', encoding='utf-8') as f:
        f.write('# LAMMPS stable stacking fault: %s \n' % label)
        f.write('# gamma = (E_fault - natoms*E0_bulk) / (nfault * Asf) \n')
        f.write('# all atoms free in x,y,z; cell relaxed to sigma_3j = 0 (zz, xz, yz) \n')

        f.write('\n%16s %16s %16s %16s %8s \n'
                % ('Asf (Ang^2)', 'a11 (Ang)', 'a22 (Ang)', 'E0_bulk (eV)', 'nfault'))
        f.write('%16.8f %16.8f %16.8f %16.8f %8d \n\n'
                % (dict_out['Asf'], dict_out['a11'], dict_out['a22'],
                   dict_out['E0bulk'], dict_out['nfault']))

        f.write('%20s: %16.8f \n' % ('gamma_ssf (mJ/m^2)', dict_out['gamma']))
        f.write('%20s: %16.8f \n' % ('dE (eV)', dict_out['dE']))
        f.write('%20s: %16s \n' % ('stable', 'yes' if dict_out['stable'] else 'no'))
        f.write('%20s: %16s \n' % ('transformed', 'yes' if dict_out['transformed'] else 'no'))
        f.write('%20s: %16.8f \n' % ('gap_change (Ang)', dict_out['gap_change']))
        f.write('%20s: %16.8f \n' % ('bulk_dev (Ang)', dict_out['bulk_dev']))

        # t1/t2 是弛豫后的真实层错矢量；dt1/dt2 是它相对理想倾斜位置的面内弛豫，
        # 非零就说明「标准法」（面内锁死）会低估这个层错的稳定性。
        f.write('\n%10s %10s %10s %10s %12s %12s \n'
                % ('t1 (-)', 't2 (-)', 'dt1 (-)', 'dt2 (-)',
                   'da33 (Ang)', 'umax (Ang)'))
        f.write('%10.4f %10.4f %10.4f %10.4f %12.6f %12.6f \n'
                % (dict_out['t1'], dict_out['t2'], dict_out['dt1'], dict_out['dt2'],
                   dict_out['da33'], dict_out['umax']))
        f.write(' \n')
    print("✅ ssf written to '%s'  gamma = %.4f mJ/m^2\n" % (save_txt_path, dict_out['gamma']))


def write_surface_output(dict_out: dict = None, label: str = 'surface',
                         save_txt_path: Union[str, Path] = './y_post_surface.txt') -> None:
    """把 :func:`post_surface` 的结果写成定宽文本表。

    Args:
        dict_out: :func:`post_surface` 的返回字典。
        label: 面名，写进表头。
        save_txt_path: 输出路径。
    """
    with open(str(save_txt_path), 'w+', encoding='utf-8') as f:
        f.write('# LAMMPS surface energy: %s \n' % label)
        f.write('# gamma = (E_slab - E_bulk) / (2*Asf),  vacuum cleaves one cell into 2 surfaces \n')

        f.write('\n%16s %16s %16s \n' % ('Asf (Ang^2)', 'E0_bulk (eV)', 'natoms'))
        f.write('%16.8f %16.8f %16d \n\n'
                % (dict_out['Asf'], dict_out['E0bulk'], dict_out['natoms']))

        f.write('%24s: %16.8f \n' % ('gamma_unrelaxed (mJ/m^2)', dict_out['gamma_unrelaxed']))
        f.write('%24s: %16.8f \n' % ('gamma_relaxed (mJ/m^2)', dict_out['gamma_relaxed']))
        f.write('%24s: %16.8f \n' % ('u3max (Ang)', dict_out['u3max']))

        f.write('\n%16s %16s \n' % ('z (Ang)', 'u3 (Ang)'))
        for zi, ui in zip(dict_out['z'], dict_out['u3']):
            f.write('%16.6f %16.6f \n' % (zi, ui))
        f.write(' \n')
    print("✅ surface written to '%s'  gamma_unrel = %.4f  gamma_rel = %.4f mJ/m^2\n"
          % (save_txt_path, dict_out['gamma_unrelaxed'], dict_out['gamma_relaxed']))


def read_scalar_line(lline: list = None, key: str = None) -> float:
    """从 ``'  name: value'`` 形式的行里取 value；找不到返回 ``np.nan``。

    Args:
        lline: 文件的全部行。
        key: 冒号左侧的标识串（做子串匹配）。

    Returns:
        float: 解析到的数值，或 ``np.nan``。
    """
    for line in lline:
        if key in line and ':' in line:
            try:
                return float(line.split(':')[-1].strip())
            except ValueError:
                continue
    return float('nan')


def read_ssf(save_txt_path: Union[str, Path] = './y_post_ssf.txt') -> Dict[str, Any]:
    """读回 :func:`write_ssf_output` 写的 txt。

    Args:
        save_txt_path: txt 路径。

    Returns:
        Dict[str, Any]: ``{'gamma': ..., 'dE': ..., 'stable': ..., 'transformed': ...}``；
        文件不存在或字段缺失时数值为 ``nan``、``stable`` 为 False。加整胞相变判据前写的
        旧 txt 没有 ``transformed`` 行，读作 False（即当时未判）。

        不返回 ``gap_change``：调用方（mymetal.post.property）要求返回值全部有限，旧 txt
        读不到它会是 nan，整族就被判失败。
    """
    path_txt = Path(save_txt_path)
    if not path_txt.is_file():
        return {'gamma': float('nan'), 'dE': float('nan'), 'stable': False, 'transformed': False}
    lline = path_txt.read_text(encoding='utf-8').splitlines()
    stable = any(line.strip().startswith('stable') and line.strip().endswith('yes')
                 for line in lline)
    transformed = any(line.strip().startswith('transformed') and line.strip().endswith('yes')
                      for line in lline)
    return {'gamma': read_scalar_line(lline, 'gamma_ssf'),
            'dE': read_scalar_line(lline, 'dE (eV)'),
            'stable': stable, 'transformed': transformed}


def read_surface(save_txt_path: Union[str, Path] = './y_post_surface.txt') -> Dict[str, Any]:
    """读回 :func:`write_surface_output` 写的 txt。

    Args:
        save_txt_path: txt 路径。

    Returns:
        Dict[str, Any]: ``{'gamma_unrelaxed': ..., 'gamma_relaxed': ...}``；
        文件不存在或字段缺失时值为 ``nan``。
    """
    path_txt = Path(save_txt_path)
    if not path_txt.is_file():
        return {'gamma_unrelaxed': float('nan'), 'gamma_relaxed': float('nan')}
    lline = path_txt.read_text(encoding='utf-8').splitlines()
    return {'gamma_unrelaxed': read_scalar_line(lline, 'gamma_unrelaxed'),
            'gamma_relaxed': read_scalar_line(lline, 'gamma_relaxed')}
