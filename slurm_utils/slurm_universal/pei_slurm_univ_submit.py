#!/usr/bin/env python3
"""
pei_slurm_univ_submit —— 通用 Slurm 作业脚本生成 / 提交 CLI。

这是 mymetal.slurm.submit.pei_slurm_univ_submit 的薄命令行封装：本文件只负责
解析 argparse 参数再原样转交给库函数，不再自带任何业务逻辑。

设计要点：
    每个命令行选项名都与库函数的形参严格同名（dest 一致），因此 main 里可以直接
    `pei_slurm_univ_submit(**vars(args))` 透传，无需手抄一份参数映射、也不会漏改。
    所有结构性 / 枚举 / None / 正整数校验都下推给库函数，这里只搬运。
"""

import argparse
import sys
import textwrap
from pathlib import Path

from mymetal.slurm.submit import pei_slurm_univ_submit

# 输出含大量 emoji（📁▶️✅❌📊🎉），在 locale 非 UTF-8 的计算节点上直接 print 会炸；
# 入口处兜底把 stdout/stderr 切到 utf-8，老 Python 没有 reconfigure 时静默跳过。
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass


# ================ 📋 常见用法预设（preset 注册表）
# 把以前散落在文件底部、彼此用同名 args 反复覆盖且从不被调用的 dict，收敛成键控注册表。
# --preset NAME 先用这里的值灌默认，命令行再显式覆盖个别字段（如 --ncores 64）。
# 只放各场景的「区别性」字段；path_root / lsubdir 等沿用 CLI 默认，不在 import 期固化 cwd。
PRESETS = {
    # zcm6-vasp-0：每个子目录一个 VASP 作业，启动器经 MY_LAUNCHER 传递
    "zcm6-vasp-0": {
        "mode": "each-subdir",
        "dir_root": Path("./y_dir"),
        "chunks": 5,
        "module_profile_type": "zcm6-vasp-0",
        "launcher_type": "srun",
        "cmd": "pei_vasp_univ_sbatch",
        "if_use_my_launcher": True,
        "partition": "amd_512",
        "nodes": 1,
        "ncores": 128,
    },
    # zcm6-n2p2-scaling-0：nnp-scaling，单核估算并行扩展性
    "zcm6-n2p2-scaling-0": {
        "mode": "each-subdir",
        "dir_root": Path("./y_dir"),
        "chunks": 5,
        "module_profile_type": "zcm6-n2p2-0",
        "launcher_type": "mpirun",
        "cmd": "nnp-scaling 10000",
        "if_use_my_launcher": False,
        "partition": "amd_512",
        "nodes": 1,
        "ncores": 1,
    },
    # zcm6-n2p2-train-0：nnp-train，16–32 核是推荐区间
    "zcm6-n2p2-train-0": {
        "mode": "each-subdir",
        "dir_root": Path("./y_dir"),
        "chunks": 5,
        "module_profile_type": "zcm6-n2p2-0",
        "launcher_type": "mpirun",
        "cmd": "nnp-train",
        "if_use_my_launcher": False,
        "partition": "amd_512",
        "nodes": 1,
        "ncores": 24,
    },
    # zcm6-lammps-0：lmp -in lmp.in
    "zcm6-lammps-0": {
        "mode": "each-subdir",
        "dir_root": Path("./y_dir"),
        "chunks": 5,
        "module_profile_type": "zcm6-lammps-0",
        "launcher_type": "mpirun",
        "cmd": "lmp -in lmp.in",
        "if_use_my_launcher": False,
        "partition": "amd_512",
        "nodes": 1,
        "ncores": 24,
    },
}
# to here ================

# 受 preset / 命令行共同提供的「必填」字段。argparse 不再用 required=True 强制——否则 --preset
# 经 set_defaults 灌进来的值满足不了 required 检查；改为合并后在 check_required 里统一兜底。
REQUIRED_FIELDS = (
    "mode", "module_profile_type", "launcher_type",
    "cmd", "partition", "nodes", "ncores",
)


def fail(msg: str):
    # 统一失败出口：打印 + 退出码 1，调用处不必各自 raise。
    print("❌ ERROR: " + msg)
    raise SystemExit(1)


def resolve_preset(name: str) -> dict:
    # 查表取预设；未知名直接 fail 并列出可用项。返回副本，避免调用方改到全局表。
    if name not in PRESETS:
        fail("未知预设: " + repr(name) + "；可用：" + ", ".join(sorted(PRESETS)))
    return dict(PRESETS[name])


def print_presets():
    # --list-presets：列出所有预设及关键字段，方便挑选。
    print("📋 可用预设（--preset NAME）：")
    for name in sorted(PRESETS):
        p = PRESETS[name]
        print("  • " + name
              + "  →  module=" + str(p.get("module_profile_type"))
              + ", launcher=" + str(p.get("launcher_type"))
              + ", cmd=" + repr(p.get("cmd"))
              + ", nodes=" + str(p.get("nodes"))
              + ", ncores=" + str(p.get("ncores")))


def show_preset(name: str):
    # --show-preset NAME：摊开单个预设的全部字段，核对后再 --preset 调用。
    preset = resolve_preset(name)
    print("📋 预设 " + name + "：")
    for key in sorted(preset):
        print("  " + key + " = " + repr(preset[key]))


def check_required(args):
    # 合并 preset + 命令行之后兜底：必填字段缺一不可（任一来源提供即可）。
    missing = [f for f in REQUIRED_FIELDS if getattr(args, f, None) is None]
    if missing:
        fail("缺少必填参数：" + ", ".join("--" + m for m in missing)
             + "；请在命令行直接指定，或用 --preset 提供（可用预设见 --list-presets）。")


def parse_bool(value) -> bool:
    # argparse 的 type 只在「显式传了值」时才被调用；裸写 --if_sbatch 走的是 const=True，
    # 不会进这里。所以这里只管把 True/False（及常见同义词）字符串解析成 bool。
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in ("true", "1", "yes", "y", "on", "t"):
        return True
    if normalized in ("false", "0", "no", "n", "off", "f"):
        return False
    raise argparse.ArgumentTypeError(
        "if_sbatch 只接受布尔值 True/False（或 1/0、yes/no），收到: " + repr(value)
    )


def build_parser() -> argparse.ArgumentParser:
    # 选项名与 pei_slurm_univ_submit 形参一一对应；枚举的合法取值只写在 help 里做提示，
    # 真正的校验交给库函数，避免在 CLI 再抄一份枚举导致两边漂移。
    parser = argparse.ArgumentParser(
        prog="pei_slurm_univ_submit",
        description="通用 Slurm 作业脚本生成 / 提交引擎（mymetal.slurm 的 CLI 封装）。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            示例（对应库里的 test_args1/2/3 三种模式）：

              # parallel：每个子目录一个独立作业，只提交不等待
              pei_slurm_univ_submit --path_root /public3/home/scg6928/mywork/test \\
                  --mode parallel --dir_root ./y_dir --chunks 5 \\
                  --module_profile_type zcm6-vasp-0 --launcher_type srun \\
                  --cmd "echo Hello, World!" --partition amd_512 --nodes 2 --ncores 16 --if_sbatch

              # each-subdir：编排作业里逐个 sbatch --wait 子作业
              pei_slurm_univ_submit --path_root /public3/home/scg6928/mywork/test \\
                  --mode each-subdir --dir_root ./y_dir --chunks 5 \\
                  --module_profile_type zcm6-vasp-0 --launcher_type srun \\
                  --cmd pei_vasp_univ_sbatch --partition amd_512 --nodes 2 --ncores 16 --if_sbatch

              # single-alloc：单次分配内顺序跑各子目录
              pei_slurm_univ_submit --path_root /public3/home/scg6928/mywork/test \\
                  --mode single-alloc --dir_root ./y_dir --chunks 5 \\
                  --module_profile_type zcm6-vasp-0 --launcher_type srun \\
                  --cmd pei_vasp_univ_sbatch --partition amd_512 --nodes 2 --ncores 16 --if_sbatch

            预设（preset）：常见场景已封装在 PRESETS 里，可一行直呼，必填项由预设提供：

              pei_slurm_univ_submit --preset zcm6-vasp-0                  # 直接用预设
              pei_slurm_univ_submit --preset zcm6-vasp-0 --ncores 64      # 预设基础上覆盖个别项
              pei_slurm_univ_submit --preset zcm6-n2p2-train-0 --if_sbatch True
              pei_slurm_univ_submit --list-presets                      # 查看所有预设
              pei_slurm_univ_submit --show-preset zcm6-vasp-0             # 查看某预设全部字段

            不加 --preset 时，--mode/--module_profile_type/--launcher_type/--cmd/
            --partition/--nodes/--ncores 为必填；命令行显式值始终覆盖预设。
            不加 --if_sbatch 时只生成脚本、不提交（dry）。
            """
        ),
    )

    # —— 预设 / 发现 ——
    # --preset 先把某场景的默认值灌进来，命令行其余 flag 再覆盖；--list/--show 只读不跑。
    parser.add_argument("--preset",
                        help="使用内置预设填充默认值（命令行可逐项覆盖）。可用：" + ", ".join(sorted(PRESETS)) + "。")
    parser.add_argument("--list-presets", action="store_true",
                        help="列出所有内置预设后退出。")
    parser.add_argument("--show-preset", metavar="NAME",
                        help="打印指定预设的全部字段后退出。")

    # —— 运行的基本模式 ——
    parser.add_argument("--path_root", type=Path, default=Path.cwd(),
                        help="根目录，必须是绝对路径；相对路径会被解析为绝对。默认当前目录。")
    parser.add_argument("--mode",
                        help="运行模式：parallel / each-subdir / single-alloc。无 --preset 时必填。")
    parser.add_argument("--dir_root", type=Path, default=Path("./dir"),
                        help="作业子目录的父目录，相对 path_root。默认 ./dir。")
    parser.add_argument("--lsubdir", nargs="*", default=None,
                        help="子目录名列表（basename，需位于 dir_root 下）；留空则自动遍历 dir_root。")
    parser.add_argument("--chunks", type=int, default=1,
                        help="把作业目录分成多少块（仅 each-subdir / single-alloc 生效）。默认 1。")

    # —— script 相关参数 ——
    parser.add_argument("--module_profile_type",
                        help="环境 module profile：zcm6-vasp-0 / zcm6-lammps-0 / zcm6-n2p2-0。无 --preset 时必填。")
    parser.add_argument("--launcher_type",
                        help="并行启动器：srun / mpirun / none。无 --preset 时必填。")
    parser.add_argument("--cmd",
                        help="每个子目录里实际执行的命令。无 --preset 时必填。")
    parser.add_argument("--if_use_my_launcher", action="store_true",
                        help="通过 MY_LAUNCHER 环境变量传递启动器，而非在脚本里写死。")

    # —— 作业资源参数 ——
    parser.add_argument("--partition", help="Slurm 分区（-p）。无 --preset 时必填。")
    parser.add_argument("--nodes", type=int, help="节点数（-N）。无 --preset 时必填。")
    parser.add_argument("--ncores", type=int, help="核数（-n）。无 --preset 时必填。")

    # —— 提交开关 ——
    # 接收布尔值：--if_sbatch True / --if_sbatch False；裸写 --if_sbatch 等价于 True；
    # 完全不写则取默认 False（只生成脚本 dry run）。
    parser.add_argument("--if_sbatch", type=parse_bool, nargs="?", const=True, default=False,
                        help="是否真正 sbatch 提交：True/False（默认 False，只生成脚本 dry run）；"
                             "裸写 --if_sbatch 等价于 True。")

    return parser


def main():
    # ====== check ======
    # 第一段只认 --preset：用 parse_known_args 容忍此刻其余必填项尚未给全。
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--preset")
    known, _ = pre.parse_known_args()

    parser = build_parser()
    if known.preset is not None:
        # 预设灌默认；下面 parse_args 时命令行显式给的同名 flag 会再覆盖它。
        parser.set_defaults(**resolve_preset(known.preset))

    args = parser.parse_args()

    # 只读动作：列出 / 查看预设，打印完即走，不进入提交流程。
    if args.list_presets:
        print_presets()
        return
    if args.show_preset is not None:
        show_preset(args.show_preset)
        return

    # path_root 必须绝对：先把相对路径解析成绝对，库函数内部仍会再 check 一次兜底。
    args.path_root = args.path_root.resolve()

    # 合并 preset + 命令行后做结构性必填兜底（取代被去掉的 argparse required=True）。
    check_required(args)
    # ====== to here ======

    # ====== main ======
    # preset/list_presets/show_preset 是 CLI 自身的控制项、不是库函数形参，透传前剔除；
    # 其余 dest 与形参同名，直接整体透传。枚举 / None / 正整数等校验仍由库函数负责。
    for key in ("preset", "list_presets", "show_preset"):
        vars(args).pop(key, None)

    pei_slurm_univ_submit(**vars(args))
    # ====== to here ======


if __name__ == "__main__":
    main()

################################### ZCM6
# 以下常见用法已收敛进顶部的 PRESETS 注册表，改用 --preset 一行直呼，例如：
#   pei_slurm_univ_submit --preset zcm6-vasp-0
#   pei_slurm_univ_submit --preset zcm6-n2p2-scaling-0
#   pei_slurm_univ_submit --preset zcm6-n2p2-train-0 --ncores 32
#   pei_slurm_univ_submit --preset zcm6-lammps-0
# 发现 / 核对：--list-presets、--show-preset NAME。要新增场景就往 PRESETS 里加一项。

###### python 通用 sbatch 包装（与本 CLI 无关，仅备忘）
# python
# #!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 1
# /public3/home/scg6928/mysoft/env/pyenv/codex/bin/python "$@"

