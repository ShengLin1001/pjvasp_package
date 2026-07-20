"""Generate and submit tool-agnostic Slurm job scripts.

Functions:
    check_wall_time: Validate an optional Slurm wall-time value.
    generate_slurm_script_base: Generate one calculation job script.
    generate_slurm_script_sequential: Generate one sequential chunk worker.
    generate_slurm_script_shared_parent: Generate one parent for chunk workers.
    get_y_dir_lsubdir: Recursively discover jobs below every y_dir directory.
    split_chunks: Split calculation directories into balanced scheduling lanes.
    pei_slurm_univ_submit: Orchestrate script generation and optional submission.
"""

import os
import re
import shlex
import textwrap
from pathlib import Path

import numpy as np

from mymetal.universal.print.print import warn, fail
from mymetal.universal.check.type import check_positive_int, check_none, check_basename, check_absolute_path

# 启动失败重试包装器：slurm_utils/slurm_universal/pei_slurm_univ_launch_retry，与 cmd
# (如 pei_vasp_univ_sbatch) 一样靠 PATH 找到。集中成常量，别在生成 shell 处散写字符串。
LAUNCH_RETRY_WRAPPER = "pei_slurm_univ_launch_retry"
# sbatch「提交」失败重试包装器：slurm_utils/slurm_universal/pei_slurm_univ_sbatch_retry，同样
# 靠 PATH 找到。它只在「提交没被 slurmctld 接住」时轮询重试（默认 99 次 ×10s），一旦作业被
# 创建（输出含 "Submitted batch job" 或退出码 0）就原样透传退出码，绝不因作业算崩而重投。
SBATCH_RETRY_WRAPPER = "pei_slurm_univ_sbatch_retry"
CHUNK_PARENT_LAYOUTS = ["auto", "shared", "per-chunk"]
WALL_TIME_PATTERN = re.compile(r"^[0-9]+(?:-[0-9]+)?(?::[0-9]+){0,2}$")

# generate_script

def check_wall_time(wall_time: str = None) -> str | None:
    """检查可选的 Slurm wall time 是否可安全写入 ``#SBATCH --time``。

    这里只检查结构，具体时间范围和 partition/QOS 上限仍交给 Slurm。支持纯分钟数、
    ``HH:MM:SS`` 和 ``D-HH:MM:SS`` 等由数字、连字符和冒号组成的标准写法。

    Args:
        wall_time: Slurm wall-time 字符串；None 表示不生成 ``--time`` 行。

    Returns:
        去除首尾空白后的 wall-time 字符串，或 None。

    Raises:
        SystemExit: 当 wall_time 为空或格式非法时退出。
    """
    if wall_time is None:
        return None

    wall_time = str(wall_time).strip()
    if not wall_time or WALL_TIME_PATTERN.fullmatch(wall_time) is None:
        fail(
            "Invalid wall time: " + repr(wall_time)
            + "; use minutes, HH:MM:SS, or D-HH:MM:SS"
        )
    return wall_time


def generate_script_header(
        partition: str,
        nodes: int,
        ncores: int,
        module_profile_type: str,
        MODULE_BLOCKS: dict[str, str],
        wall_time: str = None) -> str:
    """生成 Slurm 脚本头部和环境模块加载部分。

    Args:
        partition: Slurm 分区名称。
        nodes: 申请的节点数。
        ncores: 申请的任务数或核心数。
        module_profile_type: 环境模块配置类型，对应 MODULE_BLOCKS 中的键。
        MODULE_BLOCKS: 环境模块配置字典，键为配置类型，值为对应的 shell 代码块。
        wall_time: 可选的 Slurm 最大运行时间；None 时不生成 ``--time`` 行。

    Returns:
        生成的 Bash 脚本头部字符串。

    Raises:
        SystemExit: 当 module_profile_type 不存在于 MODULE_BLOCKS 中时退出。
    """
    wall_time = check_wall_time(wall_time)
    line_wall_time = ""
    if wall_time is not None:
        line_wall_time = "#SBATCH --time=" + wall_time + "\n"

    header = (
        "#!/bin/bash\n"
        "#SBATCH -p " + str(partition) + "\n"
        "#SBATCH -N " + str(nodes) + "\n"
        "#SBATCH -n " + str(ncores) + "\n"
        + line_wall_time
    )
    module_block = MODULE_BLOCKS.get(module_profile_type)
    if module_block is None:
        fail("Unknown module profile type: " + str(module_profile_type))
    header += module_block
    return header

def generate_launcher_command(launcher: str, cmd: str, if_use_my_launcher: bool = False,
                              launch_retry_wrapper: str = LAUNCH_RETRY_WRAPPER):
    """生成用于启动目标计算命令的 shell 代码。

    Args:
        launcher: 启动器类型，例如 "srun"、"mpirun" 或 "none"。
        cmd: 实际需要执行的计算命令。
        if_use_my_launcher: 是否通过 MY_LAUNCHER 环境变量传递启动器命令。
        launch_retry_wrapper: 包在真正 launcher 外面的启动失败重试包装器；置空则不包。

    Returns:
        生成的 shell 命令字符串。

    Raises:
        SystemExit: 当 launcher 类型未知时退出。
    """
    line_launcher = ""

    # 编译器 / MPI 抖动会让 srun / mpirun 偶发起不来（作业步创建被临时禁用等），程序根本没执行。
    # 包一层 pei_slurm_univ_launch_retry：它只在命中启动失败特征时重试（默认 99 次 × 10s），
    # 程序真跑起来了才失败的非 0（VASP 不收敛、输入错误、被杀）原样透传退出码、绝不重试，
    # 否则必崩的算例会占着整个分配反复重跑。重试次数/间隔由 PEI_LAUNCH_RETRY_MAX /
    # PEI_LAUNCH_RETRY_SLEEP 环境变量调，不必重新生成脚本。
    prefix = launch_retry_wrapper + " " if launch_retry_wrapper else ""

    # 裸变量写法两处通用：$SLURM_NTASKS 是纯数字，无需内层引号防 word splitting。
    if launcher == "srun":
        line_launcher_1 = prefix + 'srun -n $SLURM_NTASKS'
    elif launcher == "mpirun":
        line_launcher_1 = prefix + 'mpirun -np $SLURM_NTASKS'
    elif launcher == "none":
        # 没有 launcher 就没有「启动失败」可言，cmd 直接在本 shell 里跑，不包重试。
        line_launcher_1 = ""
    else:
        fail("Unknown launcher type: " + str(launcher))

    line_launcher_2 = cmd

    # 对于一些封装好的脚本，可能需要通过环境变量传递launcher类型，而不是直接在脚本里写死，例如 srun -n "$SLURM_NTASKS" cmd
    # 而对于 n2p2, vasp 一些直接调用可执行文件的，直接在脚本里写死 srun -n "$SLURM_NTASKS" cmd 可能更方便
    if if_use_my_launcher:
        # 方案A1：用双引号包裹。该脚本运行在计算节点的 Slurm 作业中，此刻 SLURM_NTASKS
        # 已由 Slurm 设好，于是在 export 这一刻就把它展开成真实核数，直接烤进 MY_LAUNCHER
        # 的值（例如 'pei_slurm_univ_launch_retry srun -n 16'）。消费端（如 pei_vasp_univ_sbatch）
        # 拿到的已是现成数字，无需再做二次展开 / eval。注意不能用单引号，否则 $SLURM_NTASKS 会被原样保存。
        # 重试包装器在这里只是 MY_LAUNCHER 的一段前缀：消费端 `$launcher "$exe"` 不加引号，
        # 词分割后自然展开成 `pei_slurm_univ_launch_retry srun -n 16 vasp_std`，消费端无需改动。
        line_launcher = ('export MY_LAUNCHER="' + line_launcher_1 + '"\n'
                         + line_launcher_2)
    else:
        # 直接写入分支：该行在作业里执行，运行时 $SLURM_NTASKS 同样会展开。
        line_launcher = line_launcher_1 + "  " + line_launcher_2
    return line_launcher

def generate_slurm_script_base(partition: str, nodes: int, ncores: int, module_profile_type: str, MODULE_BLOCKS: dict[str, str],
                        launcher: str, cmd: str, if_use_my_launcher: bool,
                        if_output: bool = True, path_save: Path = Path('./sub_slurm_univ.sh'),
                        wall_time: str = None):
    """生成单个计算目录使用的基础 Slurm 脚本。

    Args:
        partition: Slurm 分区名称。
        nodes: 申请的节点数。
        ncores: 申请的任务数或核心数。
        module_profile_type: 环境模块配置类型，对应 MODULE_BLOCKS 中的键。
        MODULE_BLOCKS: 环境模块配置字典，键为配置类型，值为对应的 shell 代码块。
        launcher: 启动器类型，例如 "srun"、"mpirun" 或 "none"。
        cmd: 在作业脚本中执行的计算命令。
        if_use_my_launcher: 是否通过 MY_LAUNCHER 环境变量传递启动器命令。
        if_output: 是否将生成的脚本写入文件。
        path_save: 生成脚本的保存路径。
        wall_time: 计算子作业的最大运行时间；None 时不生成 ``--time`` 行。

    Returns:
        生成的 Slurm 脚本内容字符串。
    """
    line_header = generate_script_header(
        partition, nodes, ncores, module_profile_type, MODULE_BLOCKS,
        wall_time=wall_time,
    )
    line_launcher = generate_launcher_command(launcher, cmd, if_use_my_launcher)
    line_myheader = (
        "# Auto-generated by pei_slurm_univ_submit (base job script). Do not edit by hand.\n"
        "# 该脚本由程序自动生成，每次运行都会被覆盖。\n"
    )
    line = line_header + "\n\n" + line_myheader + "\n\n" + line_launcher + "\n"

    if if_output:
        with open(path_save, 'w') as f:
            f.write(line)

    return line

def generate_slurm_script_sequential(partition: str, nodes: int, ncores: int, module_profile_type: str, MODULE_BLOCKS: dict[str, str],
                        launcher: str, cmd: str, if_use_my_launcher: bool,
                        group: list, mode: str,
                        if_output: bool = True, path_save: Path = Path('./sub_slurm_each_subdir.sh'),
                        wall_time: str = None):
    """生成按子目录顺序执行的 Slurm 编排脚本。

    Args:
        partition: Slurm 分区名称。
        nodes: 申请的节点数。
        ncores: 申请的任务数或核心数。
        module_profile_type: 环境模块配置类型，对应 MODULE_BLOCKS 中的键。
        MODULE_BLOCKS: 环境模块配置字典，键为配置类型，值为对应的 shell 代码块。
        launcher: 启动器类型，例如 "srun"、"mpirun" 或 "none"。
        cmd: 在每个子目录中执行的计算命令。
        if_use_my_launcher: 是否通过 MY_LAUNCHER 环境变量传递启动器命令。
        group: 当前编排脚本负责处理的子目录列表。
        mode: 执行模式，例如 "each-subdir" 或 "single-alloc"。
        if_output: 是否将生成的脚本写入文件。
        path_save: 生成脚本的保存路径。
        wall_time: 父/编排脚本的最大运行时间；None 时不生成 ``--time`` 行。

    Returns:
        生成的 Slurm 编排脚本内容字符串。
    """
    line_header = generate_script_header(
        partition, nodes, ncores, module_profile_type, MODULE_BLOCKS,
        wall_time=wall_time,
    )
    line_launcher = generate_launcher_command(launcher, cmd, if_use_my_launcher)
    line_myheader = (
        "# Auto-generated by pei_slurm_univ_submit (" + str(mode)
        + " job script). Do not edit by hand.\n"
        "# 该脚本由程序自动生成，每次运行都会被覆盖。\n"
    )
    line_loop = generate_loop(group, mode, line_launcher)
    line = line_header + "\n\n" + line_myheader + "\n\n" + line_loop

    if if_output:
        with open(path_save, 'w') as f:
            f.write(line)

    return line


def generate_slurm_script_shared_parent(
        partition: str,
        module_profile_type: str,
        MODULE_BLOCKS: dict[str, str],
        lpath_worker: list[Path],
        if_output: bool = True,
        path_save: Path = Path("./sub_slurm_each_subdir_parent.sh"),
        wall_time: str = None) -> str:
    """生成一个并发运行多个 each-subdir chunk worker 的父脚本。

    worker 仍是可独立 ``sbatch`` 的历史 chunk 脚本；shared parent 只把它们作为普通
    Bash 子进程并发启动。这样调度器只看见一个 ``-n 1`` 编排父作业，而每条 worker
    仍依次执行自己的 ``sbatch --wait`` 子作业链。

    Args:
        partition: Slurm 分区名称。
        module_profile_type: 环境模块配置类型。
        MODULE_BLOCKS: 环境模块配置字典。
        lpath_worker: 本次生成的 chunk worker 绝对路径列表。
        if_output: 是否将脚本写入文件。
        path_save: shared parent 脚本保存路径。
        wall_time: shared parent 的最大运行时间；None 时不生成 ``--time`` 行。

    Returns:
        生成的 Bash 脚本内容。

    Raises:
        SystemExit: 当 worker 列表为空或路径不是绝对路径时退出。
    """
    if not lpath_worker:
        fail("generate_slurm_script_shared_parent: worker list is empty")

    lpath_worker = [Path(path_worker) for path_worker in lpath_worker]
    for path_worker in lpath_worker:
        check_absolute_path(path_worker)

    line_header = generate_script_header(
        partition, 1, 1, module_profile_type, MODULE_BLOCKS,
        wall_time=wall_time,
    )
    line_myheader = (
        "# Auto-generated by pei_slurm_univ_submit (shared each-subdir parent). Do not edit by hand.\n"
        "# 该脚本由程序自动生成，每次运行都会被覆盖。\n"
    )
    line_worker_array = "worker_scripts=(\n"
    for path_worker in lpath_worker:
        line_worker_array += "    " + shlex.quote(str(path_worker)) + "\n"
    line_worker_array += ")\n"

    # 每个 worker 单独写日志，避免 K 条后台流同时写父 stdout；父 stdout 只保留启动信息和
    # 最终聚合。job_tag 在非 Slurm 手动检查时退回 PID，防止日志名为空。
    line_parent = (
        'start_dir="$(pwd)"\n'
        'path_log_dir=' + shlex.quote(str(path_save.parent.resolve())) + '\n'
        'job_tag="${SLURM_JOB_ID:-manual-$$}"\n'
        "worker_pids=()\n"
        "worker_names=()\n"
        "worker_logs=()\n"
        "failed_workers=()\n"
        "\n"
        "cancel_workers() {\n"
        '    echo "⚠️  父作业收到终止信号，停止本地 chunk worker；已提交的 Slurm 子作业不会在此处自动取消。" >&2\n'
        '    trap - INT TERM\n'
        '    for pid in "${worker_pids[@]}"; do\n'
        '        kill "$pid" 2>/dev/null || true\n'
        "    done\n"
        "    wait 2>/dev/null || true\n"
        "    exit 143\n"
        "}\n"
        "trap cancel_workers INT TERM\n"
        "\n"
        'for path_worker in "${worker_scripts[@]}"; do\n'
        '    worker_name="$(basename "$path_worker" .sh)"\n'
        '    path_log="$path_log_dir/slurm-${job_tag}-${worker_name}.out"\n'
        '    echo "================ ▶️ 启动 $worker_name"\n'
        '    echo "  script: $path_worker"\n'
        '    echo "  log:    $path_log"\n'
        "    (\n"
        '        cd "$start_dir" || exit 1\n'
        '        bash "$path_worker"\n'
        '    ) > "$path_log" 2>&1 &\n'
        '    worker_pids+=("$!")\n'
        '    worker_names+=("$worker_name")\n'
        '    worker_logs+=("$path_log")\n'
        "done\n"
        "\n"
        "worker_success=0\n"
        "worker_failed=0\n"
        'for ((index = 0; index < ${#worker_pids[@]}; index++)); do\n'
        '    pid="${worker_pids[$index]}"\n'
        '    worker_name="${worker_names[$index]}"\n'
        '    path_log="${worker_logs[$index]}"\n'
        '    if wait "$pid"; then\n'
        '        echo "✅ $worker_name 完成；日志: $path_log"\n'
        "        worker_success=$((worker_success + 1))\n"
        "    else\n"
        "        status=$?\n"
        '        echo "❌ $worker_name 失败 (exit $status)；日志: $path_log" >&2\n'
        "        worker_failed=$((worker_failed + 1))\n"
        '        failed_workers+=("$worker_name")\n'
        "    fi\n"
        "done\n"
        "trap - INT TERM\n"
        "\n"
        'echo ""\n'
        'echo "================ 📊 Summary (mode=each-subdir, parent=shared)"\n'
        'echo "worker_success=$worker_success    worker_failed=$worker_failed"\n'
        'if (( ${#failed_workers[@]} > 0 )); then\n'
        '    echo "❌ 失败的 chunk worker:"\n'
        '    for worker_name in "${failed_workers[@]}"; do echo "  - $worker_name"; done\n'
        "    exit 1\n"
        "fi\n"
        'echo "🎉 Done: 所有 chunk worker 均已完成！"\n'
        "exit 0\n"
    )
    line = line_header + "\n\n" + line_myheader + "\n\n" + line_worker_array + "\n" + line_parent

    if if_output:
        with open(path_save, "w") as file_out:
            file_out.write(line)

    return line


def generate_loop(group: list, mode: str, line_launcher: str):
    """生成遍历多个计算子目录的 Bash 循环代码。

    Args:
        group: 需要处理的计算子目录列表。
        mode: 执行模式，例如 "each-subdir" 或 "single-alloc"。
        line_launcher: 用于执行计算任务的 shell 命令代码块。

    Returns:
        包含目录遍历、状态统计和最终汇总的 Bash 循环字符串。
    """
    # subdir 已经是绝对路径了，在 get_lsubdir 中检查过了。
    # 用普通字符串拼接（而非 f-string）来写 bash，避免 $((...)) / ${...} 的花括号
    # 被 Python 当成替换字段；只在注入 Python 值（subdirs / mode / line_launcher）处用 +。
    lsubdir = [str(subdir) for subdir in group]
    subdirs = " ".join(shlex.quote(subdir) for subdir in lsubdir)

    # —— 计数器 + 失败/未收敛目录数组 + 目录列表 + 循环开头 ——
    line_loop = (
        'start_dir="$(pwd)"\n'
        "submit_success=0\n"
        "submit_failed=0\n"
        "submit_success_convergenced=0\n"
        "submit_success_not_convergenced=0\n"
        "failed_dirs=()\n"
        "not_convergenced_dirs=()\n"
        "lsubdir=(" + subdirs + " )\n"
        "\n"
        'for subdir in "${lsubdir[@]}"; do\n'
        '    cd "$start_dir" || exit 1\n'
        '    echo ""\n'
        '    echo "================ 📁 $subdir"\n'
        '    if ! cd "$subdir"; then\n'
        '        echo "❌ ERROR: 无法进入目录: $subdir" >&2\n'
        '        submit_failed=$((submit_failed + 1)); failed_dirs+=("$subdir")\n'
        "        continue\n"
        "    fi\n"
        '    echo "📍 当前目录: $(pwd)"\n'
    )

    if mode == "each-subdir":
        # 子作业经 sbatch --wait 提交，但外面包一层 pei_slurm_univ_sbatch_retry：提交若没被
        # slurmctld 接住（暂态繁忙 / 超时 / 通信抖动），它会轮询重试（默认 99 次 ×10s）而不是
        # 像以前那样直接跳过本子目录；重试耗尽（或命中永久性错误）仍未提交成功，才判失败并 skip。
        # 作业一旦被创建（输出含 Submitted batch job），--wait 的退出码即作业本身结果，据此判收敛，
        # 绝不因作业算崩而重投（否则会把必崩算例重跑 99 遍）。
        line_loop += (
            '    echo "▶️  ' + SBATCH_RETRY_WRAPPER + ' --wait sub_slurm_univ.sh"\n'
            '    sbatch_out="$(' + SBATCH_RETRY_WRAPPER + ' --wait sub_slurm_univ.sh 2>&1)"; status=$?\n'
            '    echo "$sbatch_out"\n'
            '    if ! grep -q "Submitted batch job" <<< "$sbatch_out"; then\n'
            '        echo "❌ 提交失败（已轮询重试仍未成功）: $subdir" >&2\n'
            '        submit_failed=$((submit_failed + 1)); failed_dirs+=("$subdir")\n'
            "        continue\n"
            "    fi\n"
            "    submit_success=$((submit_success + 1))\n"
        )
    elif mode == "single-alloc":
        # 单一分配内直接运行命令：cd 成功即视为运行成功，再按退出码判断收敛。
        line_loop += (
            '    echo "▶️  run (single-alloc)"\n'
            + textwrap.indent(line_launcher, "    ") + "\n"
            "    status=$?\n"
            "    submit_success=$((submit_success + 1))\n"
        )

    # —— 依据退出码判定收敛：pei_vasp_univ_sbatch 约定 0/10=已收敛，其余=未收敛 ——
    line_loop += (
        "    if (( status == 0 || status == 10 )); then\n"
        '        echo "✅ 已收敛 (exit $status)"\n'
        "        submit_success_convergenced=$((submit_success_convergenced + 1))\n"
        "    else\n"
        '        echo "❌ 未收敛 (exit $status): $subdir" >&2\n'
        "        submit_success_not_convergenced=$((submit_success_not_convergenced + 1))\n"
        '        not_convergenced_dirs+=("$subdir")\n'
        "    fi\n"
        "done\n"
    )

    # —— 最终 summary：四个量 + 挨个列出 not_convergenced_dirs 与 failed_dirs ——
    line_loop += (
        'cd "$start_dir" || exit 1\n'
        'echo ""\n'
        'echo "================ 📊 Summary (mode=' + mode + ')"\n'
        'echo "submit_success=$submit_success    submit_failed=$submit_failed'
        "    convergenced=$submit_success_convergenced"
        '    not_convergenced=$submit_success_not_convergenced"\n'
        "if (( ${#not_convergenced_dirs[@]} > 0 )); then\n"
        '    echo "❌ 未收敛的目录:"\n'
        '    for d in "${not_convergenced_dirs[@]}"; do echo "  - 📁 $d"; done\n'
        "fi\n"
        "if (( ${#failed_dirs[@]} > 0 )); then\n"
        '    echo "❌ 提交/运行失败的目录:"\n'
        '    for d in "${failed_dirs[@]}"; do echo "  - 📁 $d"; done\n'
        "fi\n"
        "if (( submit_failed > 0 || submit_success_not_convergenced > 0 )); then\n"
        '    echo "❌ Done: 存在失败或未收敛的计算。"\n'
        "    exit 1\n"
        "fi\n"
        'echo "🎉 Done: 全部完成且收敛！"\n'
        "exit 0\n"
    )
    return line_loop


# prepare parts
def get_y_dir_lsubdir(dir_root: Path) -> list[Path]:
    """递归发现 ``dir_root`` 下所有 ``y_dir`` 的一级计算子目录。

    搜索包含 ``dir_root`` 自身，因此既支持从工作流根目录统一发现多层 mode，
    也支持把某个 ``y_dir`` 本身直接作为搜索根。只把每个 ``y_dir`` 的一级子目录
    当作计算目录，不会把更深层的普通目录误当作独立作业。

    Args:
        dir_root: 递归搜索根目录，必须为绝对路径。

    Returns:
        按绝对路径排序的计算子目录列表。

    Raises:
        SystemExit: 当搜索根不存在、未发现 ``y_dir`` 或其中没有计算子目录时退出。
    """
    path_dir_root = Path(dir_root)
    check_absolute_path(path_dir_root)
    if not path_dir_root.is_dir():
        fail(f"dir_root {path_dir_root} is not a directory")

    ly_dir = []
    if path_dir_root.name == "y_dir":
        ly_dir.append(path_dir_root)
    ly_dir.extend(
        path_y_dir
        for path_y_dir in path_dir_root.rglob("y_dir")
        if path_y_dir.is_dir()
    )
    ly_dir = sorted(set(ly_dir))
    if not ly_dir:
        fail(
            f"No y_dir directories found recursively under dir_root {path_dir_root}"
        )

    lsubdir = []
    print(f"📁 Found {len(ly_dir)} y_dir directories under {path_dir_root}")
    for path_y_dir in ly_dir:
        ljob_dir = sorted(path_child for path_child in path_y_dir.iterdir()
                          if path_child.is_dir())
        print(f"  - {path_y_dir}: {len(ljob_dir)} job directories")
        lsubdir.extend(ljob_dir)

    lsubdir = sorted(set(lsubdir))
    if not lsubdir:
        fail(
            "No job subdirectories found below y_dir directories under "
            + str(path_dir_root)
        )

    print(f"📊 Discovered {len(lsubdir)} job directories in total")
    return lsubdir


def get_lsubdir(lsubdir: list[str] = None, dir_root: Path = None):
    """解析并检查计算子目录列表。

    Args:
        lsubdir: 可选的计算子目录 basename 过滤列表；同名目录会在所有 ``y_dir``
            中一并选中。为空时保留递归发现的全部计算目录。
        dir_root: ``y_dir`` 的递归搜索根。

    Returns:
        经过检查的绝对子目录路径列表。

    Raises:
        SystemExit: 当 dir_root 不存在、没有可用子目录，或子目录路径非法时退出。
    """
    path_dir_root = Path(dir_root)
    check_absolute_path(path_dir_root)
    ldiscovered = get_y_dir_lsubdir(path_dir_root)

    # 显式 basename 是递归结果过滤器；多个 mode 中的同名 case 必须全部保留。
    if lsubdir is None or len(lsubdir) == 0:
        lsubdir = ldiscovered
    else:
        for subdir in lsubdir:
            check_basename(subdir)
        set_requested = set(lsubdir)
        lsubdir = [
            path_subdir for path_subdir in ldiscovered
            if path_subdir.name in set_requested
        ]
        set_found = {path_subdir.name for path_subdir in lsubdir}
        lmissing = sorted(set_requested - set_found)
        if lmissing:
            fail(
                "Requested subdirectory basename(s) not found below any y_dir: "
                + ", ".join(lmissing)
            )
        print(f"📊 Selected {len(lsubdir)} job directories by --lsubdir")

    # 检查是不是绝对路径，这很重要
    for subdir in lsubdir:
        if not subdir.is_dir():
            fail(f"subdir {subdir} is not a directory")
        check_absolute_path(subdir)

    return lsubdir


def split_chunks(chunks: int, lsubdir: list[Path]) -> list[list[Path]]:
    """将计算子目录均匀划分为若干组。

    Args:
        chunks: 需要划分的组数。
        lsubdir: 待划分的计算子目录路径列表。

    Returns:
        子目录分组列表，每个元素是一组子目录路径。

    Raises:
        SystemExit: 当子目录列表为空、chunks 非法，或内部划分结果不一致时退出。
    """
    total = len(lsubdir)
    if total == 0:
        fail("split_chunks: no directories to split")
    if chunks < 1:
        fail("split_chunks: chunk count must be >= 1 (got %d)" % chunks)
    if chunks > total:
        fail("--chunks (%d) cannot exceed number of job directories (%d)" % (chunks, total))
    base, rem = divmod(total, chunks)
    lgroup = []
    offset = 0
    num_check = 0
    for c in range(1, chunks + 1):
        count = base + (1 if c <= rem else 0)
        group = lsubdir[offset:offset + count]
        lgroup.append(group)
        num_check += len(group)
        offset += count
    if np.abs(num_check - total) > 1e-10:
        fail("split_chunks: internal error: num_check %d != total %d" % (num_check, total))
    return lgroup


def check_chunk_parent_layout(
        chunk_parent_layout: str,
        mode: str,
        lchunk_parent_layout: list[str] = CHUNK_PARENT_LAYOUTS) -> str:
    """检查并解析 chunk worker 的父作业布局。

    ``auto`` 只在 each-subdir 下收拢为一个 shared 父作业；single-alloc 的每个
    chunk 必须持有自己的计算资源分配，因此始终解析为历史 ``per-chunk`` 布局。

    Args:
        chunk_parent_layout: ``auto``、``shared`` 或 ``per-chunk``。
        mode: 当前提交模式。
        lchunk_parent_layout: 允许的布局列表。

    Returns:
        解析后的 ``shared`` 或 ``per-chunk``。

    Raises:
        SystemExit: 当布局未知，或尝试为非 each-subdir 使用 shared 时退出。
    """
    if chunk_parent_layout not in lchunk_parent_layout:
        fail("Unknown chunk parent layout: " + str(chunk_parent_layout))

    if chunk_parent_layout == "auto":
        if mode == "each-subdir":
            return "shared"
        return "per-chunk"

    if chunk_parent_layout == "shared" and mode != "each-subdir":
        fail("chunk_parent_layout=shared is only supported for each-subdir mode")

    return chunk_parent_layout


# main body
def pei_slurm_univ_submit(
                    path_root: Path = None,
                    # 运行的基本模式
                    mode: str = None,
                    dir_root: Path = Path('./dir'),
                    lsubdir: str = None,
                    chunks: int = 1,
                    chunk_parent_layout: str = "auto",
                    # script 相关参数
                    module_profile_type: str = None,
                    launcher_type: str = None,
                    cmd: str = None,
                    if_use_my_launcher: bool = False,  # 如果打开，将通过定义 MY_LAUNCHER 环境变量来传递 launcher_type，而不是直接在脚本中写死例如 srun -n "$SLURM_NTASKS" cmd
                    # 作业资源参数
                    partition: str = None,
                    nodes: int = None,
                    ncores: int = None,
                    # check markers dict
                    if_sbatch: bool = False,
                    child_wall_time: str = None,
                    parent_wall_time: str = None,

                    # work_dir_name, dry_run, show_script, skip_if_file_contains

                    # global variables
                    MODULE_BLOCKS: dict[str, str] = None,
                    LAUNCHERS: list = ["srun", "mpirun", "none"],
                    MODES: list = ["parallel", "each-subdir", "single-alloc"],
                    CHUNK_PARENT_LAYOUTS: list = CHUNK_PARENT_LAYOUTS,

                    ):
    """生成并可选提交多个计算子目录的 Slurm 作业脚本。

    Args:
        path_root: 项目根目录，必须为绝对路径。
        mode: 运行模式，可选 "parallel"、"each-subdir" 或 "single-alloc"。
        dir_root: path_root 下的递归搜索根；为空列表时发现其中所有 ``y_dir``。
        lsubdir: 可选的计算子目录 basename 过滤列表；为空时汇总所有 ``y_dir`` 的一级子目录。
        chunks: 在顺序编排模式下将子目录划分的组数。
        chunk_parent_layout: chunk worker 的父作业布局；auto 在 each-subdir 下使用
            shared，在 single-alloc 下保持 per-chunk。
        module_profile_type: 环境模块配置类型，对应 MODULE_BLOCKS 中的键。
        launcher_type: 启动器类型，例如 "srun"、"mpirun" 或 "none"。
        cmd: 在每个计算目录中执行的命令。
        if_use_my_launcher: 是否通过 MY_LAUNCHER 环境变量向下游脚本传递启动器命令。
        partition: Slurm 分区名称。
        nodes: 申请的节点数。
        ncores: 申请的任务数或核心数。
        if_sbatch: 是否调用 sbatch 提交生成的脚本。
        child_wall_time: parallel / each-subdir 计算子作业的最大运行时间；None 时不在
            ``sub_slurm_univ.sh`` 中生成 ``#SBATCH --time``。single-alloc 下忽略。
        parent_wall_time: each-subdir / single-alloc 父脚本的最大运行时间；None 时不生成
            ``#SBATCH --time``。parallel 下忽略。
        MODULE_BLOCKS: 环境模块配置字典，键为配置类型，值为对应的 shell 代码块。
        LAUNCHERS: 允许使用的启动器类型列表。
        MODES: 允许使用的运行模式列表。
        CHUNK_PARENT_LAYOUTS: 允许使用的 chunk 父作业布局列表。

    Returns:
        None.

    Raises:
        SystemExit: 当必要参数缺失、路径非法、模式非法、资源参数非法或启动器类型未知时退出。

    Notes:
        值得注意的是当选择each-subdir模式时，ncores 应为单独子作业的核数，而不是总核数。父作业的核数将被自动设置为 1；而在single-alloc模式下，ncores 是整个分配的总核数。
    """
    ################################### Check
    # check if None
    for temp_value, temp_name in [(path_root, "path_root"), (mode, "mode"), (module_profile_type, "module_profile_type"),
                                  (launcher_type, "launcher_type"), (cmd, "cmd"), (partition, "partition"),
                                  (chunks, "chunks"), (nodes, "nodes"), (ncores, "ncores")]:
        check_none(temp_value, temp_name)

    # path_root是一切的根基，必须是绝对路径
    check_absolute_path(path_root)

    # check if positive int（必须接住返回值，把字符串等输入强制转成 int，
    # 否则后续 split_chunks 的 chunks < 1 / divmod 会因类型不符而报错）
    chunks = check_positive_int(chunks, "chunks")
    nodes = check_positive_int(nodes, "nodes")
    ncores = check_positive_int(ncores, "ncores")

    # check mode
    if mode not in MODES:
        fail(f"Unknown mode: {mode}")

    child_wall_time = check_wall_time(child_wall_time)
    if mode == "single-alloc" and child_wall_time is not None:
        warn("--child_wall_time ignored in single-alloc mode (no calculation child jobs)")

    parent_wall_time = check_wall_time(parent_wall_time)
    if mode == "parallel" and parent_wall_time is not None:
        warn("--parent_wall_time ignored in parallel mode (no parent job)")

    # chunks 表示并发调度流数量；父作业布局另行解析，避免把两个概念继续绑定。
    chunk_parent_layout = check_chunk_parent_layout(
        chunk_parent_layout, mode, CHUNK_PARENT_LAYOUTS
    )

    # check module_profile_type
    if module_profile_type not in list(MODULE_BLOCKS.keys()):
        fail(f"Unknown module profile type: {module_profile_type}")

    # check launcher_type
    if launcher_type not in LAUNCHERS:
        fail(f"Unknown launcher type: {launcher_type}")
    ################################### Check to here

    os.chdir(path_root)
    print(f"path_root: {path_root}")

    # 把当前路径全部转化为绝对路径
    dir_root = path_root / dir_root

    ################################### Prepare
    # get lsubdir
    lsubdir = get_lsubdir(lsubdir, dir_root) # 已经被检查了确认都存在文件夹
    lgroup = split_chunks(chunks, lsubdir)
    ################################### Prepare to here

    ################################### Main control flow
    # 不使用check_file_contain参数，这是在cmd里应该执行的，而非提前去检查

    # 提前生成所需的base slurm script
    if mode in ["parallel", "each-subdir"]:
        for subdir in lsubdir:
            path_save = subdir / "sub_slurm_univ.sh"
            generate_slurm_script_base(partition, nodes, ncores, module_profile_type, MODULE_BLOCKS,
                                        launcher_type, cmd, if_use_my_launcher,
                                        if_output=True, path_save=path_save,
                                        wall_time=child_wall_time)

    if mode == "parallel":
        if chunks > 1:
            warn("--chunks ignored in parallel mode (already one job per dir)")
        if if_sbatch:
            # parallel 只负责提交、不等待结果，因此只统计 submit_success / submit_failed。
            submit_success = 0
            submit_failed = 0
            submit_failed_dirs = []
            for subdir in lsubdir:
                path_save = subdir / "sub_slurm_univ.sh"
                print("")
                print(f"================ 📁 {subdir}")
                print(f"▶️  {SBATCH_RETRY_WRAPPER} {path_save}")
                # 这个 bash 脚本结束就会回到 path_root。sbatch 提交经 pei_slurm_univ_sbatch_retry
                # 轮询重试：提交被 slurmctld 暂态拒于门外时不立刻记失败，重试耗尽才判失败。
                rc = os.system(f"cd {subdir} && {SBATCH_RETRY_WRAPPER} {path_save}")
                exit_code = os.waitstatus_to_exitcode(rc)
                if exit_code == 0:
                    print("✅ 提交成功")
                    submit_success += 1
                else:
                    print(f"❌ 提交失败 (exit {exit_code}): {subdir}")
                    submit_failed += 1
                    submit_failed_dirs.append(subdir)
            # —— summary：两个量 + 挨个列出 submit_failed_dirs ——
            print("")
            print("================ 📊 Summary (mode=parallel)")
            print(f"submit_success={submit_success}    submit_failed={submit_failed}")
            if submit_failed_dirs:
                print("❌ 提交失败的目录:")
                for d in submit_failed_dirs:
                    print(f"  - 📁 {d}")
            print("❌ Done with submission failures." if submit_failed > 0 else "🎉 Done!")
    elif mode in ["each-subdir", "single-alloc"]:
        # each-subdir 的 chunk 脚本只负责循环 + sbatch --wait 子作业，本身不跑计算，
        # 实际计算资源由各子 base 脚本申请，因此 worker / shared parent 固定为 -N 1 -n 1。
        parent_nodes = nodes
        parent_ncores = ncores
        if mode == "each-subdir" and (nodes != 1 or ncores != 1):
            warn(f"each-subdir 编排脚本无需多资源，已将 nodes={nodes}, ncores={ncores} 修正为 1")
            parent_nodes = 1
            parent_ncores = 1
        # 编排脚本(chunk)及其 slurm-<jobid>.out 单独收进 path_root/slurm/ —— 否则会
        # 直接堆在 path_root 里，与后处理产物(y_post_*.txt 等)混在一起。
        # 从 slurm_dir 里 sbatch，默认的 slurm-*.out 也就落在这里；子作业(base script)经
        # sbatch --wait 在各自 subdir 内提交，其输出仍留在对应 subdir，不受影响。
        slurm_dir = path_root / "slurm"
        slurm_dir.mkdir(parents=True, exist_ok=True)
        # 始终保留历史 chunk 文件名：既能被 shared parent 当普通 Bash worker 调用，也能在
        # per-chunk 兼容模式或故障恢复时单独 sbatch。只把本次生成的明确列表交给 shared
        # parent，不能 glob slurm_dir，避免旧运行残留的更高编号 chunk 被误执行。
        lpath_worker = []
        for chunk_id, group in enumerate(lgroup, start=1):
            path_worker = (
                slurm_dir
                / ("sub_slurm_" + mode.replace("-", "_")
                   + "_chunk" + str(chunk_id).zfill(3) + ".sh")
            )
            generate_slurm_script_sequential(
                partition, parent_nodes, parent_ncores,
                module_profile_type, MODULE_BLOCKS,
                launcher_type, cmd, if_use_my_launcher,
                group=group, mode=mode,
                if_output=True, path_save=path_worker,
                wall_time=parent_wall_time,
            )
            lpath_worker.append(path_worker.resolve())

        lpath_submit = list(lpath_worker)
        if (mode == "each-subdir"
                and chunk_parent_layout == "shared"
                and len(lpath_worker) > 1):
            path_parent = slurm_dir / "sub_slurm_each_subdir_parent.sh"
            generate_slurm_script_shared_parent(
                partition, module_profile_type, MODULE_BLOCKS,
                lpath_worker=lpath_worker,
                if_output=True, path_save=path_parent,
                wall_time=parent_wall_time,
            )
            lpath_submit = [path_parent.resolve()]

        if if_sbatch:
            for path_submit in lpath_submit:
                print("")
                print("================ ▶️ 提交父作业")
                print("  mode: " + mode)
                print("  chunk_parent_layout: " + chunk_parent_layout)
                print("  script: " + str(path_submit))
                # 父作业提交同样经 pei_slurm_univ_sbatch_retry 轮询重试，避免提交抖动直接漏投。
                os.system(
                    "cd " + shlex.quote(str(slurm_dir))
                    + " && " + SBATCH_RETRY_WRAPPER + " " + shlex.quote(str(path_submit))
                )

    ################################### Main control flow to here
    return None
