"""Tests for the F3 guards against rewriting job scripts that a running job still executes.

2026-10-11：在同一 ``slurm/`` 里重跑提交器时，旧实现先原地截断重写全部 chunk 脚本，再由
``pei_slurm_univ_sbatch_retry`` 发现同标签父作业在队而静默跳过；在跑的 shared 父作业用
``bash <chunk 脚本>`` 从磁盘逐段读取，循环结束后从旧偏移读到新内容而报错。这里逐道核对：
A 写前查队列、B 原子替换换 inode、C 父作业去重命中按失败、E 主体整体解析。
所有 Slurm 命令都是临时目录里的假脚本，不碰真实队列。
"""

import os
import stat
import subprocess
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from mymetal.slurm.submit import (
    check_slurm_dir_idle,
    generate_slurm_script_sequential,
    generate_slurm_script_shared_parent,
    get_queue_jobs,
    pei_slurm_univ_submit,
    read_comment_tag,
    stamp_comment_tag,
    submit_parent_scripts,
    wrap_script_body,
    write_text_atomic,
)


PATH_SBATCH_RETRY = (
    Path(__file__).resolve().parents[1]
    / "slurm_utils" / "slurm_universal" / "pei_slurm_univ_sbatch_retry"
)

# 假 squeue：按 -o 格式返回测试经环境变量预置的作业表；FAKE_SQUEUE_RC 非 0 时模拟查询失败。
FAKE_SQUEUE = r"""#!/bin/bash
echo "$*" >> "${FAKE_SLURM_LOG:-/dev/null}"
if [[ "${FAKE_SQUEUE_RC:-0}" != "0" ]]; then
    echo "slurm_load_jobs error: Socket timed out on send/recv operation" >&2
    exit "$FAKE_SQUEUE_RC"
fi
for a in "$@"; do [[ "$a" == "-j" ]] && exit 0; done
fmt=""; prev=""
for a in "$@"; do [[ "$prev" == "-o" ]] && fmt="$a"; prev="$a"; done
case "$fmt" in
    "%i|%j|%Z") [[ -n "${FAKE_QUEUE_JOBS:-}" ]] && printf '%s\n' "$FAKE_QUEUE_JOBS" ;;
    "%i|%k") [[ -n "${FAKE_TAG:-}" ]] && echo "4242|${FAKE_TAG}" ;;
esac
exit 0
"""
FAKE_SBATCH = r"""#!/bin/bash
echo "$*" >> "$FAKE_SBATCH_LOG"
echo "Submitted batch job 999"
exit 0
"""
FAKE_SACCT = r"""#!/bin/bash
echo "COMPLETED|0:0"
"""
# 假提交包装器：记录 cwd 与参数；脚本名含 "fail" 时按去重命中失败退出。
FAKE_SBATCH_RETRY = r"""#!/bin/bash
echo "$PWD $*" >> "$FAKE_WRAPPER_LOG"
[[ "${@: -1}" == *fail* ]] && exit 1
exit 0
"""


def write_executable(path_save: Path, content: str) -> Path:
    """Write a fake command and make it executable."""
    path_save.write_text(content, encoding="utf-8", newline="\n")
    path_save.chmod(0o755)
    return path_save


class FakeSlurmCase(unittest.TestCase):
    """Temporary directory with fake squeue / sbatch / sacct first on PATH."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path_tmp = Path(self.tmp.name).resolve()
        self.path_bin = self.path_tmp / "fakebin"
        self.path_bin.mkdir()
        write_executable(self.path_bin / "squeue", FAKE_SQUEUE)
        write_executable(self.path_bin / "sbatch", FAKE_SBATCH)
        write_executable(self.path_bin / "sacct", FAKE_SACCT)
        write_executable(self.path_bin / "pei_slurm_univ_sbatch_retry", FAKE_SBATCH_RETRY)
        self.path_slurm_log = self.path_tmp / "squeue.log"
        self.path_sbatch_log = self.path_tmp / "sbatch.log"
        self.path_wrapper_log = self.path_tmp / "wrapper.log"
        self.dict_env = {
            "PATH": str(self.path_bin) + os.pathsep + os.environ.get("PATH", ""),
            "FAKE_SLURM_LOG": str(self.path_slurm_log),
            "FAKE_SBATCH_LOG": str(self.path_sbatch_log),
            "FAKE_WRAPPER_LOG": str(self.path_wrapper_log),
            "FAKE_QUEUE_JOBS": "",
            "FAKE_SQUEUE_RC": "0",
            "FAKE_TAG": "",
            "PEI_SUBMIT_SQUEUE_RETRY_MAX": "2",
            "PEI_SUBMIT_SQUEUE_RETRY_SLEEP": "0",
        }
        patcher = patch.dict(os.environ, self.dict_env)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.addCleanup(os.chdir, os.getcwd())

    def set_queue(self, *ljob):
        os.environ["FAKE_QUEUE_JOBS"] = "\n".join("|".join(job) for job in ljob)

    def make_job_tree(self) -> Path:
        """path_root/y_dir/{001,002} + path_root/slurm/ with an old chunk script."""
        path_root = self.path_tmp / "work"
        for name in ("001", "002"):
            (path_root / "y_dir" / name).mkdir(parents=True)
        (path_root / "slurm").mkdir()
        return path_root

    def run_dry_submit(self, path_root: Path):
        pei_slurm_univ_submit(
            path_root=path_root,
            mode="each_subdir",
            dir_root=Path("."),
            chunks=2,
            module_profile_type="none",
            launcher_type="srun",
            cmd="run-one-case",
            partition="debug",
            nodes=1,
            ncores=1,
            if_sbatch=False,
            MODULE_BLOCKS={"none": "# no modules\n"},
        )


class TestQueueGuard(FakeSlurmCase):
    """A：slurm/ 里还有作业在队时拒绝写脚本，dry run 也一样。"""

    def test_job_from_slurm_dir_is_rejected(self):
        slurm_dir = self.path_tmp / "work" / "slurm"
        self.set_queue(("43371606", "sub_slurm_each_subdir_parent.sh", str(slurm_dir)))
        with self.assertRaises(SystemExit):
            check_slurm_dir_idle(slurm_dir)

    def test_monitor_job_and_other_dirs_pass(self):
        slurm_dir = self.path_tmp / "work" / "slurm"
        self.set_queue(
            ("43371556", "sub.monitor.sh", str(slurm_dir)),
            ("43371624", "sub_slurm_each_subdir_parent.sh", str(self.path_tmp / "other" / "slurm")),
            # 计算子作业在各自 subdir 里提交，不读 slurm/ 里的脚本。
            ("43372188", "sub_slurm_univ.sh", str(slurm_dir.parent / "y_dir" / "001")),
        )
        check_slurm_dir_idle(slurm_dir)
        self.assertIn("-h -o %i|%j|%Z", self.path_slurm_log.read_text(encoding="utf-8"))

    def test_symlinked_slurm_dir_matches_resolved_workdir(self):
        path_real = self.path_tmp / "real"
        (path_real / "slurm").mkdir(parents=True)
        path_link = self.path_tmp / "link"
        path_link.symlink_to(path_real)
        # squeue 报的是解析后的真实路径（尾部斜杠也要容忍）。
        self.set_queue(("1", "sub_slurm_each_subdir_parent.sh", str(path_real / "slurm") + "/"))
        with self.assertRaises(SystemExit):
            check_slurm_dir_idle(path_link / "slurm")

    def test_squeue_failure_is_rejected_after_polled_retries(self):
        os.environ["FAKE_SQUEUE_RC"] = "1"
        with self.assertRaises(SystemExit):
            get_queue_jobs()
        # 环境变量把重试压到 2 次：失败不能被当成「队列为空」。
        self.assertEqual(len(self.path_slurm_log.read_text(encoding="utf-8").splitlines()), 2)

    def test_missing_squeue_is_rejected(self):
        path_empty = self.path_tmp / "empty_bin"
        path_empty.mkdir()
        os.environ["PATH"] = str(path_empty)
        with self.assertRaises(SystemExit):
            get_queue_jobs()

    def test_busy_dry_run_writes_nothing(self):
        path_root = self.make_job_tree()
        path_chunk = path_root / "slurm" / "sub_slurm_each_subdir_chunk001.sh"
        path_chunk.write_text("old worker\n", encoding="utf-8")
        inode_old = path_chunk.stat().st_ino
        self.set_queue(("43371606", "sub_slurm_each_subdir_parent.sh", str(path_root / "slurm")))

        with self.assertRaises(SystemExit):
            self.run_dry_submit(path_root)

        self.assertEqual(path_chunk.stat().st_ino, inode_old)
        self.assertEqual(path_chunk.read_text(encoding="utf-8"), "old worker\n")
        self.assertEqual(sorted(p.name for p in (path_root / "slurm").iterdir()),
                         ["sub_slurm_each_subdir_chunk001.sh"])
        for name in ("001", "002"):
            self.assertEqual(list((path_root / "y_dir" / name).iterdir()), [])

    def test_idle_dry_run_regenerates_with_new_inode(self):
        path_root = self.make_job_tree()
        self.set_queue(("43371556", "sub.monitor.sh", str(path_root / "slurm")))
        path_chunk = path_root / "slurm" / "sub_slurm_each_subdir_chunk001.sh"
        path_chunk2 = path_root / "slurm" / "sub_slurm_each_subdir_chunk002.sh"
        path_parent = path_root / "slurm" / "sub_slurm_each_subdir_parent.sh"

        self.run_dry_submit(path_root)
        inode_first = path_chunk.stat().st_ino
        tag_first = read_comment_tag(path_chunk)
        self.run_dry_submit(path_root)

        self.assertNotEqual(path_chunk.stat().st_ino, inode_first)
        self.assertEqual(read_comment_tag(path_chunk), tag_first)
        for path_script in (path_chunk, path_chunk2, path_parent):
            line = path_script.read_text(encoding="utf-8")
            self.assertIn("\n{\n", line)
            self.assertTrue(line.endswith("\n}; exit\n"))
            completed = subprocess.run(["bash", "-n", str(path_script)],
                                       capture_output=True, text=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(sorted(p.name for p in (path_root / "slurm").iterdir()),
                         ["sub_slurm_each_subdir_chunk001.sh", "sub_slurm_each_subdir_chunk002.sh",
                          "sub_slurm_each_subdir_parent.sh"])


class TestAtomicWrite(unittest.TestCase):
    """B：同目录临时文件 + os.replace，换 inode、保留权限、不留临时文件。"""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path_tmp = Path(self.tmp.name)

    def test_replaces_inode_and_keeps_mode(self):
        path_script = self.path_tmp / "chunk001.sh"
        path_script.write_text("old\n", encoding="utf-8")
        path_script.chmod(0o755)
        inode_old = path_script.stat().st_ino

        write_text_atomic(path_script, "#!/bin/bash\necho 新\n")

        self.assertNotEqual(path_script.stat().st_ino, inode_old)
        self.assertEqual(stat.S_IMODE(path_script.stat().st_mode), 0o755)
        self.assertEqual(path_script.read_bytes(), "#!/bin/bash\necho 新\n".encode("utf-8"))
        self.assertEqual([p.name for p in self.path_tmp.iterdir()], ["chunk001.sh"])

    def test_open_reader_keeps_old_content(self):
        # 在跑的 bash 持有旧 inode 的描述符：替换后它读到的仍是旧内容。
        path_script = self.path_tmp / "chunk001.sh"
        path_script.write_text("old content\n", encoding="utf-8")
        with open(path_script, encoding="utf-8") as f:
            write_text_atomic(path_script, "new content that is longer\n")
            self.assertEqual(f.read(), "old content\n")
        self.assertEqual(path_script.read_text(encoding="utf-8"), "new content that is longer\n")

    def test_new_file_follows_umask(self):
        path_new = self.path_tmp / "new.sh"
        umask_old = os.umask(0o022)
        try:
            write_text_atomic(path_new, "x\n")
        finally:
            os.umask(umask_old)
        self.assertEqual(stat.S_IMODE(path_new.stat().st_mode), 0o644)

    def test_failed_replace_keeps_target_and_removes_only_temp(self):
        path_script = self.path_tmp / "chunk001.sh"
        path_script.write_text("old\n", encoding="utf-8")
        with patch("mymetal.slurm.submit.os.replace", side_effect=OSError("boom")):
            with self.assertRaises(OSError):
                write_text_atomic(path_script, "new\n")
        self.assertEqual(path_script.read_text(encoding="utf-8"), "old\n")
        self.assertEqual([p.name for p in self.path_tmp.iterdir()], ["chunk001.sh"])

    def test_stamp_comment_tag_replaces_inode(self):
        path_script = self.path_tmp / "sub.sh"
        path_script.write_text("#!/bin/bash\n#SBATCH -p debug\necho hi\n", encoding="utf-8")
        path_script.chmod(0o755)
        inode_old = path_script.stat().st_ino

        tag = stamp_comment_tag(path_script)

        self.assertNotEqual(path_script.stat().st_ino, inode_old)
        self.assertEqual(stat.S_IMODE(path_script.stat().st_mode), 0o755)
        self.assertEqual(read_comment_tag(path_script), tag)


class TestRunningScriptRewrite(unittest.TestCase):
    """B + E 的行为核对：bash 执行到一半时改写脚本文件。

    原地截断重写 -> 在跑的 bash 从旧偏移读到新内容（复现 A5 事故）；
    原子替换或主体包进 ``{ ... }; exit`` -> 照原脚本跑完。
    """

    BODY = (
        'echo start\n'
        'while [ ! -e "$GO" ]; do sleep 0.05; done\n'
        'echo after\n'
    )

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path_tmp = Path(self.tmp.name)

    def run_and_rewrite(self, line_script: str, rewrite) -> subprocess.CompletedProcess:
        path_script = self.path_tmp / "worker.sh"
        path_go = self.path_tmp / "go"
        path_script.write_text(line_script, encoding="utf-8", newline="\n")
        # 新内容：前面全是注释、恰好在旧脚本循环结束的字节偏移处接一条 echo HIJACKED。
        offset = len(line_script.encode("utf-8").split(b"done\n", 1)[0]) + len(b"done\n")
        line_new = "#!/bin/bash\n" + "#" * (offset - len("#!/bin/bash\n") - 1) + "\n" + "echo HIJACKED\n"
        self.assertEqual(len(line_new.encode("utf-8").split(b"echo HIJACKED", 1)[0]), offset)

        process = subprocess.Popen(["bash", str(path_script)], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True,
                                   env={**os.environ, "GO": str(path_go)})
        try:
            self.assertEqual(process.stdout.readline(), "start\n")
            time.sleep(0.2)
            rewrite(path_script, line_new)
            path_go.write_text("", encoding="utf-8")
            stdout, stderr = process.communicate(timeout=30)
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate()
        return subprocess.CompletedProcess(process.args, process.returncode, "start\n" + stdout, stderr)

    @staticmethod
    def rewrite_in_place(path_script: Path, line_new: str):
        path_script.write_text(line_new, encoding="utf-8", newline="\n")

    def test_in_place_rewrite_hijacks_running_bash(self):
        completed = self.run_and_rewrite("#!/bin/bash\n" + self.BODY, self.rewrite_in_place)
        self.assertIn("HIJACKED", completed.stdout)
        self.assertNotIn("after", completed.stdout)

    def test_atomic_replace_keeps_running_bash_on_old_script(self):
        completed = self.run_and_rewrite("#!/bin/bash\n" + self.BODY, write_text_atomic)
        self.assertEqual(completed.stdout, "start\nafter\n")
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_wrapped_body_survives_in_place_rewrite(self):
        completed = self.run_and_rewrite("#!/bin/bash\n" + wrap_script_body(self.BODY),
                                         self.rewrite_in_place)
        self.assertEqual(completed.stdout, "start\nafter\n")
        self.assertEqual(completed.returncode, 0, completed.stderr)


class TestGeneratedScriptsWrapped(unittest.TestCase):
    """E：worker 与 shared parent 的主体都包进 ``{ ... }; exit`` 且语法正确。"""

    def test_worker_and_parent_bodies_are_wrapped(self):
        dict_args = {
            "partition": "debug",
            "module_profile_type": "none",
            "MODULE_BLOCKS": {"none": "# no modules\n"},
            "if_output": False,
        }
        line_worker = generate_slurm_script_sequential(
            **dict_args, nodes=1, ncores=1, launcher="srun", cmd="run-one-case",
            if_use_my_launcher=False, group=[Path("/jobs/001")], mode="each_subdir",
        )
        line_parent = generate_slurm_script_shared_parent(
            **dict_args, lpath_worker=[Path("/work/slurm/sub_slurm_each_subdir_chunk001.sh")],
        )
        for line in (line_worker, line_parent):
            self.assertIn("\n{\n", line)
            self.assertTrue(line.endswith("\n}; exit\n"))
            # #SBATCH 头必须留在包裹之外，sbatch 才认得。
            self.assertLess(line.index("#SBATCH"), line.index("\n{\n"))
            completed = subprocess.run(["bash", "-n"], input=line, text=True,
                                       capture_output=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)


class TestDedupHitFail(FakeSlurmCase):
    """C：顶层父作业提交命中去重标签按失败；worker 内（不带开关）仍按对账成功。"""

    TAG = "1234567890" * 10

    def setUp(self):
        super().setUp()
        self.path_case = self.path_tmp / "case"
        self.path_case.mkdir()
        (self.path_case / "sub_test.sh").write_text(
            "#!/bin/bash\n#SBATCH -p debug\n#SBATCH --comment=" + self.TAG + "\necho hi\n",
            encoding="utf-8",
        )
        for key in ("PEI_SBATCH_RETRY_SLEEP", "PEI_SBATCH_THROTTLE_SLEEP", "PEI_SBATCH_RECONCILE_SLEEP",
                    "PEI_SBATCH_WAIT_SLEEP", "PEI_SBATCH_SACCT_SLEEP"):
            os.environ[key] = "0"

    def run_retry(self, *args, if_hit: bool) -> subprocess.CompletedProcess:
        os.environ["FAKE_TAG"] = self.TAG if if_hit else ""
        return subprocess.run(["bash", str(PATH_SBATCH_RETRY), *args], cwd=self.path_case,
                              capture_output=True, text=True, timeout=120, check=False)

    def sbatch_calls(self) -> list[str]:
        if not self.path_sbatch_log.exists():
            return []
        return self.path_sbatch_log.read_text(encoding="utf-8").splitlines()

    def test_flag_and_hit_fails_without_submitting(self):
        completed = self.run_retry("-dedup_hit_fail", "sub_test.sh", if_hit=True)
        self.assertEqual(completed.returncode, 1)
        self.assertIn("-dedup_hit_fail", completed.stderr)
        self.assertNotIn("Submitted batch job", completed.stdout)
        self.assertEqual(self.sbatch_calls(), [])

    def test_flag_without_hit_submits_and_strips_flag(self):
        completed = self.run_retry("-dedup_hit_fail", "sub_test.sh", if_hit=False)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Submitted batch job 999", completed.stdout)
        self.assertEqual(self.sbatch_calls(), ["sub_test.sh"])

    def test_worker_wait_hit_still_reconciles(self):
        completed = self.run_retry("--wait", "sub_test.sh", if_hit=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Submitted batch job 4242", completed.stdout)
        self.assertEqual(self.sbatch_calls(), [])

    def test_hit_without_flag_reports_existing_job(self):
        completed = self.run_retry("sub_test.sh", if_hit=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Submitted batch job 4242", completed.stdout)
        self.assertEqual(self.sbatch_calls(), [])

    def test_submit_parent_scripts_attempts_all_then_fails(self):
        slurm_dir = self.path_tmp / "work" / "slurm"
        slurm_dir.mkdir(parents=True)
        lpath = [slurm_dir / "sub_chunk001_fail.sh", slurm_dir / "sub_chunk002.sh"]
        with self.assertRaises(SystemExit):
            submit_parent_scripts(lpath, slurm_dir)
        lcall = self.path_wrapper_log.read_text(encoding="utf-8").splitlines()
        self.assertEqual(lcall, [str(slurm_dir) + " -dedup_hit_fail " + str(path) for path in lpath])

    def test_submit_parent_scripts_returns_submitted(self):
        slurm_dir = self.path_tmp / "work" / "slurm"
        slurm_dir.mkdir(parents=True)
        lpath = [slurm_dir / "sub_slurm_each_subdir_parent.sh"]
        self.assertEqual(submit_parent_scripts(lpath, slurm_dir), lpath)


if __name__ == "__main__":
    unittest.main()
