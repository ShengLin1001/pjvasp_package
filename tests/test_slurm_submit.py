"""Tests for shared and historical per-chunk Slurm parent layouts."""

import os
import subprocess
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from mymetal.slurm.submit import (
    check_chunk_parent_layout,
    generate_slurm_script_base,
    generate_slurm_script_sequential,
    generate_slurm_script_shared_parent,
    get_lsubdir,
    get_y_dir_lsubdir,
    pei_slurm_univ_submit,
)


PATH_DISCOVERY_FIXTURE = (
    Path(__file__).resolve().parent / "fixtures" / "slurm_submit"
).resolve()


class TestJobDirectoryDiscovery(unittest.TestCase):
    """Verify recursive y_dir discovery independently of directory depth."""

    def test_recursively_collects_immediate_children_of_every_y_dir(self):
        lsubdir = get_y_dir_lsubdir(PATH_DISCOVERY_FIXTURE)

        lexpected = sorted([
            PATH_DISCOVERY_FIXTURE / "y_dir" / "root_case",
            PATH_DISCOVERY_FIXTURE / "y_dir" / "shared_case",
            PATH_DISCOVERY_FIXTURE / "mode_a" / "y_dir" / "mode_case",
            PATH_DISCOVERY_FIXTURE / "mode_a" / "y_dir" / "shared_case",
            PATH_DISCOVERY_FIXTURE / "deep" / "mode_b" / "y_dir" / "deep_case",
        ])
        self.assertEqual(lsubdir, lexpected)
        self.assertNotIn(
            PATH_DISCOVERY_FIXTURE / "deep" / "not_y_dir" / "ignored_case",
            lsubdir,
        )

    def test_y_dir_itself_can_be_the_explicit_search_root(self):
        path_y_dir = PATH_DISCOVERY_FIXTURE / "y_dir"

        self.assertEqual(
            get_lsubdir(None, path_y_dir),
            [path_y_dir / "root_case", path_y_dir / "shared_case"],
        )

    def test_explicit_lsubdir_filters_all_discovered_y_dirs(self):
        self.assertEqual(
            get_lsubdir(["shared_case"], PATH_DISCOVERY_FIXTURE),
            sorted([
                PATH_DISCOVERY_FIXTURE / "y_dir" / "shared_case",
                PATH_DISCOVERY_FIXTURE / "mode_a" / "y_dir" / "shared_case",
            ]),
        )

    def test_search_root_without_y_dir_is_rejected(self):
        with self.assertRaises(SystemExit):
            get_y_dir_lsubdir(PATH_DISCOVERY_FIXTURE / "deep" / "not_y_dir")

    def test_parent_scripts_keep_recursive_job_paths_absolute(self):
        lsubdir = get_y_dir_lsubdir(PATH_DISCOVERY_FIXTURE)
        dict_args = {
            "partition": "amd_512",
            "nodes": 1,
            "ncores": 1,
            "module_profile_type": "none",
            "MODULE_BLOCKS": {"none": "# no modules\n"},
            "launcher": "srun",
            "cmd": "run-one-case",
            "if_use_my_launcher": False,
            "group": lsubdir,
            "if_output": False,
        }

        for mode in ("each-subdir", "single-alloc"):
            with self.subTest(mode=mode):
                line = generate_slurm_script_sequential(**dict_args, mode=mode)
                self.assertIn('cd "$subdir"', line)
                for path_subdir in lsubdir:
                    self.assertIn(str(path_subdir), line)


class TestChunkParentLayout(unittest.TestCase):
    """Verify that chunks remain lanes while parent-job topology is selectable."""

    def run_submit(
            self,
            mode: str,
            chunks: int,
            chunk_parent_layout: str = "auto",
            child_wall_time: str = None,
            parent_wall_time: str = None):
        """Run the orchestration branch with filesystem and Slurm calls mocked."""
        lsubdir = [Path("/jobs/" + name) for name in ("001", "002", "003")]
        with ExitStack() as stack:
            mock_chdir = stack.enter_context(patch("mymetal.slurm.submit.os.chdir"))
            mock_mkdir = stack.enter_context(patch.object(Path, "mkdir"))
            mock_system = stack.enter_context(patch("mymetal.slurm.submit.os.system"))
            mock_get_lsubdir = stack.enter_context(
                patch("mymetal.slurm.submit.get_lsubdir", return_value=lsubdir)
            )
            mock_base = stack.enter_context(
                patch("mymetal.slurm.submit.generate_slurm_script_base")
            )
            mock_worker = stack.enter_context(
                patch("mymetal.slurm.submit.generate_slurm_script_sequential")
            )
            mock_shared = stack.enter_context(
                patch("mymetal.slurm.submit.generate_slurm_script_shared_parent")
            )

            pei_slurm_univ_submit(
                path_root=Path("/work"),
                mode=mode,
                dir_root=Path("./y_dir"),
                chunks=chunks,
                chunk_parent_layout=chunk_parent_layout,
                module_profile_type="none",
                launcher_type="srun",
                cmd="run-one-case",
                partition="amd_512",
                nodes=2,
                ncores=16,
                if_sbatch=True,
                child_wall_time=child_wall_time,
                parent_wall_time=parent_wall_time,
                MODULE_BLOCKS={"none": "# no modules\n"},
            )

        return {
            "chdir": mock_chdir,
            "mkdir": mock_mkdir,
            "system": mock_system,
            "get_lsubdir": mock_get_lsubdir,
            "base": mock_base,
            "worker": mock_worker,
            "shared": mock_shared,
        }

    def test_auto_layout_is_mode_sensitive(self):
        self.assertEqual(check_chunk_parent_layout("auto", "each-subdir"), "shared")
        self.assertEqual(check_chunk_parent_layout("auto", "single-alloc"), "per-chunk")
        self.assertEqual(check_chunk_parent_layout("auto", "parallel"), "per-chunk")

    def test_shared_layout_is_rejected_for_single_alloc(self):
        with self.assertRaises(SystemExit):
            check_chunk_parent_layout("shared", "single-alloc")

    def test_each_subdir_auto_submits_one_shared_parent(self):
        dict_mock = self.run_submit("each-subdir", chunks=3)

        self.assertEqual(dict_mock["base"].call_count, 3)
        self.assertEqual(dict_mock["worker"].call_count, 3)
        self.assertEqual(dict_mock["shared"].call_count, 1)
        self.assertEqual(dict_mock["system"].call_count, 1)

        # 子作业保留用户资源；只负责 sbatch --wait 的 chunk worker 修正为 1/1。
        self.assertEqual(dict_mock["base"].call_args_list[0].args[1:3], (2, 16))
        self.assertEqual(dict_mock["worker"].call_args_list[0].args[1:3], (1, 1))

    def test_each_subdir_per_chunk_keeps_historical_parent_count(self):
        dict_mock = self.run_submit(
            "each-subdir", chunks=3, chunk_parent_layout="per-chunk"
        )

        self.assertEqual(dict_mock["worker"].call_count, 3)
        self.assertEqual(dict_mock["shared"].call_count, 0)
        self.assertEqual(dict_mock["system"].call_count, 3)

    def test_single_alloc_chunks_remain_separate_parent_jobs(self):
        dict_mock = self.run_submit("single-alloc", chunks=3)

        self.assertEqual(dict_mock["base"].call_count, 0)
        self.assertEqual(dict_mock["worker"].call_count, 3)
        self.assertEqual(dict_mock["shared"].call_count, 0)
        self.assertEqual(dict_mock["system"].call_count, 3)
        self.assertEqual(dict_mock["worker"].call_args_list[0].args[1:3], (2, 16))

    def test_chunk_scripts_stay_in_path_root_slurm_directory(self):
        for mode in ("each-subdir", "single-alloc"):
            with self.subTest(mode=mode):
                dict_mock = self.run_submit(mode, chunks=3)
                for call_worker in dict_mock["worker"].call_args_list:
                    self.assertEqual(
                        call_worker.kwargs["path_save"].parent,
                        Path("/work/slurm"),
                    )

    def test_child_wall_time_only_reaches_calculation_child_jobs(self):
        for mode in ("parallel", "each-subdir"):
            with self.subTest(mode=mode):
                dict_mock = self.run_submit(
                    mode, chunks=1, child_wall_time="2-00:00:00"
                )
                self.assertEqual(dict_mock["base"].call_count, 3)
                self.assertEqual(
                    dict_mock["base"].call_args_list[0].kwargs["wall_time"],
                    "2-00:00:00",
                )
                if mode == "each-subdir":
                    self.assertIsNone(
                        dict_mock["worker"].call_args.kwargs["wall_time"]
                    )

        dict_mock = self.run_submit(
            "single-alloc", chunks=1, child_wall_time="2-00:00:00"
        )
        self.assertEqual(dict_mock["base"].call_count, 0)

    def test_base_script_wall_time_is_optional_and_validated(self):
        dict_args = {
            "partition": "amd_512",
            "nodes": 1,
            "ncores": 16,
            "module_profile_type": "none",
            "MODULE_BLOCKS": {"none": "# no modules\n"},
            "launcher": "srun",
            "cmd": "run-one-case",
            "if_use_my_launcher": False,
            "if_output": False,
        }

        line_default = generate_slurm_script_base(**dict_args)
        line_limited = generate_slurm_script_base(
            **dict_args, wall_time="2-00:00:00"
        )

        self.assertNotIn("#SBATCH --time=", line_default)
        self.assertIn("#SBATCH --time=2-00:00:00\n", line_limited)
        with self.assertRaises(SystemExit):
            generate_slurm_script_base(
                **dict_args, wall_time="01:00:00\n#SBATCH -p injected"
            )

    def test_parent_wall_time_only_reaches_parent_scripts(self):
        dict_mock = self.run_submit(
            "each-subdir", chunks=3, parent_wall_time="7-00:00:00"
        )
        self.assertEqual(
            dict_mock["worker"].call_args_list[0].kwargs["wall_time"],
            "7-00:00:00",
        )
        self.assertEqual(
            dict_mock["shared"].call_args.kwargs["wall_time"],
            "7-00:00:00",
        )
        self.assertIsNone(dict_mock["base"].call_args_list[0].kwargs["wall_time"])

        dict_mock = self.run_submit(
            "single-alloc", chunks=1, parent_wall_time="7-00:00:00"
        )
        self.assertEqual(
            dict_mock["worker"].call_args.kwargs["wall_time"],
            "7-00:00:00",
        )

        dict_mock = self.run_submit(
            "parallel", chunks=1, parent_wall_time="7-00:00:00"
        )
        self.assertEqual(dict_mock["worker"].call_count, 0)
        self.assertEqual(dict_mock["shared"].call_count, 0)

    def test_parent_scripts_emit_parent_wall_time(self):
        dict_args = {
            "partition": "amd_512",
            "module_profile_type": "none",
            "MODULE_BLOCKS": {"none": "# no modules\n"},
            "wall_time": "7-00:00:00",
            "if_output": False,
        }
        line_worker = generate_slurm_script_sequential(
            **dict_args,
            nodes=1,
            ncores=1,
            launcher="srun",
            cmd="run-one-case",
            if_use_my_launcher=False,
            group=[Path("/jobs/001")],
            mode="each-subdir",
        )
        line_shared = generate_slurm_script_shared_parent(
            **dict_args,
            lpath_worker=[Path("/work/slurm/chunk001.sh")],
        )

        self.assertIn("#SBATCH --time=7-00:00:00\n", line_worker)
        self.assertIn("#SBATCH --time=7-00:00:00\n", line_shared)

    def test_generated_shared_parent_is_valid_bash_and_quotes_paths(self):
        line = generate_slurm_script_shared_parent(
            partition="amd_512",
            module_profile_type="none",
            MODULE_BLOCKS={"none": "# no modules\n"},
            lpath_worker=[
                Path("/work/slurm/chunk 001.sh"),
                Path("/work/slurm/chunk002.sh"),
            ],
            if_output=False,
            path_save=Path("/work/slurm/sub_slurm_each_subdir_parent.sh"),
        )

        completed = subprocess.run(
            ["bash", "-n"], input=line, text=True, capture_output=True, check=False
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("'/work/slurm/chunk 001.sh'", line)
        self.assertIn('bash "$path_worker"', line)
        self.assertIn('if wait "$pid"', line)
        self.assertIn("worker_failed", line)
        self.assertNotIn("#SBATCH --time=", line)

    def test_shared_parent_waits_all_workers_and_aggregates_failure(self):
        job_tag = "codex-shared-parent-" + str(os.getpid())
        path_log_ok = Path("/tmp/slurm-" + job_tag + "-null.out")
        path_log_failed = Path(
            "/tmp/slurm-" + job_tag + "-codex-worker-does-not-exist.out"
        )
        line = generate_slurm_script_shared_parent(
            partition="amd_512",
            module_profile_type="none",
            MODULE_BLOCKS={"none": "# no modules\n"},
            lpath_worker=[Path("/dev/null"), Path("/codex-worker-does-not-exist")],
            if_output=False,
            path_save=Path("/tmp/sub_slurm_each_subdir_parent.sh"),
        )

        try:
            dict_env = {**os.environ, "SLURM_JOB_ID": job_tag}
            completed = subprocess.run(
                ["bash"], input=line, text=True, capture_output=True,
                env=dict_env, check=False,
            )
            self.assertEqual(completed.returncode, 1)
            self.assertIn("worker_success=1    worker_failed=1", completed.stdout)
            self.assertTrue(path_log_ok.is_file())
            self.assertTrue(path_log_failed.is_file())
        finally:
            # 两个日志路径均由本测试完整解析；逐文件清理，不做递归或通配删除。
            if path_log_ok.is_file():
                path_log_ok.unlink()
            if path_log_failed.is_file():
                path_log_failed.unlink()


if __name__ == "__main__":
    unittest.main()
