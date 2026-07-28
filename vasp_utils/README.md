# vasp_utils

当前 VASP workflow 命令统一使用 `pei_*` 文件名。完整入口、目录状态和安全边界见
Sphinx 文档：

- `docs/source/manual/vasp.rst`
- `docs/source/manual/slurm.rst`
- `docs/source/reference/scripts.rst`

仓库没有注册 console scripts。先从仓库根目录检查真实入口：

```bash
bash -n vasp_utils/vasp_universal/pei_vasp_univ_sbatch
bash -n vasp_utils/vasp_workflow_bulk/pei_vasp_run_neb
python slurm_utils/slurm_universal/pei_slurm_univ_submit.py -h
```

需要裸命令时，将实际目录加入 `PATH`。Slurm 入口默认 dry-run；先检查生成的
`sub_*.sh`、资源、module、launcher 和 `-cmd`，确认后再明确添加
`-if_sbatch`。

清理、续算与 resubmit 脚本可能改写输入或输出。必须先读脚本，并在非生产目录上
验证；不要把历史 `yin_*` 文档当作当前入口。



