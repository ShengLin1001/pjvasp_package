# DECISIONS — 设计选择记录

## D001: 今晚不 push 到远程
- 理由：未获明确授权，仅做本地 commit
- 时间：2026-08-30

## D002: 不运行 VASP/LAMMPS/n2p2/sbatch
- 理由：无人值守，安全第一
- 替代：dry-run / fixture / synthetic data
- 时间：2026-08-30

## D003: 优先补 API 覆盖而非视觉
- 理由：~30 个模块无 automodule 覆盖是最大 gap
- 时间：2026-08-30
