# general_mass.mod — 每个原子类型的质量（mass <type> <amu>）
# 运行期由 pei_lmp_run_properties 在各 workdir 重新写入：内容来自 model.py 的 self.lele
# （元素 -> ASE 原子质量），作为 runner 的第 4 个参数 mass_content 传入。
# 经 general_potential.mod 在每次「建/读结构」之后 include，故 stretch/Cij/gsfe 及其循环统一获得质量。
# 此处留作占位/文档；实际内容在作业里生成。
