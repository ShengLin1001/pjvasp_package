#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python

import argparse
from mymetal.calculate.calqm.kpoints import get_size_by_distance
from ase.io import read


# --- 1️⃣ 解析命令行参数 ---
parser = argparse.ArgumentParser(description="Generate k-points sizes by distance.")
parser.add_argument("-i", "--input", type=str, default="CONTCAR", help="Input structure file (default: CONTCAR)")
parser.add_argument("-o", "--output", type=str, default="p_get_size_by_distance.txt", help="Output text file")
parser.add_argument("--rmin", type=int, default=20, help="Minimum distance (default: 20)")
parser.add_argument("--rmax", type=int, default=100, help="Maximum distance (default: 100)")
parser.add_argument("--rstep", type=int, default=1, help="Distance step (default: 1)")
args = parser.parse_args()

# --- 2️⃣ 读取结构文件 ---
atoms = read(args.input, format="vasp")

# --- 3️⃣ 执行计算 ---
results = []
results.append("Distance Old_Size New_Size")
print("rk Old_Size(this) New_Size")
for dis in range(args.rmin, args.rmax + 1, args.rstep):
    oldsize, newsize = get_size_by_distance(atoms, rk=dis)
    results.append(f"{dis} {oldsize} {newsize}")
    print(f"{dis} {oldsize} {newsize}")

# --- 4️⃣ 保存输出 ---
with open(args.output, 'w') as f:
    f.write("\n".join(results))
print(f"\n✅ Results saved to: {args.output}")

