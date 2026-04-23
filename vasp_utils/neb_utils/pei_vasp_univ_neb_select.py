#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python

from mymetal.io.general import general_read, general_write
from pathlib import Path
import pandas as pd
import os
import sys
import argparse


parser = argparse.ArgumentParser(description='处理NEB数据文件，提取指定frame的信息')
parser.add_argument('--frame', type=int, default=None,
                    help='要提取的frame编号')

args = parser.parse_args()

if args.frame == None:
    raise ValueError("You must input this parameter --frame")

myroot = Path(os.getcwd())
path_neb_full = myroot / "p_neb_full.dat"
path_neb = myroot / "p_neb.dat"
path_nebef = myroot / "p_nebef.dat"

# select one frame or all frames
def get_neb_nebef_of_selected_frame(frame: int = None, neb_full_df: pd.DataFrame = None) -> tuple:
    dfc = neb_full_df.copy()
    df = dfc[dfc['frame'] == frame].copy()
    #print(df)

    # 1. 生成 nebdf
    nebdf = df[["id", "dist_cum(Å)", "energy_rel(eV)", "force_b(eV/Å)"]].copy()
    nebdf["image"] = nebdf["id"]

    # 调整列顺序
    nebdf = nebdf[["id", "dist_cum(Å)", "energy_rel(eV)", "force_b(eV/Å)", "image"]]

    # 2. 生成 nebefdf
    nebefdf = pd.DataFrame({
        "id": df["id"],
        "force_max(eV/Å)": 0.0,   # 或 np.nan
        "energy(eV)": df["energy_abs(eV)"],
        "energy_rel(eV)": df["energy_rel(eV)"],
    })

    print("nebdf:")
    print(nebdf)

    print("\nnebefdf:")
    print(nebefdf)

    return nebdf, nebefdf

df = general_read(path_neb_full)
nebdf, nebefdf = get_neb_nebef_of_selected_frame(frame = args.frame, neb_full_df = df)

general_write(path_neb, dfc = nebdf, if_write_col_num=True)
general_write(path_nebef, dfc = nebefdf, if_write_col_num=True)
#print(df)
