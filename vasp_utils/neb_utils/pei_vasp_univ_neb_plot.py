#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python

from mymetal.io.general import general_read
from pathlib import Path
import os
from mymetal.universal.plot.workflow import my_plot_neb_full
import argparse

parser = argparse.ArgumentParser(description='绘制NEB完整数据图，支持切片选择特定frames')
parser.add_argument('--start', type=int, default=None,
                    help='切片起始位置（包含），支持负数，默认None表示从开头')
parser.add_argument('--end', type=int, default=None,
                    help='切片结束位置（不包含），支持负数，默认None表示到结尾')

args = parser.parse_args()

myroot = Path(os.getcwd())
path_neb_full = myroot / "p_neb_full.dat"
path_neb = myroot / "p_neb.dat"
path_nebef = myroot / "p_nebef.dat"
savefig_path = myroot / 'p_post_neb_full_selected.pdf'

neb_full_df = general_read(path_neb_full)

my_plot_neb_full(neb_full_df, savefig_path=savefig_path, slice_start=args.start, slice_end=args.end)

