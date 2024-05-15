from mpi4py import MPI
import numpy as np
from time import time
import os

def parallel_sum():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # 创建一个随机数组
    np.random.seed(42 + rank)
    data = np.random.rand(1000000)

    # 使用 OpenMP 获取并行部分的运行时间
    start_time = time()

    # 每个进程计算部分和
    local_sum = np.sum(data)

    # 使用 MPI 合并计算结果
    total_sum = comm.reduce(local_sum, op=MPI.SUM, root=0)

    # 计时结束
    end_time = time()
    if rank == 0:
        print(f"Total sum is {total_sum}, computed in {end_time - start_time} seconds.")

# 运行函数
if __name__ == "__main__":
    parallel_sum()
