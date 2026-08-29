############################ MD heating ladder ############################
# 单个相的升温 MD：0 K 弛豫 -> 逐档升温 -> 每档采样若干帧（含每原子力与每原子能量）。
# 由 md_heating_template.in 在每个相里 include，进入本文件前必须已经：
#   - include general_init.mod（etol/ftol/maxiter/maxeval、boundary p p p）
#   - include stretch_model.mod（建好该相的盒子与原子，并定义 ${lat} 与 ${lat_hex}）
#   - replicate 到目标超胞
#   - include general_potential.mod（pair_style/pair_coeff/mass、compute peatom1/pe1）
#   - include md_params.mod（${phase} ${t_start} ${t_step} ${n_temp} ${n_equil}
#                            ${n_prod} ${n_dump} ${dt} ${t_damp} ${p_damp}）
#   - include md_seed.mod（${md_seed}，由 pei_lmp_run_md_heating 在作业里现取，不固定）
#
# 产物（落在作业目录，即 cwd）：
#   y_md_<phase>.lammpstrj    逐帧 id type x y z fx fy fz c_peatom1（sort id，定宽科学计数）
#   y_md_<phase>_temp.txt     每档一行：k  T_set  step_begin  step_end
#
# 能量**从 dump 的 c_peatom1 列逐帧求和**得到，而不是另开一个 fix print 文件再按 step 对齐：
# dump 会在 run 起点（步号恰为 n_dump 倍数时）写一帧，而 fix print 的 end_of_step 不在
# 起点触发，两者必然差一帧；把能量放进 dump 就从根上没有对齐问题。温度则按
# y_md_<phase>_temp.txt 的 step 区间归属，同样不依赖逐帧精确匹配。

# ===== 0 K 弛豫：让升温从**该势函数自己的**平衡晶格出发，而不是 md_params.mod 里的猜测值 =====
# box/relax 的自由度与 stretch_full_relax.mod 同款：立方相 iso 各向同性；
# 六方框架（${lat_hex}==1，含 HCP 与以 (111) 为底的 FCC）只把 x/y 耦合（couple xy），
# c 轴独立弛豫，保留 c/a 自由度。判据用 stretch_model.mod 定义的 ${lat_hex} 而不是
# 硬编码 lat 编号：再加一个六方相时这里不必改。
if "${lat_hex} == 1" then &
    "fix br all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0" &
else &
    "fix br all box/relax iso 0.0 fixedpoint 0 0 0"

min_style cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
unfix br

print "================ MD heating: phase ${phase}, ${n_temp} temperature(s) from ${t_start} K, step ${t_step} K"
print "MD seed = ${md_seed} (not pinned: pei_lmp_run_md_heating draws a fresh one per job run)"

# ===== 采样文件：每次跑都从头写 =====
# dump_modify append yes 是为了让**同一相的所有温度档**追加进同一个文件；
# 若上一轮留下同名残骸不先删，新一轮会接在旧帧后面，帧计数与温度归属全部错位。
shell rm -f y_md_${phase}.lammpstrj
shell rm -f y_md_${phase}_temp.txt

timestep ${dt}
reset_timestep 0

# 初速度只在最低温给一次；之后靠 NPT 恒温器把体系一档一档拉上去 —— 这才叫「升温」。
# 每档重新 velocity create 会抹掉上一档的构型历史，退化成 n_temp 个互相独立的定温模拟。
velocity all create ${t_start} ${md_seed} mom yes rot yes dist gaussian

print "k T_set step_begin step_end" append y_md_${phase}_temp.txt

variable k loop ${n_temp}
label temp_loop

variable t_now equal ${t_start}+(${k}-1)*${t_step}
print "================ ${phase}: temperature ${k}/${n_temp} = ${t_now} K"

# NPT 零压：立方相 iso；六方框架三个方向独立、x/y 耦合，让 c/a 随温度自由变化。
# 盒子是 prism（六方框架有 xy 倾斜），未被控压的倾斜量由 LAMMPS 随盒边等比缩放。
if "${lat_hex} == 1" then &
    "fix md all npt temp ${t_now} ${t_now} ${t_damp} x 0.0 0.0 ${p_damp} y 0.0 0.0 ${p_damp} z 0.0 0.0 ${p_damp} couple xy" &
else &
    "fix md all npt temp ${t_now} ${t_now} ${t_damp} iso 0.0 0.0 ${p_damp}"

# 平衡段不采样：恒温器要把体系从上一档温度拉到本档，这段构型不具代表性。
run ${n_equil}

variable step_begin equal $(step)
dump md all custom ${n_dump} y_md_${phase}.lammpstrj id type x y z fx fy fz c_peatom1
dump_modify md sort id append yes &
    format line "%8d %4d %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e"

run ${n_prod}

undump md
unfix md

# 该档采样覆盖的 step 区间；后处理据此给每一帧标上 T_set（区间归属，不做逐帧精确匹配）
print "${k} ${t_now} ${step_begin} $(step)" append y_md_${phase}_temp.txt

next k
jump SELF temp_loop
