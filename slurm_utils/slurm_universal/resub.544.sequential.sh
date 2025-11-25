#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128         
source /public3/soft/modules/module.sh
module load mpi/intel/17.0.7-thc
export PATH=/public3/home/scg6928/mysoft/vasp/vasp/544-yin/vasp.5.4.4.pl2/bin:$PATH

echo "🚀 ==> Resubmit unfinished jobs under y_dir sequentially"

for dir in ./y_dir/*/; do
    # ||, or operator to skip if cd fails
    cd "$dir" || continue  
    echo "================ 📁 $dir"
    echo "📍 Current dir: $(pwd)"
    

    tag_resub=F
    # Check if OUTCAR exists
    if [[ -f "OUTCAR" ]]; then
        # if "reached required accuracy" is found in OUTCAR, skip resubmission
        if grep -q "reached required accuracy" OUTCAR; then
            echo "✅ Calculation completed successfully, skipping..."
        # if not found, resubmit the job
        else
            echo "❌ Calculation not completed, resubmitting..."
            tag_resub=T
        fi
    else
        # if OUTCAR does not exist, resubmit the job
        echo "🔍 No OUTCAR found, resubmitting job..."
        tag_resub=T
    fi

    if [[ $tag_resub == "T" ]]; then
        #pei_vasp_univ_find_and_change -ediff 1e-10
        cp CONTCAR POSCAR
        srun vasp_std
    fi
    
    # Now, we have the OUTCAR, and do double check
    if grep -q "reached required accuracy" OUTCAR; then
        echo "✅ Calculation completed successfully after resubmission."
    else
        echo "❌ Calculation still not completed after resubmission."
    fi

    cd - > /dev/null
    echo ""  # 空行分隔不同任务输出
done

echo "🎉 Done!"

