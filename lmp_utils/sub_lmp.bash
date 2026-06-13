#!/bin/bash

myroot=$(pwd)
pair_style="eam"
pair_coeff="${myroot}/potential/Au_u3.eam"
python_path="/home/louis/mysoft/env/pyenv/mydft/bin/python"
lmp_template_path="/home/louis/mywork/github/myrepo/PJ-advanced-reasearch/20250505_machine_learning/workflow_for_n2p2_test/template"

rm -rf y_stretch
rm -rf y_Cij_energy
rm -rf y_gsfe

############ Found equilibrium lattice constant ############
for lat in fcc bcc hcp; do
    workdir="${myroot}/y_stretch/${lat}/y_dir"
    dumpdir="${myroot}/y_stretch/${lat}/dump"
    mkdir -p ${workdir}
    cd ${workdir}

    for file in \
        "${lmp_template_path}/stretch.mod" \
        "${lmp_template_path}/stretch_template.in" \
        "${lmp_template_path}/stretch_model.mod" \
        "${lmp_template_path}/stretch_full_relax.mod"\
        "${lmp_template_path}/stretch_constrained_relax.mod"\
        "${lmp_template_path}/general_init.mod" \
        "${lmp_template_path}/general_potential.mod" \
        "${lmp_template_path}/general_output.mod" \
        "${lmp_template_path}/general_structural_info.mod";
    do
        cp "${file}" .
    done

    # bcc:3, fcc:2, hcp:1
    case "${lat}" in
        bcc)
            latnum=3
            ini_aa=3.20
            ;;
        fcc)
            latnum=2
            ini_aa=4.08
            ;;
        hcp)
            latnum=1
            ini_aa=2.88
            ;;
        *)
            echo "Unknown lattice: ${lat}"
            exit 1
            ;;
    esac

    sed -i "s/lat_template/${latnum}/g" stretch_model.mod
    sed -i "s/aa_template/${ini_aa}/g"  stretch_model.mod
    sed -i "s/pair_style_template/${pair_style}/g" general_potential.mod
    sed -i "s|pair_coeff_template|${pair_coeff}|g" general_potential.mod
    lmp -in stretch_template.in


    cd $myroot

done

${python_path} post/stretch.py


############ Cij_energy ############
for lat in fcc bcc hcp; do
    workdir="${myroot}/y_Cij_energy/${lat}/y_dir"
    dumpdir="${myroot}/y_Cij_energy/${lat}/dump"
    mkdir -p ${workdir}
    cd ${workdir}

    # read equilibrium lattice constant
    stretch_output_dir="${myroot}/y_stretch/${lat}"
    a0=`sed -n '21,21p' $stretch_output_dir/p_post_stretch.txt | awk '{printf "%.16f", $1}'`
    c0=`sed -n '21,21p' $stretch_output_dir/p_post_stretch.txt | awk '{printf "%.16f", $3}'`

    for file in \
        "${lmp_template_path}/stretch_model.mod" \
        "${lmp_template_path}/stretch_constrained_relax.mod"\
        "${lmp_template_path}/general_init.mod" \
        "${lmp_template_path}/general_potential.mod" \
        "${lmp_template_path}/general_structural_info.mod" \
        "${lmp_template_path}/Cij_energy_template.in"\
        "${lmp_template_path}/Cij_energy.mod" \
        "${lmp_template_path}/Cij_energy_output.mod";
    do
        cp "${file}" .
    done


    # bcc:3, fcc:2, hcp:1
    case "${lat}" in
        bcc) latnum=3 ;;
        fcc) latnum=2 ;;
        hcp) latnum=1 ;;
        *) echo "Unknown lattice: ${lat}" ; exit 1 ;;
    esac

    sed -i "s/lat_template/${latnum}/g" stretch_model.mod
    sed -i "s/aa_template/${a0}/g"      stretch_model.mod
    sed -i "s/pair_style_template/${pair_style}/g" general_potential.mod
    sed -i "s|pair_coeff_template|${pair_coeff}|g" general_potential.mod
    lmp -in Cij_energy_template.in

    cd $myroot

done

${python_path} post/Cij_energy.py

############ gsfe ############

# HCP: HCP_basal, HCP_prism1w, HCP_pyr1w, HCP_pyr2
# FCC: FCC_111, FCC_100

for lat in fcc hcp; do

    # bcc:3, fcc:2, hcp:1
    case "${lat}" in
        bcc)
            latnum=3
            ;;
        fcc)
            latnum=2
            lgsfe_type=("FCC_111" "FCC_100")
            ;;
        hcp)
            latnum=1
            lgsfe_type=("HCP_basal" "HCP_prism1w" "HCP_pyr1w" "HCP_pyr2")
            ;;
        *)
            echo "Unknown lattice: ${lat}"
            exit 1
            ;;
    esac

    # read equilibrium lattice constant
    stretch_output_dir="${myroot}/y_stretch/${lat}"
    a0=`sed -n '21,21p' $stretch_output_dir/p_post_stretch.txt | awk '{printf "%.16f", $1}'`
    c0=`sed -n '21,21p' $stretch_output_dir/p_post_stretch.txt | awk '{printf "%.16f", $3}'`

    for gsfe_type in "${lgsfe_type[@]}"; do
        workdir="${myroot}/y_gsfe/${lat}/${gsfe_type}/y_dir"
        dumpdir="${myroot}/y_gsfe/${lat}/${gsfe_type}/dump"
        mkdir -p ${workdir}
        cd ${workdir}

        for file in \
            "${lmp_template_path}/gsfe_model.py" \
            "${lmp_template_path}/gsfe_template.in" \
            "${lmp_template_path}/gsfe.mod"\
            "${lmp_template_path}/general_init.mod" \
            "${lmp_template_path}/general_potential.mod" \
            "${lmp_template_path}/general_output.mod" ;
        do
            cp "${file}" .
        done

        sed -i "s|gsfe_type_template|${gsfe_type}|g"     gsfe_model.py
        sed -i "s/aa_template/${a0}/g"                   gsfe_model.py
        sed -i "s/cc_template/${c0}/g"                   gsfe_model.py
        sed -i "s|pair_style_template|${pair_style}|g" general_potential.mod
        sed -i "s|pair_coeff_template|${pair_coeff}|g" general_potential.mod
        sed -i "s|python_path_template|${python_path}|g" gsfe_template.in
        sed -i "s|gsfe_type_template|${gsfe_type}|g"     gsfe_template.in
        lmp -in gsfe_template.in

        cd $myroot
    done

done

${python_path} post/gsfe.py