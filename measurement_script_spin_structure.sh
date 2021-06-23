#!/bin/bash

# *-----------------------------------*
# | shell script for spin structure   |
# |         simulations               |
# |                                   |
# |  please check the paths to the    |
# |       files before running        |
# |                                   |
# |    installation of Mayavi is      |
# |    required for the 3D plots      |
# |     (see description below)       |
# *-----------------------------------*


measurement="spin_structure"
echo $measurement

# directory for outputs
#output_dir="/Users/tobiaskohler/PhD/APB-Paper/Monte-Carlo_Simulation/New_simulations_different_sizes/D11/D11_brute_dip/10T_noAPB/"
#output_dir="/Users/tobiaskohler/PhD/thesis/Simulations/spin_structure_sims_for_thesis/D8_sphere/APB_5T_Fd3m/"

#output_dir="/Users/tobiaskohler/PhD/thesis/Simulations/spin_structure_sims_for_thesis/D11_sphere/Vacancy_ordering/noAPB_5T_P43212_0d33Fe2/"
output_dir="/Users/tobiaskohler/PhD/thesis/Simulations/spin_structure_sims_for_thesis/D6_sphere/no_APB_3T/"
#output_dir="/Users/tobiaskohler/PhD/thesis/Simulations/spin_structure_sims_for_thesis/D11_sphere/Vacancy_ordering/noAPB_5T_P43212_0d83Fe2Fe3Fe4/"
#output_dir="/Users/tobiaskohler/PhD/thesis/Simulations/spin_structure_sims_for_thesis/D11_sphere/Vacancy_ordering/noAPB_5T_Fd3m_0d89Feoct/"

#--------------------------
# particle settings
#--------------------------
outer_loop=20
# inner loop determines how many processes are started in parallel
# do not use too many for larger particles --> the system might crash!
inner_loop=1
# total number of particles to be calculated
num_particles=$(($outer_loop*$inner_loop))

# cif-file for crystal structure (only Vesta cif supported)
# if needed import the cif-file in Vesta and save it again under a different
# name. This will generate the cif-file in the right format.

#ciffile="/Users/tobiaskohler/Desktop/MCS3P/cif-files/Fd-3m_perfect.cif"
ciffile="/Users/tobiaskohler/Desktop/MCS3P/cif-files/P4_32_12_perfect.cif"

# atom labels from cif-file used for occupancies (need to be comma separated)
#atomlabels="Fe(oct),Fe(tet)"
#occupancies="0.89,1.0"
atomlabels="Fe1,Fe2,Fe3,Fe4"
#occupancies="1.0,0.33,1.0,1.0"
occupancies="1.0,0.83,0.83,0.83"

# unitcell parameters in Angstrom
#latticepars="8.35,8.35,8.35"
latticepars="8.3965,8.3965,8.3965"

# particle shape (either "Sphere" or "Cube")
shape="Sphere"

# particle size in unit cells
particle_size=6

# particle orientation angles
alpha=0.0
beta=0.0
gamma=-45.0

# particle with or without APB
APB=false

# effective anisotropy constant: K[kJ/m^3] / N_atoms/particle_volume
anisotropy_constant=3.25e-25

# exchange constants: Fe_TT, Fe_TO, Fe_OO in K
exchange_constants="-21.0,-8.6,-28.1"
#exchange_constants="-42.0,-17.2,-56.2"
# exchange contant accross APB in K
APB_constant=-106.28
#APB_constant=-212.56

# seed for particle generation
seed=20210503

#--------------------------
# field settings
#--------------------------
measurement_field=3.0

#--------------------------
# temperature settings
#--------------------------
temperature=0.01


#--------------------------
# Monte-Carlo settings
#--------------------------
# Number of Monte Carlo steps
# (statistical number of trial moves on one spin in the structure)
steps=5000

# sigma parameter for opening of gaussian cone in trial move
sigma=0.03

# setting for method of dipole interaction calculation
dipole_interactions="brute_force"
#dipole_interactions="None"
#dipole_interactions="macrocell_method"

# macrocell size only used if macrocell method is selected
macrocell_size=0.2

#----------------------------------
# generate crystal structure files
#----------------------------------
# change path to "generate_nanoparticle.py" if necessary
#
if [ "$APB" = true ] ; then
    #python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --occ_oct=1.0 --occ_tet=1.0 --diameter=$particle_size --seed=$seed --APB --template
    
    python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars --APB --template
else
    #python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --occ_oct=1.0 --occ_tet=1.0 --diameter=$particle_size --seed=$seed --template
    
    python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars --template
fi


for (( c=1; c<=$num_particles; c++ )) 
do
    if [ "$APB" = true ] ; then
        #python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --occ_oct=0.88 --occ_tet=1.0 --diameter=$particle_size --number=$c --seed=$seed --APB
        
        python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --number=$c --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars --APB
    else
        #python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --occ_oct=0.88 --occ_tet=1.0 --diameter=$particle_size --number=$c --seed=$seed
        python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --number=$c --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars
    fi
done

#--------------------------
# generate JSON files
#--------------------------
# these files are used as input parameter files for the simulation program
# change path to "generate_json_files.py" if necessary
#
if [ "$APB" == true ] ; then
    python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_json_files.py --measurement=$measurement --output=$output_dir --structure_path=$output_dir --num_particles=$num_particles --particle_size=$particle_size --alpha=$alpha --beta=$beta --gamma=$gamma --meas_field=$measurement_field --temperature=$temperature --steps=$steps --dipole=$dipole_interactions --macrocell_size=$macrocell_size --sigma=$sigma --latticepar=$latticepars --APB --anisotropy_constant=$anisotropy_constant --exchangeconstants=$exchange_constants --APB_constant=$APB_constant
else
    python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_json_files.py --measurement=$measurement --output=$output_dir --structure_path=$output_dir --num_particles=$num_particles --particle_size=$particle_size --alpha=$alpha --beta=$beta --gamma=$gamma --meas_field=$measurement_field --temperature=$temperature --steps=$steps --dipole=$dipole_interactions --macrocell_size=$macrocell_size --sigma=$sigma --latticepar=$latticepars --anisotropy_constant=$anisotropy_constant --exchangeconstants=$exchange_constants
fi

#--------------------------
# Run Simulation
#--------------------------
# change path to the compiled simulation program if necessary
#
trap "kill 0" EXIT

if [ "$APB" == true ] ; then
    /Users/tobiaskohler/Library/Developer/Xcode/DerivedData/MCS3P-fdtoextkqpexyiceuusyyotrcfyb/Build/Products/Release/MCS3P ${output_dir}/spin_structure_APB_D${particle_size}_template.json 
else
    /Users/tobiaskohler/Library/Developer/Xcode/DerivedData/MCS3P-fdtoextkqpexyiceuusyyotrcfyb/Build/Products/Release/MCS3P ${output_dir}/spin_structure_noAPB_D${particle_size}_template.json 
fi

count=0
for (( i=1; i<=$outer_loop; i++ )); do
    for ((c=1; c<=$inner_loop; c++ )); do
        ((count=count+1))
        if [ "$APB" == true ] ; then
            /Users/tobiaskohler/Library/Developer/Xcode/DerivedData/MCS3P-fdtoextkqpexyiceuusyyotrcfyb/Build/Products/Release/MCS3P ${output_dir}/${count}_spin_structure_APB_D${particle_size}.json &
        else
            /Users/tobiaskohler/Library/Developer/Xcode/DerivedData/MCS3P-fdtoextkqpexyiceuusyyotrcfyb/Build/Products/Release/MCS3P ${output_dir}/${count}_spin_structure_noAPB_D${particle_size}.json &
        fi
    done
    wait
done
wait


#---------------------------------
# Calculate average spin structure
#---------------------------------

t=${temperature%.*}
mf=${measurement_field%.*}

if [ "$APB" == true ] ; then
    template_filename=D${particle_size}_structure_APB_template_spin_structure_${t}K_${steps}MCS_${mf}T_dipNone.txt
else
    template_filename=D${particle_size}_structure_noAPB_template_spin_structure_${t}K_${steps}MCS_${mf}T_dipNone.txt
fi

echo $template_filename

python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/average_spin_structure.py --input_dir=$output_dir --template_file=$template_filename --particle_size=$particle_size --num_particles=$num_particles --dipole=$dipole_interactions

if [ "$APB" == true ] ; then
    averaged_file=D${particle_size}_structure_APB_averaged_spin_structure_${t}K_${steps}MCS_${mf}T_dip${dipole_interactions}.txt
else
    averaged_file=D${particle_size}_structure_noAPB_averaged_spin_structure_${t}K_${steps}MCS_${mf}T_dip${dipole_interactions}.txt
fi

#---------------------------------
# Generate polar plot
#---------------------------------
python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/polar_plot.py --input_dir=$output_dir --template_file=$template_filename --averaged_file=$averaged_file --particle_size=$particle_size --num_particles=$num_particles --dipole=$dipole_interactions

#---------------------------------
# Generate 3D image with Mayavi
#---------------------------------
# Mayavi has to be installed 
# use: pip install mayavi
# and: pip install PyQt5
# a conda environment with python version 3.7 is required
# create the environment with: conda create --name py37 python=3.7

source activate py37
python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/3dplot_MC_Mayavi.py --input_dir=$output_dir --file=$averaged_file
conda deactivate







