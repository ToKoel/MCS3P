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

SCRIPT_DIR=$(pwd)

# directories and filepaths
output_dir="/absolute/path/to/outfolder/"

PROGRAM_FILES_DIR="/absolute/path/to/program_files"
CIF_DIR="/absolute/path/to/cif/files"

generate_nanoparticle=$PROGRAM_FILES_DIR"/MCS3P/generate_nanoparticle.py"
generate_json_files=$PROGRAM_FILES_DIR"/MCS3P/generate_json_files.py"
average_spin_structure=$PROGRAM_FILES_DIR"/MCS3P/average_spin_structure.py"
generate_polar_plot=$PROGRAM_FILES_DIR"/MCS3P/polar_plot.py"
plot_3d=$PROGRAM_FILES_DIR"/MCS3P/3dplot_MC_Mayavi.py"
MCS3P=$PROGRAM_FILES_DIR"/bin/MCS3P/Build/Products/Debug/MCS3P"

# setup folder structure
cd $output_dir

date=$(date '+%Y-%m-%d')
OUT_DIR="experiment-$date-$RANDOM"

if [ -d "$OUT_DIR" ]; then
    echo "$OUT_DIR already exists."
    exit 1
else
    mkdir $OUT_DIR
fi

cd $OUT_DIR

STRUCT_DIR="crystal_structure_files"
SPIN_DIR="spin_structure_files"
JSON_DIR="json_files"
ANALYSIS_DIR="analysis"

mkdir $STRUCT_DIR
mkdir $SPIN_DIR
mkdir $JSON_DIR
mkdir $ANALYSIS_DIR

STRUCT_DIR_ABS=$output_dir$OUT_DIR/$STRUCT_DIR
SPIN_DIR_ABS=$output_dir$OUT_DIR/$SPIN_DIR
JSON_DIR_ABS=$output_dir$OUT_DIR/$JSON_DIR
ANALYSIS_DIR_ABS=$output_dir$OUT_DIR/$ANALYSIS_DIR

#--------------------------
# particle settings
#--------------------------
outer_loop=4
# inner loop determines how many processes are started in parallel
# do not use too many for larger particles --> the system might crash!
inner_loop=1
# total number of particles to be calculated: inner_loop * outer_loop
num_particles=$(($outer_loop*$inner_loop))

# cif-file for crystal structure (only Vesta cif supported)
# if needed import the cif-file in Vesta and save it again under a different
# name. This will generate the cif-file in the right format.
ciffile=$CIF_DIR"/P4_32_12_perfect.cif"

# atom labels from cif-file used for occupancies (need to be comma separated)
atomlabels="Fe1,Fe2,Fe3,Fe4"
occupancies="1.0,0.83,0.83,0.83"

# unitcell parameters in Angstrom
latticepars="8.3965,8.3965,8.3965"

# particle shape (either "Sphere" or "Cube")
shape="Sphere"

# particle size in unit cells
particle_size=5

# particle orientation angles
alpha=0.0
beta=0.0
gamma=-45.0

# particle with or without APB
APB=true

# effective anisotropy constant: K[kJ/m^3] / N_atoms/particle_volume
anisotropy_constant=3.25e-25

# exchange constants: Fe_TT, Fe_OO, Fe_TO in K
exchange_constants="-21.0,-8.6,-28.1"

# exchange constant accross APB in K
APB_constant=-106.28

# seed for particle generation (determines random placement of vacancies)
seed=20210503

#--------------------------
# field settings
#--------------------------
measurement_field=5.0

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
#dipole_interactions="brute_force"
dipole_interactions="None"
#dipole_interactions="macrocell_method"

# macrocell size only used if macrocell method is selected
macrocell_size=0.2

#----------------------------------
# generate crystal structure files
#----------------------------------
if [ "$APB" = true ] ; then
    python ${generate_nanoparticle} --output=$STRUCT_DIR_ABS --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars --numAPBs=1 --template
else
    python ${generate_nanoparticle} --output=$STRUCT_DIR_ABS --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars --template
fi


for (( c=1; c<=$num_particles; c++ )) 
do
    if [ "$APB" = true ] ; then
        python ${generate_nanoparticle} --output=$STRUCT_DIR_ABS --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --number=$c --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars --numAPBs=1
    else
        python ${generate_nanoparticle} --output=$STRUCT_DIR_ABS --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --number=$c --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars
    fi
done

#--------------------------
# generate JSON files
#--------------------------
if [ "$APB" == true ] ; then
    python $generate_json_files --measurement=$measurement --output=$SPIN_DIR_ABS --structure_file_path=$STRUCT_DIR_ABS --json_file_path=$JSON_DIR_ABS --num_particles=$num_particles --particle_size=$particle_size --alpha=$alpha --beta=$beta --gamma=$gamma --meas_field=$measurement_field --temperature=$temperature --steps=$steps --dipole=$dipole_interactions --macrocell_size=$macrocell_size --sigma=$sigma --latticepar=$latticepars --APB --anisotropy_constant=$anisotropy_constant --exchangeconstants=$exchange_constants --APB_constant=$APB_constant --number=0 --seed=$seed
else
    python $generate_json_files --measurement=$measurement --output=$SPIN_DIR_ABS --structure_file_path=$STRUCT_DIR_ABS --json_file_path=$JSON_DIR_ABS --num_particles=$num_particles --particle_size=$particle_size --alpha=$alpha --beta=$beta --gamma=$gamma --meas_field=$measurement_field --temperature=$temperature --steps=$steps --dipole=$dipole_interactions --macrocell_size=$macrocell_size --sigma=$sigma --latticepar=$latticepars --anisotropy_constant=$anisotropy_constant --exchangeconstants=$exchange_constants --number=0 --seed=$seed
fi


for (( c=1; c<=$num_particles; c++ ))
do
if [ "$APB" == true ] ; then
    python $generate_json_files --measurement=$measurement --output=$SPIN_DIR_ABS --structure_file_path=$STRUCT_DIR_ABS --json_file_path=$JSON_DIR_ABS --num_particles=$num_particles --particle_size=$particle_size --alpha=$alpha --beta=$beta --gamma=$gamma --meas_field=$measurement_field --temperature=$temperature --steps=$steps --dipole=$dipole_interactions --macrocell_size=$macrocell_size --sigma=$sigma --latticepar=$latticepars --APB --anisotropy_constant=$anisotropy_constant --exchangeconstants=$exchange_constants --APB_constant=$APB_constant --number=$c --seed=$seed
else
    python $generate_json_files --measurement=$measurement --output=$SPIN_DIR_ABS --structure_file_path=$STRUCT_DIR_ABS --json_file_path=$JSON_DIR_ABS --num_particles=$num_particles --particle_size=$particle_size --alpha=$alpha --beta=$beta --gamma=$gamma --meas_field=$measurement_field --temperature=$temperature --steps=$steps --dipole=$dipole_interactions --macrocell_size=$macrocell_size --sigma=$sigma --latticepar=$latticepars --anisotropy_constant=$anisotropy_constant --exchangeconstants=$exchange_constants --number=$c --seed=$seed
fi
done

#--------------------------
# Run Simulation
#--------------------------
trap "kill 0" EXIT

# template simulation
$MCS3P $JSON_DIR_ABS/settings_0

count=0
for (( i=1; i<=$outer_loop; i++ )); do
    for ((c=1; c<=$inner_loop; c++ )); do
        ((count=count+1))
        if [ "$APB" == true ] ; then
            $MCS3P $JSON_DIR_ABS/settings_$count &
        else
            $MCS3P $JSON_DIR_ABS/settings_$count &
        fi
    done
    wait
done
wait

cd $SCRIPT_DIR


##---------------------------------
## Calculate average spin structure
##---------------------------------
t=${temperature%.*}
mf=${measurement_field%.*}
stepsf=${steps%.*}

cd $SPIN_DIR_ABS
template_filename=$(find . -maxdepth 1 -name "*template*" | sed "s|^\./||")
cd $SCRIPT_DIR

python $average_spin_structure --input_dir=$SPIN_DIR_ABS --output_dir=$ANALYSIS_DIR_ABS --template_file=$template_filename --particle_size=$particle_size

cd $ANALYSIS_DIR_ABS
averaged_file=$(find . -maxdepth 1 -name "*averaged*" | sed "s|^\./||")
cd $SCRIPT_DIR

##---------------------------------
## Generate polar plot
##---------------------------------
python $generate_polar_plot --input_dir=$SPIN_DIR_ABS --output_dir=$ANALYSIS_DIR_ABS --template_file=$template_filename --averaged_file=$averaged_file --particle_size=$particle_size --num_particles=$num_particles --dipole=$dipole_interactions

##---------------------------------
## Generate 3D image with Mayavi
##---------------------------------
## Mayavi has to be installed
## a conda environment with python version 3.7 is required
## create the environment with: conda create --name py37 python=3.7
## use: pip install mayavi
## and: pip install PyQt5
#
source activate py37
python $plot_3d --input_dir=$ANALYSIS_DIR_ABS --file=$averaged_file
conda deactivate

