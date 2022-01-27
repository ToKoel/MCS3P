#!/bin/bash

# *-----------------------------------*
# |     shell script for M vs. T      |
# |           simulations             |
# |                                   |
# |  please check the paths to the    |
# |       files before running        |
# |                                   |
# |    installation of Mayavi is      |
# |    required for the 3D plots      |
# |     (see description below)       |
# *-----------------------------------*


measurement="MvT"
echo $measurement

# directory for outputs
output_dir="output_dir/"


#--------------------------
# particle settings
#--------------------------
outer_loop=1
# inner loop determines how many processes are started in parallel
# do not use too many for larger particles --> the system might crash!
inner_loop=1
# total number of particles to be calculated
num_particles=$(($outer_loop*$inner_loop))

# cif-file for crystal structure (only Vesta cif supported)
# if needed import the cif-file in Vesta and save it again under a different
# name. This will generate the cif-file in the right format.
ciffile="./MCS3P/cif-files/Fd-3m_perfect.cif"

# atom labels from cif-file used for occupancies (need to be comma separated)
atomlabels="Fe(oct),Fe(tet)"
occupancies="0.83,1.0"
latticepars="8.3965,8.3965,8.3965"

# particle shape (either "Sphere" or "Cube")
shape="Sphere"

# particle size (diameter) in unit cells
particle_size=6

# particle with or without APB
APB=false

# seed for particle generation
seed=20210503

# particle starting orientations
particle_orientations=9999

# magnetocrystalline anisotropy
anisotropy_constant=3.25e-25

# exchange constants: Fe_TT, Fe_OO, Fe_TO
exchange_constants="-21.0,-8.6,-28.1"
#exchange_constants="-4.9,28.0,-129.4"
#exchange_constants="-9.8,56.0,258.8"
# exchange contant accross APB in K
APB_constant=-106.28

#--------------------------
# field settings
#--------------------------
measurement_field=0.005
cooling_field=0.0

#--------------------------
# temperature settings
#--------------------------
temperature_upper=30.0
temperature_lower=0.001
temperature_step=0.5

#--------------------------
# Measurement settings
#--------------------------
# number of complete Monte Carlo steps = steps x totalNumAtoms
relaxation_steps=500.0
averaging_steps=0

# sigma parameter for opening of gaussian cone in trial move
sigma=0.03

# setting for method of dipole interaction calculation
#dipole_interactions="brute_force"
dipole_interactions="None"
#dipole_interactions="macrocell_method"

# macrocell size only used if macrocell method is selected
macrocell_size=0.2

# select the measurements
ZFC=true
FC=true


#----------------------------------
# generate crystal structure files
#----------------------------------
# change path to "generate_nanoparticle.py" if necessary
#
for (( c=1; c<=$num_particles; c++ )) 
do
    if [ "$APB" = true ] ; then      
        python ./MCS3P/generate_nanoparticle.py --output=$output_dir --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --number=$c --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars --APB
    else       
        python ./MCS3P/generate_nanoparticle.py --output=$output_dir --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --number=$c --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars
    fi
done

#--------------------------
# generate JSON files
#--------------------------
# these files are used as input parameter files for the simulation program
# change path to "generate_json_files.py" if necessary
#
if [ "$APB" == true ] ; then
    if [ "$ZFC" == true ] ; then
        if [ "$FC" == true ] ; then
            python ./MCS3P/generate_json_files.py --measurement=$measurement --output=$output_dir --structure_path=$output_dir --num_particles=$num_particles --particle_size=$particle_size --meas_field=$measurement_field --dipole=$dipole_interactions --steps=$relaxation_steps --av_steps=$averaging_steps --num_or=$particle_orientations --cool_field=$cooling_field --Tupper=$temperature_upper --Tlower=$temperature_lower --T_step=$temperature_step --anisotropy_constant=$anisotropy_constant --sigma=$sigma --ZFC --FC --APB --latticepar=$latticepars --exchangeconstants=$exchange_constants --APB_constant=$APB_constant
        else
            python ./MCS3P/generate_json_files.py --measurement=$measurement --output=$output_dir --structure_path=$output_dir --num_particles=$num_particles --particle_size=$particle_size --meas_field=$measurement_field --dipole=$dipole_interactions --steps=$relaxation_steps --av_steps=$averaging_steps --num_or=$particle_orientations --cool_field=$cooling_field --Tupper=$temperature_upper --Tlower=$temperature_lower --T_step=$temperature_step --anisotropy_constant=$anisotropy_constant --sigma=$sigma --ZFC --APB --latticepar=$latticepars --exchangeconstants=$exchange_constants --APB_constant=$APB_constant
        fi
    fi
else
    if [ "$ZFC" == true ] ; then
        if [ "$FC" == true ] ; then
            python ./MCS3P/generate_json_files.py --measurement=$measurement --output=$output_dir --structure_path=$output_dir --num_particles=$num_particles --particle_size=$particle_size --meas_field=$measurement_field --dipole=$dipole_interactions --steps=$relaxation_steps --av_steps=$averaging_steps --num_or=$particle_orientations --cool_field=$cooling_field --Tupper=$temperature_upper --Tlower=$temperature_lower --T_step=$temperature_step --anisotropy_constant=$anisotropy_constant --sigma=$sigma --ZFC --FC --latticepar=$latticepars --exchangeconstants=$exchange_constants --APB_constant=$APB_constant
        else
            python ./MCS3P/generate_json_files.py --measurement=$measurement --output=$output_dir --structure_path=$output_dir --num_particles=$num_particles --particle_size=$particle_size --meas_field=$measurement_field --dipole=$dipole_interactions --steps=$relaxation_steps --av_steps=$averaging_steps --num_or=$particle_orientations --cool_field=$cooling_field --Tupper=$temperature_upper --Tlower=$temperature_lower --T_step=$temperature_step --anisotropy_constant=$anisotropy_constant --sigma=$sigma --ZFC --latticepar=$latticepars --exchangeconstants=$exchange_constants --APB_constant=$APB_constant
        fi
    fi
fi

#--------------------------
# Run Simulation
#--------------------------
# change path to the compiled simulation program if necessary
#
trap "kill 0" EXIT

count=0
for (( i=1; i<=$outer_loop; i++ )); do
    for ((c=1; c<=$inner_loop; c++ )); do
        ((count=count+1))
        if [ "$APB" == true ] ; then
            build_dir/MCS3P ${output_dir}/${count}_MvsT_simulation_APB_D${particle_size}.json &
        else
            build_dir/MCS3P ${output_dir}/${count}_MvsT_simulation_noAPB_D${particle_size}.json &
        fi
    done
    wait
done
wait
