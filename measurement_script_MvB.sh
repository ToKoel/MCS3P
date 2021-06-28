#!/bin/bash

# *-----------------------------------*
# | shell script for hysteresis loop  |
# |         simulations               |
# |                                   |
# |  please check the paths to the    |
# |       files before running        |
# |                                   |
# *-----------------------------------*


measurement="MvB"
echo $measurement

# directory for outputs
output_dir="/Users/tobiaskohler/PhD/thesis/Simulations/Hysteresis_sims/P43212_random_vacancies/D6/"


#--------------------------
# particle settings
#--------------------------
outer_loop=1
# inner loop determines how many processes are started in parallel
# do not use too many for larger particles --> the system might crash!
inner_loop=2
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

# particle orientations per configuration used for averaging
particle_orientations=20

# particle with or without APB
APB=false

# seed for particle generation
seed=20210503

anisotropy_constant=3.25e-25

# exchange constants: Fe_TT, Fe_TO, Fe_OO
exchange_constants="-21.0,-8.6,-28.1"
# exchange contant accross APB in K
APB_constant=-106.28

#--------------------------
# field (in tesla) and temperature (in Kelvin) settings
#--------------------------
# settings used for the field sweep
lower_field_limit=-1.5
upper_field_limit=1.5
field_step=0.1

# settings for zero field/field cooling
# starting temp. and final temp. have to be different
cooling_field=0.0
starting_temperature=301.0
temperature_step=1.0
final_temperature=5.0

#--------------------------
# Monte-Carlo settings
#--------------------------
# Number of Monte Carlo steps per temperature or field step
steps=8000
cooling_steps=1

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
# change path to "generate_nanoparticle.py" if necessary
#
#for (( c=1; c<=$num_particles; c++ ))
#do
#    if [ "$APB" = true ] ; then
 #       python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --occ_oct=0.88 --occ_tet=1.0 --diameter=$particle_size --number=$c --seed=$seed --APB
#    else
#        python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --occ_oct=0.88 --occ_tet=1.0 --diameter=$particle_size --number=$c --seed=$seed
#    fi
#done

for (( c=1; c<=$num_particles; c++ ))
do
    if [ "$APB" = true ] ; then
        python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_nanoparticle.py --output=$output_dir --atomlabels=$atomlabels --occupancies=$occupancies --diameter=$particle_size --number=$c --seed=$seed --ciffile=$ciffile --shape=$shape --latticepar=$latticepars --APB
    else
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
    python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_json_files.py --measurement=$measurement --output=$output_dir --structure_path=$output_dir --num_particles=$num_particles --particle_size=$particle_size --steps=$steps --dipole=$dipole_interactions --macrocell_size=$macrocell_size --num_or=$particle_orientations --Bupper=$upper_field_limit --Blower=$lower_field_limit --Bstep=$field_step --start_temperature=$starting_temperature --cooling_steps=$cooling_steps --temperature=$final_temperature --T_step=$temperature_step --cool_field=$cooling_field --anisotropy_constant=$anisotropy_constant --APB --sigma=$sigma --latticepar=$latticepars --exchangeconstants=$exchange_constants --APB_constant=$APB_constant
else
    python /Users/tobiaskohler/Desktop/MCS3P/MCS3P/generate_json_files.py --measurement=$measurement --output=$output_dir --structure_path=$output_dir --num_particles=$num_particles --particle_size=$particle_size --steps=$steps --dipole=$dipole_interactions --macrocell_size=$macrocell_size --num_or=$particle_orientations --Bupper=$upper_field_limit --Blower=$lower_field_limit --Bstep=$field_step --start_temperature=$starting_temperature --cooling_steps=$cooling_steps --temperature=$final_temperature --T_step=$temperature_step --cool_field=$cooling_field --anisotropy_constant=$anisotropy_constant --sigma=$sigma --latticepar=$latticepars --exchangeconstants=$exchange_constants
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
            /Users/tobiaskohler/Library/Developer/Xcode/DerivedData/MCS3P-fdtoextkqpexyiceuusyyotrcfyb/Build/Products/Release/MCS3P ${output_dir}/${count}_MvsB_simulation_APB_D${particle_size}.json &
        else
            /Users/tobiaskohler/Library/Developer/Xcode/DerivedData/MCS3P-fdtoextkqpexyiceuusyyotrcfyb/Build/Products/Release/MCS3P ${output_dir}/${count}_MvsB_simulation_noAPB_D${particle_size}.json &
        fi
    done
    wait
done
wait

