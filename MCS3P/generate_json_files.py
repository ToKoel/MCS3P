#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 16:43:01 2021

@author: tobiaskohler
"""

import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--measurement', type=str)
    parser.add_argument('--output', type=str)
    parser.add_argument('--structure_path', type=str)
    parser.add_argument('--num_particles', type=int)
    parser.add_argument('--particle_size', type=int)
    parser.add_argument('--meas_field', type=float)
    parser.add_argument('--temperature', type=float)
    parser.add_argument('--dipole', type=str)
    parser.add_argument('--steps', type=float)
    parser.add_argument('--av_steps', type=int)
    parser.add_argument('--num_or', type=int)
    parser.add_argument('--cool_field', type=float)
    parser.add_argument('--Tupper', type=float)
    parser.add_argument('--Tlower', type=float)
    parser.add_argument('--T_step', type=float)
    parser.add_argument('--alpha', type=float)
    parser.add_argument('--beta', type=float)
    parser.add_argument('--gamma', type=float)
    parser.add_argument('--macrocell_size', type=float)
    parser.add_argument('--anisotropy_constant', type=float)
    parser.add_argument('--ZFC', action='store_true')
    parser.add_argument('--FC', action='store_true')
    parser.add_argument('--APB', action='store_true')
    parser.add_argument('--Bupper', type=float)
    parser.add_argument('--Blower', type=float)
    parser.add_argument('--Bstep', type=float)
    parser.add_argument('--start_temperature', type=float)
    parser.add_argument('--cooling_steps', type=float)
    parser.add_argument('--sigma', type=float)
    parser.add_argument('--latticepar', type=str)
    parser.add_argument('--exchangeconstants', type=str)
    parser.add_argument('--APB_constant', type=float)
    
    args = parser.parse_args()
    
    measurement = args.measurement
    #print("Measurement: ",measurement)
    save_path = args.output
    #print("Output path: ",save_path)
    structure_path = args.structure_path
    #print("Structure path: ",structure_path)
    num_particles = args.num_particles
    #print("Number of particles: ",num_particles)
    particle_size = args.particle_size
    #print("Particle diameter [UC]: ",particle_size)
    measurement_field = args.meas_field
    #print("Magnetic field: ",measurement_field)
    temperature = args.temperature
    #print("Temperature: ",temperature)
    dipole = args.dipole
    #print("Dipole interactions: ",dipole)
    steps = args.steps
    #print("Steps: ",steps)
    averaging_steps = args.av_steps
    #print("Averaging steps: ",averaging_steps)
    num_orientations = args.num_or
    #print("Number of particle orientations: ",num_orientations)
    cooling_field = args.cool_field
    #print("Cooling field: ",cooling_field)
    T_upper = args.Tupper
    #print("Upper limit T: ",T_upper)
    T_lower = args.Tlower
    #print("Lower limit T: ",T_lower)
    T_stepsize = args.T_step
    #print("Temperature step size: ",T_stepsize)
    ZFC = args.ZFC
    #print("ZFC: ",ZFC)
    FC = args.FC
    #print("FC: ",FC)
    APB = args.APB
    #print("APB: ",APB)
    alpha = args.alpha
    #print("alpha: ",alpha)
    beta = args.beta
    #print("beta: ",beta)
    gamma = args.gamma
    #print("gamma: ",gamma)
    macrocell_size = args.macrocell_size
    #print("macrocell size: ",macrocell_size)
    anisotropy_constant= args.anisotropy_constant
    #print("anisotropy constant: ", anisotropy_constant)
    sigma = args.sigma
    
    exchange_constants = (args.exchangeconstants).split(",")
    Fe_TT = float(exchange_constants[0])
    Fe_OO = float(exchange_constants[1])
    Fe_TO = float(exchange_constants[2])
    APB_constant = args.APB_constant
    
    latticepars = (args.latticepar).split(",")
    lattice_a = float(latticepars[0])
    lattice_b = float(latticepars[1])
    lattice_c = float(latticepars[2])
    
    Bupper = args.Bupper
    Blower = args.Blower
    Bstep = args.Bstep
    start_T = args.start_temperature
    cool_steps= args.cooling_steps
    
    
    if measurement == "MvB":
        print("\tStarting M vs B simulation\n")
        print("Output path: ",save_path)
        print("Number of particles: ",num_particles)
        print("Number of orientations per particle: ", num_orientations)
        print("Measurement field: ",measurement_field)
        print("Cooling field: ", cooling_field)
        print("Field range: %.1f - %.1f T (%.2f T steps)"%(Blower, Bupper, Bstep))
        print("Monte Carlo steps per temperature step: ", steps)
        print("Averaging Monte Carlo steps: ", averaging_steps)
        print("Sigma (gaussian cone): ", sigma)
        print(" ")
        for i in range(1,num_particles+1):
            if APB:
                file_name = "%d_MvsB_simulation_APB_D%d.json"%(i, particle_size)
                structure_file = structure_path + "D%d_structure_APB_%d"%(particle_size,i)
            else:
                file_name = "%d_MvsB_simulation_noAPB_D%d.json"%(i, particle_size)
                structure_file = structure_path + "D%d_structure_noAPB_%d"%(particle_size,i)
                
            completeName = save_path+file_name
            
            f = open(completeName, "w+")
            f.write("{\n")
            f.write('"Measurement": "M vs B",\n') 
            f.write('"output_dir": "%s",\n' %(save_path))
            f.write('"structure_file": "%s",\n'%(structure_file))
            f.write('"dipole_interactions": "%s",\n'%(dipole))
            f.write('"steps": %d,\n'%(steps)) 
            f.write('"num_orientations": %d,\n'%(num_orientations))
            f.write('"temperature": %.4f,\n'%(temperature))
            f.write('"B_upper": %.4f,\n'%(Bupper))
            f.write('"B_lower": %.4f,\n'%(Blower))
            f.write('"B_step": %.4f,\n'%(Bstep))
            f.write('"cooling_field": %.4f,\n'%(cooling_field))
            f.write('"cooling_steps": %d,\n'%(cool_steps))
            f.write('"start_temperature": %.4f,\n'%(start_T))
            f.write('"temperature_step": %.4f,\n'%(T_stepsize))
            f.write('"FeTT": %.2f,\n'%(Fe_TT))
            f.write('"FeOO": %.2f,\n'%(Fe_OO))
            f.write('"FeTO": %.2f,\n'%(Fe_TO))
            f.write('"FeOO_APB": %.2f,\n'%(APB_constant))
            f.write('"anisotropy constant": %.2e,\n'%(anisotropy_constant))
            f.write('"lattice_a": %.4f,\n'%(lattice_a))
            f.write('"lattice_b": %.4f,\n'%(lattice_b))
            f.write('"lattice_c": %.4f,\n'%(lattice_c))
            f.write('"sigma": %.2f,\n'%(sigma))
            f.write('"particle center": %.2f\n}'%(particle_size/2.0+0.5))
            f.close() 
    
    if measurement == "MvT":
        print("\tStarting M vs T simulation\n")
        print("Output path: ",save_path)
        print("Number of particles: ",num_particles)
        print("Number of orientations per particle: ", num_orientations)
        print("Measurement field: ",measurement_field)
        print("Cooling field: ", cooling_field)
        print("Temperature range: %.4f - %.4f K (%.4f K steps)"%(T_lower, T_upper, T_stepsize))
        print("Monte Carlo steps per temperature step: ", steps)
        print("Averaging Monte Carlo steps: ", averaging_steps)
        print("Sigma (gaussian cone): ", sigma)
        Meas_modes ="Measurements: "
        if ZFC:
            Meas_modes += "ZFC"
        if FC:
            Meas_modes += ", FC"
        print(Meas_modes)
        print(" ")
        
        for i in range(1,num_particles+1):
            if APB:
                file_name = "%d_MvsT_simulation_APB_D%d.json"%(i, particle_size)
                structure_file = structure_path + "D%d_structure_APB_%d"%(particle_size,i)
            else:
                file_name = "%d_MvsT_simulation_noAPB_D%d.json"%(i, particle_size)
                structure_file = structure_path + "D%d_structure_noAPB_%d"%(particle_size,i)
                
            completeName = save_path+file_name
            
            f = open(completeName, "w+")
            f.write("{\n")
            f.write('"Measurement": "M vs T",\n') 
            f.write('"output_dir": "%s",\n' %(save_path))
            f.write('"structure_file": "%s",\n'%(structure_file))
            f.write('"dipole_interactions": "%s",\n'%(dipole))
            f.write('"steps": %.2f,\n'%(steps))
            f.write('"averaging steps": %d,\n'%(averaging_steps))
            f.write('"num_orientations": %d,\n'%(num_orientations))
            f.write('"measurement_field": %.4f,\n'%(measurement_field))
            f.write('"cooling_field": %.4f,\n'%(cooling_field))
            f.write('"TUpperLimit": %.4f,\n'%(T_upper))
            f.write('"TLowerLimit": %.4f,\n'%(T_lower))
            f.write('"TstepSize": %.4f,\n'%(T_stepsize))
            f.write('"FeTT": %.2f,\n'%(Fe_TT))
            f.write('"FeOO": %.2f,\n'%(Fe_OO))
            f.write('"FeTO": %.2f,\n'%(Fe_TO))
            f.write('"FeOO_APB": %.2f,\n'%(APB_constant))
            f.write('"anisotropy constant": %.2e,\n'%(anisotropy_constant))
            if ZFC:
                f.write('"ZFC": true,\n')
            else:
                f.write('"ZFC": false,\n')
            if FC:
                f.write('"FC": true,\n')
            else: 
                f.write('"FC": false,\n')
            f.write('"particle center": %.2f,\n'%(particle_size/2.0+0.5))
            f.write('"lattice_a": %.4f,\n'%(lattice_a))
            f.write('"lattice_b": %.4f,\n'%(lattice_b))
            f.write('"lattice_c": %.4f,\n'%(lattice_c))
            f.write('"sigma": %.2f\n}'%(sigma))
            f.close() 

    if measurement == "spin_structure":
        print("\tStarting spin structure simulation\n")
        print("Output path: ",save_path)
        print("Number of particles: ",num_particles)
        print("Particle orientation angles: %.1f°, %.1f°, %.1f°"%(alpha, beta, gamma))
        print("Measurement field: ",measurement_field)
        print("Temperature: ", temperature)
        print("Monte Carlo steps: ", steps)
        print("Sigma (gaussian cone): ", sigma)
        print("Dipole interactions: ", dipole)
        print(" ")
        
        if APB:
            file_name = "spin_structure_APB_D%d_template.json"%(particle_size)
            structure_file = structure_path + "D%d_structure_APB_template"%(particle_size)
        else: 
            file_name = "spin_structure_noAPB_D%d_template.json"%(particle_size)
            structure_file = structure_path + "D%d_structure_noAPB_template"%(particle_size)
            
        completeName = save_path+file_name
             
        f = open(completeName, "w+")
        f.write("{\n")
        f.write('"Measurement": "spin structure",\n')
        f.write('"output_dir": "%s",\n' %(save_path))
        f.write('"structure_file": "%s",\n'%(structure_file))
        f.write('"dipole_interactions": "None",\n')
        f.write('"steps": %d,\n'%(steps))
        f.write('"magnetic field": %.4f,\n'%(measurement_field))
        f.write('"temperature": %.4f,\n'%(temperature))
        f.write('"FeTT": %.2f,\n'%(Fe_TT))
        f.write('"FeOO": %.2f,\n'%(Fe_OO))
        f.write('"FeTO": %.2f,\n'%(Fe_TO))
        f.write('"FeOO_APB": %.2f,\n'%(APB_constant))
        f.write('"anisotropy constant": %.2e,\n'%(anisotropy_constant))
        #f.write('"anisotropy constant": %.2e,\n'%(anisotropy_constant))
        f.write('"alpha": %.1f,\n'%(alpha))
        f.write('"beta": %.1f,\n'%(beta))
        f.write('"gamma": %.1f,\n'%(gamma))
        f.write('"macrocell_size": %.2f,\n'%(macrocell_size))
        f.write('"particle center": %.2f,\n'%(particle_size/2.0+0.5))
        f.write('"lattice_a": %.4f,\n'%(lattice_a))
        f.write('"lattice_b": %.4f,\n'%(lattice_b))
        f.write('"lattice_c": %.4f,\n'%(lattice_c))
        f.write('"sigma": %.2f\n}'%(sigma))
        f.close()      
        
        for i in range(1,num_particles+1):
            if APB:
                file_name = "%d_spin_structure_APB_D%d.json"%(i, particle_size)
                structure_file = structure_path + "D%d_structure_APB_%d"%(particle_size,i)
            else: 
                file_name = "%d_spin_structure_noAPB_D%d.json"%(i, particle_size)
                structure_file = structure_path + "D%d_structure_noAPB_%d"%(particle_size,i)
            
            completeName = save_path+file_name
             
            f = open(completeName, "w+")
            f.write("{\n")
            f.write('"Measurement": "spin structure",\n')
            f.write('"output_dir": "%s",\n' %(save_path))
            f.write('"structure_file": "%s",\n'%(structure_file))
            f.write('"dipole_interactions": "%s",\n'%(dipole))
            f.write('"steps": %d,\n'%(steps))
            f.write('"magnetic field": %.4f,\n'%(measurement_field))
            f.write('"temperature": %.4f,\n'%(temperature))
            f.write('"FeTT": %.2f,\n'%(Fe_TT))
            f.write('"FeOO": %.2f,\n'%(Fe_OO))
            f.write('"FeTO": %.2f,\n'%(Fe_TO))
            f.write('"FeOO_APB": %.2f,\n'%(APB_constant))
            f.write('"anisotropy constant": %.2e,\n'%(anisotropy_constant))
            f.write('"alpha": %.1f,\n'%(alpha))
            f.write('"beta": %.1f,\n'%(beta))
            f.write('"gamma": %.1f,\n'%(gamma))
            f.write('"macrocell_size": %.2f,\n'%(macrocell_size))
            f.write('"particle center": %.2f,\n'%(particle_size/2.0+0.5))
            f.write('"lattice_a": %.4f,\n'%(lattice_a))
            f.write('"lattice_b": %.4f,\n'%(lattice_b))
            f.write('"lattice_c": %.4f,\n'%(lattice_c))
            f.write('"sigma": %.2f\n}'%(sigma))
            f.close() 
            
    
if __name__ == "__main__":
    main()

