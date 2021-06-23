#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 22:11:48 2021

@author: tobiaskohler
"""
import numpy as np
import scipy.stats as stats
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', type=str)
parser.add_argument('--template_file', type=str)
parser.add_argument('--particle_size', type=int)
parser.add_argument('--num_particles', type=int)
parser.add_argument('--dipole', type=str)
      
args = parser.parse_args()

template_file = args.template_file
particle_size = args.particle_size
num_particles = args.num_particles
dipole = args.dipole
center= particle_size/2.0+0.5
input_dir = args.input_dir

files = []
for i in range(1,num_particles+1):
    fname = template_file.replace("template","%d"%(i))
    fname = fname.replace("None", dipole)
    files.append(input_dir+fname)

outfile = input_dir+template_file.replace("template", "averaged")
outfile = outfile.replace("None", dipole)
    
def averaging_spins(files, template_file, filename, center):
    xt,yt,zt,tt,at = np.loadtxt(input_dir+template_file,delimiter=',', usecols=(0,1,2,6,7), unpack = True)      

    ut = np.zeros_like(xt)
    vt = np.zeros_like(xt)
    wt = np.zeros_like(xt)
    
    coord_spin_t = list(zip(xt,yt,zt,ut,vt,wt))
    dict1 = {(x,y,z): (u,v,w) for x,y,z,u,v,w in coord_spin_t}
        
    for file in files:
        
        x,y,z,u,v,w,t,a = np.loadtxt(file,delimiter=',', unpack=True)
        print(len(x))
 
        coord_spin = list(zip(x,y,z,u,v,w))
        dict2 = {(x,y,z): (u,v,w) for x,y,z,u,v,w in coord_spin}
        
        matchxyz = set(dict1) & set(dict2)
        
        for xyz in matchxyz:
            a = list(dict1[xyz])
            b = list(dict2[xyz])
            for n in range(3):
                a[n] += b[n]
            dict1[xyz] = a
     
    uf,vf,wf = [],[],[]
    xf,yf,zf = [],[],[]
    for key in dict1:
        x,y,z = key
        u,v,w = dict1[key]
        uf.append(u)
        vf.append(v)
        wf.append(w)
        xf.append(x)
        yf.append(y)
        zf.append(z)
    
    uf,vf,wf = np.array(uf), np.array(vf), np.array(wf)
    print(len(ut))
    print(len(uf))
    print(len(xf))    
    
    lengths = np.sqrt(uf**2+vf**2+wf**2)
    uf /= lengths
    vf /= lengths
    wf /= lengths
    np.savetxt(filename, np.column_stack((xf,yf,zf,uf,vf,wf,tt,at)), fmt='%.5f')
     
averaging_spins(files, template_file, outfile, center)


    
    
    