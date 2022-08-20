#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 22:11:48 2021

@author: tobiaskohler
"""
import numpy as np
import scipy.stats as stats
import argparse
import os


parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', type=str)
parser.add_argument('--output_dir', type=str)
parser.add_argument('--template_file', type=str)
parser.add_argument('--particle_size', type=int)

args = parser.parse_args()

def averaging_spins(args):
    files = os.listdir(args.input_dir)
    files = [file for file in files if not "template" in file]
    outfile = os.path.join(args.output_dir, args.template_file.replace("template", "averaged"))
    
    xt,yt,zt,at = np.loadtxt(os.path.join(args.input_dir, args.template_file),
                                delimiter=',', usecols=(0,1,2,7), unpack = True)
    position = np.loadtxt(os.path.join(args.input_dir, args.template_file),
                                delimiter=',', usecols=(6), unpack = True, dtype=str)

    ut = np.zeros_like(xt)
    vt = np.zeros_like(xt)
    wt = np.zeros_like(xt)

    coord_spin_t = list(zip(xt,yt,zt,ut,vt,wt))
    dict1 = {(x,y,z): (u,v,w) for x,y,z,u,v,w in coord_spin_t}

    for file in files:

        x,y,z,u,v,w,a = np.loadtxt(os.path.join(args.input_dir,file),delimiter=',', unpack=True, usecols=(0,1,2,3,4,5,7))

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
    lengths = np.sqrt(uf**2+vf**2+wf**2)
    uf /= lengths
    vf /= lengths
    wf /= lengths
    
    uf = np.nan_to_num(uf, nan=0.0)
    vf = np.nan_to_num(vf, nan=0.0)
    wf = np.nan_to_num(wf, nan=0.0)
    
    tt = [1 if pos == " tetrahedral" else 0 for pos in position]
    np.savetxt(outfile, np.column_stack((xf,yf,zf,uf,vf,wf,tt,at)), fmt='%.5f')

averaging_spins(args)
    
    
