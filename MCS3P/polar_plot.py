#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 22:22:03 2021

@author: tobiaskohler
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', type=str)
parser.add_argument('--output_dir', type=str)
parser.add_argument('--template_file', type=str)
parser.add_argument('--averaged_file', type=str)
parser.add_argument('--particle_size', type=int)
parser.add_argument('--num_particles', type=int)
parser.add_argument('--dipole', type=str)
      
args = parser.parse_args()

template_file = args.template_file
averaged_file = args.averaged_file
particle_size = args.particle_size
num_particles = args.num_particles
dipole = args.dipole
center= particle_size/2.0+0.5
input_dir = args.input_dir

files = os.listdir(args.input_dir)
files = [os.path.join(args.input_dir, file) for file in files if not "template" in file]

def get_polar(vector):
    phi = np.arctan2(vector[2], vector[1])
    theta = np.arctan(np.sqrt(vector[1]**2+vector[2]**2)/vector[0])*180/np.pi
    return phi,theta

def spin_calcs(file, delimiter,center, averaged):
    if not averaged:
        x,y,z,u,v,w,a = np.loadtxt(file ,delimiter=delimiter, unpack=True, usecols=(0,1,2,3,4,5,7))
        t = np.loadtxt(file,delimiter=delimiter, unpack=True, usecols=(6), dtype=str)
        t= [1 if pos == " tetrahedral" else 0 for pos in t]
    else:
        x,y,z,u,v,w,t,a = np.loadtxt(file ,delimiter=delimiter, unpack=True, usecols=(0,1,2,3,4,5,6,7))
    
    xn = u/np.sqrt(u**2+v**2+w**2)
    yn = v/np.sqrt(u**2+v**2+w**2)
    zn = w/np.sqrt(u**2+v**2+w**2)
    
    u = np.nan_to_num(xn, nan=0.0)
    v = np.nan_to_num(yn, nan=0.0)
    w = np.nan_to_num(zn, nan=0.0)

    theta, theta_apb = [],[]
    phi, phi_apb = [],[]
    position = []
    
    for i in range(len(u)):
        if u[i] != 0.0:
            if a[i] == 1:
                theta_apb.append(np.arctan(np.sqrt(v[i]**2+w[i]**2)/u[i])*180/np.pi)
                phi_apb.append(np.arctan2(w[i], v[i]))
        
            position.append(int(t[i]))
            theta.append(np.arctan(np.sqrt(v[i]**2+w[i]**2)/u[i])*180/np.pi)
            phi.append(np.arctan2(w[i], v[i]))
        
    total_spin = (u.sum(),v.sum(),w.sum())
    
    phi_total,theta_total = get_polar(total_spin)
     
    right = np.array((0.0,0.0,0.0))
    left = np.array((0.0,0.0,0.0))
    
    
    for n,i in enumerate(x):
        if y[n] > center+0.1:
            right += np.array([u[n],v[n],w[n]])
        elif y[n] < center-0.1:
            left += np.array([u[n],v[n],w[n]])
    
    phi_left,theta_left = get_polar(left)
    phi_right,theta_right = get_polar(right)
    print(left, right)
    
    count_in_80 = 0
    count_in_85 = 0
    for i in theta:
        if 90-abs(i) >= 80:
            count_in_80 += 1
        if 90-abs(i) >= 85:
            count_in_85 += 1
    
    return position,theta,phi,theta_apb,phi_apb, phi_total, theta_total,phi_left,theta_left,phi_right,theta_right,xn,yn,zn


def plot_polar(ax,fig,file, tet, octa,apb,delimiter,center, averaged):
    position,theta,phi,theta_apb,phi_apb, phi_total, theta_total,phi_left,theta_left,phi_right,theta_right,x,y,z = spin_calcs(file,delimiter,center,averaged)
    
    if octa:
        alpha_o = 0.6
    else:
        alpha_o = 0.0
        
    if tet:
        alpha_t = 0.6
    else:
        alpha_t = 0.0
        
    if apb:
        alpha_a = 0.9
    else:
        alpha_a = 0.0
    
        
    colors = [(0.1804,0.4235,0.72157),
          (0.2471,0.9333,0.90196,alpha_t),
          (0.988,0.2667,0.2706,alpha_o),
          (0.039,0.094,0.1608,alpha_a)]
    
    colormap = np.array([colors[2],colors[1]])
    
    ax.set_rlim(bottom=90,top=70)
    ax.set_rticks([90,80,70])
    ax.set_yticklabels(["90°","80°","70°"], fontsize=24)
    ax.set_xticklabels(["0°","45°","90°","135°","180°","225°","270°","315°"],fontsize=24)

    bbox = dict(boxstyle="round,pad=0.0", ec=(1,1,1,0.0), fc="white", alpha=0.6)
    plt.setp(ax.get_ymajorticklabels(), bbox=bbox)
    ax.scatter(phi,90-abs(np.array(theta)),color=colormap[position], marker='.', lw=0,s=40)
    ax.scatter(phi_apb,90-abs(np.array(theta_apb)),color=colors[3],marker='.',lw=0,s=40)
    
    ax2 = fig.add_axes([0.27,0.7,0.24,0.24],polar=True)
    bbox = dict(boxstyle="round,pad=0.0", ec=(1,1,1,0.0), fc="white", alpha=0.6)
    ax2.set_rlim(bottom=90,top=70)
    ax2.set_rticks([80])
    ax2.set_yticklabels(["80°"],fontsize=14)
    ax2.set_xticklabels([])
    ax2.set_axisbelow(True)
    plt.setp(ax2.get_ymajorticklabels(), bbox=bbox)
    ax2.scatter(phi_left,90-abs(np.array(theta_left)),color=(0.9,0.3,0.3,1.0),marker='<',lw=0,s=40)
    ax2.scatter(phi_right,90-abs(np.array(theta_right)),color=(0.3,0.3,0.9,1.0),marker='>',lw=0,s=40)
    ax2.scatter(phi_total,90-abs(np.array(theta_total)),color=(0.1,0.1,0.1,1.0),marker='.',lw=0,s=100)
    
    
fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)    
for i in files:
    plot_polar(ax,fig,i, False, False,False, ",",center, False)
    
plot_polar(ax,fig,os.path.join(args.output_dir, averaged_file), True, True,True, " ",center, True)
plt.savefig(os.path.join(args.output_dir, "polar_plot.png"), dpi=300)
#plt.show()









