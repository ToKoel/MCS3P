#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 11:57:04 2021

@author: tobiaskohler

Running of this script requires MayaVi installed in python 3.7 environment.
The script has to be executed in this environment.
"""

import numpy as np
from mayavi import mlab
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', type=str)
parser.add_argument('--file', type=str)
parser.add_argument('--delimiter', type=str)

args = parser.parse_args()

file = args.file
input_dir = args.input_dir

if args.delimiter == None:
    delimiter=" "
else:
    delimiter=args.delimiter

sfile = input_dir+file

file_r = os.path.splitext(sfile)[0]

def plot_mayavi(spins, size, x,y,z,u,v,w,t,a):
    
    # setting up the data arrays
    # octahedral iron positions
    x_oct = np.array([xi for xi,ti,ai in zip(x,t,a) if ti == 0 and ai == 0])
    y_oct = np.array([yi for yi,ti,ai in zip(y,t,a) if ti == 0 and ai == 0])
    z_oct = np.array([zi for zi,ti,ai in zip(z,t,a) if ti == 0 and ai == 0])
    u_oct = np.array([ui for ui,ti,ai in zip(u,t,a) if ti == 0 and ai == 0])
    v_oct = np.array([vi for vi,ti,ai in zip(v,t,a) if ti == 0 and ai == 0])
    w_oct = np.array([wi for wi,ti,ai in zip(w,t,a) if ti == 0 and ai == 0])
    
    # tetrahedral iron positions
    x_tet = np.array([xi for xi,ti,ai in zip(x,t,a) if ti == 1 and ai == 0])
    y_tet = np.array([yi for yi,ti,ai in zip(y,t,a) if ti == 1 and ai == 0])
    z_tet = np.array([zi for zi,ti,ai in zip(z,t,a) if ti == 1 and ai == 0])
    u_tet = np.array([ui for ui,ti,ai in zip(u,t,a) if ti == 1 and ai == 0])
    v_tet = np.array([vi for vi,ti,ai in zip(v,t,a) if ti == 1 and ai == 0])
    w_tet = np.array([wi for wi,ti,ai in zip(w,t,a) if ti == 1 and ai == 0])
    
    # APB positions
    x_apb = np.array([xi for xi,ai in zip(x,a) if ai == 1])
    y_apb = np.array([yi for yi,ai in zip(y,a) if ai == 1])
    z_apb = np.array([zi for zi,ai in zip(z,a) if ai == 1])
    u_apb = np.array([ui for ui,ai in zip(u,a) if ai == 1])
    v_apb = np.array([vi for vi,ai in zip(v,a) if ai == 1])
    w_apb = np.array([wi for wi,ai in zip(w,a) if ai == 1])

    # plotting balls for the atomic positions and vectors for the spins with transparency
    # and color coding according to the spin orientation (Sx (u) component)
    
    transp = 0.4 # transparency setting for the balls
    sf = 0.2 # scale factor for spin vectors
    
    # APB atoms
    if x_apb.size != 0:
        apb_atom = mlab.points3d(x_apb, y_apb, z_apb, u_apb,
                                 scale_mode='none', scale_factor=size, resolution=30, figure=fig)
        apb_atom.glyph.color_mode = 'color_by_scalar'
        apb_atom.module_manager.scalar_lut_manager.data_range = np.array([-1., 1.])
        apb_atom.module_manager.scalar_lut_manager.lut.alpha_range = np.array([transp, transp])
        
        apb_spin = mlab.quiver3d(x_apb,y_apb,z_apb,
                                 u_apb,v_apb,w_apb,scalars=np.array(u_apb),mode='arrow',
                                 scale_factor=sf, resolution=30, figure=fig)
        apb_spin.glyph.color_mode = 'color_by_scalar'
        apb_spin.module_manager.scalar_lut_manager.data_range = np.array([-1., 1.])
    
    # tetrahedral positions
    tet_atom = mlab.points3d(x_tet,y_tet,z_tet,u_tet,
                             scale_mode='none', scale_factor=size, resolution=30, figure=fig)
    tet_atom.glyph.color_mode = 'color_by_scalar'
    tet_atom.module_manager.scalar_lut_manager.data_range = np.array([-1., 1.])
    tet_atom.module_manager.scalar_lut_manager.lut.alpha_range = np.array([transp, transp])
    
    tet_spin = mlab.quiver3d(x_tet,y_tet,z_tet,
                             u_tet,v_tet, w_tet,
                             scalars=np.array(u_tet),mode='arrow',
                             scale_factor=sf, resolution=30, figure=fig)
    tet_spin.glyph.color_mode = 'color_by_scalar'
    tet_spin.module_manager.scalar_lut_manager.data_range = np.array([-1., 1.])
    
    # octahedral positions
    oct_atom = mlab.points3d(x_oct,y_oct,z_oct,u_oct,
                             scale_mode='none', scale_factor=size, resolution=30, figure=fig)
    oct_atom.glyph.color_mode = 'color_by_scalar'
    oct_atom.module_manager.scalar_lut_manager.data_range = np.array([-1., 1.])
    oct_atom.module_manager.scalar_lut_manager.lut.alpha_range = np.array([transp, transp])
    
    oct_spin = mlab.quiver3d(x_oct,y_oct,z_oct,
                             u_oct,v_oct,w_oct,
                             scalars=np.array(u_oct),mode='arrow',
                             scale_factor=sf, resolution=30, figure=fig)
    oct_spin.glyph.color_mode = 'color_by_scalar'
    oct_spin.module_manager.scalar_lut_manager.data_range = np.array([-1., 1.])
    
    # Frame for APB location
    p= 4.5
    xp = np.array([0.5,0.5,7,7,0.5])
    yp = np.array([p,p,p,p,p])
    zp = np.array([0.5,7,7,0.5,0.5])
    #mlab.plot3d(xp,yp,zp, tube_radius=0.01,color=(0.3,0.3,0.3))
    
    # coordinate axis
    mlab.quiver3d(0.5,0,1.0,1.5,0,0,color=(0.4,0.4,0.4),mode='arrow', resolution=40)
    mlab.quiver3d(0.5,0,1.0,0,1.5,0,color=(0.4,0.4,0.4),mode='arrow', resolution=40)
    mlab.quiver3d(0.5,0,1.0,0,0,1.5,color=(0.4,0.4,0.4),mode='arrow', resolution=40)

# setting up the figure (size, lighting, camera position etc.)
fig = mlab.figure(size=(1000,1000), bgcolor=(1,1,1))

# in case a colorbar is needed
# get the current lut manager
#lut_manager = mlab.colorbar(orientation='vertical')
# fix the range
#lut_manager.data_range = (-1, 1)

# these lines are needed for the lights to work properly, not sure what they do though..
from pyface.api import GUI
_gui = GUI()

while fig.scene.light_manager is None:
    _gui.process_events()

fig.scene.light_manager.lights[0].elevation = 32.37
fig.scene.light_manager.lights[0].azimuth = -5.0
fig.scene.light_manager.lights[0].intensity = 0.9322

fig.scene.light_manager.lights[1].elevation = -34.29
fig.scene.light_manager.lights[1].azimuth = -10.0
fig.scene.light_manager.lights[1].intensity = 0.1582

fig.scene.light_manager.lights[2].elevation = -41.43
fig.scene.light_manager.lights[2].azimuth = 15.48
fig.scene.light_manager.lights[2].intensity = 0.4915

fig.scene.light_manager.lights[3].elevation = 0.0
fig.scene.light_manager.lights[3].azimuth = 0.0
fig.scene.light_manager.lights[3].intensity = 0.0


x,y,z,u,v,w,t,a = np.loadtxt(sfile, usecols=(0,1,2,3,4,5,6,7), delimiter=delimiter, unpack=True)
plot_mayavi(True,0.15,x,y,z,u,v,w,t,a)

mlab.gcf().scene.parallel_projection = True

mlab.view(azimuth=90,elevation=180, distance=20)
#mlab.view(azimuth=90,elevation=160, distance=20)
#mlab.view(azimuth=10,elevation=10, distance=20)

fig.scene.anti_aliasing_frames = 20
mlab.draw()
fig.scene.camera.zoom(1.2)

mlab.savefig(file_r+".png", size=(2000,2000))
mlab.show()

