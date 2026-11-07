# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:43:49 2023

@author: la506
"""

import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy
from magpylib.magnet import Cuboid, CylinderSegment
import itertools
from scipy.spatial.transform import Rotation as R
import pandas as pd
import cProfile
import sys
from multiprocessing.pool import ThreadPool
from multiprocessing import Pool, freeze_support
from os import getpid
import time
import pickle
import scipy
from loguru import logger
from magpylib_material_response.demag import apply_demag
from magpylib_material_response.meshing import mesh_all
# import addcopyfighandler

#%%


# function defines magnet positions in space. this function can be replaced to incorporate more complex shapes/locations
def define_magnet_positions(xmin,
                            xmax,
                            ymin,
                            ymax,
                            zmin,
                            zmax,
                            spacing,
                            angleforcurvature,
                            plotting=0,
                            magnetlength=25.4/2,
                            tofile='positions.xlsx'):
    
    print('generating all locations for possible magnet locations')

    
    total_positions = np.zeros((1,3)) 

    # assuming spacing is the same in x y and z.
    xincrement = spacing
    yincrement = spacing
    zincrement = spacing
    
    x_vect_odd = np.reshape(np.round(np.arange(xmin, xmax, xincrement),2),(-1,1))
    x_vect_even = np.zeros((x_vect_odd.shape[0]-1,1))
    y_vect = np.round(np.arange(ymin, ymax, yincrement),2)
    z_vect = np.round(np.arange(zmin, zmax, zincrement),2)
    
    layers = y_vect.shape[0]
    
    for ii in range(x_vect_odd.shape[0]-1):
        # print(ii)
        x_vect_even[ii] = (x_vect_odd[ii]+x_vect_odd[ii+1])/2
    
    # here I need to implment different layers
    for layernum in range(layers):
        
        ycurr = y_vect[layernum] # change this for each iteration of y

        if (layernum%2==0): # check if odd layer (starting at index 0)
            points_for_magnets = np.zeros((x_vect_odd.shape[0]*1*z_vect.shape[0],3))
            aa=0
            for ii in range(x_vect_odd.shape[0]):
                for kk in range(z_vect.shape[0]):
                    # print(aa)
                    points_for_magnets[aa,0] = x_vect_odd[ii]
                    points_for_magnets[aa,1] = ycurr
                    points_for_magnets[aa,2] = z_vect[kk]
                    aa=aa+1
                    
            dmin = np.sqrt(zincrement**2 + yincrement**2)
            ballrad = ycurr # setting the radius to the y height of this iteration
            
            desiredphaseinc = np.arcsin(dmin/(2*ballrad))*180/np.pi #in degrees
            # desiredphaseinc = 45 #in degrees

            phase_vect = np.arange(desiredphaseinc, angleforcurvature, desiredphaseinc)
            z_curved_vect = np.round(zmin - ballrad*np.sin(phase_vect*np.pi/180),2)
            y_curved_vect = np.round(ballrad*np.cos(phase_vect*np.pi/180),2)
            
            points_for_curved_section = np.zeros((x_vect_odd.shape[0]*y_curved_vect.shape[0],3))
            aa=0
            for ii in range(x_vect_odd.shape[0]):
                for jj in range(y_curved_vect.shape[0]):
                    # print(aa)
                    points_for_curved_section[aa,0] = x_vect_odd[ii]
                    points_for_curved_section[aa,1] = y_curved_vect[jj]
                    points_for_curved_section[aa,2] = z_curved_vect[jj]
                    aa=aa+1
            
            phase_vect_mirrored = np.arange(desiredphaseinc, angleforcurvature, desiredphaseinc)
            z_curved_vect_mirrored = np.round(z_vect[-1] + ballrad*np.sin(phase_vect_mirrored*np.pi/180),2)
            y_curved_vect_mirrored = np.round(ballrad*np.cos(phase_vect_mirrored*np.pi/180),2)
            
            points_for_curved_section_mirrored = np.zeros((x_vect_odd.shape[0]*y_curved_vect.shape[0],3))
            aa=0
            for ii in range(x_vect_odd.shape[0]):
                for jj in range(y_curved_vect.shape[0]):
                    # print(aa)
                    points_for_curved_section_mirrored[aa,0] = x_vect_odd[ii]
                    points_for_curved_section_mirrored[aa,1] = y_curved_vect_mirrored[jj]
                    points_for_curved_section_mirrored[aa,2] = z_curved_vect_mirrored[jj]
                    aa=aa+1
            
            # points_for_magnets = np.concatenate((points_for_magnets, points_for_curved_section, points_for_curved_section_mirrored), axis=0)
        else:
            points_for_magnets = np.zeros((x_vect_even.shape[0]*1*z_vect.shape[0],3))
            aa=0
            for ii in range(x_vect_even.shape[0]):
                for kk in range(z_vect.shape[0]):
                    # print(aa)
                    points_for_magnets[aa,0] = x_vect_even[ii]
                    points_for_magnets[aa,1] = ycurr
                    points_for_magnets[aa,2] = z_vect[kk]
                    aa=aa+1
                    
            dmin = np.sqrt(zincrement**2 + yincrement**2)
            ballrad = ycurr # setting the radius to the y height of this iteration
            
            desiredphaseinc = np.arcsin(dmin/(2*ballrad))*180/np.pi #in degrees
            # desiredphaseinc = 45 #in degrees

            phase_vect = np.arange(desiredphaseinc, angleforcurvature, desiredphaseinc)
            z_curved_vect = np.round(zmin - ballrad*np.sin(phase_vect*np.pi/180),2)
            y_curved_vect = np.round(ballrad*np.cos(phase_vect*np.pi/180),2)
            
            points_for_curved_section = np.zeros((x_vect_even.shape[0]*y_curved_vect.shape[0],3))
            aa=0
            for ii in range(x_vect_even.shape[0]):
                for jj in range(y_curved_vect.shape[0]):
                    # print(aa)
                    points_for_curved_section[aa,0] = x_vect_even[ii]
                    points_for_curved_section[aa,1] = y_curved_vect[jj]
                    points_for_curved_section[aa,2] = z_curved_vect[jj]
                    aa=aa+1
            
            phase_vect_mirrored = np.arange(desiredphaseinc, angleforcurvature, desiredphaseinc)
            z_curved_vect_mirrored = np.round(z_vect[-1] + ballrad*np.sin(phase_vect_mirrored*np.pi/180),2)
            y_curved_vect_mirrored = np.round(ballrad*np.cos(phase_vect_mirrored*np.pi/180),2)
            
            points_for_curved_section_mirrored = np.zeros((x_vect_even.shape[0]*y_curved_vect.shape[0],3))
            aa=0
            for ii in range(x_vect_even.shape[0]):
                for jj in range(y_curved_vect.shape[0]):
                    # print(aa)
                    points_for_curved_section_mirrored[aa,0] = x_vect_even[ii]
                    points_for_curved_section_mirrored[aa,1] = y_curved_vect_mirrored[jj]
                    points_for_curved_section_mirrored[aa,2] = z_curved_vect_mirrored[jj]
                    aa=aa+1
            
        points_for_magnets = np.concatenate((points_for_magnets, points_for_curved_section, points_for_curved_section_mirrored), axis=0)

        total_positions = np.concatenate((total_positions,points_for_magnets), axis=0) 
    
    total_positions = np.delete(total_positions, 0, axis=0)
    position_ledger = pd.DataFrame(total_positions)
    position_ledger.columns =['X-pos', 'Y-pos', 'Z-pos']
    position_ledger['X-rot'] = np.zeros((total_positions.shape[0],1))
    position_ledger['Y-rot'] = np.zeros((total_positions.shape[0],1))
    position_ledger['Z-rot'] = np.zeros((total_positions.shape[0],1))
    position_ledger['Searched'] = np.zeros((total_positions.shape[0],1))
    position_ledger['CostValue'] = np.ones((total_positions.shape[0],1))*1e6
    position_ledger['Used'] = np.zeros((total_positions.shape[0],1))
    position_ledger['Placement_index'] = np.zeros((total_positions.shape[0],1))
    position_ledger['Bmag'] = np.zeros((total_positions.shape[0],1))
    position_ledger['Magnet_length'] = np.ones((total_positions.shape[0],1))*magnetlength
    position_ledger['Tag'] = 'planar'

# #LA validate this...taking advantage of symmetry thus only taking into account points that are negative X and Z

#     indx_rows =position_ledger.loc[(position_ledger["X-pos"] <=0) & (position_ledger["Z-pos"] <= 0)].index
#     position_ledger=position_ledger.drop(indx_rows)
    
    if(plotting==1):
        fig = plt.figure()
        ax = plt.axes(projection='3d')           
        ax.scatter3D(total_positions[:,0], total_positions[:,1], total_positions[:,2], c=1000*np.ones((total_positions.shape[0],1)));
        
    if(tofile!=None):
        position_ledger.to_excel(tofile, index=False)
        
    print('total magnet positions: ' +str(position_ledger.shape[0]) )
    return position_ledger


#generate list magnet positions in a ring with radius r and z position.
def generate_ring_of_magnets(r,z,cube_side_length,number_mags_aximuthal,theta_offset,tag):
    
        # theta_delta = np.mean(theta[0:2])
    point_list = []
    
    mintheta = 360/number_mags_aximuthal
    theta = np.linspace(0,360-mintheta,number_mags_aximuthal)*np.pi/180 + theta_offset
    
    if(mintheta < 2*np.arcsin(((np.sqrt(3)*cube_side_length)/2)/r)*180/np.pi):
        sys.exit('angle between two different magnets too small, decrease number of magnets azimuthally')
        
    for jj in theta:
        # print(jj)
        x= r*np.cos(jj)
        y= r*np.sin(jj)
        rotx = 0
        roty = 0
        rotz = 2*jj*180/np.pi
        point_list.append([x,y,z,rotx,roty,rotz,0,0,0,0,0,cube_side_length,tag])
    
    return point_list



# creates points for halbach simulations. rmin is the minimum magnet radius, h is the height, layers is the number of concentric layers, cube_size_length is the length of the cube, all in mm. plotting is if you want to plot the points in 3D space and to file is to save the information into ledger format in excel
def define_magnet_positions_halbach(rmin,
                                    number_mags_z,
                                    dist_z,
                                    number_mags_aximuthal,
                                    cube_side_length, 
                                    tag, 
                                    plotting=True,
                                    tofile='halbach_positions.xlsx'):

    # rmin=150
    # des_margin=0
    layers = number_mags_aximuthal.shape[0]
    # cube_side_length = 25.4
    # d = round(np.sqrt(3)*cube_side_length,4) 
       
    # mintheta = 2*np.arcsin((d/2)/rmin)*180/np.pi
    d = dist_z
            
    zloc_odd = np.round(np.r_[-number_mags_z*d/2+d:number_mags_z*d/2:d],4)
    zloc_even = np.round(np.r_[-number_mags_z*d/2+d/2+d:number_mags_z*d/2-d/2:d],4)
    
    # theta_delta = np.mean(theta[0:2])
    point_list = []
    
    for layer in range(layers):
        mintheta = 360/number_mags_aximuthal[layer]
        
        if(mintheta < 2*np.arcsin(((np.sqrt(3)*cube_side_length)/2)/rmin)*180/np.pi):
            sys.exit('angle between two different magnets too small, decrease number of magnets azimuthally')
            
        r = rmin+layer*d #each layer, the radius increases by the minimum distance increase of d
        
        if(np.mod(layer,2)==0):
            for z in zloc_odd:
                print(z)
                tmplist = generate_ring_of_magnets(r,z,cube_side_length,int(number_mags_aximuthal[layer]),0,tag+'_layer_'+str(layer))
                point_list = point_list + tmplist
        else:
            for z in zloc_even:
                print(z)
                tmplist = generate_ring_of_magnets(r,z,cube_side_length,int(number_mags_aximuthal[layer]),0,tag+'_layer_'+str(layer))
                point_list = point_list + tmplist
                      
    ledger = pd.DataFrame(point_list,columns =['X-pos', 'Y-pos', 'Z-pos','X-rot','Y-rot','Z-rot','Searched','CostValue','Used','Placement_index','Bmag','Magnet_length','Tag'])


    cols_to_round = ['X-pos', 'Y-pos', 'Z-pos','X-rot','Y-rot','Z-rot','Searched','CostValue','Used','Placement_index','Bmag','Magnet_length']
    ledger[cols_to_round] = ledger[cols_to_round].round(3)

    if(plotting==True):
        fig = plt.figure()
        ax = plt.axes(projection='3d')        
        large_mags = ledger.loc[ledger['Magnet_length'] == cube_side_length].to_numpy()
        
        ax.scatter3D(large_mags[:,0], large_mags[:,1], large_mags[:,2], c=1000*np.ones((large_mags.shape[0],1)));

        
    if(tofile!=None):
        ledger.to_excel(tofile, index=False)
        
    print('total magnet positions: ' +str(ledger.shape[0]) )
    return ledger


def add_endrings(ledger,number_endrings,mags_per_endring,cube_side_length,tofile):
    if(np.mod(number_endrings,2) != 0):
        print('uneven endrings...adding one')
        number_endrings = number_endrings+1

    tmp=ledger.copy()
    tmp['R']=np.round(np.sqrt(tmp['Y-pos'].to_numpy()**2+tmp['X-pos'].to_numpy()**2),3)
    d = round(np.sqrt(3)*cube_side_length,4) 
    r = tmp['R'].max()+d
    z_s=np.sort(ledger['Z-pos'].unique())
    z_one_end = z_s[:int(number_endrings):2]
    z_second_end = z_s[-1:-int(number_endrings):-2]
    z_arr  = np.concatenate((z_one_end, z_second_end), axis=0)
    
    point_list=[]
    for z in z_arr:
        endring = generate_ring_of_magnets(r,z,cube_side_length,mags_per_endring,0,'shims')
        point_list = point_list + endring

    
    tmp = pd.DataFrame(point_list,columns =['X-pos', 'Y-pos', 'Z-pos','X-rot','Y-rot','Z-rot','Searched','CostValue','Used','Placement_index','Bmag','Magnet_length','Tag'])
    
    
    cols_to_round = ['X-pos', 'Y-pos', 'Z-pos','X-rot','Y-rot','Z-rot','Searched','CostValue','Used','Placement_index','Bmag','Magnet_length']
    tmp[cols_to_round] = tmp[cols_to_round].round(3)
    
    ledg = pd.concat([ledger, tmp])
    ledg = ledg.reset_index(drop=True) #needed correction since rows are dropped

    if(tofile!=None):
        ledg.to_excel(tofile, index=False)
            
    return ledg

#function that takes the ledger as input and plots in 3D the points in the ledger
def plot_3D_points_ledger(ledger):
    fig = plt.figure()
    ax = plt.axes(projection='3d')           
    ax.scatter3D(ledger['X-pos'].to_numpy(), ledger['Y-pos'].to_numpy(), ledger['Z-pos'].to_numpy(),        
                 c=1000*np.ones((ledger.shape[0],1)));
    
#this function laods the magnet positions from an xslx file to a ledger format
def load_magnet_positions(filename,mag_constant=[0,1270,0]):
        
       ledger = pd.read_excel(filename, 
                  dtype={'X-pos': np.float64, 
                         'Y-pos': np.float64, 
                         'Z-pos': np.float64, 
                         'X-rot': np.float64, 
                         'Y-rot': np.float64,
                         'Z-rot': np.float64,
                         'Searched': np.float64,
                         'CostValue': np.float64,
                         'Used': np.float64,
                         'Placement_index': np.float64,
                         'Bmag': np.float64,
                         'Magnet_length': np.float64,
                         'Tag': str})   
       
       ledger.columns = ['X-pos','Y-pos','Z-pos','X-rot','Y-rot','Z-rot','Searched','CostValue','Used','Placement_index','Bmag','Magnet_length','Tag']
       
       mags = magpy.Collection(style_label='magnets')
       cubelist = [None] * ledger.shape[0] #multiplied times two since at each iteration a maximum of 2 magnets can be placed...

       ledger = ledger.reset_index(drop=True) #needed correction since rows are dropped

       aa=0
       for idx_ledger, ledger_row in ledger.iterrows():
           # print(idx_ledger)
           mag_length = ledger.at[idx_ledger,'Magnet_length']
           cubelist[aa] = magpy.magnet.Cuboid(magnetization=mag_constant, dimension=(mag_length,mag_length,mag_length))
           cubelist[aa].position= (ledger.at[idx_ledger,'X-pos'],
                           ledger.at[idx_ledger,'Y-pos'],
                           ledger.at[idx_ledger,'Z-pos'])
                    
           cubelist[aa].orientation = R.from_euler('zyx', [ledger.at[idx_ledger,'Z-rot'],
                                                           ledger.at[idx_ledger,'Y-rot'],
                                                           ledger.at[idx_ledger,'X-rot']], degrees=True)
           
           mags.add(cubelist[aa])
           aa=aa+1
           
       return ledger, mags

def define_sensor_points_on_sphere(num_pts,r,base_coord):
    indices = np.arange(0, num_pts, dtype=float) + 0.5
    
    phi = np.arccos(1 - 2*indices/num_pts)
    theta = np.pi * (1 + 5**0.5) * indices
    x, y, z = r*np.cos(theta) * np.sin(phi), r*np.sin(theta) * np.sin(phi), r*np.cos(phi);
    
    # pp.figure().add_subplot(111, projection='3d').scatter(x, y, z);
    # pp.show()

    arr = np.stack((x,y,z),axis=1)
    arr = np.vstack([arr, np.array([0,0,0])])
    
    arr[:,0] = arr[:,0]+base_coord[0]
    arr[:,1] = arr[:,1]+base_coord[1]
    arr[:,2] = arr[:,2]+base_coord[2]

    sensor = magpy.Sensor(position=arr,style_size=2)
    return sensor

# sensor points on the surface of a cylinder along z (radius r, length h): rings on the side wall plus sunflower-spread points on both end caps.
# the field extremes of a source-free region lie on its boundary, so the surface is enough to evaluate homogeneity of the whole cylinder
def define_sensor_points_on_cylinder(num_pts,r,h,base_coord):
    side_area = 2*np.pi*r*h
    cap_area = np.pi*r**2
    n_side = int(num_pts*side_area/(side_area+2*cap_area))
    n_cap = (num_pts-n_side)//2

    n_rings = max(2, int(round(np.sqrt(n_side*h/(2*np.pi*r)))))
    n_per_ring = max(4, n_side//n_rings)
    theta = np.linspace(0, 2*np.pi, n_per_ring, endpoint=False)
    zs = np.linspace(-h/2, h/2, n_rings)
    side = np.array([(r*np.cos(t+0.5*k*np.pi/n_per_ring), r*np.sin(t+0.5*k*np.pi/n_per_ring), z) for k, z in enumerate(zs) for t in theta])

    idx = np.arange(n_cap) + 0.5
    rr = r*np.sqrt(idx/n_cap)
    tt = np.pi*(1 + 5**0.5)*idx
    cap = np.stack((rr*np.cos(tt), rr*np.sin(tt), np.zeros(n_cap)), axis=1)
    arr = np.vstack([side, cap + [0,0,h/2], cap - [0,0,h/2], np.array([0,0,0])])
    arr = arr + np.array(base_coord)

    sensor = magpy.Sensor(position=arr,style_size=2)
    return sensor

def define_sensor_points_on_filled_sphere(num_pts,r,num_r_points,base_coord):
    r_points = np.linspace(0.1,r,num_r_points)

    for rr in r_points:
        # print(rr)
        num_pts_each_rad = np.ceil(num_pts*rr / r) #normalizing number of points for each radius based on the ratio of the radius
        indices = np.arange(0, num_pts_each_rad, dtype=float) + 0.5
        
        phi = np.arccos(1 - 2*indices/num_pts_each_rad)
        theta = np.pi * (1 + 5**0.5) * indices
        x, y, z = rr*np.cos(theta) * np.sin(phi), rr*np.sin(theta) * np.sin(phi), rr*np.cos(phi);
        
        arr = np.stack((x,y,z),axis=1)

        # pp.figure().add_subplot(111, projection='3d').scatter(x, y, z);
        # pp.show()
        # print(arr)
        if(rr==0.1):
            arr_full = arr
        else:
            arr_full = np.vstack([arr_full, arr])
    
    arr_full[:,0] = arr_full[:,0]+base_coord[0]
    arr_full[:,1] = arr_full[:,1]+base_coord[1]
    arr_full[:,2] = arr_full[:,2]+base_coord[2]

    sensor = magpy.Sensor(position=arr_full,style_size=2)
    return sensor
#%%  
# function takes the magnet position ledger and rotation space ledger and compuates all rotations for all spaces and returns an aggregated ledger for the search                      
def setup_permutationmatrix(active_magnet_ledger,rot_search_space):
    aa=0
    pdrot_search_space = rot_search_space.to_numpy()
    temparr= np.zeros([active_magnet_ledger.shape[0]*rot_search_space.shape[0],11])
    for ii in range(active_magnet_ledger.shape[0]):
        xpos = active_magnet_ledger['X-pos'].iloc[ii]
        ypos = active_magnet_ledger['Y-pos'].iloc[ii]
        zpos = active_magnet_ledger['Z-pos'].iloc[ii] 
        
        for jj in range(pdrot_search_space.shape[0]):
    
            temparr[aa,0] = xpos
            temparr[aa,1] = ypos
            temparr[aa,2] = zpos       
            
            temparr[aa,3] = pdrot_search_space[jj,0]
            temparr[aa,4] = pdrot_search_space[jj,1]
            temparr[aa,5] = pdrot_search_space[jj,2]
            aa=aa+1
        if(ii%200==0):
            print('setting permutation matrix: running index:' + str(ii) + ' out of:' + str(active_magnet_ledger.shape[0]))
            
    
    multiplicative_ledger = pd.DataFrame(temparr,columns =['X-pos', 'Y-pos', 'Z-pos','X-rot','Y-rot','Z-rot','Searched','CostValue','Used','Placement_index','Bmag','Magnet_length','Tag'])
    #**
    multiplicative_ledger['Magnet_length'] = active_magnet_ledger['Magnet_length'].iloc[0]
    multiplicative_ledger['Tag'] = active_magnet_ledger['Tag'].iloc[0]

    return multiplicative_ledger


# this function computes the homogeneity and mean B0 field given a set of sensor locations
def cost_func(B,component=0):
    # Bcomponent=0 for X, Bcomponent=1 for Y and Bcomponent=2 for Z
    Bcomp_flat = B[:,component]   
    # Bcomp_flat = Bcomp.flatten()
    meanfield = np.mean(B[:,component]) 
    eta = 1e6*np.abs((np.max(Bcomp_flat)-np.min(Bcomp_flat))/(meanfield))  
    return eta, meanfield

# this function computes the homogeneity and mean B0 field given a set of sensor locations
def cost_func_one_dim(B):
    # Bcomponent=0 for X, Bcomponent=1 for Y and Bcomponent=2 for Z
    Bcomp_flat = B
    # Bcomp_flat = Bcomp.flatten()
    meanfield = B[-1] # last element is the field at the center, defined by another function.
    eta = 1e6*np.abs((np.max(Bcomp_flat)-np.min(Bcomp_flat))/(meanfield))  
    return eta, meanfield


# function used to compute the B field with parallelization option...used to speed up computation
# TO DO - implement the acceptance of collection of sensors and collection of magnets.

def task_Bcalc(posx,posy,posz,rotx,roty,rotz,idx,maxidx,col_sensors,magnets,mag_constant,mag_cube_dimension,Bcomponent):
     
      
    magnet = magpy.magnet.Cuboid(magnetization=mag_constant, dimension=mag_cube_dimension)
    magnet.position= (posx,posy,posz)
    magnet.orientation = R.from_euler('zyx', [rotz,roty,rotx], degrees=True)
    
    magnets.add(magnet)
    # # computing B
    B = col_sensors.getB(magnets)
    # # computing cost function for only one magnet
    
    eta, meanfield = cost_func(B,Bcomponent)

    # if( idx%500 == 0):
    #     print('task calcB running index:' + str(idx) + ' out of:' + str(maxidx))
    #cleanup
    magnets.remove(magnet)

    return eta, meanfield,idx

#computes all possible rotation possibilities for each magnet
def generate_local_rotation_combinations(num_angles_x,
                                         num_angles_y,
                                         num_angles_z,
                                         min_angle_x,
                                         max_angle_x,
                                         min_angle_y,
                                         max_angle_y,
                                         min_angle_z,
                                         max_angle_z):
    
    print('generating all magnet possible rotation permutations')

    rot_angles_x = np.linspace(min_angle_x,max_angle_x,num_angles_x)
    rot_angles_y = np.linspace(min_angle_y,max_angle_y,num_angles_y)
    rot_angles_z = np.linspace(min_angle_z,max_angle_z,num_angles_z)

    rotation_combinations = np.array(list(itertools.product(rot_angles_x, rot_angles_y, rot_angles_z)))  
    
    df_rotation_combinations = pd.DataFrame(rotation_combinations)
    df_rotation_combinations.columns =['X-rot','Y-rot','Z-rot']
    return df_rotation_combinations

#once the results are obtained for all the rotations, this function picks the best rotation for each location in space since magnets cannot overlap in space.
def consoledate_BmagMax(active_magnet_ledger,results_ledger,rot_search_space):
    
    updated_ledger_np = np.zeros(active_magnet_ledger.shape)
    tmp = results_ledger.to_numpy()
    
    Bs = results_ledger['Bmag'].to_numpy()
    Bs = np.reshape(Bs,(active_magnet_ledger.shape[0],rot_search_space.shape[0]))
    indxmax = np.argmax(Bs,axis=1)
    
    for ii in range(updated_ledger_np.shape[1]):
        # print(ii)
        coldata = tmp[:,ii]
        coldata = np.reshape(coldata,(active_magnet_ledger.shape[0],rot_search_space.shape[0]))
        opt_col_data = coldata[np.arange(Bs.shape[0]),list(indxmax)]
        updated_ledger_np[:,ii] = opt_col_data
#**
    res_ledg = pd.DataFrame(updated_ledger_np,columns =['X-pos', 'Y-pos', 'Z-pos','X-rot','Y-rot','Z-rot','Searched','CostValue','Used','Placement_index','Bmag','Magnet_length','Tag'])

    return res_ledg

# function that computes the mean B0 field, given the number of magnets specified
def calc_mean_field(numberofmagnets, res_ledg):
    df = res_ledg['Bmag'].sort_values(ascending=False).iloc[0:numberofmagnets]
    sumBfield = df.sum()
    return df, sumBfield

# plots the field strngth per number of magnets iteratively. this function can be used to decided how many magnets are to be used in the design.
def plot_mean_field_vs_magnet_number(ledger,maxmagnets):
    
    nummags =  np.zeros((maxmagnets,1))
    fieldstrength =  np.zeros((maxmagnets,1))
    
    for ii in range(1,maxmagnets):
        nummags[ii]=ii
        tt, fieldstrength[ii] = calc_mean_field(ii, ledger)
        
    plt.figure()
    plt.plot(nummags,fieldstrength)

    
# function used to visualize in 3D the magnet locations alongside the position of the sensors. 
def plot_magnets3D(ledger,localsensors=None,mag_constant=[0,1270,0]):
    
    magnets = magpy.Collection(style_label='magnets')
    cubelist = [None] * ledger.shape[0] #multiplied times two since at each iteration a maximum of 2 magnets can be placed...

    ledger = ledger.reset_index(drop=True) #needed correction since rows are dropped

    aa=0
    for idx_ledger, ledger_row in ledger.iterrows():
        print(idx_ledger)
#**
        maglength = ledger.at[idx_ledger,'Magnet_length']
        
        cubelist[aa] = magpy.magnet.Cuboid(magnetization=mag_constant, dimension=(maglength,maglength,maglength))
        cubelist[aa].position= (ledger.at[idx_ledger,'X-pos'],
                        ledger.at[idx_ledger,'Y-pos'],
                        ledger.at[idx_ledger,'Z-pos'])
        
        cubelist[aa].orientation = R.from_euler('zyx', [ledger.at[idx_ledger,'Z-rot'],ledger.at[idx_ledger,'Y-rot'],ledger.at[idx_ledger,'X-rot']], degrees=True)
        
        magnets.add(cubelist[aa])
        aa=aa+1
    if(localsensors==None):
        magpy.show(magnets)
    else:
        magpy.show(magnets, localsensors)

#simulate and plot the B field from an arbitrary number of magnets. uses the ledger as an input to define which magnets to plot. It returns the homogeneity, meanB0 and magnet collection to be used by other functions (e.g. when optimizing homogeneity)
def simulate_ledger(mags, 
                     col_sensors,
                     mag_constant, 
                     ledger,
                     cubexi,
                     meshelements_percube, 
                     extrinsicrot=True, 
                     bool_runmeshingfordemag=False,
                     backgroundB = None,
                     bool_plotting2D=False):
    
    Bcomponent=np.argmax(mag_constant)
    col_magnet = mags.copy()
    df = ledger.copy()

    cubelist = []
    
    for df_ledger, df_row in df.iterrows():
        # print(df_ledger)
        # print(df_row)
#**
        cube = magpy.magnet.Cuboid(magnetization=mag_constant, dimension=ledger_dimension(df_row))
        # print(df_row)
        cube.position= (df_row['X-pos'],
                        df_row['Y-pos'],
                        df_row['Z-pos'])
        cube.xi = cubexi # µr=1 +cube.xi
        # cube.orientation = R.from_rotvec((df_row['X-rot'],
        #                                   df_row['Y-rot'],
        #                                   df_row['Z-rot']), 
        #                                   degrees=True)
        if(extrinsicrot==True): #extrinsic=global
            cube.orientation = R.from_euler('zyx', [df_row['Z-rot'],df_row['Y-rot'],df_row['X-rot']], degrees=True) #global rot
        else:
            cube.orientation = R.from_euler('ZYX', [df_row['Z-rot'],df_row['Y-rot'],df_row['X-rot']], degrees=True)
        

        cubelist.append(cube)
        
 #iterate over elements of magnet list       
    for i in cubelist:
        col_magnet.add(i)

    
    if(bool_plotting2D==True):
        magpy.show(col_magnet, col_sensors)
    
    # B = col_sensors.getB(col_magnet)
    # print(B)

    if(bool_runmeshingfordemag==True):
        coll_meshed = mesh_all(col_magnet, target_elems=meshelements_percube, per_child_elems=True, min_elems=1)    
        coll_demag = apply_demag(coll_meshed,style={"label": f"Coll_demag ({len(coll_meshed.sources_all):3d} cells)"},)
        B = magpy.getB(coll_demag, col_sensors)
        # print(B)
    else:
        B = col_sensors.getB(col_magnet)

    if backgroundB is not None:
        B = B + backgroundB
    
    eta, meanB0 = cost_func(B,Bcomponent)
    
    
    if(bool_plotting2D==True):
        # create an observer grid in the xz-symmetry plane, around 0
        ts = np.linspace(-10, 10, 50) 
        grid = np.array([[(x,0,z) for x in ts] for z in ts])
        B_grid = col_magnet.getB(grid)
    
        plt.figure()
        plt.imshow(B_grid[:,:,Bcomponent], extent=[-10, 10, -10, 10])
        plt.colorbar()
        tmpB = B_grid[:,:,Bcomponent]
        tmpB = tmpB[-1]
        match Bcomponent:
            case 0:
                strcomponent = 'x'
            case 1:
                strcomponent = 'y'            
            case 2:
                strcomponent = 'z'
                
        plt.title('Magnitude of B'+strcomponent + ' mean:' + str(np.round(np.mean(tmpB),5)) + ' std:' + str(np.round(np.std(tmpB),6)))         
        plt.show()

    return eta, meanB0, col_magnet, B


def movefromoneledgertoanother(target_ledger,source_ledger,column,condition):
    
    cond = source_ledger[column] ==condition
    rows = source_ledger.loc[cond, :]
    target_ledger = pd.concat([target_ledger, pd.DataFrame(rows)])
    target_ledger.to_excel('interim_used_ledger.xlsx', index=False)
    source_ledger.drop(rows.index, inplace=True)
    
    return source_ledger, target_ledger

def write_array_to_ledger_format(arr,arrcolumns=['X-pos', 'Y-pos', 'Z-pos','X-rot','Y-rot','Z-rot','Searched','CostValue','Used','Placement_index','Bmag','Magnet_length'],fl='ledger.xlsx'):
    
    ledger = pd.DataFrame(arr,columns=arrcolumns)

    ledger.to_excel(fl, index=False)
    return ledger

#function that does the mirror operation for a vector. can do this across 'yz', 'xy' or 'xz' planes
def mirror_op(arr,plane):
    rotx=arr[0]
    roty=arr[1]
    rotz=arr[2]
    
    RotOp =  R.from_euler('zyx',[rotz,rotx,roty],degrees=True).as_quat()
    # print(RotOp)
    mirrorplane = np.array([0,0,0,0])
    if(plane=='yz'):
        mirrorplane = np.array([RotOp[0], -RotOp[1],-RotOp[2],RotOp[3]])
    if(plane=='xy'):
        mirrorplane = np.array([-RotOp[0], -RotOp[1],RotOp[2],RotOp[3]])
    if(plane=='xz'):
        mirrorplane = np.array([-RotOp[0], RotOp[1], -RotOp[2],RotOp[3]])
        
    RRR= R.from_quat(mirrorplane)
    # print(mirrorplane)
    return RRR

def populate_full_geom_from_quad(ledger_quad,ledger_zeros_):
    # need to implement 3 mirrors
    
    # first across x
    ledger_x_mirror = ledger_quad.copy()
    ledger_x_mirror['X-pos'] = -ledger_quad['X-pos']
    
    tmp = ledger_x_mirror[['X-rot','Y-rot','Z-rot']].to_numpy()
    
    for ii in range(tmp.shape[0]):
        Rmirrored = mirror_op(tmp[ii,:],'yz').as_euler('xyz',degrees=True)
        tmp[ii,:] = Rmirrored
    
    ledger_x_mirror[['X-rot','Y-rot','Z-rot']]=tmp
        
    ledger_mirror1 = pd.concat([ledger_quad, ledger_x_mirror])
    ledger_mirror1=ledger_mirror1.drop_duplicates()
    ledger_mirror1 = ledger_mirror1.reset_index(drop=True) #needed correction since rows are dropped

    #  across y

    ledger_y_mirror = ledger_mirror1.copy()
    ledger_y_mirror['Y-pos'] = -ledger_mirror1['Y-pos']
    
    tmp = ledger_y_mirror[['X-rot','Y-rot','Z-rot']].to_numpy()
    
    for ii in range(tmp.shape[0]):
        Rmirrored = mirror_op(tmp[ii,:],'xz').as_euler('xyz',degrees=True)
        tmp[ii,:] = Rmirrored
    ledger_y_mirror[['X-rot','Y-rot','Z-rot']]=tmp

        
    ledger_mirror2 = pd.concat([ledger_mirror1, ledger_y_mirror])
    ledger_mirror2=ledger_mirror2.drop_duplicates()
    ledger_mirror2 = ledger_mirror2.reset_index(drop=True) #needed correction since rows are dropped

    
    # #  across z

    ledger_z_mirror = ledger_mirror2.copy()
    ledger_z_mirror['Z-pos'] = -ledger_mirror2['Z-pos']
    
    tmp = ledger_z_mirror[['X-rot','Y-rot','Z-rot']].to_numpy()
    
    for ii in range(tmp.shape[0]):
        Rmirrored = mirror_op(tmp[ii,:],'xy').as_euler('xyz',degrees=True)
        tmp[ii,:] = Rmirrored
    
    ledger_z_mirror[['X-rot','Y-rot','Z-rot']]=tmp
    
    ledger_mirror3 = pd.concat([ledger_mirror2, ledger_z_mirror])
    ledger_mirror3=ledger_mirror3.drop_duplicates()
    ledger_mirror3 = ledger_mirror3.reset_index(drop=True) #needed correction since rows are dropped

    new_ledg = pd.concat([ledger_mirror3, ledger_zeros_])
    new_ledg = new_ledg.reset_index(drop=True) 
    
    rotx = new_ledg['X-rot'].to_numpy()
    new_ledg['X-rot'] = np.mod(rotx, 360)

    roty = new_ledg['Y-rot'].to_numpy()
    new_ledg['Y-rot'] = np.mod(roty, 360)
    
    rotz = new_ledg['Z-rot'].to_numpy()
    new_ledg['Z-rot'] = np.mod(rotz, 360)
    
    return new_ledg

def magpylib_angles_to_opera(rotx,roty,rotz):
    Rotobj = R.from_euler('xyz', [rotx,roty,rotz], degrees=True)
    ZYZrot = Rotobj.as_euler('ZYZ',degrees=True)
    return ZYZrot

#uses an input of the rotation in Z, Y, and Z (intrinsic) and returns a 3x1 angles in degrees corresponding to the global xyz rotations used in magpylib

def opera_angles_to_magpylib(rotz,roty,rotzp):
    Rotobj = R.from_euler('ZYZ', [rotz,roty,rotzp], degrees=True)
    xyzrot = Rotobj.as_euler('xyz',degrees=True)
    return xyzrot


def add_colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar


# extracts and returns a dictionary of the data
# Bfield - the magnetic field extracted
# xmin - minimum x coordinate
# xmax - maximum x coordinate
# ymin - minimum y coordinate
# ymax - maximum y coordinate
# zmin - minimum z coordinate
# zmax - maximum z coordinate
# numberpoints_per_ax -  number of points in each axis. e.g 50 will generate a 50x50x50x3 matrix
# filename - name of filename being saved.
def extract_3Dfields(mags,xmin=-100,xmax=100,ymin=-100,ymax=100,zmin=-100, zmax=100, numberpoints_per_ax = 51,filename='output.pkl',plotting=True,Bcomponent=0):
    
    xs = np.linspace(xmin, xmax, numberpoints_per_ax)
    ys = np.linspace(ymin, ymax, numberpoints_per_ax)
    zs = np.linspace(zmin,zmax,numberpoints_per_ax)
    
    x_mesh, y_mesh, z_mesh=np.meshgrid(xs,ys,zs)
    coordinates = [(a2, b2, c2,) for a, b, c in zip(x_mesh, y_mesh, z_mesh) for a1, b1, c1 in zip(a, b, c) for a2, b2, c2 in zip(a1, b1, c1)]

        
    B = np.zeros((xs.shape[0],ys.shape[0],zs.shape[0],3))
    
    for kk in range(zs.shape[0]):
        print('extracting slice:' + str(kk+1) + ' out of:' + str(zs.shape[0]))

        grid_xy = np.array([[(x,y,zs[kk]) for x in xs] for y in ys])

        B[:,:,kk,:] = mags.getB(grid_xy)
                    
    data = {'Bfield': B,'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'zmin': zmin, 'zmax': zmax, 'numberpoints_per_ax': numberpoints_per_ax, 'filename': filename, 'coordinates':coordinates}
    
    if(filename != None):
        fileObj = open(filename, 'wb')
        pickle.dump(data,fileObj)
        fileObj.close()
    
    if(plotting==True):
        B_mT = B
        
        # minmin = np.min(Bcomponent.flatten())
        # maxmax = np.max(Bcomponent.flatten())
        meanmean = np.mean(B_mT[:,:,:,Bcomponent].flatten())
        stdstd = np.std(B_mT[:,:,:,Bcomponent].flatten())

        fig, (ax1, ax2, ax3) = plt.subplots(3)
        im1=ax1.imshow(B_mT[:,:,int(zs.shape[0]/2),Bcomponent], vmin=meanmean-2*stdstd, vmax=meanmean+2*stdstd, cmap='jet', aspect='auto', extent=[xmin, xmax, ymin, ymax])
        add_colorbar(im1)
        ax1.title.set_text('xy plane')
        
        
        im2=ax2.imshow(B_mT[:,int(ys.shape[0]/2),:,Bcomponent], vmin=meanmean-2*stdstd, vmax=meanmean+2*stdstd, cmap='jet', aspect='auto', extent=[xmin, xmax, zmin, zmax])
        add_colorbar(im2)
        ax2.title.set_text('xz plane')
        im3=ax3.imshow(B_mT[int(xs.shape[0]/2),:,:,Bcomponent], vmin=meanmean-2*stdstd, vmax=meanmean+2*stdstd, cmap='jet', aspect='auto', extent=[ymin, ymax, zmin, zmax])
        add_colorbar(im3)
        ax3.title.set_text('yz plane')
        
        fig.suptitle('B field (mT). Mean: ' + str(round(meanmean,3)) + ', std: '  + str(round(stdstd,3)) + '\n Mean B field (MHz): ' + str(round(42.580000*meanmean/1000,3)) + ' std (Hz)' + str(round(42580000*stdstd/1000,3)))
        fig.show()
        
        
    return data


   
def data_to_matlab(picklefilename):
    data = pickle.load(open(picklefilename, "rb"))
    scipy.io.savemat(picklefilename[:-3]+'mat', mdict=data)
    
def data_from_matlab(picklefilename):
    arr = scipy.io.loadmat(picklefilename)
    return arr

#this assumes rotations in X,Y and Z
def get_max_magnets_per_radius(r,magnet_length):
    d = 2*np.arcsin(((np.sqrt(3)*magnet_length)/2)/r)*180/np.pi #result in degrees
    n=np.floor(360/d)
    return n

# remember B field is inputted in mT, not T
def plot_3D_field(B_mT,Bcomponent=0,xmin=-100,xmax=100,ymin=-100,ymax=100,zmin=-100, zmax=100):
    (nx,ny,nz) = B_mT.shape
    # minmin = np.min(Bcomponent.flatten())
    # maxmax = np.max(Bcomponent.flatten())
    meanmean = np.mean(B_mT[:,:,:,Bcomponent].flatten())
    stdstd = np.std(B_mT[:,:,:,Bcomponent].flatten())

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    im1=ax1.imshow(B_mT[:,:,int(nz/2),Bcomponent], vmin=meanmean-2*stdstd, vmax=meanmean+2*stdstd, cmap='jet', aspect='auto', extent=[xmin, xmax, ymin, ymax])
    add_colorbar(im1)
    ax1.title.set_text('xy plane')
    
    
    im2=ax2.imshow(B_mT[:,int(ny/2),:,Bcomponent], vmin=meanmean-2*stdstd, vmax=meanmean+2*stdstd, cmap='jet', aspect='auto', extent=[xmin, xmax, zmin, zmax])
    add_colorbar(im2)
    ax2.title.set_text('xz plane')
    im3=ax3.imshow(B_mT[int(nx/2),:,:,Bcomponent], vmin=meanmean-2*stdstd, vmax=meanmean+2*stdstd, cmap='jet', aspect='auto', extent=[ymin, ymax, zmin, zmax])
    add_colorbar(im3)
    ax3.title.set_text('yz plane')

    
    fig.suptitle('B field (mT). Mean: ' + str(round(meanmean,3)) + ', std: '  + str(round(stdstd,3)) + '\n Mean B field (MHz): ' + str(round(42.580000*meanmean/1000,3)) + ' std (Hz)' + str(round(42580000*stdstd/1000,3)))
    fig.show()
    

#%% demagnetization check

# iHc (kA/m) of the grade series in magnet_properties/CY-Mag-NdFeB.pdf, with the max working temperature (C). N50/N52 are only 875 kA/m.
NDFEB_IHC_GRADES = [('N', 955, 80), ('M', 1114, 100), ('H', 1353, 120), ('SH', 1592, 150), ('UH', 1989, 180), ('EH', 2387, 200), ('AH', 2787, 220)]

# magnet dimension from a ledger row: a cube of Magnet_length, or a x a x Magnet_height (height along local z) if that column is set
def ledger_dimension(df_row):
    a = df_row['Magnet_length']
    h = df_row.get('Magnet_height', np.nan)
    return (a, a, a if pd.isna(h) else h)

# builds a magpylib collection from a ledger, the same way simulate_ledger does (no demag meshing)
def build_magnet_collection(ledger, mag_constant=[1270,0,0], extrinsicrot=True):
    col_magnet = magpy.Collection(style_label='magnets')
    for df_ledger, df_row in ledger.iterrows():
        cube = magpy.magnet.Cuboid(magnetization=mag_constant, dimension=ledger_dimension(df_row))
        cube.position = (df_row['X-pos'], df_row['Y-pos'], df_row['Z-pos'])
        if(extrinsicrot==True):
            cube.orientation = R.from_euler('zyx', [df_row['Z-rot'],df_row['Y-rot'],df_row['X-rot']], degrees=True)
        else:
            cube.orientation = R.from_euler('ZYX', [df_row['Z-rot'],df_row['Y-rot'],df_row['X-rot']], degrees=True)
        col_magnet.add(cube)
    return col_magnet

# computes the reverse (demagnetizing) H field each magnet sees: field from all other magnets plus its own self-demagnetization,
# projected on the magnet's magnetization direction. positive = opposing the magnetization. units kA/m (magpylib 4: mm, mT -> H in kA/m).
# H is sampled on a points_per_axis^3 grid spanning edge_fraction of the cube; the H of an ideal sharp-edged cube diverges (logarithmically)
# at its edges, so the sampled max grows slowly as edge_fraction -> 1. the center value is the robust bulk number.
def compute_demag_fields(ledger, mag_constant=[1270,0,0], extrinsicrot=True, points_per_axis=5, edge_fraction=0.9):
    ledger = ledger.reset_index(drop=True)
    col_magnet = build_magnet_collection(ledger, mag_constant, extrinsicrot)
    mag_dir = np.array(mag_constant, dtype=float)/np.linalg.norm(mag_constant)

    u = np.linspace(-0.5, 0.5, points_per_axis)*edge_fraction
    local_grid = np.array(list(itertools.product(u, u, u)))
    local_grid = np.vstack([np.zeros((1,3)), local_grid]) # first point is always the cube center

    H_center = np.zeros(len(ledger))
    H_max = np.zeros(len(ledger))
    for idx, cube in enumerate(col_magnet.children):
        pts = cube.orientation.apply(local_grid*cube.dimension) + cube.position
        H = magpy.getH(col_magnet, pts)
        H_reverse = -H @ cube.orientation.apply(mag_dir)
        H_center[idx] = H_reverse[0]
        H_max[idx] = np.max(H_reverse)

    result = ledger.copy()
    result['Hrev_center_kAm'] = H_center
    result['Hrev_max_kAm'] = H_max
    return result

# minimum iHc (kA/m, at 20 C datasheet conditions) so that the worst reverse field stays below the knee of the intrinsic curve at max_temp_C.
# knee_fraction: Hk/iHc squareness (~0.9 for sintered NdFeB). ihc_temp_coeff: %/C, -0.6 from the CY-Mag datasheet.
def required_ihc(H_reverse_kAm, max_temp_C=40, safety_factor=1.2, knee_fraction=0.9, ihc_temp_coeff=-0.6):
    ihc_at_temp = H_reverse_kAm*safety_factor/knee_fraction
    derating = 1 + ihc_temp_coeff/100*max(max_temp_C-20, 0)
    return ihc_at_temp/derating

def suggest_grade(ihc_needed_kAm, max_temp_C=40):
    for series, ihc, tmax in NDFEB_IHC_GRADES:
        if ihc >= ihc_needed_kAm and tmax >= max_temp_C:
            return series
    return None


#%% tolerance (as-built) homogeneity

# fits the field on the sensor points with the given shim order and returns the residual.
# order 0: nothing removed but the mean, 1: linear gradients (free shim with gradient offsets), 2: + second-order harmonics
def shim_residual(B, pos, order=1):
    x, y, z = pos[:,0], pos[:,1], pos[:,2]
    terms = [np.ones_like(x)]
    if order >= 1:
        terms += [x, y, z]
    if order >= 2:
        terms += [x*x-y*y, 2*z*z-x*x-y*y, x*y, x*z, y*z]
    A = np.stack(terms, axis=1)
    coef = np.linalg.lstsq(A, B, rcond=None)[0]
    return B - A @ coef

# homogeneity of an as-built magnet with random errors on every magnet: Br (fraction), position (mm) and angle of the magnetization (deg).
# the field change is linear in small errors, so each magnet's sensitivities are computed once (6 field evaluations) and n_samples
# random builds are drawn from them. returns (nominal ppm, B0 in mT, percentile ppm of the as-built magnet after shimming with shim_order)
def tolerance_homogeneity(ledger, sensor_pos, mag_constant=[1320,0,0], sigma_Br=0.01, sigma_pos=0.1, sigma_ang_deg=1.0,
                          n_samples=100, percentile=90, shim_order=1, extrinsicrot=True, seed=0):
    col = build_magnet_collection(ledger, mag_constant, extrinsicrot)
    cubes = col.children
    Br = np.linalg.norm(mag_constant)
    comp = int(np.argmax(np.abs(mag_constant)))

    B = magpy.getB(cubes, sensor_pos, sumup=False)[..., comp]                    # (n_mag, n_sens) nominal per-magnet field
    for c in cubes: c.magnetization = (0, Br, 0)
    B_rot1 = magpy.getB(cubes, sensor_pos, sumup=False)[..., comp]               # magnetization tipped towards local y
    for c in cubes: c.magnetization = (0, 0, Br)
    B_rot2 = magpy.getB(cubes, sensor_pos, sumup=False)[..., comp]               # magnetization tipped towards local z
    for c in cubes: c.magnetization = tuple(mag_constant)
    h = 0.01
    grad = [(magpy.getB(cubes, sensor_pos - h*e, sumup=False)[..., comp] - B)/h for e in np.eye(3)]  # moving a magnet by +d = moving sensors by -d

    # sensitivity matrix: rows are error sources scaled by their sigma, columns are sensors
    G = np.concatenate([B*sigma_Br, B_rot1*np.deg2rad(sigma_ang_deg), B_rot2*np.deg2rad(sigma_ang_deg)] + [g*sigma_pos for g in grad], axis=0)
    B_nom = B.sum(axis=0)
    rng = np.random.default_rng(seed)
    samples = B_nom[None, :] + rng.standard_normal((n_samples, G.shape[0])) @ G   # (n_samples, n_sens)

    meanB0 = np.mean(B_nom)
    eta_nom = 1e6*np.abs((B_nom.max()-B_nom.min())/meanB0)
    res = shim_residual(samples.T, sensor_pos, shim_order).T
    eta_built = 1e6*np.abs((res.max(axis=1)-res.min(axis=1))/meanB0)
    return eta_nom, meanB0, np.percentile(eta_built, percentile)
