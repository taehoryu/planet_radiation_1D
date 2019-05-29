#!/usr/bin/env python
# -*- coding: utf-8 -*-



#Nov6 : Centering the position + add vorticity slice plot

import yt
from yt.mods import *
from yt.visualization.base_plot_types import get_multi_plot
from mpl_toolkits.mplot3d import Axes3D

import glob
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib
import matplotlib.colorbar as cb
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
matplotlib.use('Agg')
import matplotlib.ticker as ticker
import os
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from yt.units import kpc
#from yt import derived_field
from matplotlib import pylab
from yt.units import dimensions
np.set_printoptions(threshold=np.nan)
matplotlib.rcParams.update({'font.size': 15})

i=0


#print 'Listing all profile/dat files'
#profilefilelist = sorted(glob.glob('plt*'), key=lambda name: int(name[3:]))
profilefilelist = sorted(glob.glob('smallplt*'), key=lambda name: int(name[8:]))
#raw_input('Press ENTER to continue...')
print profilefilelist
number_plt=0
max_time=0
for i in profilefilelist:
    #  if np.mod(j,every_nthfile)==0:
  ds=yt.load(i)
  #ds.add_fields(["pressure"])
  #print  ds.field_list
  number_plt=number_plt+1
  max_time=np.maximum(max_time,float(ds.current_time))
  if number_plt==1:
    ds.index
    for j in sorted(ds.field_list):
      if j[1]=="rad":
        radiation_field=1
        print "rad"
    my_ray_x = ds.ortho_ray(0,(0,0))
    srt_x=np.argsort(my_ray_x['x'])
    my_ray_y = ds.ortho_ray(1,(0,0))
    srt_y=np.argsort(my_ray_y['y'])
    my_ray_z = ds.ortho_ray(2,(0,0))
    srt_z=np.argsort(my_ray_z['z'])
    print(srt_x,srt_y,srt_z)
    if(int(len(np.array(my_ray_z['z'][srt_z])))>1):
      Dimension=3
      print "1. Dimension= ",Dimension
    elif (int(len(np.array(my_ray_y['y'][srt_y])))>1):
      Dimension=2
      print "2. Dimension= ",Dimension
    else:
      Dimension=1
      print "3. Dimension= ",Dimension
    if Dimension==1:
      number_cell_width=float(len(np.array(my_ray_x['x'][srt_x])))
      number_cell_z=number_cell_width

      y_range=[np.min(np.array(my_ray_x['x'][srt_x])),np.amax(np.array(my_ray_x['x'][srt_x]))]
      T_range=[np.min(np.array(my_ray_x['Temp'][srt_x])),np.amax(np.array(my_ray_x['Temp'][srt_x]))]
      rho_range=[np.min(np.array(my_ray_x['density'][srt_x])),np.amax(np.array(my_ray_x['density'][srt_x]))]
      P_range=[np.min(np.array(my_ray_x['pressure'][srt_x])),np.amax(np.array(my_ray_x['pressure'][srt_x]))]
      index_1mbar=numpy.argmin(numpy.abs(np.array(my_ray_x['pressure'][srt_x]) - 1.e3))
      index_1bar=numpy.argmin(numpy.abs(np.array(my_ray_x['pressure'][srt_x]) - 1.e6))
      z_1mbar=np.array(my_ray_x['x'][srt_x])[index_1mbar]
      z_1bar=np.array(my_ray_x['x'][srt_x])[index_1bar]
      delta_x=(y_range[1]-y_range[0])/number_cell_z
    elif Dimension==2:
      number_cell_z=float(len(np.array(my_ray_y['y'][srt_y])))
      number_cell_width=float(len(np.array(my_ray_x['x'][srt_x])))
      x_range=[np.min(np.array(my_ray_x['x'][srt_x])),np.amax(np.array(my_ray_x['x'][srt_x]))]
      
      y_range=[np.min(np.array(my_ray_y['y'][srt_y])),np.amax(np.array(my_ray_y['y'][srt_y]))]
      T_range=[np.min(np.array(my_ray_y['Temp'][srt_y])),np.amax(np.array(my_ray_y['Temp'][srt_y]))]
      rho_range=[np.min(np.array(my_ray_y['density'][srt_y])),np.amax(np.array(my_ray_y['density'][srt_y]))]
      P_range=[np.min(np.array(my_ray_y['pressure'][srt_y])),np.amax(np.array(my_ray_y['pressure'][srt_y]))]
      index_1mbar=numpy.argmin(numpy.abs(np.array(my_ray_y['pressure'][srt_y]) - 1.e3))
      index_1bar=numpy.argmin(numpy.abs(np.array(my_ray_y['pressure'][srt_y]) - 1.e6))
      z_1mbar=np.array(my_ray_y['y'][srt_y])[index_1mbar]
      z_1bar=np.array(my_ray_y['y'][srt_y])[index_1bar]
      delta_x=(y_range[1]-y_range[0])/number_cell_z
    elif Dimension==3:
      number_cell_height=float(len(np.array(my_ray_z['z'][srt_z])))
      number_cell_width=float(len(np.array(my_ray_y['y'][srt_y])))
      print(number_cell_height,number_cell_width)
      print(np.array(my_ray_x['x'][srt_x]))
      number_cell_x=float(len(np.array(my_ray_x['x'][srt_x])))
      number_cell_z=float(len(np.array(my_ray_z['z'][srt_z])))
      number_cell_y=float(len(np.array(my_ray_y['y'][srt_y])))
      x_range=[np.min(np.array(my_ray_z['x'][srt_z])),np.amax(np.array(my_ray_x['x'][srt_x]))]
      y_range=[np.min(np.array(my_ray_y['y'][srt_y])),np.amax(np.array(my_ray_y['y'][srt_y]))]
      z_range=[np.min(np.array(my_ray_z['z'][srt_y])),np.amax(np.array(my_ray_z['z'][srt_z]))]
      T_range=[np.min(np.array(my_ray_z['Temp'][srt_z])),np.amax(np.array(my_ray_z['Temp'][srt_z]))]
      rho_range=[np.min(np.array(my_ray_z['density'][srt_z])),np.amax(np.array(my_ray_z['density'][srt_z]))]
     # ds.add_fields(["pressure"])
      press=3.55319e7*np.multiply(np.array(my_ray_z['Temp'][srt_z]),np.array(my_ray_z['density'][srt_z]))
      P_range=[np.min(press),np.amax(press)]
      index_1mbar=numpy.argmin(numpy.abs(press - 1.e3))
      index_50mbar=numpy.argmin(numpy.abs(press - 50.e3))
      index_1bar=numpy.argmin(numpy.abs(press - 1.e6))
      index_5mbar=numpy.argmin(numpy.abs(press - 5.e3))
      z_1mbar=np.array(my_ray_z['z'][srt_z])[index_1mbar]
      z_1bar=np.array(my_ray_z['z'][srt_z])[index_1bar]
      z_50mbar=np.array(my_ray_z['z'][srt_z])[index_50mbar]
      z_5mbar=np.array(my_ray_z['z'][srt_z])[index_5mbar]
      delta_x=(z_range[1]+z_range[0])/number_cell_z
      T_1mbar=np.array(my_ray_z['Temp'][srt_z])[index_1mbar]
      if(T_1mbar>=2500.0):
        teddy=1734.0
      else:
        teddy=1343.0
#P_range=[x for x in P_range] #[bar]
#x_range=[x for x in x_range] #[km]
#y_range=[x for x in y_range] #[km]
number_slice_x=number_cell_y
choose_label=int(raw_input('Labels? (0:time or else:set manually) ' ))

maximumtime_array=raw_input('Set Maximum time? (Y or N)  ' )
if maximumtime_array=="Y" or maximumtime_array=="y" or maximumtime_array=="yes":
  max_time2=float(raw_input('maximum time in s : '))
  for i in profilefilelist:
    #  if np.mod(j,every_nthfile)==0:
    ds=load(i)
    if(max_time2<=float(ds.current_time)):
      max_time=float(ds.current_time)
      break
print "The maximum time is set at t= ", max_time, "s"


time_difference_T_P = int(raw_input("time_difference_T_P :"))

print "preliminary step"
print "----------------ATM profile--------------"
print "DIMENSION                : ",Dimension
print "height at P=1mbar   [km] : ",  z_1mbar/1e5
print "height at P= 1bar   [km] : ",  z_1bar/1e5
print "x range             [km] : ",x_range[0]/1e5,x_range[1]/1e5
if Dimension==2:
  print "y range           [km] : ",y_range[0]/1e5,y_range[1]/1e5
if Dimension==3:
  print "y range           [km] : ",y_range[0]/1e5,y_range[1]/1e5
  print "z range           [km] : ",z_range[0]/1e5,z_range[1]/1e5

print "cell number (height)     : ",number_cell_z
print "cell number (width)      : ",number_cell_y
print "T range              [T] : ",T_range
print "P range            [bar] : ",P_range[0]/1e6,P_range[1]/1e6
print "rho range       [g/cm^3] : ",rho_range[0],rho_range[1]

#print "Time difference between lines:",max_time/numberoflines
number_slice_x = int(raw_input("numberofslices :"))
#time_difference_T_P=max_time/numberoflines
old_time=0
time_difference=100.0
old_time=-time_difference
j=0
l=0
jj=0
l=0
R_gas=3.553188696490115*1e7
for i in profilefilelist:

  pf=load(i)

  time=float(pf.current_time)
  if(time>max_time):
    break
  time_tau=str(int(time/teddy))
  
  time=int(time-(time % 10.0))
  filename_time=str(time)
  filename=filename_time[:5]+'s'
  print(filename)

  if ((time>=float(jj+1)*time_difference_T_P or time==max_time) and time>0.0):
    jj+=1
    position=[0.33,0.66]
    position_y=(y_range[0]+position[0]*(y_range[1]+y_range[0]),\
                y_range[0]+position[1]*(y_range[1]+y_range[0]))
    slicevector=[0,position_y[0],z_1mbar*0.55]

    ps= SlicePlot(pf,'y','Temp',origin="native",center=slicevector,width=(z_1mbar*0.9,(y_range[0]+y_range[1])*0.5))
    pg0_frb=ps.data_source.to_frb((y_range[0]+y_range[1])*0.44,int(number_cell_width)*2,slicevector)
    center = ((x_range[0]+x_range[1])*0.5,(y_range[0]+y_range[1])*0.5)
    y=[x_range[0]-delta_x*0.5,x_range[1]+delta_x*0.5]-center[0]
    x=[0.0,z_1mbar*1.1]
    y=y/1e8
  
    extent = np.min(x)/1e5, np.max(x)/1e5, np.max(y)*1e3/2.0, np.min(y)*1e3/2.0
    ps_array=np.array(pg0_frb['Temp'])

    def format_func(value, tick_number):
    # find number of multiples of pi/2

      if value == 10000.0:
        return r"$10^{-3}$"
      elif value == 8000.0:
        return "0.006"#r"$6\times10^{-3}$"
      elif value == 7543.0:
        return "0.01"#r"$4\times10^{-2}$"

      elif value == 6000.0:
        return "0.04"#r"$4\times10^{-2}$"
      elif value == 5000.0:
        return "0.1"#r"$3\times10^{-1}$"

      elif value == 4000.0:
        return "0.3"#r"$3\times10^{-1}$"
      elif value == 2600.0:
        return r"$1.0$"

      elif value == 2000.0:
        return r"$1.6$"
      elif value == 170.0:
        return r"$10$"
      elif value == 0.0:
        return r"$11$"



    def format_func2(value, tick_number):
    # find number of multiples of pi/2
      return int(value/1000.0)
     # if value == 6000.0:
     #   return r"$6$"
     # elif value == 4000.0:
     #   return "4"#r"$6\times10^{-3}$"
     # elif value == 7543.0:
     #   return "0.01"#r"$4\times10^{-2}$"

     # elif value == 6000.0:
     #   return "0.04"#r"$4\times10^{-2}$"
     # elif value == 5000.0:
     #   return "0.1"#r"$3\times10^{-1}$"

     # elif value == 4000.0:
     #   return "0.3"#r"$3\times10^{-1}$"
     # elif value == 2600.0:
     #   return r"$1.0$"

     # elif value == 2000.0:
     #   return r"$1.6$"
     # elif value == 170.0:
     #   return r"$10$"
     # elif value == 0.0:
     #   return r"$11$"


    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    
    #plt.imshow(ps_array,origin="navtive",cmap='jet',extent=extent)
        #plt.clim(2800,3500)
    #fig = plt.figure()
    #ax = fig.add_axes ([0.92,0.1,0.8,0.02])
    #ax1 = fig.add_subplot(111)
    im=ax.imshow(ps_array,origin="navtive",cmap='jet',extent=extent, vmin=2700, vmax=3800)
    cbaxes=fig.add_axes([0.84,0.175,0.03,0.75])#(x-position,y-position of left-bottom,width,height)
    cbar=fig.colorbar(im,ticks=[2800,3000,3300,3600],cax=cbaxes)#,orientation='horizontal')
    cbar.set_label(r'$T$ [K]',rotation=270, labelpad=15,fontsize=15)
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation=270)
    cbar.set_clim(2700,3800)
   # cbar.set_clim(3800,2700)
    cbar.ax.invert_yaxis()
    #colormap_r = ListedColormap(colormap.colors[::-1])
    ax.set_xticks([170,2600,5000,7543.0,10000])

    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func2))
    scaleheight=R_gas*T_1mbar/1000.0/1e5

    plt.setp(ax.get_xticklabels(), rotation=270,fontsize=12)
    plt.setp(ax.get_yticklabels(), rotation=270,fontsize=12)
    ax.set_xlabel(r'$P$ [bar]', labelpad=-2,fontsize=15)
    ax.set_ylabel(r'$x$ [$10^{3}$ km]',rotation=270, labelpad=12,fontsize=15)
    ax.yaxis.set_label_coords(-0.15,0.5)
#ax11 = fig.add_subplot(111, aspect='equal')
    p1=patches.FancyArrowPatch(
      (6500,-5e3),
      (6500+2*scaleheight,-5e3),
      arrowstyle='|-|',edgecolor="white",
      mutation_scale=5,linewidth=5.0)
    ax.add_patch(p1)
    p=patches.FancyArrowPatch(
        (6500,-5e3),
        (6500+2*scaleheight,-5e3),
        arrowstyle='|-|',edgecolor="magenta",
        mutation_scale=5,linewidth=2.0)
    ax.add_patch(p)
#    ax.text(7050,-4.e3,r"2$H$",fontsize=20,rotation=270,color="white",fontweight='bold')
    ax.text(7050,-4.e3,r"2$H$",fontsize=20,rotation=270,color="magenta")
    plt.tight_layout()
    plt.gcf().subplots_adjust(left=0.13)
    plt.gcf().subplots_adjust(right=0.83)
    #plt.gcf().subplots_adjust(top=0.149)
    plt.gcf().subplots_adjust(bottom=0.15)
    print("DRAWING")
    fig.savefig('plot_Temp_'+filename+'.pdf')
    plt.close()
