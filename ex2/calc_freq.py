#!/usr/bin/env python

'''takes a directory with Gaussian output files for H2O or H2S as a commandline argument, extracts the energy as a function of bond length and angle, fits a quadratic function around the minimum and calculates vibrational frequencies and plots the data.'''


import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import sys

def extract(path):
   '''function for extracting bond lengths, angles and energies from all files in a given directory. Assumes directory contains only Gaussian single point output files with a name of the form $molecule$.r$length$theta$angle$.out'''

   lens, angs, ens = [], [], []
   try:
      for filename in os.listdir(path):   #look at all files in directory
         lens.append(float(filename.split('theta')[0].split('.r')[1]))
         angs.append(float(filename.split('theta')[1].split('.out')[0])) 
         with open(str(path+filename), 'r') as f:
            for line in f:
               if 'SCF Done:' in line:   #this line contains the final energy as the fifth word
                  ens.append(float(line.split()[4]))
   except FileNotFoundError:
      print('ERROR: Directory not found\nUse: python calc_freq.py *name_of_directory*\n')
      exit()
   except IndexError:
      print('ERROR: File format incorrect\nPlease specify directory containing only Gaussian single point calculation output files for either H2S or H2O\n')
      exit()

   #Find unique bond lengths and angles as well as the corresponding indices in the original array using numpy.unique 
   x_vals, x_idx = np.unique(lens, return_inverse=True)
   y_vals, y_idx = np.unique(angs, return_inverse=True)
   xs, ys = np.meshgrid(x_vals, y_vals)   #construct meshgrids with all r and theta values
   vals_array = np.empty(x_vals.shape + y_vals.shape)   #construct empty grid for inserting energies
   vals_array[x_idx, y_idx] = ens   #for each original pair of (r,ang), add the energy to the corresponding position on the grid of r and ang given by the uique values (as determined by the indices)
   zs=vals_array.T-np.amin(vals_array)   #consider the lowest energy as the zero-point energy and subtract from all values
   mol=os.listdir(path)[0].split('.r')[0]   #determine which molecule we're looking at

   print('finished extracting data\n')
   return xs, ys, zs, mol

def reg(xs, ys, zs, mol):
   '''Performs quadratic regression on data around minimum. xs, ys are meshgrids, zs an array of data. Fits x and y separately and returns the coefficients of the squared terms as well as the corresponding vibrational frequencies'''

   xmin, ymin = [int(x) for x in np.where(zs==np.amin(zs))] #find coordinates of minimum energy
   #extract the 5 energies surrounding the minimum along both the x and y axes and the corresponding bond lengths and angles. These are used for regression.
   xreg, zregx = [a[xmin, ymin-2:ymin+3] for a in (xs, zs)]
   yreg, zregy = [a[xmin-2:xmin+3, ymin] for a in (ys, zs)]

   p_r=np.polyfit(xreg-xreg[2], zregx, 2)   #units of Eh*angstrom^-2
   p_ang=np.polyfit(yreg-yreg[2], zregy, 2)  #units of Eh*degrees^-2

   k_r=2*p_r[0] * 4.35974*10**(-18) * (10**10)**2   #units of J*m^-2
   k_ang=2*p_ang[0] * 4.35974*10**(-18) * (180/np.pi)**2   #units of J*radian^-2

   m_u = 0.001/(6.022*10**23)
   #use formulae given in handout to calculate frequencies
   nu_r = (k_r/(2*m_u))**0.5 / (2*np.pi)
   nu_ang = (k_ang/( 0.5*m_u * (xreg[2] * 10**(-10))**2 ))**0.5 / (2*np.pi)

   print('Stretching frequency: '+str(int(nu_r/10**10))+'x10^10 s^-1')
   print('Bending frequency: '+str(int(nu_ang/10**10))+'x10^10 s^-1')
   return p_r[0], p_ang[0], nu_r, nu_ang

def plot(xs, ys, zs, mol, p_r, p_ang, nu_r, nu_ang, open=True):
   '''takes xs, ys: meshgrids of x and y values to be plotted. zs: array of z values to be plotted. mol: name of molecule to be used as title. Creates and opens colour coded 2D surface of z-values.'''
   print('\nplotting data\n')

   fig=plt.figure()
   ax=fig.gca(projection='3d')
   #some features of the plot will depend on whether we're looking at H2O or H2S. Adjust POV (elev), maximum of z-axis (zmax) and location of text box (height) accordingly.
   if mol[-1]=='O': zmax, elev, height = np.amax(zs), 40, 0.35
   else: zmax, elev, height = 0.32, 50, 0.25
   xmin, ymin = [int(x) for x in np.where(zs==np.amin(zs))]   #get coordinates of minimum energy point
   zs[zs>zmax]=np.nan   #only plot values lower than our cutoff
   surf=ax.plot_surface(xs, ys, zs, cstride=1, rstride=4, cmap=cm.coolwarm, linewidth=0, vmin=0, vmax=zmax)
   #surf=ax.plot_surface(xs, ys, zs, linewidth=0, c=zs)

   #Generate data corresponding to our idealized harmonic vibrations
   xreg=xs[xmin, ymin-5:ymin+6]
   zregx=[p_r*(x-xreg[5])**2 for x in xreg]
   yreg=ys[xmin-18:xmin+19, ymin]
   zregy=[p_ang*(y-yreg[18])**2 for y in yreg]

   ax.plot(xreg, [ys[xmin,0] for x in xreg], zregx, c='r')   #plot harmonic data
   ax.plot([xs[0,ymin] for y in yreg], yreg, zregy, c='r')

   #Adjust parameters to make plot look nice
   ax.view_init(elev=elev, azim=75)
   ax.set_zlim(0, zmax)
   ax.set_xlabel('r / angstrom')
   ax.set_ylabel('$\Theta$ / degrees')
   ax.set_zlabel('E - $\mathrm{E_{min}}$ / $\mathrm{E_h}$')
   ax.set_title(mol)

   #Generate and add text box with regression data and vibrational frequencies
   textstr="$E = E0 +$ "+str(np.round(p_r,2))+" $(r-r_{eq})^2$\n      $+ $"+str(np.round(p_ang,7))+" $(\Theta-\Theta_{eq})^2$\n"
   textstr += "$\\nu_1 = $"+str(int(nu_r/10**10))+"$x10^{10} \, s^{-1}$\n"
   textstr += "$\\nu_2 = $"+str(int(nu_ang/10**10))+" $x10^{10} \, s^{-1}$"
   ax.text(0.3, 80, height, textstr, fontsize=8, bbox=dict(facecolor='grey', alpha=0.5), horizontalalignment='right')

   fig.savefig(mol+'.pdf', bbox_inches='tight')
   if open: plt.show()


if __name__ == '__main__':

   if len(sys.argv) > 1: path = sys.argv[1]+str('/')
   else:
      print('\nERROR: No Directory specified')
      print('Use: python calc_freq.py *name_of_directory*\n')
      quit()

   print('\nExtracting and analyzing data from', path, '\n')

   xs, ys, zs, mol = extract(path)   #Parse files in target directory
   p_r, p_ang, nu_r, nu_ang = reg(xs, ys, zs, mol)   #Perform linear regression and calculate vibrational frequencies
   plot(xs, ys, zs, mol, p_r, p_ang, nu_r, nu_ang)   #Plot the data

   print('\nTerminating program\n')

