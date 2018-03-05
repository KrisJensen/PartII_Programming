#!/usr/bin/env python

from huckel import *

print('\n\nYou are about to use a Huckel solver! Fasten your seatbelts and prepare for a wild ride')

def runPolyCalc(N, lin):
   '''creates an instance of the polyene class, solves the Huckel problem and prints the energies. N is number of atoms, 'lin' is True for linear polyene, 'False' for cyclic'''
   p=Polyene(N, lin)
   p.constructHuckel()
   p.solveHuckel()
   p.sort()   #sort eigenvalues and count degeneracies

   print("{:<10}".format('Energy')+"{:<10}".format('Degeneracy'))
   for i in range(len(p.degs)):
      print("{:<10}".format(str(np.round(p.ens[i],4)))+"{:<10}".format(p.degs[i]))   #print energy and degeneracy

def runSolidCalc(typ):
   '''creates an instance of the 'Solids' class for the platonic solids, solves the Huckel problem and prints the energies. 'typ' specifies which solid the calculation is run for'''
   p=Solid(typ)
   p.constructHuckel()
   p.solveHuckel()
   p.sort()   #sort eigenvalues and count degeneracies

   print("{:<10}".format('Energy')+"{:<10}".format('Degeneracy'))
   for i in range(len(p.degs)):
      print("{:<10}".format(str(np.round(p.ens[i],4)))+"{:<10}".format(p.degs[i]))   #print energy and degeneracy


def pols(lin=True):
   '''function for running a series of calculations on polyenes. lin is True for linear polyenes, False for cyclic'''
   try: N=int(input('\nHow many carbon atoms would you like in your molecule? '))
   except:
      print('Sorry, that was not a valid input. We currently only support integer inputs')   #if the input cannot be converted to an integer, something is wrong and we ask for a different input
      pols(lin)
      return

   if N==0:   #Catch Alex trying to break the program by running a calculation with 0 carbon atoms
      print('Sorry, we cannot construct a molecule with 0 carbon atoms, please try a different number')
      pols(lin)
      return

   print('\n- - - - - - - - - - - - - - - - - - - - -\n')
   print('Running calculation for N='+str(N)+'\n')
   runPolyCalc(N, lin)   #function creates a polyene, calculates energies and prints results
   print('\n- - - - - - - - - - - - - - - - - - - - -\n')
   new=str(input('Would you like to run a new calculation on the same type of polyene? (y/n) '))
   if new == 'y': pols(lin)   #give the user the option to run another calculation on the same type of molecule.
   else: return

def plats():
   '''function for running a series of calculations on the platonic solids'''
   print("\nyou can pick between a tetrahedron (tet), octahedron (oct), cube (cube), dodecahedron (dod), icosahedron (icos) or Buckminsterfullerene (b)\n")

   typ=str(input('What would you like to look at? '))
   avails=[['tet','oct','cube','dod','icos','b'],['tetrahedron','octahedron','cube','dodecahedron','icosahedron','Buckminsterfullerene']]
   try: assert typ in avails[0]
   except AssertionError:   #if the input is not one of the platonic solids implemented in the class, print an error message and try again
      print('Sorry, that was not a valid input. Please read the instructions and try again')
      plats()
      return

   print('\n- - - - - - - - - - - - - - - - - - - - -\n')
   print('Running calculation on '+avails[1][avails[0].index(typ)]+'\n')
   runSolidCalc(typ)   #function creates a Solid, calculates energies and prints results
   print('\n- - - - - - - - - - - - - - - - - - - - -\n')
   new=str(input('Would you like to run a new calculation on the Platonic solids or Buckminsterfullerene? (y/n) '))
   if new == 'y': plats()   #give the user the option to run another calculation on a platonic solid
   else: return None


new='y'
while new == 'y':   #outer loop in which calculations are run.

   print("\n\nType l to look at linear polyenes, type c to look at circular polyenes, or type p to look at Platonic solids and Buckminsterfullerene\n")
   struc=str(input('Which structure would you like to look at? '))

   if struc == 'l':
      print('\nYou have picked the linear polyenes!')
      print('- - - - - - - - - - - - - - - - - - - - -\n')
      pols(True)   #Environment for calculating energies of linear polyenes

   elif struc == 'c':
      print('\n- - - - - - - - - - - - - - - - - - - - -\n')
      print('You have picked the circular polyenes!')
      pols(False)   #Environment for calculating energies of cyclic polyenes
   
   elif struc == 'p':
      print('\n- - - - - - - - - - - - - - - - - - - - -\n')
      print('You have picked the Platonic solids!')
      plats()   #Environment for calculating energies of platonic solids

   else: print("\nWe are afraid that was not a viable option. Please type l, c or p\n")   #complain if the user does not pick an available option

   new=str(input('\nWould you like to run another calculation? (y/n) '))   #provides the user with the options of running multiple calculations. If anything else than 'y' is picked, the program is terminated.

print('\n- - - - - - - - - - - - - - - - - - - - -\n')
print('Thank you for shopping with us, we hope to see you again')
print('\n- - - - - - - - - - - - - - - - - - - - -\n')
