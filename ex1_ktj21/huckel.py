#!/usr/bin/env python

import numpy as np

class Molecule():
   '''Class that contains a general huckel solver. Polyenes and platonic solids will be implemented as subclasses as they both need to find solutions of a Huckel matrix'''

   def __init__(self):
      '''initializes class with empty arrays corresponding to eigenvalues and eigenvectors'''
      self.H=np.zeros()
      self.eig_vecs=np.array()
      self.eig_vals=np.array()

   def solveHuckel(self):
      '''Class method that uses numpy to diagonalize the Hamiltonian for the molecule; assumes this has already been built. This could have just been included in the 'sort' class method which sorts and counts degeneracies, but I figured it was nice to separate the two steps as one could imagine other uses of the 'raw' eigenvalues and eigenvectors which can then also be implemented as a separate method'''
      self.eig_vals, self.eig_vecs = np.linalg.eig(self.H)
      return self.eig_vals, self.eig_vecs

   def sort(self):
      '''takes the eigenvalues and eigenvectors generated for a Hamiltonian, sorts them according to decreasing eigenvalues and counts degeneracies. Returns lists of energies and corresponding degeneracies'''

      idx = self.eig_vals.argsort()[::-1]   #Generate index list for sorting eig vals and eig vecs
      vals = self.eig_vals[idx]   #create list of sorted eig vals
      self.orbs = self.eig_vecs[:,idx]   #create list of sorted eig vecs
      ens=[]
      degs=[]
      for e in vals:   #Count degeneracies
         try:
            if ens[-1] == np.round(np.real(e),8): degs[-1] += 1   #if eigenvalue identical to previous to 8th decimal, consider degenerate
            else:
               ens.append(np.round(np.real(e),8))   #if not, add to list of unique eigenvalues and start new degeneracy counter
               degs.append(1)
         except IndexError:   #Catch the IndexError raised when the first eigenvalue is considered
               ens.append(np.round(np.real(e),8))
               degs.append(1)
      self.ens=ens
      self.degs=degs


class Polyene(Molecule):
   '''Class for generating Huckel matrices for linear and cyclic polyenes. Subclass of molecule; inherits functions for diagonalizing Huckel matrices and sorting eigenvalues'''

   def  __init__(self, N, linear=True):
      try: assert type(N) == int   #make sure the size of the molecule is actually given by an integer. Otherwise complain.
      except AssertionError: print('sorry, N has to be an integer')

      if N > 5000: print('WARNING: for large N, apparent degeneracies are likely to occur as a result of very close spacings of energy levels')

      self.H=np.zeros((N,N))
      self.lin=linear   #We define the molecule as circular or linear as we initiate an instance of the class
      self.N=N

   def constructHuckel(self,a=0,b=1):
      '''Class method to generate a Huckel matrix'''
      for i in range(self.N):
         self.H[i,i]=a   #set all diagonal elements to alpha (default 0)
         if i!= 0:
            self.H[i,i-1]=b   #set all interactions between adjacent atoms to beta (default 1)
            self.H[i-1,i]=b
      if not self.lin:   #for cyclic polyenes, we add an interaction between atoms 1 and N
            self.H[0,self.N - 1]=b
            self.H[self.N - 1,0]=b
      return self.H

class Solid(Molecule):
   '''Class for generating Huckel matrices for the platonic solids. The adjacencies for atoms in the platonic solids have been hardcoded and are given by tuples of tuples of indices of adjacent atoms'''
   def __init__(self, type='tet'):   #specify the type of platonic solid as an instance of the class is initiated
      self.type=type
      self.adj={}
      self.adj['tet']=((2,3,4),(1,3,4),(1,2,4),(1,2,3))   #Adjacencies for tetrahedron
      self.adj['cube']=((2,4,6),(1,3,7),(2,4,8),(1,3,5),(6,8,4),(5,7,1),(6,8,2),(3,5,7))  #Adjacencies for cube
      self.adj['dod']=((2,5,8),(1,3,10),(2,4,12),(3,5,14),(4,1,6),(5,7,15),(6,8,17),(7,9,1),(8,10,18),(9,11,2),(10,12,19),(11,13,3),(12,14,20),(13,15,4),(14,16,6),(15,17,20),(16,18,7),(17,19,9),(18,20,11),(19,16,13))   #Adjacencies for dodecahedron
      self.adj['oct']=((2,3,4,5),(1,3,5,6),(1,2,4,6),(3,5,1,6),(2,4,1,6),(2,3,4,5))   #Adjacencies for octahedron
      self.adj['icos']=((2,3,4,5,6),(1,3,6,7,8),(1,2,4,8,9),(1,3,5,9,10),(1,4,6,10,11),(1,5,2,11,7),(2,6,8,11,12),(2,3,7,9,12),(3,4,8,10,12),(4,5,9,11,12),(5,6,10,7,12),(7,8,9,10,11))   #Adjacencies for icosahedron
      self.adj['b']=((58,59,60),(56,57,60),(55,57,59),(53,54,60),(52,54,58),(50,51,59),(51,55,49),(47,48,57),(46,48,56),(44,45,54),(43,45,53),(41,42,52),(39,40,51),(38,40,42),(37,38,39),(35,36,49),(33,34,48),(32,34,36),(31,32,33),(30,31,43),(28,29,44),(26,27,38),(25,27,32),(25,26,29),(23,24,31),(22,24,28),(22,23,37),(21,26,42),(21,24,30),(20,29,45),(19,20,25),(18,19,23),(17,19,46),(17,18,35),(16,34,47),(16,18,37),(15,27,36),(14,15,22),(13,15,49),(13,14,41),(12,40,50),(12,14,28),(11,20,46),(10,21,52),(10,11,30),(9,33,43),(8,35,55),(8,9,17),(7,16,39),(6,41,58),(6,7,13),(5,12,44),(4,11,56),(4,5,10),(3,7,47),(2,9,53),(2,3,8),(1,5,50),(1,3,6),(1,2,4))   #Adjacencies for the Bucky ball

   def constructHuckel(self,a=0,b=1):
      '''Class method for generating a Huckel matrix'''
      self.H=np.zeros((len(self.adj[self.type]),len(self.adj[self.type])))   #Initiate Hamiltonian
      for i, r in enumerate(self.adj[self.type]):
         self.H[i,i]=a   #Set all diagonal elements to alpha (default 0)
         for j in r:
            self.H[i,j-1]=b   #Set all interactions between adjacent carbons to beta (default 1)
      return self.H


