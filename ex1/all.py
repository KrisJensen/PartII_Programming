#!/usr/bin/env python

import numpy
from huckel import *



lins=[]
circs=[]

for i in range(2,11):

   print 'new N: ', i

   lins.append(Polyene(i))
   circs.append(Polyene(i, False))
   
   print 'lin'

   lins[-1].constructHuckel()
   lins[-1].solveHuckel()

   lins[-1].sort()

   print 'circ'

   circs[-1].constructHuckel()
   circs[-1].solveHuckel()
   circs[-1].sort()


solids=['tet','cub','dod']
pols=[]
for s in solids:
   p=Solid(s)
   p.constructHuckel()
   p.solveHuckel()
   p.sort()
   pols.append(p)


sols=[lins,circs]

for s in sols:
 
   if s==lins: print('Linear Polyenes\n')
   elif s==circs: print('\n\nCircular Polyenes\n')
 

   for n in s:
      print("\nN = "+str(n.N))
      print("{:<10}".format('Energy')+"{:<10}".format('Degeneracy'))

      for i in range(len(n.degs)):
         print("{:<10}".format(str(n.ens[i]))+"{:<10}".format(n.degs[i]))

for p in pols:
   print("\n"+str(p.type))
   print("{:<10}".format('Energy')+"{:<10}".format('Degeneracy'))

   for i in range(len(p.degs)):
      print("{:<10}".format(str(p.ens[i]))+"{:<10}".format(p.degs[i]))
 

