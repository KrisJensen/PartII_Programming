'''script for calculating the equilibrium geometry of a system of 7 particles interacting via a given pair poential.
Use: python simulation.py *potential (LJ/M1/M2)* (OPTIONAL:*d0* *nrand* *stepsize* *cutoff* *nmax*)
See readme for further details'''

import numpy as np
import sys
import random

class particle():
   '''class for manipulating and moving particles'''
   def __init__(self, position):
      '''initialize particle with position=[x, y, z]'''
      self.pos=np.array(position)

   def getCoords(self):
      return self.pos

   def move(self, dr):
      '''move particle by an amount dr=[dx, dy, dz]'''
      self.pos += np.array(dr)

   def getDist(self, other):
      '''returns the geometric distance between self and another particle'''
      return np.linalg.norm(self.pos-other.pos)
      
   def __str__(self):
      return str(tuple(self.pos))


def calcPot(p1, p2, pot):
   '''Calculate the potential energy resulting from interaction of particles 1 and 2 given a potential (LJ=Lennard-Jones, M1=Morse req/sigma=1, M2=Morse req/sigma=2)'''
   r=p1.getDist(p2)
   if pot=='LJ': u=4*(1/r**12-1/r**6)
   elif pot=='M1': u=(1-np.exp(-r+1))**2
   elif pot=='M2': u=(1-np.exp(-r+2))**2
   else:
      print('Potential not recognized, please try again')
      print('Use: python simulation.py *potential (LJ/M1/M2)* (OPTIONAL:*d0* *nrand* *stepsize* *cutoff* *nmax*)\n')
      exit()
   return u


class system():
   '''class for constructing and manipulating a system of particles. Contains methods for calculating the potential energy of the system and minimizing this energy'''

   def __init__(self, potential):
      '''initialize system with 0 particles specifying the potential that is used to hold the particles together'''
      self.particles=[]
      self.nParticles=0
      self.pot=potential
      self.gradList=[np.inf]

   def addParticle(self, position):
      '''add new particle at position=[x, y, z]'''
      self.particles.append(particle(position))
      self.nParticles += 1

   def calcPot(self):
      '''return the total potential energy of the system'''
      u = 0
      for p1 in range(self.nParticles):
         for p2 in range(p1+1, self.nParticles):   #consider interactions between all pairs of particles
            u += calcPot(self.particles[p1], self.particles[p2], self.pot)
      return u

   def calcPartGrad(self, pi, delta):
      '''calculate the gradient of the potential energy with respect to the coordinates [x, y, z] of a single particle i. return [grad(ix), grad(iy), grad(iz)]. delta specifies the displacement used for the numerical calculation'''
      U0=self.calcPot()   #original potential
      self.particles[pi].move([delta, 0, 0])   #move x and calculate new potential
      Ux=self.calcPot()
      self.particles[pi].move([-delta, delta, 0])   #move y and calculate new potential (moving x back)
      Uy=self.calcPot()
      self.particles[pi].move([0, -delta, delta])   #move z and calculate new potential (moving y back)
      Uz=self.calcPot()
      self.particles[pi].move([0, 0, -delta])   #restore original position
      
      grads=[]
      for i, U in enumerate([Ux, Uy, Uz]):
         grads.append((U-U0)/(2*delta))   #calculate gradient with respect to coordinate q as (Uq-U0)/(2delta)
      return grads

   def calcGrad(self, delta):
      '''calculate the gradient of Epot with respect to all coordinates and return as [[grad(p1], [grad(p2)]...]'''
      grads=[]
      for p in range(self.nParticles):
         grads.append(self.calcPartGrad(p, delta))   #calculate gradient with respect to coordinates of particle p and do this for all particles

      gradList=[]
      for grad in grads:
         for g in grad: gradList.append(g)
      self.gradList=gradList   #maintain an updated list of gradients that can be used for assessing convergence

      return grads

   def takeStep(self, delta):
      '''take a step in the direction of steepest descent given by gradients*delta'''
      grads=np.array(self.calcGrad(delta))  #calculate gradients
      for i, grad in enumerate(grads):
         if np.linalg.norm(grad)*delta > 1.0: grad *= 1/(np.linalg.norm(grad)*delta)  #we impose a maximum stepsize to avoid particles flying apart when using random starting conditions
         self.particles[i].move(grad*delta*(-1))   #move each particle i by [dx, dy, dz] = [dU/dxi*delta, dU/dyi*delta, dU/dzi*delta]
   
   def getCoords(self):
      '''return list of coordinates: [1x, 1y, 1z, 2x, 2y...]'''
      coords=[]
      for particle in self.particles:
         for coordinate in particle.getCoords():
            coords.append(coordinate)
      return np.array(coords)

   def getDists(self):
      '''returns an array with distances between particles i, j'''
      dists=np.zeros((self.nParticles, self.nParticles))
      for i, p0 in enumerate(self.particles):
         for j, p1 in enumerate(self.particles):
            dists[i, j]=p0.getDist(p1)
      return dists

   def getGradNorm(self):
      '''returns the 2 norm of the gradient vector for the system'''
      return np.linalg.norm(np.array(self.gradList))

   def __str__(self):
      return str(np.round(self.getDists(), 3))

def optimizeSys(sys, step, cutoff, nmax):
   '''given a system sys, returns the system with optimized coordinates and convergence parameter. step: delta used at each iteration and used for gradient calculations. cutoff: norm of change in position matrix at which the calculation is considered to have converged. nmax: maximum number of iterations.'''

   n = 0
   d0, distnorm = sys.getDists(), cutoff+1

   while distnorm > cutoff and n < nmax:   #interate until converges or maximum iterations exceeded
      n += 1
      sys.takeStep(step)
      d1=sys.getDists()
      distnorm = np.linalg.norm(d1-d0)

      if n%1000 == 0:   #print some information every 1000 iterations to monitor progress
         print('\nn=', n, '\nDistance matrix=\n', sys)
         print('Potential:', sys.calcPot())
         print('Gradient norm:', sys.getGradNorm())
         print('Distance norm', distnorm, '\n')

      d0=d1
         
   return sys, distnorm

if __name__ == '__main__':

   try:   #initialize the variables used in the optimization if given as arguments
      pot=sys.argv[1]
      defaults=[1.0, 0, 0.0001, 10**(-10), 500000] #default [d0, nrand, step, cutoff, nIter]
      if len(sys.argv)==2: d0, nrand, step, cutoff, nmax = defaults[0:5]
      elif len(sys.argv)==3: d0, nrand, step, cutoff, nmax = float(sys.argv[2]), int(defaults[1]), float(defaults[2]), float(defaults[3]), float(defaults[4])
      elif len(sys.argv)==4:
         d0, nrand, step, cutoff, nmax = float(sys.argv[2]), int(sys.argv[3]), defaults[2], defaults[3], defaults[4]
      elif len(sys.argv)==5: d0, nrand, step, cutoff, nmax = float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), defaults[3], defaults[4]
      elif len(sys.argv)==6: d0, nrand, step, cutoff, nmax = float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), defaults[4]
      else: d0, nrand, step, cutoff, nmax = float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])

   except:
         print('\nError in parsing arguments, please try again')
         print('Use: python simulation.py *potential (LJ/M1/M2)* (OPTIONAL:*d0* *nrand* *stepsize* *cutoff* *nmax*)\n')
         exit()

   if not nrand==0:   #run algorithm with random starting conditions
      print('\nRunning simulation for', nrand, 'random initial conditions')
      sysses=[]
      distnorms=[]
      ens=[]
      for guesses in range(nrand):
         sys=system(pot)
         for i in range(7): sys.addParticle([d0*random.random(), d0*random.random(), d0*random.random()])   #add 7 particles at random locations in box of side length d0
         print('\nInitial system'+str(guesses+1)+':')
         print(sys)
         print('Potential energy:', sys.calcPot())
         sys, distnorm = optimizeSys(sys, step, cutoff, nmax)   #optimize by steepest descent from initial guess
         sysses.append(sys)   #save converged syste
         distnorms.append(distnorm)

      for num, sys in enumerate(sysses):   #print a summary of the results at the end
         gradNorm, potE = sys.getGradNorm(), sys.calcPot()
         ens.append(potE)
         print('\nSystem'+str(num+1))
         print('Calculation converged to a gradient norm of', gradNorm, 'and distance norm of', distnorms[num])
         print('Final distance matrix:') 
         print(sys)
         print('Potential energy:', potE, '\n')

      emin=min(ens)   #print some information about the supposedly global minimum (the lowest energy converged system)
      imin = ens.index(emin)
      minsys  = sysses[imin]
      print('\nBest system is system'+str(imin+1))
      print(minsys)
      print('Potential energy:', emin, '\n')
 
   else:
      init=[]   #build an array of initial positions; we will just use corners of a with sidelength d0 to make life simple. Different starting configurations may converge to different local minima
      for i in range(2):
         for j in range(2):
            for k in range(2):
               init.append([float(d0*i), float(d0*j), float(d0*k)])
      sys=system(pot)
      for i in range(7): sys.addParticle(init[i])   #add 7 particles at corners of a cube

      print('\nInitial system:')
      print(sys)
      print('Potential energy:', sys.calcPot())

      sys, distnorm = optimizeSys(sys, step, cutoff, nmax)   #optimize system

      gradNorm, potE = sys.getGradNorm(), sys.calcPot()   #print some information about the converged system
      print('\nCalculation converged to a gradient norm of', gradNorm, 'and distance norm of', distnorm)
      print('Final distance matrix:')
      print(sys)
      print('Potential energy:', potE, '\n')

