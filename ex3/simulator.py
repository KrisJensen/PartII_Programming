'''A script that is used to run one of three different kinetics simulations as specified in the input. 'fold' simulates folding a protein in the presence of increasing urea concentrations. 'oregonator' simulates the oregonator reaction for 90 seconds. 'diffusion' simulates a set of cells in which the oregonator reaction runs as it diffuses along the cells.
Use: python simulation.py *simulation type* (#cells D)'''

import numpy as np
import copy
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pickle
import sys
import time

class cell():
   '''cell class used for keeping track of reactions and concentrations during simulations'''
   def __init__(self):
      self.reactants={}
      self.reactions={}
      self.time=0.0

   def addReactant(self, name, conc):
      '''name is name of reactant, conc is a float giving the concentration in M'''
      self.reactants[name]=conc

   def addReaction(self, label, reactants, products, rate):
      '''reactants/products are lists of reactants/products (given as strings). rate is a float giving the rate constant'''
      self.reactions[label]=[reactants, products, rate]

   def step(self, t):
      '''This method advances the reaction taking place in the cell by time t. First calculates the change in concentration of each species and then adds this change to current values'''

      deltas={}   #dict of how much the concentration of each species changes
      for reactant in self.reactants: deltas[reactant]=0.0

      #for each reaction, calculate the change in reactants and products in time t
      for reaction in self.reactions.keys():
         rxn=self.reactions[reaction]   #rxn is a list of reactants, products, rate constant

         delta=rxn[2]*t   #for the change in reactant/product we multiply rate constant with time step and concentration of each of the reactants
         for reactant in rxn[0]:
            delta *= self.reactants[reactant]   #multiply change by current concentration

         for reactant in rxn[0]: deltas[reactant] -= delta   #calculate change in reactants
         for product in rxn[1]: deltas[product] += delta   #calculate change in products

      for species in deltas.keys():
         self.reactants[species] += deltas[species]   #update concentrations
      self.time += t   #update elapsed time

   def returnReactants(self):
      '''returns array of current reactant concentrations sorted alphabetically.'''
      return np.array([self.reactants[k] for k in sorted(self.reactants.keys())])

   def folder(self, U=0.0, D=1.0, I=0.0, N=0.0):
      '''Adds the appropriate species and reactions for initializing a protein folding simulation'''
      self.addReactant('D', D)
      self.addReactant('I', I)
      self.addReactant('N', N)
      self.addReaction('3.3', ['D'], ['I'], 26000*np.exp(-1.68*U))
      self.addReaction('3.4', ['I'], ['D'], 0.06*np.exp(0.95*U))
      self.addReaction('3.5', ['I'], ['N'], 730*np.exp(-1.72*U))
      self.addReaction('3.6', ['N'], ['I'], 0.00075*np.exp(1.20*U))

   def Oregonator(self, A=0.06, B=0.06, P=0.0, Q=0.0, X=10**(-9.8), Y=10**(-6.52), Z=10**(-7.32)):
      '''Adds the appropriate species and reactions for initializing an Oregonator simulation'''
      self.addReactant('A', A)
      self.addReactant('B', B)
      self.addReactant('P', P)
      self.addReactant('Q', Q)
      self.addReactant('X', X)
      self.addReactant('Y', Y)
      self.addReactant('Z', Z)
      self.addReaction('1', ['A','Y'], ['X','P'], 1.34)
      self.addReaction('2', ['X','Y'], ['P'], 1.6*10**9)
      self.addReaction('3', ['B','X'], ['X','X','Z'], 8*10**3)
      self.addReaction('4', ['X','X'], ['Q'], 4*10**7)
      self.addReaction('5', ['Z'], ['Y'], 1)

   def __str__(self):
      '''prints a list of reactants and concentrations'''
      string=""
      for r in self.reactants.keys():
         string += str(r)+': '+str(self.reactants[r])+'\n'
      return string

class cellLine():
   '''Class for modelling spatial diffusion of reactions occurring in cells. Based on the diffusion equation: c(x, t) = c0/sqrt(4 pi D t) exp(-x^2/ (4 D t)). Uses coordinates 0:1 for the cells. This might have been done better as a monte carlo diffusion simulation'''

   def __init__(self, cellnumber, timestep, D):
      ''''cellnumber' specifies the number of cells used for modelling diffusion (integer). 'timestep' specifies the timestep used for propagating the simulation (in seconds). 'D' specifies the diffusion constant used for the simulation. Initially all reactants are contained at a point (first cell) and and they are subsequently allowed to diffuse to more distant cells as the reaction proceeds'''
      self.time = 0
      self.timestep = timestep
      self.cellnumber=cellnumber
      self.D=D

      self.cells=[cell()]
      self.cells[0].Oregonator()   #initializes first cell with default concentrations
      for i in range(cellnumber-1):   #add empty cells as specified by 'cellnumber'
         self.cells.append(cell())
         self.cells[-1].Oregonator(0, 0, 0, 0, 0, 0, 0)

      self.coeffs=[]
      #calculate diffusion coefficients for cells a given distance apart as these are contant. Normalize length so final cell is at 1.0. To reduce computational times, diffusion is only considered to occur between cells less than 3 cells apart.
      for i in range(4):
            self.coeffs.append( (4*timestep*np.pi*D)**(-0.5)/(cellnumber-1) * np.exp( - i**2 / ((cellnumber-1)**2*timestep*4*D) ))

      #print(self.coeffs)
      #for c in self.cells: print(c)

   def step(self):
      '''propagates the simulation by one timestep. This is done by first allowing the reactions in every cell to proceed by self.timestep seconds, then allowing diffusion to occur for self.timestep seconds'''

      for cell in self.cells: cell.step(self.timestep)   #let the reaction proceed in each cell

      newconcs=np.array([np.array([0.0 for x in range(7)]) for x in range(self.cellnumber)])   #initialize array for storing new concentrations with the structure [ cell1[[A], [B]..], cell2[A], [B]...] ]

      for orig in range(self.cellnumber):   #for each cell, calculate the amount of reactants that has diffused from this cell to any other cell
         coeff_sum=0.0
         deltas=np.array([np.array([0.0 for x in range(7)]) for x in range(self.cellnumber)])   #array for storing the diffused reactants to each other cell

         for adj in range(orig-3, orig+4):   #diffusion to cells up to three cells away
            if adj < 0 or adj >= len(self.cells): continue
            coeff = self.coeffs[np.abs(orig-adj)]   #retrieve pre-calculated prefactor for diffusion
            coeff_sum += coeff
            deltas[adj] = self.cells[orig].returnReactants() * coeff  #this is the change in concentrations due to diffusion from cell orig to cell adj
         newconcs += deltas / coeff_sum  #normalize changes so no material is lost (this is necessary because we're using a small and finite array)

      for i, cell in enumerate(self.cells):
         for j, species in enumerate(['A','B','P','Q','X','Y','Z']): cell.reactants[species]=newconcs[i][j]   #update concentrations in each cell

      self.time += self.timestep


def diffusionSim(cells, t, tmax, D):
   '''Runs a diffusion simulation for the oregonator using 'cells' cells and 'steps' steps of time t.
   returns an array of time values and a list (data) of lists corresponding to different cells, each with 7 lists corresponding to concentrations of the species at each time step ([ cell1[ [A0, A1, ...], [B0, B1, ...], ...], cell2[ [A0, A1, ...], [B0, B1, ...], ...], cell... ])'''

   cl = cellLine(cells, t, D)
   data = [[[] for x in range(7)] for x in range(cells)]

   ts=[0]
   for i, cell in enumerate(cl.cells):
      for j, conc in enumerate(list(cell.returnReactants())):
         data[i][j].append(conc)   #add initial concentrations for t=0

   n = 0   #keep counter of iterations
   for time in np.linspace(t, tmax, tmax/t):
      n += 1
      cl.step()   #propagate the simulation by one timestep

      if n % 10000 == 0:
         ts.append(cl.time)
         for i, cell in enumerate(cl.cells):
            for j, conc in enumerate(list(cell.returnReactants())):
               data[i][j].append(conc)   #update concentrations

      if n % 10000 == 0:   #occasionally print the progress so we know what's going on
         print('time:', cl.time)
         for cell in cl.cells: print(cell)

   return ts, data

def plotDif(ts, data, D):
   '''takes datastructure as specified in 'diffusionSim' and returns plots of concentrations vs time for each cell'''

   print('plotting data')

   for i in range(len(data)):   #make new plot for each cell (named according to cell number and value for diffusion constant)

      name='D'+str(D)+'cell'+str(i)
      ax=plt.gca()
      for j, d in enumerate(data[i]):
         plt.semilogy(ts, d, style[j], label=leg2[j])

      lgd = ax.legend(loc='lower right')
      plt.ylabel('Concentration / M')
      plt.xlabel('Time / s')
      plt.title('Oregonator Diffusion Simulation Cell #'+str(i))
      plt.xlim([min(ts), max(ts)])
      plt.savefig(name+'.pdf', bbox_inches='tight')
      plt.close()


def simOne(sim, t, cutoff=10**(-11), tmax=np.inf, U=0.0):
   '''Runs a single simulation with a given cutoff in terms of change in square norm (cutoff) or maximum time (tmax) for completion. U species concentration of urea for protein folding simulation. 'sim' specifies whether the simulation is for folding of a protein ('folder'; returns final concentration) or an oregonator ('Oregonator'; returns all concentrations over the course of the simulation'''

   Cell=cell()
   if sim=='folder': Cell.folder(U)
   else: Cell.Oregonator()

   c1 = Cell.returnReactants()   #array of concentrations of species
   dif = cutoff+1
   n=0   #set limit to prevent simulation from running indefinitely
   concs=[[],[],[],[],[],[],[]]   #store concentrations for oregonator
   ts=[]   #store times for oregonator

   while dif > cutoff and Cell.time < tmax and n < 100000000:   #if norm change in concentrations is greater than cutoff, the simulation has run for long enough, or we have run too many iterations, we quit the simulation
      n += 1
      if sim=='Oregonator' and n%1000==0:   #if we're running an Oregonator simulation, we store concentrations every 1000th timepoint
         ts.append(Cell.time)
         for i, c in enumerate(list(c1)): concs[i].append(c)

      Cell.step(t)   #advance reaction by one timestep
      c0=copy.deepcopy(c1)  
      c1=Cell.returnReactants()   #update concentrations
      dif=np.linalg.norm(c1-c0)   #calculate norm of change in concentrations
      if n%50000==0:
         print('\n', n, '   t:', np.round(Cell.time, 5), '\nc0:', c0, '\nc1:', c1, '\ndif:', dif)
   if sim=='folder': return c1   #return final concentration if folding a protein
   else: return ts, concs   #return all concentrations if running oregonator


def simProtMany(umin=0.0, umax=10, steps=101, t=0.00005):
   '''runs a series of protein folding simulations with urea concentrations between umin and umax. The number of concentrations is given by the 'steps' argument. t specifies the timestep used. returns tuple of ([c(urea)], [[c(D)],[c(I)],[c(N)]])'''
   xs=[]
   concs=[[],[],[]]   
   for urea in np.linspace(umin, umax, steps):
     xs.append(urea)
     newconc = list(simOne('folder', t, U=urea))   #run a simulation with the given urea concentration
     print('U:', urea, 'concs:', newconc) #print some output to monitor progress
     for i in range(3):
        concs[i].append(newconc[i])   #add concentrations from simulation to lists of concentrations
   return xs, concs

def plot(xs, data, legends, log=False, name='noname'):
   '''data is list of lists of data to be plotted against xs. log specifies semilog y axis, name specifies name of file'''
   print('plotting data')
   for i, d in enumerate(data):
      #for each list in 'data', plot against xs
      if log: plt.semilogy(xs, d, style[i], label=legends[i])
      else: plt.plot(xs, d, style[i], label=legends[i])
   ax = plt.gca()

   if name=='Protein':
      lgd = ax.legend()
      plt.xlabel('[Urea] / M')
      plt.ylabel('Equilibrium Fraction')
      plt.title('Protein Folding')
      plt.ylim([0,1])
   else:
      lgd = ax.legend(loc='lower right')
      plt.ylabel('Concentration / M')
      plt.xlabel('Time / s')
      plt.title('Oregonator Time Simulation')
      plt.ylim([10**(-11), 10**(-1)])

   plt.xlim([min(xs), max(xs)])
   plt.savefig(name+'.pdf', bbox_inches='tight')
   plt.close()


if __name__ == '__main__':

   #lets just hard-code some legends and colours for plots
   leg1=['D', 'I', 'N']
   leg2=['A','B','P','Q','X','Y','Z']
   style=['--r','--g','--b', '--c', '--m', '--k', '--y']

   if len(sys.argv) < 2:
      print('\nERROR: Calculation type not specified.\n\nUse: python simulator.py *Calculation type (fold/oregonator/diffusion)* (OPTIONAL: *#cells* *D*)\n')
      exit()

   if sys.argv[1]=='fold':   #Run simulation for protein folding
      xs, concs = simProtMany(t=0.00005)
      plot(xs, concs, leg1, name='Protein')

   if sys.argv[1]=='oregonator': #Run Oregonator timecourse
      ts, concs=simOne('Oregonator', 0.000001, cutoff=0.0, tmax=90.0)
      plot(ts, concs, leg2, log=True, name='Oregonator')

   if sys.argv[1]=='diffusion':

      if len(sys.argv) > 2:
         try:
            cells = int(sys.argv[2])
            if len(sys.argv) > 3:   D = float(sys.argv[3])
            else: D = 500
         except ValueError:
            print('\nERROR: wrong specification of #cells or D. Must be int and int/float respectively.\n\nUse: python simulator.py *Calculation type (fold/oregonator/diffusion)* (OPTIONAL: *#cells* *D*)\n')
            exit()
      else: cells, D = 5, 500

      if D==0:
         print('cannot run a diffusion simulation without diffusion; try an oregonator simulation instead')
         exit()


      print('running diffusion simulation with', cells, 'cells and D =', D)
      ts, data = diffusionSim(cells, 0.000001, 120, D)
      plotDif(ts, data, int(D))

   else: print('\nERROR: Calculation type not specified.\n\nUse: python simulator.py *Calculation type (fold/oregonator/diffusion)* (OPTIONAL: *#cells* *D*)\n')
 
