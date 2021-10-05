from dolfin import Expression, Constant

# Inputs
# Mesh File
caseId = 'TestFlow10m_'
meshPath = 'Mesh/'
meshFile = 'pipeMesh_tri_138k' 
resPath = 'res/'+caseId+'/' 
replaceGeometry = False
# Choose between "Pure Difusion" OR "Flow and Difusion"
simulationType = "Flow and Difusion" # "Pure Difusion" # 
numCores = 10    # If num Cores == 1, Runs in Series
                # Parallel computing is indicated for "Flow and Difusion" simulations. 
                # "Pure Difusion" simulations ignores this parameter and runs in series.

# Gravitational Acceleration ( (0.0,0.0) if neglegted )
g = Constant((0.0,0.0))               # [m/s2]

# Fluid Properties
# Density            
rho_Displaced = 1000        # [kg/m3]
rho_Displacer = 1000      # [kg/m3]

# Viscosity
mu_Displaced = 1e-1        # [Pa.s]
mu_Displacer = 1e-3      # [Pa.s]

# Diffusion Coefficient
# No Mixture
D_NoMixture = 1e-11;       # [m²/s]
# Test small mixture
D_SmallMixture = 1e-9;     # [m²/s]

# Fluid Properties set for Simulation
# Fluid 1: Displaced
Fluid1Tag = 0   #  Inside Chamber (do not change)
rho_1 = rho_Displaced       # [kg/m3]
mu_1 = mu_Displaced         # [Pa.s]

# Fluid 2: Displacer
Fluid2Tag = 1   # Outside chamber (do not change)
rho_2 = rho_Displacer        # [kg/m3] 
mu_2 = mu_Displacer    # [Pa.s]

# Diffusion of Fluid 2 in Fluid 1--------------------------------------------
D = D_NoMixture        # [m²/s]

# Mesh Elements
# Velocity
velocityElementfamily = 'Lagrange'
velocityElementOrder = 2
# Pressure
pressureElementfamily = 'Lagrange'
pressureElementOrder = 1
# Concentration
concentrationElementfamily = 'Lagrange'
concentrationElementOrder = 1

# Solver Parameters
absTol = 1e-10          # absolute tolerance: residual value
relTol = 1e-12          # relative tolerance: change with respect to previous
maxIter = 15            # Maximum iterations for non-linear solver
nlinSolver = 'newton'   # Non-Linear Solver(Coupled Pressure/Velocity)
linSolver = 'mumps'     # Linear Solver(Concentration)
alpha = 0.75            # Flow relaxation coefficient
alphaC = 0.8            # Mass Transport relaxation coefficient

# Simulation Time
tEnd =  100 # [s]
# Auto-time step controls 
# Passing equal dt0=dtMin=dtMax, auto-timestep is turned off.
dt0 = 1e-6              # [s]
dtMin = 1e-8            # [s]
dtMax = 1               # [s]
# If iterations for convergence <= minDTIter => increase(DOUBLE) time step if possible
# If iterations for convergence > minDTIter + deltaDTIter => reduce time step(HALF) if possible
minDtIter = 3           # Positive integer 
deltaDtIter = 1         # Positive integer 
# Saving Time Step
dtSav = dtMax*1       # [s]
# Multiple Saving Time Steps
# tChangeSave = 10        # [s]
# dtSav = {0:       1e-3,
#          1:       1,
#         tChangeSave: dtMax*3}

# Pressure Boundary Conditions
psi2Pa = 6894.76
P0 = 12000           # [psi]
## Pressure Difference
Pout = 12000         # [psi]
Pin = Pout + 0.1/psi2Pa  # [psi] 25 bbl/min
# Pin = 12000.0           # [psi] No Flow
# Pa Conversion
P0 = P0; Pout = Pout; Pin = Pin*psi2Pa

## Fluids Concentration Difference
Cin = 0.0 # [s : mol/m³]
Cout = 1000.0 # [s : mol/m³]

# Concentraton
xChange = 0.5
C0 = Expression('( (CMax-CMin) / (1 + exp( IntIncl*(-x[2]+x0) ) ) ) + CMin',degree=2,IntIncl = 20,x0=xChange,CMax=Cout,CMin=Cin)

# Constant Initial Conditions
U0x = 0.0                       # [m/s]
U0y = 0.0                       # [m/s]
U0z = 0.0                       # [m/s]

vAxialInlet =  0.363
vRadialInlet = 0.0              # [m/s]
vInlet =  vAxialInlet
# vInlet =  Expression('( 0.0, 0.0, ( (vxMax) / (1 + exp( IntIncl*(-t+t0) ) ) ) )',degree=2,t=0,IntIncl = 1,t0=10,vxMax=vAxialInlet)