from dolfin import Expression, Constant

# Inputs
# Mesh File
caseId = 'SimplePipeTest_parallel_16cores_VIn_29k' 
# caseId = 'PureDifusion_HClBrine_229kEl_varBC_12h_12h'
meshPath = 'Mesh/'
meshFile = 'pipeMesh_tri_29k'
resPath = 'res/'+caseId+'/' 
replaceGeometry = False
# Choose between "Pure Difusion" OR "Flow and Difusion"
simulationType = "Flow and Difusion" # "Pure Difusion" # 
numCores = 16    # If num Cores == 1, Runs in Series
                # Parallel computing is indicated for "Flow and Difusion" simulations. 
                # "Pure Difusion" simulations ignores this parameter and runs in series.

# Gravitational Acceleration ( (0.0,0.0) if neglegted )
g = Constant((0.0,0.0,0.0))               # [m/s2]

# Fluid Properties measured at LMMP - PUC-Rio
# Density            
rho_Displacer = 1000      # [kg/m3]
rho_Displaced = 1000      # [kg/m3]

# Viscosity
mu_Displacer = 1.8        # [Pa.s]
mu_Displaced = 8.1        # [Pa.s]


# Diffusion Coefficient
# Fluid 1 in Fluid 2
D_NoMixture = 1.87e-9;     # [m²/s]

# Fluid Properties set for Simulation
# Fluid 1: Inside the pipe
Fluid1Tag = 0   #  Inside the pipe
rho_1 = rho_Displaced       # [kg/m3]
mu_1 = mu_Displaced         # [Pa.s]

# Fluid 2: Entering the Pipe
Fluid2Tag = 1   # Entering the pipe
rho_2 = rho_Displacer        # [kg/m3] 
mu_2 = mu_Displacer          # [Pa.s]

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
alpha = 0.9             # Flow relaxation coefficient
alphaC = 0.7            # Mass Transport relaxation coefficient

# Simulation Time
tEnd = 0.5        # [s]
# Auto-time step controls 
# Passing equal dt0=dtMin=dtMax, auto-timestep is turned off.
dt0 = 1e-6              # [s]
dtMin = 1e-8         # [s]
dtMax = 5e-2        # [s]
# If iterations for convergence <= minDTIter => increase(DOUBLE) time step if possible
# If iterations for convergence > minDTIter + deltaDTIter => reduce time step(HALF) if possible
minDtIter = 1.0           # Positive integer 
deltaDtIter = 3.0         # Positive integer 
# Saving Time Step
# dtSav = dtMax*1     # [s]
tChangeSave = 10 # [s]
dtSav =    {0: 1e-6,
            1e-3: 1e-3}

# Pressure Boundary Conditions
P0 = 0           # [Pa]
## Pressure Difference
Pout = P0         # [Pa]
Pin = Pout + 225.04*10  # [Pa] 25 bbl/min
Cin = 0.0
Cout = 1000.0
# Constant condition for Inlet Concentration
# Cin = pHin
# Variable condition for Inlet Concentration

# Constant Initial Conditions
U0x = 0.0                       # [m/s]
U0y = 0.0                       # [m/s]
U0z = 0.0                       # [m/s]
C0 = Expression('( (CMax-CMin) / (1 + exp( IntIncl*(-x[2]+x0) ) ) ) + CMin',degree=2,IntIncl = 20,x0=0.5,CMax=Cout,CMin=Cin)
vAxialInlet = 0.363           # [m/s]