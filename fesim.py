#!/bin/python3

# Libraries Import
import matplotlib.pyplot as plt
import numpy as np
import inputs as inputs
import dolfin as df
import os, sys
import timeit
import datetime
from mpi4py import MPI

# Class Definitions
class SavingData():
    def __init__(self):
        self.ParaViewFilenames = [];                    self.ParaViewTitles = []
        self.ParaViewFilenames.append("velocity");      self.ParaViewTitles.append('Velocity (m/s)')
        self.ParaViewFilenames.append("pressure");      self.ParaViewTitles.append('Pressure (Pa)')
        self.ParaViewFilenames.append("concentration"); self.ParaViewTitles.append('Mass Fraction (Fluid Tags)')
        # self.ParaViewFilenames.append("condutivity");    self.ParaViewTitles.append('Condutivity (-)')
        # self.ParaViewFilenames.append("reynolds");       self.ParaViewTitles.append('Re (-)')
        
        #Output Files
        self.paraFullPath = inputs.resPath
        self.paraFiles = savingPreparation(self.paraFullPath,self.ParaViewFilenames)

# Function Definitions
# Variable Boundary Condition
def updateBC(bcInputs, conditionalVariable):


    return bc

# Read subdomains from .msh file
def readDomains(inPath,inFile):
    # Read .msh File
    fid = open(inPath+inFile+'.msh', 'r')
    # Initialize variables
    found = 0
    finished = 0
    physicalNames = {}
    # Loop througn .msh lines
    for line in fid:
        if '$EndPhysicalNames' in line:
            finished == 1
            break
        elif '$PhysicalNames' in line:
            found = 1
        elif found==1 and finished == 0:
            word=line.split()
            if len(word)==3:
                physicalNames[word[2][1:len(word[2])-1]] = int(word[1])

    return physicalNames


#%% 3D->2D Strain Rate Tensor 
def DD(u):
    #Cartesian
    D = df.sym(df.grad(u))
    
    # Cylindrical Coordinates(RZ)
    # D = df.sym(df.as_tensor([[u[0].dx(0), 0, u[0].dx(1)],
                            # [0, u[0]/x[0], 0],
                            # [u[1].dx(0), 0, u[1].dx(1)]]))
    return D
# Stress Tensor
def TT(u, p, mu):
    #Cartesian
    T = 2*mu*DD(u) - p*df.Identity(len(u))

    #Cylindrical(RZ)
    # T = 2*mu*DD(u,x) - p*Identity(3)
    return T

# Cylindrical coordinates: Spatial Integration Domains dx(bulk) and ds(walls)
def dxCyl(x):
    return 2*np.pi*x[0]
def dsCyl(x):
    return 2*np.pi*x[0]

# Calculate Auto-time step
def autoTimestep(no_iterations,dt,inputs,increment=2):

    # Check if 
    if no_iterations <= inputs.minDtIter:
        # Increase timestep if possible
        dt = min(increment*dt, inputs.dtMax)

    elif no_iterations > inputs.minDtIter + inputs.deltaDtIter:
        # reduce timestep if necessary
        dt = max(inputs.dtMin, (1/increment)*dt)

    else:
        # Keep the timestep - Passing equal dtMin==dtMax, auto-timestep is turned off.
        dt = dt
    
    return dt

# ProgressBar
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '=', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()



## Fluids smooth Interface
def initialInterface(C,inputs):
    c0 = df.Function(C)
    # c0.assign(df.project(inputs.smoothstep,C))
    c0.assign(df.project(inputs.C0,C))
    return c0

# Calculates Fluid properties by Mesh Cell
def assignFluidProperties(inputs,C,c): #,u,t,x):
    
    # Mixture Density
    # rho = df.Constant(inputs.rho)
    rho = df.project(df.Expression('rho_1*c + rho_2*(1-c)',degree=1,\
                            c=c,rho_1=inputs.rho_1,rho_2=inputs.rho_2),C)

    # Mixture Viscosity (Can be replaced by more complex rheological model)
    # mu = df.Constant(inputs.mu)
    mu = df.project(df.Expression('mu_1*c + mu_2*(1-c)',degree=1,\
                            c=c,mu_1=inputs.mu_1,mu_2=inputs.mu_2),C)
    
    # Diffusivity
    D = df.Constant(inputs.D)
    # Reynolds Calculation
    # Dcar = 2*(inputs.ROut-inputs.RIn)
    # Re = project((rho*np.abs(u[1])*Dcar)/mu,C)
    
    return rho, mu, D
    
# Paraview Results
def savingPreparation(paraFullPath,ParaViewFilenames):
    from dolfin import File,XDMFFile
    
    # Create New Result Files (Paraview Files: .pvd)
    fileNames = []
    paraFiles = []
    for i in range(0,len(ParaViewFilenames)): # ParaViewFilenames defined in Problem Imputs
        if inputs.numCores <=1:
            fileNames.append(ParaViewFilenames[i]+".pvd")
            paraFiles.append(File(paraFullPath+fileNames[i]))
        else:

            fileNames.append(ParaViewFilenames[i]+".xdmf")
            paraFiles.append(XDMFFile(paraFullPath+fileNames[i]))

    return paraFiles

# Save Current Results
def saveResults(results,paraFiles,t,inputs,ParaViewFilenames,ParaViewTitles):
    # Save Data for Paraview
    for i in range(0,len(results)):
        
        if inputs.numCores <=1:
            # Rename output Titles
            results[i].rename(ParaViewFilenames[i],ParaViewTitles[i])
            paraFiles[i] << results[i]
        else:
            results[i].rename(ParaViewFilenames[i],ParaViewTitles[i])
            paraFiles[i].write(results[i],t)
    # Save time conditions
    if type(inputs.dtSav) == float:
        dtSav = inputs.dtSav
    else: 
        for key in inputs.dtSav.keys():
            if t >= float(key):
                dtSav = inputs.dtSav[key]
                break
            else:
                pass
    

    return t + dtSav

def main(mpi_comm,inputs,savingData):

    # Disable dolfin log
    df.set_log_active(False)

    # Set lastTry flag to False
    lastTry = False
    # Set next saving time
    if type(inputs.dtSav) == float:
        nextSavTime = inputs.dtSav
    else:
        nextSavTime = inputs.dtSav[0]

    #  Mesh Reading
    # Load Subdomains
    Subdomains = readDomains(inputs.meshPath,inputs.meshFile)
    print(Subdomains)

    # if inputs.numCores <= 1:
        # # Option 2 - Gmsh Serial XML Converted Files
        # meshObj = df.Mesh(inputs.meshPath+inputs.meshFile +'.xml')
    
        # # Initialize boundaries (inlet, outlet, etc...)
        # boundaries = df.MeshFunction('size_t',meshObj,inputs.meshPath+inputs.meshFile + "_facet_region.xml")
        
        # # Initialize subdomain (fluid)
        # markers = df.MeshFunction('size_t',meshObj,inputs.meshPath + inputs.meshFile + '_physical_region.xml')
    # else:
        # Option 3 - Parallel H5 Converted Files
    meshObj = df.Mesh()
    hdf = df.HDF5File(meshObj.mpi_comm(), inputs.meshPath + inputs.meshFile+".h5", "r")
    hdf.read(meshObj, "/mesh", False)
    
    # Initialize subdomain (fluid)
    markers = df.MeshFunction("size_t",meshObj,meshObj.topology().dim())
    hdf.read(markers, "/subdomains")

    # Initialize boundaries (inlet, outlet, etc...)
    boundaries = df.MeshFunction("size_t",meshObj,meshObj.topology().dim()-1)
    hdf.read(boundaries, "/boundaries")


    # Open new figure
    # plt.figure(figsize=(20, 10), dpi=180, facecolor='w', edgecolor='k')
    
    # plot Mesh and save image
    # plt.figure(figsize=(20, 10), dpi=180, facecolor='w', edgecolor='k')
    # df.plot(meshObj)
    # plt.savefig(inputs.resPath+'Mesh.png')

    # Get Element Shape from mesh: Triangle, etc...
    elementShape = meshObj.ufl_cell()

    # Set Mesh Elements
    # Mixed VelocityXPressure Element
    Uel = df.VectorElement(inputs.velocityElementfamily, elementShape, inputs.velocityElementOrder) # Velocity vector field
    Pel = df.FiniteElement(inputs.pressureElementfamily, elementShape, inputs.pressureElementOrder) # Pressure field
    UPel = df.MixedElement([Uel,Pel])
    # Concentration Element
    Cel = df.FiniteElement(inputs.concentrationElementfamily, elementShape, inputs.concentrationElementOrder) # Concentration field

    # Define any measure associated with domain and subdomains
    dx = df.Measure('dx', domain=meshObj)
    ds = df.Measure('ds', domain=meshObj, subdomain_data=boundaries)

    # Vectors Normal to the Mesh
    n = df.FacetNormal(meshObj) # Normal vector to mesh facets
    
    # Function Spaces: 
    ## Flow: Mixed Function Space for Pressure and Velocity
    W = df.FunctionSpace(meshObj,UPel)
    ### Split into Velocity and Pressure
    (U, P) = W.split()
    
    ### Trial and Test function(s)
    dw = df.TrialFunction(W)       # non-linear solver particularity
    (v, q) = df.TestFunctions(W)   # v -> u and q -> p
    w = df.Function(W)             # current time step
    (u, p) = (df.as_vector((w[0], w[1], w[2])), w[3])
    w1 = df.Function(W)

    w0 = df.Function(W)            # previous time step
    (u0, p0) = w0.leaf_node().split()

    # ## Diffusion: Function Space for Concentration
    C = df.FunctionSpace(meshObj,Cel)
    # ### Trial and Test function(s)
    c = df.TrialFunction(C)         # current time step
    l = df.TestFunction(C)          
    c_1 = df.Function(C)            # result functions
    
    # Uniform Boundary Conditions
    
    ## Dirichlet 
    ### No-slip wall velocity conditions
    # bcU1 = df.DirichletBC(U,df.Constant((0.0,0.0)),boundaries,Subdomains['BottomWall'])
    # bcU2 = df.DirichletBC(U,df.Constant((0.0,0.0)),boundaries,Subdomains['TopWall'])
    # bcU3 = df.DirichletBC(U,df.Constant((0.0,0.0)),boundaries,Subdomains['InnerWalls'])
    bcU1 = df.DirichletBC(U,df.Constant((0.0, 0.0, 0.0)),boundaries,Subdomains['wall'])
    # bcU4 = df.DirichletBC(U,df.Constant((inputs.vAxialInlet,inputs.vRadialInlet)),boundaries,Subdomains['Inlets'])
    ### Inlet velocity conditions
    # try:
    #     inputs.vInlet.t = 0
    #     bcU4 = df.DirichletBC(U,inputs.vInlet,boundaries,Subdomains['Inlets'])
    #     varVBC = True
    # except:
    #     bcU4 = df.DirichletBC(U,df.Constant((inputs.vAxialInlet,inputs.vRadialInlet)),boundaries,Subdomains['Inlets'])
    #     varVBC = False

    bcU = [bcU1]
    # bcU3 = df.DirichletBC(U,df.Constant((inputs.vAxialInlet,inputs.vRadialInlet)),boundaries,Subdomains['Inlets'])
    # bcU = [bcU1,bcU2,bcU3]

    ### Concentration inlet conditions
    if type(inputs.Cin) == float:
        bcC1 = df.DirichletBC(C,df.Constant(inputs.Cin),boundaries,Subdomains['inlet'])
        bcC = [bcC1]
    else:
        inputs.Cin.t = 0
        bcC1 = df.DirichletBC(C,inputs.Cin,boundaries,Subdomains['inlet'])
        # varBC = True

    # Initial Conditions
    ## Initial Time Step   
    t = 0 
    dt = inputs.dt0
    ## Initial Flow Conditions
    e_u0 = df.Constant((inputs.U0x, inputs.U0y, inputs.U0z))
    e_p0 = df.Constant(inputs.P0)
    u0 = df.interpolate(e_u0, W.sub(0).collapse())
    p0 = df.interpolate(e_p0, W.sub(1).collapse())
    df.assign(w0, [u0, p0])

    ## Initial Interface Conditions
    c_0 = initialInterface(C,inputs) # previous time step
    
    # Time execution time
    start = timeit.default_timer()

    while t <= inputs.tEnd and dt>0.0:

        # Run only once
        if rank == 0:
            sys.stdout = inputs.console
            
            if inputs.numCores <=1 or inputs.simulationType ==  "Pure Difusion":
                # Plot Progress Bar
                printProgressBar (t, inputs.tEnd)
            else:
                print(' ')

            sys.stdout = open(os.getcwd()+'/res/'+inputs.caseId+'/Execution.txt', 'a')
        
        # Initialize results Array
        results = []

        # Previous Time Step
        (u0, p0) = w0.leaf_node().split()
        
        # Update Fluid Properties
        rho,mu,D = assignFluidProperties(inputs,C,c_0)
        # Constant Properties
        # rho,mu,D = inputs.rho, inputs.mu, inputs.D
        
        # Time Step 
        Dt = df.Constant(dt)
        
        if rank == 0:
            sys.stdout = open(os.getcwd()+'/res/'+inputs.caseId+'/Execution.txt', 'a')
            print('Time:{:.3f}s and dt:{:.5f}s - Simulation at {:.2f}%'.format(t,dt,t/inputs.tEnd))
        
        # Check simulationType to solve flow or solve only concentration.
        if inputs.simulationType == "Pure Difusion":
            converged = True
            no_iterations = 1
            w=w0
            w1.assign(w)
        else:
            # Flow Problem
            ## Transient Flow Bilinear Form:
            ### Linear Momentum Conservation: F = a - L
            #         # Time dependent term                      # Inertia Term                                     # Cauchy Stresses Term  
            a1 = rho*df.dot((u-u0)/Dt,v)*dx() + inputs.alpha *(rho*df.dot(df.dot(u ,df.grad(u) ),v) + df.inner(TT(u,p,mu),DD(v)))*dx() + \
                                                (1-inputs.alpha)*(rho*df.dot(df.dot(u0,df.grad(u0)),v) + df.inner(TT(u0,p0,mu),DD(v)))*dx()  # Relaxation 
            # Pressure Boundary Conditions: Naturally Stated
                    # Inlet Pressure                                    # Outlet Pressure                                      # Gravity
            L1 = - (inputs.Pin)*df.dot(n,v)*ds(Subdomains['inlet']) - (inputs.Pout)*df.dot(n,v)*ds(Subdomains['outlet']) + df.inner(rho*inputs.g,v)*dx()
                    # Inlet Flowrate 
            # L1 = 0 - (inputs.Pout)*df.dot(n,v)*ds(Subdomains['Outlets']) #+ df.inner(rho*inputs.g,v)*dx()
                    # Cavity
            # L1 = 0 - (inputs.Pout)*df.dot(n,v)*ds(Subdomains['Outlets']) #+ df.inner(rho*inputs.g,v)*dx()
            
            ## Mass Conservation(Continuity)
            a2 = (q*df.div(u))*dx()
            L2 = 0
            
            ## Complete Weak Form
            F = (a1 + a2) - (L1 + L2)
            
            ## Jacobian Matrix Calculation
            J = df.derivative(F,w,dw)

            # if varVBC:
            #     inputs.vInlet.t = 0
            #     bcU4 = df.DirichletBC(U,inputs.vInlet,boundaries,Subdomains['Inlets'])

            # bcU[3] = bcU4
            ##########   Numerical Solver Properties
            # Problem and Solver definitions
            problemU = df.NonlinearVariationalProblem(F,w,bcU,J)
            solverU = df.NonlinearVariationalSolver(problemU)
            # # Solver Parameters
            prmU = solverU.parameters
            # #info(prmU,True)  #get full info on the parameters
            prmU['nonlinear_solver'] = inputs.nlinSolver
            prmU['newton_solver']['absolute_tolerance'] = inputs.absTol
            prmU['newton_solver']['relative_tolerance'] = inputs.relTol
            prmU['newton_solver']['maximum_iterations'] = inputs.maxIter
            prmU['newton_solver']['linear_solver'] = inputs.linSolver

            # Solve Flow Problem        
            try:
                (no_iterations,converged) = solverU.solve()
                if rank == 0:
                    print('Solved Flow - Iterations:{:.0f}'.format(no_iterations))
                
                # Update previous time step velocity & pressure values
                w1.assign(w)
                
            except:
                no_iterations = inputs.maxIter
                if rank == 0:
                    print('Could''t Solve Flow = Iterations:{:.0f}'.format(no_iterations))
                converged=False
                # w=w0
                

        if converged:
            (u1, p1) = w1.leaf_node().split()
            
            # Concentration Equation
            ## Cartesian
                            # Transient Term   #                  # Advection Term                                                                              # Diffusion Term                            
            Fc = df.inner((rho*c - rho*c_0)/Dt,l)*dx() + inputs.alphaC *(df.inner(u1,(df.grad(c * rho))*l) + df.dot(u1, df.grad(c )) *l + c * df.div(u1)*l + (D)*df.dot(df.grad(c) , df.grad(l)))*dx() +\
                                                    (1- inputs.alphaC)*(df.inner(u1,(df.grad(c_0*rho))*l)+ df.dot(u1, df.grad(c_0))*l + c_0*df.div(u1)*l + (D)*df.dot(df.grad(c_0), df.grad(l)))*dx() # Relaxation
            # Cylindrical
            #             # Transient Term   #                                                   # Advection Term                                       # Diffusion Term                            
            # Fc = inner((rho_i_t*c_i - rho_i_t0*c_i0)/Dt,l)*dx() + alphaC*(inner(rho_i_t *u1,(grad(c_i ))*l) +  c_i *rho_i_t*div2d(u1,x)*l + (D)*dot(grad(c_i ),grad(l)))*dxCil(x)*dx() +\
            #                                                     (1 - alphaC)*(inner(rho_i_t0*u1,(grad(c_i0))*l) + c_i0 *rho_i_t0*div2d(u1,x)*l + (D)*dot(grad(c_i0),grad(l)))*dxCil(x)*dx() # Relaxation
            ac, Lc = df.lhs(Fc), df.rhs(Fc)
            
            # Variable Boundary Condition
            # if varBC:
            #     # Variable Inlet Concentration Condition
            #     inputs.Cin.t = t
            #     bcCin = df.DirichletBC(C,inputs.Cin,boundaries,Subdomains['Inlets'])
            #     for key in inputs.pHin.keys():
            #         if t >= key:
            #             pHin = inputs.pHin[key]
            #             break
            #         else:
            #             pass
                
                        
            #     bcC = [bcCin]

            # Solve Mass Transport Problem
            # Variable Boundary Condition
            # if varBC:
            #     if rank == 0:
            #         print('pHin: {:.2f} '.format(pHin))
            # df.solve(ac == Lc, c_1, bcC)
            # if rank == 0:
            #     print('Solved Mass Transport')
            
            # Calculate pH from Ion Concentration (mol/m³)
            pHExpression = df.Expression('-log10(Cion/1000)',Cion=c_1, degree=2)
            

            # Save Paraview Files
            if t == 0:
                if rank == 0:
                    print('------ Saving Initial Conditions -------')
                ## Paraview Results
                results.append(u0)       #0
                results.append(p0)       #1
                results.append(c_0)      #2
                # results.append(Re)     #3
            
                nextSavTime0 = saveResults(results,savingData.paraFiles,t,inputs,\
                                        savingData.ParaViewFilenames,savingData.ParaViewTitles)    
                
            elif t >= nextSavTime or t==inputs.tEnd:
                if rank == 0:
                    print('---------------- Saving ----------------')
                # 
                # Append and save Results
                ## Paraview Results
                results.append(u1)       #0
                results.append(p1)       #1
                results.append(c_1)      #2
                # results.append(Re)     #3
                
                nextSavTime = saveResults(results,savingData.paraFiles,t,inputs,\
                                            savingData.ParaViewFilenames,savingData.ParaViewTitles)
            
            t = t+dt
            # Update previous time step concentration, pressure and velocity values
            c_0.assign(c_1)
            w0.assign(w)

        # Auto Time Steps
        dt = autoTimestep(no_iterations,dt,inputs)
        if dt <= inputs.dtMin and lastTry:
            if rank == 0:
                print('Did not converge: Minimum timestep achieved')
            # df.end()
            break
        elif dt <= inputs.dtMin:
            # Set lastTry flag to True
            lastTry = True
        
        # Adjust for last time step end with tEnd
        if t+dt>inputs.tEnd:
            dt = inputs.tEnd-t

    #####################  Post Processing
    # Update execution End Time
    stop = timeit.default_timer()
    total_time = stop - start
    
    return total_time
    

    # Test plots
    # U:
    # plt.figure(figsize=(20, 10), dpi=180, facecolor='w', edgecolor='k')
    # df.plot(u1)
    # plt.savefig(inputs.resPath+'Velocity.png')
    # # P:
    # plt.figure(figsize=(20, 10), dpi=180, facecolor='w', edgecolor='k')
    # df.plot(p1)
    # plt.savefig(inputs.resPath+'Pressure.png')
    # # C:
    # plt.figure(figsize=(20, 10), dpi=180, facecolor='w', edgecolor='k')
    # df.plot(c_1)
    # plt.savefig(inputs.resPath+'Concentration.png')
    
if __name__ == '__main__':
    
    mpi_comm = MPI.COMM_WORLD
    rank = mpi_comm.Get_rank()

    # Run only once
    if rank == 0:
        # Change output to file
        inputs.console = sys.stdout
        print('Starting Simulation - Log is in Execution.txt file inside created folder.')
        sys.stdout = open(os.getcwd()+'/res/'+inputs.caseId+'/Execution.txt', 'w')
        # Store Simulation Execution Info
        print('--------------------------------------------------------------------------')
        print('Simulation Execution: '+str(datetime.datetime.now()))
        print('Case: '+inputs.caseId)
        print('Source: '+os.path.basename(os.getcwd()))
        print('--------------------------------------------------------------------------')
        sys.sdtout = inputs.console

    # Output Variables
    savingData = SavingData()
    
    # Run main in parallel
    total_time = main(mpi_comm,inputs,savingData)

    # Execution time calculation.
    mins, secs = divmod(total_time, 60)
    hours, mins = divmod(mins, 60)
    
    # Run only once
    if rank == 0:
        # Change output to file
        print('--------------------------------------------------------------------------')
        print("Total running time: %dh:%dmin:%ds \n" % (hours, mins, secs))
        print('Finished')
        print('--------------------------------------------------------------------------')
        sys.sdtout = inputs.console

