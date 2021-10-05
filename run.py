#!/bin/python3
"""
Created on Tue Jan 05 15:35:36 2021
@author: Sergio Ribeiro
@description: calls fesim.py execution
"""
from meshConversion import os, inputs, datetime, meshConversion
from shutil import copyfile
import errno



# Solve Command Compilation
def solveCommand(inputs):
    if inputs.numCores <=1 or inputs.simulationType == "Pure Difusion":
        # return 'python3 fesim.py>>res/'+inputs.caseId+'/Execution.txt'
        return 'python3 fesim.py'
    else:
        # return 'mpirun -np '+str(inputs.numCores)+' --allow-run-as-root python3 fesim.py>>res/'+inputs.caseId+'/Execution.txt'
        return 'mpirun -np '+str(inputs.numCores)+' --allow-run-as-root python3 fesim.py'


def main():
    try:
        os.makedirs('res/'+inputs.caseId,exist_ok=True)
    except OSError as e:
        print(e)
        pass

    try:
        meshConversion(inputs)
    except Exception as e:
        print('Unable to generate mesh files:'+e)       
    
    
    copyfile(os.path.abspath('inputs.py'),os.path.abspath('res/'+inputs.caseId+'/inputs.py'))

    # Skip Line
    print('\n')

    # Call Simulation and Append Execution.txt with Simulation Messages
    os.system(solveCommand(inputs))

if __name__ == '__main__':
    # Execute Main
    main()


