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

# Convert .xml files to .HDF5
def xml2hdf5(inPath,inFile,outPath,outFile):
    mesh = df.Mesh(inPath+inFile+".xml")                                                  
    subdomains = df.MeshFunction("size_t", mesh, inPath+inFile+"_physical_region.xml")    
    boundaries = df.MeshFunction("size_t", mesh, inPath+inFile+"_facet_region.xml")
    hdf = df.HDF5File(mesh.mpi_comm(), outPath+outFile+".h5", "w")
    hdf.write(mesh, "/mesh")
    hdf.write(subdomains, "/subdomains")
    hdf.write(boundaries, "/boundaries")

# Convert .msh file to .xml files
def msh2xml(inPath,inFile,outPath,outFile):
    cmd = 'dolfin-convert '+inPath+inFile+'.msh '+outPath+outFile+'.xml'
    print(cmd)
    os.system(cmd)

def meshConversion(inputs):

    allXMLExist = ( os.path.isfile(inputs.meshPath + inputs.meshFile+'_facet_region.xml') and \
                    os.path.isfile(inputs.meshPath + inputs.meshFile+'_physical_region.xml') and \
                    os.path.isfile(inputs.meshPath + inputs.meshFile+'.xml'))
    
    allHDF5Exist = os.path.isfile(inputs.meshPath + inputs.meshFile+'.h5')

    # Recreate XML and H5 Files if needed
    if not allXMLExist or not allHDF5Exist or inputs.replaceGeometry:
        print('Regenerating mesh files.')

        # Convert .msh into xml
        msh2xml(inputs.meshPath,inputs.meshFile,\
                inputs.meshPath,inputs.meshFile)
        
        if inputs.numCores >1:
            # Convert .xml into hdf5
            xml2hdf5(inputs.meshPath,inputs.meshFile,\
                    inputs.meshPath,inputs.meshFile)          

        return True
    else:
        return False

