#!/usr/bin/python3

import numpy
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.cm as cm
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def read_vtk_file(filename, vtuFile=True):
    ### 1.1 Read data
    if vtuFile: # read vtu file
        reader = vtk.vtkXMLUnstructuredGridReader()
    else: # read vtk file
        reader = vtk.vtkUnstructuredGridReader()

    gridreader = vtk.vtkXMLUnstructuredGridReader()
    gridreader.SetFileName(filename)
    gridreader.Update()
    vtkOut = gridreader.GetOutput() # get an vtkUnstructuredGridReader
    return vtkOut

def get_points(gridOutput):
    vtkData = gridOutput.GetPoints().GetData()
    coords = numpy.array([vtkData.GetTuple3(x)
                          for x in range(vtkData.GetNumberOfTuples())])
    return coords

def get_connectivity_matrix(gridOutput):
    cell_connectivity_matrix = []
    for i in range(gridOutput.GetNumberOfCells()):
        n_ids = gridOutput.GetCell(i).GetPointIds().GetNumberOfIds()
        assert gridOutput.GetCell(i).GetNumberOfPoints() == 3
        cell_connectivity_matrix.append(
            [gridOutput.GetCell(i).GetPointIds().GetId(j)
             for j in range(n_ids)] )
    cell_connectivity_matrix = numpy.array(cell_connectivity_matrix,
                                       dtype=numpy.float)
    return cell_connectivity_matrix

def get_vtk_data(data_id=0):
    # Get the first scalar in my vtk file
    scalar_vtk_array = gridOutput.GetPointData().GetArray(data_id)
    # Convert to numpy
    scalar_numpy_array = vtk_to_numpy(scalar_vtk_array)
    return scalar_numpy_array

def plot_mesh(coords, cell_connectivity):
    plt.triplot(coords[:, 0], coords[:, 1], triangles=cell_connectivity)
    plt.gcf().set_size_inches(16, 8)
    plt.gca().set_aspect('equal')
    plt.show()

def contour_plot(coords, cell_connectivity, data):
    # create an unstructured triangular grid instance
    nodes_x = coords[:, 0]
    nodes_y = coords[:, 1]
    triangulation = tri.Triangulation(nodes_x, nodes_y, cell_connectivity)
    # plt.tricontour(triangulation, data)
    plt.tricontourf(triangulation, data, cmap=cm.terrain)
    plt.show()

filename = "/tmp/testI000000.vtu"
gridOutput = read_vtk_file(filename)
coords = get_points(gridOutput)
cell_connectivity = get_connectivity_matrix(gridOutput)
data = get_vtk_data()

if False:
    plot_mesh(coords, cell_connectivity)

contour_plot(coords, cell_connectivity, data)
