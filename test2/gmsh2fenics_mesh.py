#
# Convert GMSH mesh to FEniCS xml mesh.
#
# We use meshio, https://pypi.org/project/meshio/
# Installing with conda or pip:
#```shell
#   conda install pip
#   pip install meshio
#```

import sys
import meshio

if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} <GMSH input file> <FEniCS output file>")

input_mesh_file = sys.argv[1]
fenics_mesh_file = sys.argv[2]

mesh = meshio.read(input_mesh_file)
mesh.write(fenics_mesh_file)
