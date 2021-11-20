MESHDIR=../mesh
MESHFILE=spain_finer
# MESHFILE=spain_fine
# MESHFILE=spain_coarse
# MESHFILE=spain_coarser

all: test2.png DistribucionIndividuos.html

mesh: $(MESHDIR)/$(MESHFILE).msh
	dolfin-convert $< $@

$(MESHFILE).xml: $(MESHDIR)/$(MESHFILE).msh
	python3 gmsh2fenics_mesh.py $< $@

test2.png: test2.py
	python3 $<

DistribucionIndividuos.html: DistribucionIndividuos.ipynb
	jupyter-nbconvert --exec DistribucionIndividuos.ipynb
