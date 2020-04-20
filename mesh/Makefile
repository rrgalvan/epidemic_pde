COARSER_SCALE=0.8
COARSE_SCALE=0.5
FINE_SCALE=0.2
FINER_SCALE=0.075

TARGETS=spain_coarser.msh spain_coarse.msh spain_fine.msh spain_finer.msh

all: $(TARGETS)

spain_coarser.msh: natural_spain.geo
	gmsh -2 $< -clscale $(COARSER_SCALE) -o $@

spain_coarse.msh: natural_spain.geo
	gmsh -2 $< -clscale $(COARSE_SCALE) -o $@

spain_fine.msh: natural_spain.geo
	gmsh -2 $< -clscale $(FINE_SCALE) -o $@

spain_finer.msh: natural_spain.geo
	gmsh -2 $< -clscale $(FINER_SCALE) -o $@

clean:
	rm -f $(TARGETS)
