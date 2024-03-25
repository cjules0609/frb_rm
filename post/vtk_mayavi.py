from mayavi import mlab
from mayavi.modules.surface import Surface

dir = "/Volumes/T7Shield/ncfa/athena/out/"
file_name = "frb.block0.out1.00000.vtk"  # minimal example vtk file

# create a new figure, grab the engine that's created with it
fig = mlab.figure()
engine = mlab.get_engine()

# open the vtk file, let mayavi figure it all out
vtk_file_reader = engine.open(dir + file_name)

# plot surface corresponding to the data
surface = Surface()
engine.add_filter(surface, vtk_file_reader)

# block until figure is closed
mlab.show()