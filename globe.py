import numpy as np
from mayavi import mlab
from tvtk.api import tvtk # python wrappers for the C++ vtk ecosystem

def auto_sphere(image_file):
    # create a figure window (and scene)
    fig = mlab.figure(size=(600, 600))

    # load and map the texture
    img = tvtk.JPEGReader()
    img.file_name = image_file
    texture = tvtk.Texture(input_connection=img.output_port) #, interpolate=1
    # (interpolate for a less raster appearance when zoomed in)

    # use a TexturedSphereSource, a.k.a. getting our hands dirty
    R = 1
    Nrad = 180

    # create the sphere source with a given radius and angular resolution
    sphere = tvtk.TexturedSphereSource(radius=R, theta_resolution=Nrad,
                                       phi_resolution=Nrad)

    # assemble rest of the pipeline, assign texture    
    sphere_mapper = tvtk.PolyDataMapper(input_connection=sphere.output_port)
    sphere_actor = tvtk.Actor(mapper=sphere_mapper, texture=texture)
    fig.scene.add_actor(sphere_actor)
    return fig
    
def overlay_noaa_nowcast():
	lats,lons = np.mgrid[-90:90:512j,0:360:1024j]
	r = (6371.2+300.)/6371.2
	th = np.radians(90.-lats)
	ph = np.radians(lons)
	x = r*np.sin(th)*np.cos(ph)
	y = r*np.sin(th)*np.sin(ph)
	z = r*np.cos(th)
	aur_prob = np.genfromtxt('aurora-nowcast-map.txt').reshape((512,1024))

	mesh = mlab.mesh(x,y,z,scalars=aur_prob,transparent=True,vmin=1.,vmax=np.nanmax(aur_prob))
	mesh.parent.scalar_lut_manager.lut.alpha_range = np.array([0.,.8])
	mlab.axes()

if __name__ == "__main__":
    image_file = 'blue_marble.jpg'
    fig = auto_sphere(image_file)
    overlay_noaa_nowcast()
    mlab.show()