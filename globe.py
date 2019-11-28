import pdb
import tempfile,os
import numpy as np
import requests
from mayavi import mlab
from tvtk.api import tvtk # python wrappers for the C++ vtk ecosystem

def map_image_to_unit_sphere_actor(image_file):
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
    return sphere_actor

def lut_set_transparent_below_range(lut):
    """Set data below vmin to a transparent color"""
    lut.use_below_range_color=1
    lut.below_range_color=np.array([0.,0.,0.,0.])

def lut_set_alpha_range(lut,amin=0.,amax=1.):
    """following
    https://docs.enthought.com/mayavi/mayavi/auto/example_custom_colormap.html
    with option of range of alphas other than 0 to 1
    """
    lut_table = lut.table.to_array()
    n_colors = lut_table.shape[0]
    # docs say look up table (lut.table) needs integers
    iamin = np.floor(255*amin)
    iamax = np.floor(255*amax)
    lut_table[:,-1] = np.floor(np.linspace(iamin,iamax,n_colors))
    lut.table = lut_table

def legend(scalar_lut_manager,label=None):
    scalar_lut_manager.show_scalar_bar = True
    scalar_lut_manager.show_legend = True
    scalar_lut_manager.data_name = label if label is not None else u''

def geosurf(glats,glons,geodata,alt=300.,
                vmin=None,vmax=None,cmap='viridis',
                colorbar=False,label=None):
    """Draw gridded data in geographic coordinates as a surface (mesh) """
    Re = 6371.2
    r = (Re+alt)/Re
    th = np.radians(90.-glats)
    ph = np.radians(glons)
    x = r*np.sin(th)*np.cos(ph)
    y = r*np.sin(th)*np.sin(ph)
    z = r*np.cos(th)
    mesh = mlab.mesh(x,y,z,
                        scalars=geodata,
                        vmin=vmin,
                        vmax=vmax,
                        colormap=cmap)
    lut_set_transparent_below_range(mesh.parent.scalar_lut_manager.lut)
    lut_set_alpha_range(mesh.parent.scalar_lut_manager.lut,amin=.5,amax=1.)
    if colorbar:
        legend(mesh.parent.scalar_lut_manager,label=label)
    mlab.draw()
    #mlab.process_ui_events()

def download_blue_marble():
    """from https://visibleearth.nasa.gov/collection/1484/blue-marble"""
    imgfn = os.path.join(tempfile.gettempdir(),'blue_marble.jpg')
    if not os.path.exists(imgfn):
        print('Downloading NASA blue marble image...')
        url = ('https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73909/'
              +'world.topo.bathy.200412.3x5400x2700.jpg')
        r = requests.get(url, allow_redirects=True)
        open(imgfn,'wb').write(r.content)
    return imgfn

def download_auroral_nowcast():
    nowcastfn = os.path.join(tempfile.gettempdir(),'auroral_nowcast.txt')
    url = 'https://services.swpc.noaa.gov/text/aurora-nowcast-map.txt'
    r = requests.get(url,allow_redirects=True)
    open(nowcastfn,'w').write(r.content)
    return nowcastfn

def visual_test_geosurf():
    """Make a plot of the current NOAA auroral probability on the NASA
    blue marble image
    """
    imagefn = download_blue_marble()
    nowcasttxtfn = download_auroral_nowcast()

    # create a figure window (and scene)
    fig = mlab.figure(size=(600, 600))
    sphere_actor = map_image_to_unit_sphere_actor(imagefn)
    fig.scene.add_actor(sphere_actor)
    
    glats,glons = np.mgrid[-90:90:512j,0:360:1024j]
    aur_prob = np.genfromtxt(nowcasttxtfn).reshape((512,1024))
    geosurf(glats,glons,aur_prob,
                label='Probability of Aurora',
                cmap='viridis',vmin=1.,vmax=15.,
                colorbar=True)
    mlab.show()

if __name__ == "__main__":
    visual_test_geosurf()
