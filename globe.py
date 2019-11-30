import pdb
import tempfile,os,datetime
import numpy as np
import matplotlib.pyplot as plt
import requests
from mayavi import mlab
from tvtk.api import tvtk # python wrappers for the C++ vtk ecosystem

def solar_position_approx(dt,degrees=False):
    """
    This is from github.com/lkilcommons/geospacepy-lite, in
    astrodynamics.py. This was written by Liam Kilcommons
    From C.T. Russell, (1971) "Geophysical Coordinate Transformations", 
    Cosmic. Electrodyn. 2, 184-196
    ...
    G.D. Mead (private communication) has written a simple subroutine to\
    calculate the position of the Sun in GEI coordinates. It is accurate 
    for years 1901 through 2099, to within 0.006 deg. The input is the 
    year, day of year and seconds of the day in UT. The output is 
    Greenwich Mean Sideral Time in degrees, the ecliptic longitude, 
    apparent right ascension and declination of the Sun in degrees. 
    The listing of this program follows. We note that the cartesian 
    coordinates of the vector from the Earth to the Sun are:

      X = cos(SRASN) cos(SDEC)
      Y = sin(SRASN) cos(SDEC)
      Z = sin(SDEC)
      SUBROUTINE SUN(IYR, IDAY, SECS, GST, SLONG, SRASN, SDEC)
    C PROGRAM TO CALCULATE SIDEREAL TIME AND POSITION OF THE SUN. 
    C GOOD FOR YEARS 1901 THROUGH 2099. ACCURACY 0.006 DEGREE.
    C INPUT IS IYR, IDAY (INTEGERS), AND SECS, DEFINING UN. TIME. 
    C OUTPUT IS GREENWICH MEAN SIDEREAL TIME (GST) IN DEGREES,
    C LONGITUDE ALONG ECLIPTIC (SLONG), AND APPARENT RIGHT ASCENSION
    C AND DECLINATION (SRASN, SDEC) OF THE SUN, ALL IN DEGREES
      DATA RAD /57.29578/ 
      DOUBLE PRECISION DJ, FDAY 
      IF(IYR. LT. 1901. OR. IYR. GT. 2099) RETURN
      FDAY = SECS/86400
      DJ = 365* (IYR-1900) + (IYR-1901)/4 + IDAY + FDAY -0.5D0 
      T = DJ / 36525 
      VL = DMOD (279.696678 + 0.9856473354*DJ, 360.D0) 
      GST = DMOD (279.690983 + 0.9856473354*DJ + 360.*FDAY + 180., 360.D0)
      G = DMOD (358.475845 + 0.985600267*DJ, 360.D0) / RAD 
      SLONG = VL + (1.91946 -0.004789*T)*SIN(G) + 0.020094*SIN (2.*G) 
      OBLIQ = (23.45229 -0.0130125*T) / RAD 
      SLP = (SLONG -0.005686) / RAD 
      SIND = SIN (OBLIQ)*SIN (SLP) 
      COSD = SQRT(1.-SIND**2)
      SDEC = RAD * ATAN (SIND/COSD) 
      SRASN = 180. -RAD*ATAN2
      (COTAN (OBLIQ)*SIND/COSD, -COS (SLP)/COSD) 
      RETURN 
      END
    """
    iyear = dt.year
    iday = dt.timetuple().tm_yday
    secs = dt.hour*3600.+dt.minute*60.+dt.second
    fday = secs/86400.
    dj = 365*(iyear-1900)+(iyear-1901)/4 + iday + fday - .5
    t = dj/36525.
    vl = np.mod(279.696678 + 0.9856473354*dj, 360)
    gst = np.mod(279.690983 + 0.9856473354*dj + 360.*fday + 180., 360.)
    g = np.mod(358.475845 + 0.985600267*dj, 360.) * np.pi/180.
    slong = vl + (1.91946 -0.004789*t)*np.sin(g) + 0.020094*np.sin(2.*g) 
    obliq = (23.45229 -0.0130125*t) * np.pi/180. 
    slp = (slong - 0.005686) * np.pi/180.
    sin_d = np.sin(obliq)*np.sin(slp)
    cos_d = np.sqrt(1-sin_d**2)
    sdec = np.arctan(sin_d/cos_d)
    sransn = np.pi - np.arctan2(1/np.tan(obliq)*sin_d/cos_d,
                                -1*np.cos(slp)/cos_d)
    #Since GST is in degrees convert declination and right ascension
    if degrees:
        sdec = sdec * 180./np.pi
        sransn = sransn * 180./np.pi
        return gst,sdec,sransn
    else:
        gst = np.radians(gst)
        return gst,sdec,sransn

def solar_zenith_angle(dt,lats,lons,degrees=True):
    """
    This is from github.com/lkilcommons/geospacepy-lite, in
    astrodynamics.py. This was written by Liam Kilcommons

    Finds solar zenith angle using Russel solar position
    """
    lam = np.radians(lats)
    phi = np.radians(lons)
    gst,sdec,sra = solar_position_approx(dt)
    #Calculate hour angle
    sha = sra - (gst+phi)
    cossza = np.sin(lam)*np.sin(sdec) + np.cos(lam)*np.cos(sdec)*np.cos(sha)
    if degrees:
        return np.degrees(np.arccos(cossza))
    else:
        return np.arccos(cossza)

def manual_sphere(image_file):
    # caveat 1: flip the input image along its first axis
    img = plt.imread(image_file) # shape (N,M,3), flip along first dim
    outfile = image_file.replace('.jpg', '_flipped.jpg')
    # flip output along first dim to get right chirality of the mapping
    img = img[::-1,...]
    ny = img.shape[1]
    recentered_img = img.copy()
    recentered_img[:,:ny/2,:]=img[:,ny/2:,:]
    recentered_img[:,ny/2:,:]=img[:,:ny/2,:]
    plt.imsave(outfile, recentered_img)
    image_file = outfile  # work with the flipped file from now on

    # parameters for the sphere
    R = 1 # radius of the sphere
    Nrad = 180 # points along theta and phi
    phi = np.linspace(0, 2 * np.pi, Nrad)  # shape (Nrad,)
    theta = np.linspace(0, np.pi, Nrad)    # shape (Nrad,)
    phigrid,thetagrid = np.meshgrid(phi, theta) # shapes (Nrad, Nrad)

    # compute actual points on the sphere
    x = R * np.sin(thetagrid) * np.cos(phigrid)
    y = R * np.sin(thetagrid) * np.sin(phigrid)
    z = R * np.cos(thetagrid)

    # create meshed sphere
    mesh = mlab.mesh(x,y,z)
    mesh.actor.actor.mapper.scalar_visibility = False
    mesh.actor.enable_texture = True  # probably redundant assigning the texture later

    # load the (flipped) image for texturing
    img = tvtk.JPEGReader(file_name=image_file)
    texture = tvtk.Texture(input_connection=img.output_port, interpolate=0, repeat=0)
    mesh.actor.actor.texture = texture

    # tell mayavi that the mapping from points to pixels happens via a sphere
    mesh.actor.tcoord_generator_mode = 'sphere' # map is already given for a spherical mapping
    cylinder_mapper = mesh.actor.tcoord_generator
    # caveat 2: if prevent_seam is 1 (default), half the image is used to map half the sphere
    cylinder_mapper.prevent_seam = 0 # use 360 degrees, might cause seam but no fake data
    #cylinder_mapper.center = np.array([0,0,0])  # set non-trivial center for the mapping sphere if necessary

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
                vmin=None,vmax=None,cmap='viridis',reverse_cmap=False,
                amin=None,amax=None,alpha=None,
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
    
    if amin is not None and amax is not None:
        lut_set_alpha_range(mesh.parent.scalar_lut_manager.lut,
                                                amin=amin,amax=amax)

    mesh.parent.scalar_lut_manager.reverse_lut = reverse_cmap 
    if colorbar:
        legend(mesh.parent.scalar_lut_manager,label=label)
    if alpha is not None:
        mesh.parent.scalar_lut_manager.lut.alpha=alpha
    mlab.draw()
    return mesh

def set_camera_position(x,y,z):
    """Set the camera to particular location"""
    fig = mlab.gcf()
    scene = fig.scene
    #Set camera to equator
    camera = scene.camera
    camera.position=np.array([5.,0.,5.])
    #camera.view_angle=0.
    camera.compute_view_plane_normal()
    scene.render()

def lights_to_solar_position(dt):
    """Modify the Mayavi light locations to be a the sun's right ascension
    and declination (i.e. pointed at the subsolar point)
    """
    gst,sdec,sransn = solar_position_approx(dt)
    print('Solar Declination:',sdec)
    print('Solar Right Ascension:',sransn)
    #pdb.set_trace()
    
    #Vtk mode is one light
    fig = mlab.gcf()
    scene = fig.scene
    scene.light_manager.light_mode = "vtk"
    for light in scene.light_manager.lights:
        light.elevation = 0.
        light.azimuth = srans
        light.elevation = 90.-sdec
    
def darkness_geosurf(dt):
    """Draw a surface over the dark portion of the planet,
    i.e. regions with solar zenith angle > 90
    """
    print(dt)
    glats,glons = np.mgrid[-90:90:512j,0:360:1024j]
    gst,sdec,sransn = solar_position_approx(dt)
    
    sza = solar_zenith_angle(dt,glats,glons)
    dark = (sza-90.)/90.+1.
    #gist_yarg is a white to black colormap
    mesh = geosurf(glats,glons,dark,alt=50.,vmin=1.,vmax=1.01,alpha=.2,cmap='bone',reverse_cmap=True)
    lut = mesh.parent.scalar_lut_manager.lut
    lut.use_above_range_color=1
    lut.above_range_color=np.array([0.,0.,0.,.75])
    mlab.draw()

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
    #sphere_actor = map_image_to_unit_sphere_actor(imagefn)
    #fig.scene.add_actor(sphere_actor)
    manual_sphere(imagefn)

    #dt = datetime.datetime(2019,11,30,6)
    dt = datetime.datetime.utcnow()
    darkness_geosurf(dt)
    set_camera_position(5.,0.,5.)

    glats,glons = np.mgrid[-90:90:512j,0:360:1024j]
    aur_prob = np.genfromtxt(nowcasttxtfn).reshape((512,1024))
    geosurf(glats,glons,aur_prob,
                label='Probability of Aurora',
                cmap='viridis',vmin=1.,vmax=15.,amin=.5,amax=1.,
                colorbar=True)

    mlab.show()

if __name__ == "__main__":
    visual_test_geosurf()
