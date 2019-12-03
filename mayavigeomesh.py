import pdb,time
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

    print('Greenwich Sidereal Time:',gst)    
    print('Solar Declination:',sdec)
    print('Solar Right Ascension:',sransn)

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

def get_blue_marble_texture_imgfn():
    """Flips the image inside-out so that it has the right
    'chirality' for use as a texture map in Mayavi, and also
    divides the image in half at the (approximate) prime meridian, exchanging
    the two halves of the image so that the prime meridan is at the edges of the
    image, which means Mayavi will wrap it around the sphere such that the 
    X axis of the figure is aligned with the prime meridian, and Z axis
    is in the direction of the North pole. This makes figure
    coordinates a rough appoximation of 
    Earth-Centered Earth-Fixed coordinates
    """
    imgfn = download_blue_marble() # shape (N,M,3), flip along first dim
    outfile = image_file.replace('.jpg', '_flipped.jpg')
    # flip output along first dim to get right chirality of the mapping
    flipped_img = img.copy()
    flipped_img = img[::-1,...]

    ny = img.shape[1]
    recentered_img = img.copy()
    recentered_img[:,:ny/2,:]=img[:,ny/2:,:]
    recentered_img[:,ny/2:,:]=img[:,:ny/2,:]
    
    plt.imsave(outfile, recentered_img)
    return outfile

def map_image_to_unit_sphere_actor(image_file):
    """This is a more low-level way of accomplishing a unit sphere
    textured with an image in Mayavi. Prefer image_textured_unit_sphere
    because it generates a mlab object which can be manipulated through
    the Mayavi GUI
    This code is largely borrowed from:
    https://stackoverflow.com/questions/53074908/map-an-image-onto-a-sphere-and-plot-3d-trajectories
    """

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

def image_textured_unit_sphere_mesh(image_file):
    """Generates a spherical Mayavi mesh and overlays an jpg image
    as a texture. The image is assumed to be a rectangular map in an
    projection which maps appropriately onto the meshace of a sphere
    without distortion (e.g. the NASA Blue Marble image). The image
    is also assumed to completely cover the sphere without repeats or overlap
    This code is largely borrowed from:
    https://stackoverflow.com/questions/53074908/map-an-image-onto-a-sphere-and-plot-3d-trajectories
    """

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
    return mesh

def blue_marble_unit_sphere_mesh():
    """Generates a unit-sphere mayavi mesh textured with 
    a NASA blue marble map
    """
    return image_textured_unit_mesh(get_blue_marble_texture_imgfn())

def lut_set_out_of_range_transparent(lut,above=False,below=False):
    """Set data below vmin to a transparent color"""
    if not above and not below:
        raise ValueError(('Out of range transparency called but above or below'
                            'kwargs not specified.'))
    if below:
        lut.use_below_range_color=1
        lut.below_range_color=np.array([0.,0.,0.,0.])
    if above:
        lut.use_above_range_color=1
        lut.above_range_color=np.array([0.,0.,0.,0.])

def lut_set_alpha_range(lut,amin=0.,amax=1.,diverging=False):
    """
    Modify a colormap (look up table, lut in Mayavi language)
    have a transparency gradient.

    'diverging' boolean kwarg will make the transparency minimum at the
    middle of the colormap and maximum at the ends

    See also:
    https://docs.enthought.com/mayavi/mayavi/auto/example_custom_colormap.html
    """
    lut_table = lut.table.to_array()
    n_colors = lut_table.shape[0]
    # docs say look up table (lut.table) needs integers
    iamin = np.floor(255*amin)
    iamax = np.floor(255*amax)
    if not diverging:
        alphas = np.floor(np.linspace(iamin,iamax,n_colors))
    else:
        half_range_alphas = np.floor(np.linspace(iamin,iamax,n_colors/2))
        alphas = np.concatenate([half_range_alphas.flip(),half_range_alphas]) 
        
    lut_table[:,-1] = alphas
    lut.table = lut_table

def legend(scalar_lut_manager,label=None):
    """
    Show a colorbar
    """
    scalar_lut_manager.show_scalar_bar = True
    scalar_lut_manager.show_legend = True
    scalar_lut_manager.data_name = label if label is not None else u''

def geo2cart(glats,glons,alt=110.):
    Re = 6371.2
    r = (Re+alt)/Re
    th = np.radians(90.-glats)
    ph = np.radians(glons)
    x = r*np.sin(th)*np.cos(ph)
    y = r*np.sin(th)*np.sin(ph)
    z = r*np.cos(th)
    return x,y,z

def update_geomesh(mesh,glats,glons,geodata,alt=300.):
    x,y,z = geo2cart(glats,glons,alt)
    mesh.mlab_source.trait_set(x=x,y=y,z=z,scalars=geodata)
    mlab.draw()

def geomesh(glats,glons,geodata,alt=300.,
                vmin=None,vmax=None,cmap='viridis',
                amin=None,amax=None,diverging_tranparency=False,alpha=None,
                colorbar=False,label=None,
                transparent_below_range=False):
    """Draw gridded data in geographic coordinates as a meshace (mesh) """

    if cmap[-2:] == '_r':
        reverse_cmap = True
        cmap = cmap.replace('_r','')
    else:
        reverse_cmap = False

    x,y,z = geo2cart(glats,glons,alt)
    mesh = mlab.mesh(x,y,z,
                        scalars=geodata,
                        vmin=vmin,
                        vmax=vmax,
                        colormap=cmap)
    
    mesh.parent.scalar_lut_manager.reverse_lut = reverse_cmap

    if transparent_below_range:
        lut_set_transparent_below_range(mesh.parent.scalar_lut_manager.lut)
    
    if amin is not None and amax is not None:
        lut_set_alpha_range(mesh.parent.scalar_lut_manager.lut,
                                                amin=amin,amax=amax,
                                                diverging=diverging_tranparency)
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
    camera.position=np.array([x,y,z])
    print('Camera at [{},{},{}]'.format(x,y,z))
    #camera.view_angle=0.
    camera.compute_view_plane_normal()
    scene.render()

def camera_to_sun_sync(dt,elevation=50.,distance=5.,azi_lag=90.):
    """Positions the camera theta_lag degrees behind
    a fixed point relative to the solar position"""
    gst,sdec,sransn = solar_position_approx(dt)
    azimuth = -1*np.degrees(gst)-azi_lag
    # th = np.radians(90.-elevation)
    # phi = -1*gst
    # xf = np.cos(th)*np.cos(phi)
    # yf = np.cos(th)*np.sin(phi)
    # zf = np.sin(th) 
    mlab.view(azimuth=azimuth,
                elevation=elevation,
                distance=distance)

def lights_to_solar_position(dt):
    """Modify the Mayavi light locations to be a the sun's right ascension
    and declination (i.e. pointed at the subsolar point)
    """
    gst,sdec,sransn = solar_position_approx(dt)
    #pdb.set_trace()
    
    #Vtk mode is one light
    fig = mlab.gcf()
    scene = fig.scene
    scene.light_manager.light_mode = "vtk"
    for light in scene.light_manager.lights:
        light.elevation = 0.
        light.azimuth = srans
        light.elevation = 90.-sdec
    
def calculate_dark(dt):
    glats,glons = np.mgrid[-90:90:512j,0:360:1024j]
    gst,sdec,sransn = solar_position_approx(dt)
    sza = solar_zenith_angle(dt,glats,glons)
    dark = (sza-90.)/90.+1.
    return glats,glons,dark

def darkness_geomesh(dt):
    """Draw a meshace over the dark portion of the planet,
    i.e. regions with solar zenith angle > 90
    """
    glats,glons,dark = calculate_dark(dt)
    #gist_yarg is a white to black colormap
    mesh = geomesh(glats,glons,dark,alt=50.,
                    vmin=1.,vmax=1.01,alpha=.2,transparent_below_range=True,
                    cmap='bone_r')
    lut = mesh.parent.scalar_lut_manager.lut
    lut.use_above_range_color=1
    lut.above_range_color=np.array([0.,0.,0.,.75])
    mlab.draw()
    return mesh

def download_auroral_nowcast():
    nowcastfn = os.path.join(tempfile.gettempdir(),'auroral_nowcast.txt')
    url = 'https://services.swpc.noaa.gov/text/aurora-nowcast-map.txt'
    r = requests.get(url,allow_redirects=True)
    open(nowcastfn,'w').write(r.content)
    return nowcastfn

def nowcast_auroral_probability():
    nowcasttxtfn = download_auroral_nowcast()
    glats,glons = np.mgrid[-90:90:512j,0:360:1024j]
    aur_prob = np.genfromtxt(nowcasttxtfn).reshape((512,1024))
    return glats,glons,aur_prob

def visual_test_geomesh():
    """Make a plot of the current NOAA auroral probability on the NASA
    blue marble image
    """
    
    # create a figure window (and scene)
    fig = mlab.figure(size=(600, 600))
    #fig.scene.scene.background_color=[0,0,0]
    #sphere_actor = map_image_to_unit_sphere_actor(imagefn)
    #fig.scene.add_actor(sphere_actor)
 
    earth_mesh = blue_marble_unit_sphere_mesh()

    dt = datetime.datetime.utcnow()
    dark_mesh = darkness_geomesh(dt)
    camera_to_sun_sync(dt,distance=3.)

    glats,glons,aur_prob = nowcast_auroral_probability()
    aur_mesh = geomesh(glats,glons,aur_prob,
                label='Probability of Aurora',
                cmap='viridis',vmin=1.,vmax=15.,amin=.5,amax=1.,
                colorbar=True,transparent_below_range=True)

    ttext = mlab.text(.05,.05,dt.strftime('%m/%d/%y %H:%M'),line_width=.75)
    ttext = mlab.text(.05,.02,'NOAA SWPC Aurora Nowcast',line_width=.75)

    mlab.savefig('geomesh_test.png')

    mlab.show()

def visual_test_update_geomesh():
    """Make a plot of the terminator and sync the camera to
    the solar position to mimic earth rotation
    """
    fig = mlab.figure(size=(600, 600))
    
    earth_mesh = blue_marble_unit_sphere_mesh()

    dt = datetime.datetime.utcnow()
    dark_mesh = darkness_geomesh(dt)
    set_camera_position(5.,0.,5.)

    for i in range(24):
        dt = datetime.datetime(2019,11,30,i)
        glats,glons,dark = calculate_dark(dt)
        camera_to_sun_sync(dt)
        update_geomesh(dark_mesh,glats,glons,dark,alt=50.)
        mlab.savefig('earth_{}.png'.format(i))

    mlab.show()

if __name__ == "__main__":
    visual_test_geomesh()
    #visual_test_update_geomesh()