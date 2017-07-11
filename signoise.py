from __future__ import print_function
import flipper
from flipper import liteMap
from astropy.io import fits
import orphics.tools.io as io
from scipy.ndimage.interpolation import zoom
import numpy as np
from orphics.tools.stats import bin2D
from orphics.tools.io import Plotter,dictFromSection,listFromConfig
from orphics.analysis import flatMaps as fmaps
import numpy as np
import sys, os
from skimage.transform import resize
import orphics.tools.stats as stats

# For Mat on BNL
# map_location = "/astro/astronfs01/workarea/msyriac/act/cmb_maps/c7v5_release/s2/ACTPol_148_D56_PA1_S2_1way_I_src_free.fits"
# cat_location = "/astro/astronfs01/workarea/msyriac/act/ACTPol.fits"
# out_dir = os.environ['WWW']

# For Teva
#map_location = '../ACTdata/ACTPol_148_D56_PA1_S2_1way_I_src_free.fits'
map_location = '../ACTdata/s15_pa2_s15_pa2_day_s15_pa3_150_fullex_I_map.fits'
#for ACT catalog
# cat_location = '../ACTdata/ACTPol.fits'
#for SDSS catalog
cat_location='../ACTdata/ACTPol_BOSS-N_BJF.fits'
# cat_location='../SDSSdata/asu.fit'
out_dir = "./"

hdulist=fits.open(cat_location)
catalog= hdulist[1].data

lmap=liteMap.liteMapFromFits(map_location)
#io.highResPlot2d(lmap.data,out_dir+"newbigmap.png")
RAs=[]
DECs=[]
widthStampArcminute = 60.
pixScale = 0.5
Nrands = 10000

lmap.info()

# select non-noisy part by trial and error
lmap = lmap.selectSubMap(245,125,-4.,21.)
#io.highResPlot2d(lmap.data,out_dir+"cutmap.png")
lmap.info()
#sys.exit()

ra_range = [lmap.x1-360.,lmap.x0] 
dec_range = [lmap.y0,lmap.y1]
randomstack=np.load('newrandomstamp.npy')
profiles=[]
#lmap.data = lmap.data - lmap.data.mean()


#defines ra and dec for point sources
RAps=[]
DECps=[]
file1=open('../ACTdata/boss_north_s15_f150_pa2_20170402 2.csv','rb')
psources=np.genfromtxt(file1)
for i in range(0,2130):
    n=i*12
    RAps.append(psources[2+n])
    DECps.append(psources[3+n])

def stack_on_map(lite_map,width_stamp_arcminute,pix_scale,ra_range,dec_range,catalog=None,n_random_points=None):
    width_stamp_degrees = width_stamp_arcminute /60.
    Np = np.int(width_stamp_arcminute/pix_scale+0.5)
    pad = np.int(Np/2+0.5)
    print ("Expected width in pixels = ", Np)
    lmap = lite_map
    stack=0
    N=0

    if catalog is not None:
        looprange = range(0,len(catalog))
        assert n_random_points is None
        random = False
    else:
        assert n_random_points is not None
        assert len(ra_range)==2
        assert len(dec_range)==2
        looprange = range(0,n_random_points)
        random = True
    print(looprange)
    for i in looprange:
        banana=True
        mass=catalog[i][10]
        if random:
            ra = np.random.uniform(*ra_range)
            dec =np.random.uniform(*dec_range)
        if random==False:
            ra=catalog[i][1] #1 for ACT catalog 2 for SDSS
            dec=catalog[i][2] #2 for ACT catalog 3 for SDSS
        for j in range(0,2130):
            distance=np.sqrt((ra-RAps[j])**2+(dec-DECps[j])**2)
            crit=0.25
            if distance<crit:
                banana=False
                print('too close')
        ix, iy = lmap.skyToPix(ra,dec)
        if ix>=pad and ix<lmap.Nx-pad and iy>=pad and iy<lmap.Ny-pad and banana==True and mass>8:
            print(i)
            #print(ra,dec)
            smap = lmap.selectSubMap(ra-width_stamp_degrees/2.,ra+width_stamp_degrees/2.,dec-width_stamp_degrees/2.,dec+width_stamp_degrees/2.)
            #print (smap.data.shape)
            #cutout = zoom(smap.data.copy(),zoom=(float(Np)/smap.data.shape[0],float(Np)/smap.data.shape[1])
            cutout = resize(smap.data.copy(),output_shape=(Np,Np))
            xMap,yMap,modRMap,xx,yy = fmaps.getRealAttributes(smap)
            dt = pix_scale
            arcmax = 20.
            thetaRange = np.arange(0.,arcmax,dt)
            breal = bin2D(modRMap*180.*60./np.pi,thetaRange)
            a=breal.bin(cutout)
            profiles.append(a)
            io.quickPlot2d(cutout,str(i)+"cutout.png")
            #print (cutout.shape)
            stack = stack + cutout
            N=N+1
        else:
            print ("skip")
    stack=stack/N-randomstack
    #print(stack.shape())
    #print(smap.data.shape)
    # print(stack)
    print(N)
    stats.getStats(profiles)
    #np.save('newrandomstamp',stack)
    # io.quickPlot2d(stack,out_dir+"newACTstack.png")
    # dt = pix_scale
    # arcmax = 20.
    # thetaRange = np.arange(0.,arcmax,dt)
    # breal = bin2D(modRMap*180.*60./np.pi,thetaRange)
    # cents,recons = breal.bin(stack)
    # pl = Plotter(labelX='Distance from Center (arcminutes)',labelY='Temperature Fluctuation ($\mu K$)', ftsize=10)
    # pl.add(cents,recons)
    # pl._ax.axhline(y=0.,ls="--",alpha=0.5)
    # pl.done(out_dir+"newACTprofiles.png")
    return(stack, cents, recons)

#stack, cents, recons = stack_on_map(lmap,widthStampArcminute,pixScale,ra_range=ra_range,dec_range=dec_range,n_random_points=Nrands)
print(type(lmap))
print(type(widthStampArcminute))
print(type(pixScale))
print(type(ra_range))
print(type(dec_range))
print(type(catalog))
print(len(catalog))
stack, cents, recons = stack_on_map(lmap,widthStampArcminute,pixScale,ra_range=ra_range,dec_range=dec_range,catalog=catalog)


