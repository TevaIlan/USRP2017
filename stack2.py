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

ps_and_lowmass_avoidance=True #set to false to not avoid point sources and include all masses
wholes=True #set to True to make while hole style map
stat=True #set to True to do statistical analysis

# For Mat on BNL
# map_location = "/astro/astronfs01/workarea/msyriac/act/cmb_maps/c7v5_release/s2/ACTPol_148_D56_PA1_S2_1way_I_src_free.fits"
# cat_location = "/astro/astronfs01/workarea/msyriac/act/ACTPol.fits"
# out_dir = os.environ['WWW']

# For Teva
#map_location = '../ACTdata/ACTPol_148_D56_PA1_S2_1way_I_src_free.fits'
map_location = '../ACTdata/s15_pa2_s15_pa2_day_s15_pa3_150_fullex_I_map.fits'

#cat_location='../ACTdata/ACTPol_BOSS-N_BJF.fits'
#cat_location='../redmapper_dr8_public_v6.3_catalog.fits'
cat_location='../SDSSdata/asu.fit'
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
#lmap = lmap.selectSubMap(348,41,-8.3,5.3) #for D56 map
lmap = lmap.selectSubMap(245,125,-4.,21.) #for fullex map
#io.highResPlot2d(lmap.data,out_dir+"cutmap.png")
lmap.info()
#sys.exit()

ra_range = [lmap.x1-360.,lmap.x0] 
dec_range = [lmap.y0,lmap.y1]
# randomstack=np.load('D56randomstamp.npy')
randomstack=np.load('fullexrandomstamp.npy')
#lmap.data = lmap.data - lmap.data.mean()

#make empty arrays to use later if we want to make white hole style map
ixs=[]
iys=[]

#make empty list of cutouts to later use for statistical analysis
cutouts=[]


#eliminate undesireable objects (e.g. low mass or close to point source)
def eliminate_objects(minmass,masscolumn,catalog):
    goodobjects=[]
    #defines ra and dec for point sources
    RAps=[]
    DECps=[]
    file1=open('../ACTdata/boss_north_s15_f150_pa2_20170402 2.csv','rb')
    psources=np.genfromtxt(file1)
    lprange=range(0,len(catalog))
    for i in range(0,2130):
        n=i*12
        RAps.append(psources[2+n])
        DECps.append(psources[3+n])
    for i in lprange:
        ra=catalog[i][2] #1 for ACT catalog 2 for SDSS and redmapper
        dec=catalog[i][3] #2 for ACT catalog 3 for SDSS and redmapper
        proximity=True
        mass=catalog[i][masscolumn] #10 for SDSS catalog anyway
        for j in range(0,2130):
            distance=np.sqrt((ra-RAps[j])**2+(dec-DECps[j])**2)
            crit=0.25
            if distance<crit:
                proximity=False
        if proximity==True and mass>minmass:
            goodobjects.append(i)
    return(goodobjects)

if ps_and_lowmass_avoidance==True:
    goodobjects=eliminate_objects(minmass=8,masscolumn=10,catalog=catalog)
else:
    goodobjects=range(0,len(catalog))
                
def stack_on_map(lite_map,width_stamp_arcminute,pix_scale,ra_range,dec_range,catalog=None,n_random_points=None): 
    from skimage.transform import resize
    import orphics.tools.stats as stats

    width_stamp_degrees = width_stamp_arcminute /60.
    Np = np.int(width_stamp_arcminute/pix_scale+0.5)
    pad = np.int(Np/2+0.5)
    print ("Expected width in pixels = ", Np)

    lmap = lite_map
    stack=0
    N=0

    if catalog is not None:
        looprange = goodobjects 
        print(looprange)
        assert n_random_points is None
        random = False
    else:
        assert n_random_points is not None
        assert len(ra_range)==2
        assert len(dec_range)==2
        looprange = range(0,n_random_points)
        random = True
        
    for i in looprange:
        if random:
            ra = np.random.uniform(*ra_range)
            dec =np.random.uniform(*dec_range)
        if random==False:
            ra=catalog[i][1] #1 for ACT catalog 2 for SDSS and redmapper
            dec=catalog[i][2] #2 for ACT catalog 3 for SDSS and redmapper
        ix, iy = lmap.skyToPix(ra,dec)
        if ix>=pad and ix<lmap.Nx-pad and iy>=pad and iy<lmap.Ny-pad:
            print(i)
            smap = lmap.selectSubMap(ra-width_stamp_degrees/2.,ra+width_stamp_degrees/2.,dec-width_stamp_degrees/2.,dec+width_stamp_degrees/2.)
            #cutout = zoom(smap.data.copy(),zoom=(float(Np)/smap.data.shape[0],float(Np)/smap.data.shape[1]))
            cutout = resize(smap.data.copy(),output_shape=(Np,Np))
            cutouts.append(cutout-randomstack)
            stack = stack + cutout
            xMap,yMap,modRMap,xx,yy = fmaps.getRealAttributes(smap)
            N=N+1.
            ixs.append(ix)
            iys.append(iy)
        else:
            print ("skip")
    stack=stack/N-randomstack
    print(N)
    if catalog is not None:
        io.quickPlot2d(stack,out_dir+"stack.png")
    else:
        np.save('randomstamp',stack)
    
    dt = pix_scale
    arcmax = 20.
    thetaRange = np.arange(0.,arcmax,dt)
    breal = stats.bin2D(modRMap*180.*60./np.pi,thetaRange)
    cents,recons = breal.bin(stack)
    pl = Plotter(labelX='Distance from Center (arcminutes)',labelY='Temperature Fluctuation ($\mu K$)', ftsize=10)
    pl.add(cents,recons)
    pl._ax.axhline(y=0.,ls="--",alpha=0.5)
    pl.done(out_dir+"profiles.png")
    return stack, cents, recons

def whiteholes(xs,ys,map):
        mask = maskLiteMap(lmap,iys,ixs,holeArc=5,holeFrac=0.6)
        io.highResPlot2d(mask.data,"mask.png")

def stat_analysis(cutouts,binsize,arcmax,cents):
    profiles=[]
    for cutout in cutouts:
        thetaRange = np.arange(0.,arcmax,binsize)
        breali = bin2D(modRMap*180.*60./np.pi,thetaRange)
        a=breali.bin(cutout)[1]
        profiles.append(a)
    statistics=stats.getStats(profiles)
    mean=statistics['mean'] 
    error=statistics['errmean']
    covmat=statistics['cov']
    corrcoef=stats.cov2corr(covmat)
    io.quickPlot2d(corrcoef,'corrcoef.png')
    pl = Plotter(labelX='Distance from Center (arcminutes)',labelY='Temperature Fluctuation ($\mu K$)', ftsize=10)
    pl.add(cents,mean)
    pl.addErr(cents,mean,yerr=error)
    pl._ax.axhline(y=0.,ls="--",alpha=0.5)
    pl.done(out_dir+"error.png")


#stack, cents, recons = stack_on_map(lmap,widthStampArcminute,pixScale,ra_range=ra_range,dec_range=dec_range,n_random_points=Nrands)
stack, cents, recons = stack_on_map(lmap,widthStampArcminute,pixScale,ra_range=ra_range,dec_range=dec_range,catalog=catalog)

if wholes==True:
    whiteholes(xs=ixs,ys=iys,map=lmap)

if stat==True:
    stat_analysis(cutouts=cutouts,binsize=pixScale,arcmax=20.,cents=cents)