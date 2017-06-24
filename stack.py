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

# For Mat on BNL
map_location = "/astro/astronfs01/workarea/msyriac/act/cmb_maps/c7v5_release/s2/ACTPol_148_D56_PA1_S2_1way_I_src_free.fits"
cat_location = "/astro/astronfs01/workarea/msyriac/act/ACTPol.fits"
out_dir = os.environ['WWW']

# For Teva
#map_location = '../ACTdata/ACTPol_148_D56_PA1_S2_1way_I.fits'
#cat_location = '../ACTdata/ACTPol.fits'
#out_dir = "./"

hdulist=fits.open(cat_location)
catalog= hdulist[1].data
lmap=liteMap.liteMapFromFits(map_location)
#io.highResPlot2d(lmap.data,out_dir+"bigmap.png")
RAs=[]
DECs=[]
widthStampArcminute = 60.
pixScale = 0.5
Nrands = 1000

lmap.info()

# select non-noisy part by trial and error
lmap = lmap.selectSubMap(348,41,-8.3,5.3)
io.highResPlot2d(lmap.data,out_dir+"cutmap.png")
lmap.info()
#sys.exit()

ra_range = [lmap.x1-360.,lmap.x0] 
dec_range = [lmap.y0,lmap.y1]


#lmap.data = lmap.data - lmap.data.mean()



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
        
	for i in looprange:
		if random:
                        ra = np.random.uniform(*ra_range)
                        dec =np.random.uniform(*dec_range)
		else:
                        ra=catalog[i][1]
			dec=catalog[i][2]
		ix, iy = lmap.skyToPix(ra,dec)
		if ix>=pad and ix<lmap.Nx-pad and iy>=pad and iy<lmap.Ny-pad:
			print(i)
			#print(ra,dec)
			smap = lmap.selectSubMap(ra-width_stamp_degrees/2.,ra+width_stamp_degrees/2.,dec-width_stamp_degrees/2.,dec+width_stamp_degrees/2.)
                        #print (smap.data.shape)
                	#cutout = zoom(smap.data.copy(),zoom=(float(Np)/smap.data.shape[0],float(Np)/smap.data.shape[1]))
                	cutout = resize(smap.data.copy(),output_shape=(Np,Np))
                	#print (cutout.shape)
			stack = stack + cutout
			xMap,yMap,modRMap,xx,yy = fmaps.getRealAttributes(smap)
			N=N+1.
                else:
                        print ("skip")
	stack=stack/N
	#print(stack.shape())
	#print(smap.data.shape)
	print(stack)
	print(N)
	io.quickPlot2d(stack,out_dir+"stackrandom.png")

	dt = pix_scale
	arcmax = 20.
	thetaRange = np.arange(0.,arcmax,dt)
	breal = bin2D(modRMap*180.*60./np.pi,thetaRange)
	cents,recons = breal.bin(stack)
	pl = Plotter(labelX='Distance from Center (arcminutes)',labelY='Temperature Fluctuation ($\mu K$)', ftsize=10)
	pl.add(cents,recons)
	pl._ax.axhline(y=0.,ls="--",alpha=0.5)
	pl.done(out_dir+"randomprofiles.png")
        return stack, cents, recons

#stack, cents, recons = stack_on_map(lmap,widthStampArcminute,pixScale,ra_range=ra_range,dec_range=dec_range,n_random_points=Nrands)
stack, cents, recons = stack_on_map(lmap,widthStampArcminute,pixScale,ra_range=ra_range,dec_range=dec_range,catalog=catalog)


