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

hdulist=fits.open('../ACTdata/ACTPol.fits')
catalog=hdulist[1].data
map=liteMap.liteMapFromFits('../ACTdata/ACTPol_148_D56_PA1_S2_1way_I.fits')
RAs=[]
DECs=[]
widthStampArcminute = 20.
widthStampDegrees = widthStampArcminute /60.

pixScale = 0.5
Np = np.int(widthStampArcminute/pixScale+0.5)
print ("Expected width in pixels = ", Np)


stack=0
N=0
for i in range(0,len(catalog)):
	ra=catalog[i][1]
	dec=catalog[i][2]
	ix, iy = map.skyToPix(ra,dec)
	# print ix, iy ; assert ix>=0 ; assert ix<map.Nx ; assert iy>=0 ; assert iy<map.Ny when this is uncommented we get an AssertionError which I guess means that some of the values are outside the map
	# if 53.983333<=ra and ra<=338.533333:
	if ix>=0 and ix<map.Nx and iy>=0 and iy<map.Ny:
		smap = map.selectSubMap(ra-widthStampDegrees/2.,ra+widthStampDegrees/2.,dec-widthStampDegrees/2.,dec+widthStampDegrees/2.)
                cutout = zoom(smap.data.copy(),zoom=(float(Np)/smap.data.shape[0],float(Np)/smap.data.shape[1]))
                print (cutout.shape)
		stack = stack + cutout
		xMap,yMap,modRMap,xx,yy = fmaps.getRealAttributes(smap)
		N=N+1.
stack=stack/N
#print(stack.shape())
print(smap.data.shape)
print(stack)
print(N)
io.quickPlot2d(stack,"stack.png")

dt = pixScale
arcmax = 20.
thetaRange = np.arange(0.,arcmax,dt)
breal = bin2D(modRMap*180.*60./np.pi,thetaRange)
cents,recons = breal.bin(stack)

pl = Plotter()
pl.add(cents,recons)
pl._ax.axhline(y=0.,ls="--",alpha=0.5)
pl.done("profiles.png")

