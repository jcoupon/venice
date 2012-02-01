import venicepy
import pylab
import numpy as N
import time

def grid(bounds=(35,-5,36,-4),step=.001):
    """ construct a uniform grid of test points. """
    ra = []
    dec = []

    for y in N.arange(bounds[1],bounds[3],step):
        x = N.arange(bounds[0],bounds[2],step)
        ra.append(x)
        dec.append(y*N.ones(len(x)))

    return N.concatenate(ra),N.concatenate(dec)


# Load the DS9 region file
M = venicepy.Mask("CFHTLS_W1_i_T0005_Jean.reg")

# generate some test points
ra,dec = grid()

print "testing points:",len(ra),'size (MB)',ra.nbytes*2/1024.**2
# find which points fall inside a polygon
t0 = time.time()
inside = M.check_point(ra,dec)
t = time.time()-t0
print "venice time:",t,"how many per sec:",len(ra)/t

ii = inside>0

# plot the results
pylab.plot(ra[ii],dec[ii],"b,")


""" test random catalogue """
pylab.figure()
bounds = N.array([35.,-5.,36.,-4.], dtype='d')
ra,dec = M.random_cat(bounds, 1e4, "binside")

pylab.plot(ra,dec,"b.")
print len(ra)

ra,dec = M.random_cat(bounds, 1e4, "inside")

pylab.plot(ra,dec,"y.")

pylab.show()
