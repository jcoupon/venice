import pylab,numpy as N
import venicepy

ra = N.random.uniform(35,36,1e4)
dec = N.random.uniform(-5,-4,1e4)

M = venicepy.Mask("CFHTLS_W1_i_T0005_Jean.reg")

inside = M.check_point(ra,dec)


print len(ra),len(ra[inside])

pylab.plot(ra[inside],dec[inside],",")
pylab.show()
