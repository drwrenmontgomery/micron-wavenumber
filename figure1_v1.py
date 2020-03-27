from matplotlib import transforms as mtransforms
import astropy.units as un
import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost

class Micron2WNTransform(mtransforms.Transform):
    input_dims = 1
    output_dims =1
    is_separable = False
    has_inverse = True

    def __init__(self):
        mtransforms.Transform.__init__(self)

    def transform_non_affine(self, fr):
        return (fr*un.micron).to(un.k, equivalencies=un.spectral()).value

    def inverted(self):
        return WN2MicronTransform()

class WN2MicronTransform(Micron2WNTransform):
    input_dims = 1
    output_dims = 1
    is_separable = False
    has_inverse = True

    def __init__(self):
        mtransforms.Transform.__init__(self)

    def transform_non_affine(self, wl):
        return (wl*un.k).to(un.micron, equivalencies=un.spectral()).value

    def inverted(self):
        return Micron2WNTransform()



aux_trans = mtransforms.BlendedGenericTransform(Micron2WNTransform(), mtransforms.IdentityTransform())

fig = plt.figure(1)

for n in range(1,4,1):
	ax_mn = SubplotHost(fig, 1, 3, n)

	fig.add_subplot(ax_mn)
	ax_mn.set_xlabel('Wavelength (micron)')
#	x_micron=np.array([15,10,5,3])
#	ax_mn.set_xticks(x_micron)
#	ax_mn.set_xlim(5,15)
	ax_mn.set_ylim(0,6)
	ax_mn.set_xscale('log')

	test_spectrum=np.genfromtxt('sample_spectrum.dpt', delimiter=',')
	xvals = 10000/test_spectrum[:,0]
	data= test_spectrum[:,1]

	ax_wn = ax_mn.twin(aux_trans)
	ax_wn.set_viewlim_mode("transform")

	ax_wn.set_xlabel("Wavenumber, cm$^{-1}$")
	x_wn=np.array([1600, 1400, 1200, 1000, 800])
	ax_wn.set_xlim(1700,700)
	ax_wn.set_xticks(x_wn)

	ax_wn.axis["right"].toggle(ticklabels=False)
	ax_wn.set_xscale('linear')

	ax_mn.semilogx(xvals, data)
#	ax_wn.invert_xaxis()

plt.draw()
plt.show()

