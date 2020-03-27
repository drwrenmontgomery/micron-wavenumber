from matplotlib import transforms as mtransforms
import astropy.units as un
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import fnmatch
import os
import re
import sys

def nrmlze(data, scalefactor):
    return scalefactor+(data[:,1]-min(data[:,1]))/(max(data[:,1])-min(data[:,1]))

class WN2MicronTransform(mtransforms.Transform):
    input_dims = 1
    output_dims =1
    is_separable = False
    has_inverse = True

    def __init__(self):
        mtransforms.Transform.__init__(self)

    def transform_non_affine(self, fr):
        return (fr*un.k).to(un.micron, equivalencies=un.spectral()).value

    def inverted(self):
        return Micron2WNTransform()

class Micron2WNTransform(WN2MicronTransform):
    input_dims = 1
    output_dims = 1
    is_separable = False
    has_inverse = True

    def __init__(self):
        mtransforms.Transform.__init__(self)

    def transform_non_affine(self, wl):
        return (wl*un.micron).to(un.k, equivalencies=un.spectral()).value

    def inverted(self):
        return WN2MicronTransform()



aux_trans = mtransforms.BlendedGenericTransform(WN2MicronTransform(), mtransforms.IdentityTransform())

fig = plt.figure(1)
fig.set_size_inches(7,4)

offset=0

ax_wn = SubplotHost(fig, 1, 1, 1)

fig.add_subplot(ax_wn)

ax_mn = ax_wn.twin(aux_trans)
ax_mn.set_viewlim_mode("transform")
#	x_micron=np.array([12, 10, 8, 6])
#	ax_mn.set_xticks(x_micron)
ax_mn.xaxis.tick_bottom()
ax_mn.tick_params(axis='x', direction = 'in', labelsize=11)
ax_mn.tick_params(axis='y', right=False, labelright=False)

#	ax_wn.set_xlim(675, 1700)
ax_wn.invert_xaxis()
#	x_wn=np.array([800, 1000, 1200, 1400, 1600])
#	ax_wn.set_xticks(x_wn)
ax_wn.xaxis.tick_top()
ax_wn.tick_params(axis='x', direction = 'in', labelsize=11)
ax_wn.tick_params(axis='y', left=False, labelleft=False)
ax_wn.xaxis.set_label_position('top')
ax_mn.xaxis.tick_bottom()
ax_mn.xaxis.set_label_position('bottom')

print('load data')
datat=np.loadtxt('dpt_files/Chrysene/p1_08_chry_v2.dpt', delimiter = ',')
data2=nrmlze(datat, 0)

ax_wn.plot(datat[:,0], data2 + offset, 'k-', lw=0.5, linestyle='-')

offset=offset+1.15

ax_wn.set_ylim([-.5, offset+0.5])
ax_wn.set_ylabel('intensity, a.u.', fontsize=11)
ax_mn.set_xlabel('wavelength, $\mu$m', fontsize=11)
ax_wn.set_xlabel('wavenumber, cm$^{-1}$', fontsize=11)

#fig.suptitle('wavenumber, cm$^{-1}$', fontsize=11)

#plt.subplots_adjust(top=0.88,
#bottom=0.11,
#left=0.035,
#right=0.965,
#hspace=0.2,
#wspace=0.12)

plt.draw()
plt.show()

