#!/usr/bin/env python2.7
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import healpy as hp
import pickle
import sys
from matplotlib.pyplot import *
import os

#basedir = "$SCRATCH/ground_map/coadd"

def clip(array, mask):
	array[np.logical_not(mask)] = np.nan
        minmap = np.nanpercentile(array[mask], 1)
        maxmap = np.nanpercentile(array[mask], 99)
        return np.minimum(np.maximum(array, minmap), maxmap)

#Load the file
print('Loading the file.')
#os.chdir(basedir)
f = h5py.File(sys.argv[1], 'r')
mapinfo = pickle.loads(f.attrs['mapinfo'])

print('Instantiating arrays')
npix = hp.nside2npix(mapinfo.nside)
map_I = np.empty((npix))
map_Q = np.empty((npix))
map_U = np.empty((npix))
map_w0 = np.zeros((npix))
map_w4 = np.zeros((npix))
map_I[:] = np.nan
map_Q[:] = np.nan
map_U[:] = np.nan

#file I/O
print('Creating Maps')
obspix = mapinfo.obspix
I = f['I']
Q = f['Q']
U = f['U']
w0 = f['w0']
w4 = f['w4']

print len(I)
print len(Q)
print len(U)
print len(w0)
print len(w4)


try:
	len_obs = obspix.shape[0]
	print len_obs
	for i in range(0, len_obs):
		assert obspix[i] == i
        print('Simple healpix mapping')
	map_I[:len_obs] = I
        map_Q[:len_obs] = Q
        map_U[:len_obs] = U
        map_w0[:len_obs] = w0
        map_w4[:len_obs] = w4

except AssertionError:
	print('Complex healpix mapping')
	i_arr = ([i for i in range(0, w0.shape[0]) if w0[i] > 0])
	map_I[obspix[i_arr[:]]] = I[i_arr[:]]
	map_Q[obspix[i_arr[:]]] = Q[i_arr[:]]
	map_U[obspix[i_arr[:]]] = U[i_arr[:]]
	map_w0[obspix[i_arr[:]]] = w0[i_arr[:]]
	map_w4[obspix[i_arr[:]]] = w4[i_arr[:]]

#Create the projection
print('Projections')

def vec3pix(x, y, z):
	return hp.pixelfunc.vec2pix(mapinfo.nside, x, y, z)

#Azimuth, $^circ$ range
center_az = 180.
az_range= 70. #half width of the projection
min_az = center_az - az_range
max_az = center_az + az_range
min_el = 20.
max_el = 80.

proj = hp.projector.CartesianProj(xsize=800, ysize=800, rot=[center_az,0, 0], lonra=[-az_range, az_range], latra=[min_el, max_el])
plot_I = np.fliplr(proj.projmap(map_I, vec3pix))*10**3
plot_Q = np.fliplr(proj.projmap(map_Q, vec3pix))*10**3
plot_U = np.fliplr(proj.projmap(map_U, vec3pix))*10**3
plot_w0 = np.fliplr(proj.projmap(map_w0, vec3pix))*10**3
plot_w4 = np.fliplr(proj.projmap(map_w4, vec3pix))*10**3
plot_mag = np.sqrt(plot_Q**2 + plot_U**2)
plot_arg = np.angle(plot_Q + 1j* plot_U)

#make masks
plot_w0_mask = plot_w0 > 0.	#0.001 * np.max(plot_w0)
plot_w4_mask = plot_w4 > 0.	#001 * np.max(plot_w4)

figure()
imshow(clip(plot_I * 10 **3, plot_w0_mask), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template I, $\mu K_{CMB}$')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
savefig('I_gt_%s.png'%sys.argv[2])

figure()
imshow(clip(plot_Q * 10**3, plot_w4_mask), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template Q, $\mu K_{CMB}$')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
savefig('Q_gt_%s.png'%sys.argv[2])

figure()
imshow(clip(plot_U * 10**3, plot_w4_mask), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template U, $\mu K_{CMB}$')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
savefig('U_gt_%s.png'%sys.argv[2])

figure()
imshow(clip(plot_mag*10**3, plot_w4_mask), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Az-El polarized power $\sqrt{Q^2 + U^2}$, $\mu K_{CMB}$')
xlabel('Azimuth, $^\circ$')
ylabel('Elevation, $^\circ$')
gca().grid(True)
colorbar()
savefig('Mag_gt_%s.png'%sys.argv[2])

figure()
imshow(plot_arg, origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template Polarized Phase')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
savefig('Arg_gt_%s.png'%sys.argv[2])

figure()
imshow(plot_w0, origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template 0f weight')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
savefig('w0_gt_%s.png'%sys.argv[2])

figure()
imshow(plot_w4, origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template 4f weight')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
savefig('w4_gt_%s.png'%sys.argv[2])

#show()
