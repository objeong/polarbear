#!opt/polarbear/bin/python
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import healpy as hp
import pickle
import sys
from matplotlib.pyplot import*

def clip(array, mask):
	array[np.logical_not(mask)] = np.nan
        minmap = np.nanpercentile(array[mask], 1)
        maxmap = np.nanpercentile(array[mask], 99)
        return np.minimum(np.maximum(array, minmap), maxmap)

#Load the files
print('Loading the files')
f1 = h5py.File(sys.argv[1], 'r')
mapinfo1 = pickle.loads(f1.attrs['mapinfo'])

f2 = h5py.File(sys.argv[2], 'r')
mapinfo2 = pickle.loads(f2.attrs['mapinfo'])

print('Instantiating arrays')
npix1 = hp.nside2npix(mapinfo1.nside)
map_I_1 = np.empty((npix1))
map_Q_1 = np.empty((npix1))
map_U_1 = np.empty((npix1))
map_w0_1 = np.zeros((npix1))
map_w4_1 = np.zeros((npix1))
map_I_1[:] = np.nan
map_Q_1[:] = np.nan
map_U_1[:] = np.nan

print('Instantiating arrays')
npix2 = hp.nside2npix(mapinfo2.nside)
map_I_2 = np.empty((npix2))
map_Q_2 = np.empty((npix2))
map_U_2 = np.empty((npix2))
map_w0_2 = np.zeros((npix2))
map_w4_2 = np.zeros((npix2))
map_I_2[:] = np.nan
map_Q_2[:] = np.nan
map_U_2[:] = np.nan

#file I/O
print('Creating Maps')
obspix1 = mapinfo1.obspix
I_1 = f1['I']
Q_1 = f1['Q']
U_1 = f1['U']
w0_1 = f1['w0']
w4_1 = f1['w4']

print len(I_1)
print len(Q_1)
print len(U_1)
print len(w0_1)
print len(w4_1)

obspix2 = mapinfo2.obspix
I_2 = f2['I']
Q_2 = f2['Q']
U_2 = f2['U']
w0_2 = f2['w0']
w4_2 = f2['w4']

print len(I_2)
print len(Q_2)
print len(U_2)
print len(w0_2)
print len(w4_2)

try:
	len_obs1 = obspix1.shape[0]
	print len_obs1
	for i in range(0, len_obs1):
		assert obspix1[i] == i
        print('Simple healpix mapping')
	map_I_1[:len_obs1] = I_1
        map_Q_1[:len_obs1] = Q_1
        map_U_1[:len_obs1] = U_1
        map_w0_1[:len_obs1] = w0_1
        map_w4_1[:len_obs1] = w4_1

except AssertionError:
	print('Complex healpix mapping')
	i_arr1 = ([i for i in range(0, w0_1.shape[0]) if w0_1[i] > 0])
	map_I_1[obspix[i_arr1[:]]] = I_1[i_arr1[:]]
	map_Q_1[obspix[i_arr1[:]]] = Q_1[i_arr1[:]]
	map_U_1[obspix[i_arr1[:]]] = U_1[i_arr1[:]]
	map_w0_1[obspix[i_arr1[:]]] = w0_1[i_arr1[:]]
	map_w4_1[obspix[i_arr1[:]]] = w4_1[i_arr1[:]]

try:
	len_obs2 = obspix2.shape[0]
	print len_obs2
	for i in range(0, len_obs2):
		assert obspix2[i] == i
        print('Simple healpix mapping')
	map_I_2[:len_obs2] = I_2
        map_Q_2[:len_obs2] = Q_2
        map_U_2[:len_obs2] = U_2
        map_w0_2[:len_obs2] = w0_2
        map_w4_2[:len_obs2] = w4_2

except AssertionError:
	print('Complex healpix mapping')
	i_arr2 = ([i for i in range(0, w0_2.shape[0]) if w0_2[i] > 0])
	map_I_2[obspix[i_arr2[:]]] = I_2[i_arr2[:]]
	map_Q_2[obspix[i_arr2[:]]] = Q_2[i_arr2[:]]
	map_U_2[obspix[i_arr2[:]]] = U_2[i_arr2[:]]
	map_w0_2[obspix[i_arr2[:]]] = w0_2[i_arr2[:]]
	map_w4_2[obspix[i_arr2[:]]] = w4_2[i_arr2[:]]

#Create the projection
print('Projections')

def vec3pix1(x, y, z):
	return hp.pixelfunc.vec2pix(mapinfo1.nside, x, y, z)

def vec3pix2(x, y, z):
    return hp.pixelfunc.vec2pix(mapinfo2.nside, x, y, z)

#Azimuth, $^circ$ range
center_az = 180.
az_range= 70. #half width of the projection
min_az = center_az - az_range
max_az = center_az + az_range
min_el = 20.
max_el = 80.

proj = hp.projector.CartesianProj(xsize=800, ysize=800, rot=[center_az,0, 0], lonra=[-az_range, az_range], latra=[min_el, max_el])
plot_I_1 = np.fliplr(proj.projmap(map_I_1, vec3pix1))*10**3
plot_Q_1 = np.fliplr(proj.projmap(map_Q_1, vec3pix1))*10**3
plot_U_1 = np.fliplr(proj.projmap(map_U_1, vec3pix1))*10**3
plot_w0_1 = np.fliplr(proj.projmap(map_w0_1, vec3pix1))*10**3
plot_w4_1 = np.fliplr(proj.projmap(map_w4_1, vec3pix1))*10**3
plot_mag_1 = np.sqrt(plot_Q_1**2 + plot_U_1**2)
plot_arg_1 = np.angle(plot_Q_1 + 1j* plot_U_1)

plot_I_2 = np.fliplr(proj.projmap(map_I_2, vec3pix2))*10**3
plot_Q_2 = np.fliplr(proj.projmap(map_Q_2, vec3pix2))*10**3
plot_U_2 = np.fliplr(proj.projmap(map_U_2, vec3pix2))*10**3
plot_w0_2 = np.fliplr(proj.projmap(map_w0_2, vec3pix2))*10**3
plot_w4_2 = np.fliplr(proj.projmap(map_w4_2, vec3pix2))*10**3
plot_mag_2 = np.sqrt(plot_Q_2**2 + plot_U_2**2)
plot_arg_2 = np.angle(plot_Q_2 + 1j* plot_U_2)



#make masks
plot_w0_mask_1 = plot_w0_1 > 0.	#0.001 * np.max(plot_w0)
plot_w4_mask_1 = plot_w4_1 > 0.	#001 * np.max(plot_w4)
plot_w0_mask_2 = plot_w0_2 > 0.	#0.001 * np.max(plot_w0)
plot_w4_mask_2 = plot_w4_2 > 0.	#001 * np.max(plot_w4)

figure()
imshow(clip((plot_I_1-plot_I_2) * 10 **3, plot_w0_mask), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template I difference, $\mu K_{CMB}$')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
#clim(vmin=-1900, vmax=2000)
savefig('dI_gt_%s.png'%sys.argv[3])

figure()
imshow(clip((plot_Q_1-plot_Q_2) * 10**3, plot_w4_mask), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template Q difference, $\mu K_{CMB}$')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
#clim(vmin=-120, vmax=300)
savefig('dQ_gt_%s.png'%sys.argv[3])

figure()
imshow(clip((plot_U_1-plot_U_2) * 10**3, plot_w4_mask), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template U difference, $\mu K_{CMB}$')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
#clim(vmin=-260, vmax=420)
savefig('dU_gt_%s.png'%sys.argv[3])

figure()
imshow(clip((plot_mag_1-plot_mag_2)*10**3, plot_w4_mask), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Az-El polarized power $\sqrt{Q^2 + U^2}$ difference, $\mu K_{CMB}$')
xlabel('Azimuth, $^\circ$')
ylabel('Elevation, $^\circ$')
gca().grid(True)
colorbar()
#clim(vmin=0, vmax=600)
savefig('dMag_gt_%s.png'%sys.argv[3])

figure()
imshow((plot_arg_1-plot_arg_2), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template Polarized Phase difference')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
#clim(vmin=-2.8, vmax=2.8)
savefig('dArg_gt_%s.png'%sys.argv[3])

figure()
imshow((plot_w0_1-plot_w0_2), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template 0f weight difference')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
#clim(vmin=0, vmax=3.7e14)
savefig('dw0_gt_%s.png'%sys.argv[3])

figure()
imshow((plot_w4_1-plot_w4_2), origin='lower', interpolation='nearest', extent=[min_az, max_az, min_el, max_el], aspect='auto')
title('Ground template 4f weight difference')
xlabel('Azimuth, $^circ$')
ylabel('Elevation, $^circ$')
colorbar()
#clim(vmin=0, vmax=6.7e15)
savefig('dw4_gt_%s.png'%sys.argv[3])

show()

