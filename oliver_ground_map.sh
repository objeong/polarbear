##########################################################################################
#
# Script given to Oliver to start investigating ground pickup systematics
# Author: Neil Goeckner-Wald
# ngoecknerwald@berkeley.edu
#
##########################################################################################

basedir="$SCRATCH/ground_map"

hdf5loc="/data/pb1/neil/hdf5_ces_dmds/"
mpiargs="--parmode mpi_edison -n 1 --mppwidth 32 --nodes 1 --ppn 32 \
	--queue regular --walltime 00:59:00"

#Mapmaking files
gain="/scratch/ngoecknerwald/largepatch2/gain"
psdloc="/scratch/ngoecknerwald/largepatch2/low_freq_epsilon"
filterfolder="/scratch/ngoecknerwald/largepatch2/filterdb_zeta"
nullfolder="/scratch/ngoecknerwald/largepatch2/null_test_flags_zeta"
dataselection="/scratch/ngoecknerwald/largepatch2/full_flag_zeta_with_cescut"

#Free parameters
firstgreg="20140725_000000"
lastgreg="20161008_000000" 

#NOTE for Oliver and Kolen,
# different from the fiducial pipeline in that we have added the --groundmap argument to keep only the ground synchronous signal
# and that we have changed the pointing arguments to give us a healpix map in ground coordinates
filterargs="--detrend 2 --removeleakage_pca --leak_fmax 0.4\
	--nofft --poly 1 --polpoly 1 --commonmode --commonmode_deweight"
projectionargs='--ground_coord --width 180 180 90 90 --healpix_offset Ground --healpix 256'
mapsplits="--mapsplits FIRST_SECOND LR_SUBSCAN GAIN_BY_CES QU_PIXEL LHS_RHS \
			MOON_HORIZON PWV SUN_HORIZON SUN_DIST MOON_DIST LEAK_BY_BOLO AMP_2F_BY_BOLO AMP_4F_BY_BOLO \
			STACK_Q_FKNEE STACK_U_FKNEE"

#Real data null map call function
function run_null_mapmaker() {

	python $ABPATH/whwp/process_field_whwp.py mapbolo_whwp -o $basedir -f $hdf5loc \
	 --firstgreg $firstgreg --lastgreg $lastgreg --psd_folder $psdloc --gain_folder $gain \
	 $projectionargs --singleblock $filterargs $mpiargs --alt_flag_folder $dataselection \
	 --filterfolder $filterfolder --null_folder $nullfolder $mapsplits
}


#run_null_mapmaker

#first stage coadding
#real data
#python $ABPATH/whwp/process_field_whwp.py coadd_hdf5_whwp -f $basedir \
#--name realmap --mppwidth 22 --nodes 1 --ppn 22 --walltime 00:20:00 --queue regular \
#--parmode mpi_edison -n 1 --save_separate $mapsplits


#final coadd for the null test splits, running on real data
python $ABPATH/whwp/coadd_hdf5_whwp.py --field $basedir --name realmap \
--coaddfinally --save_separate --load_separate $mapsplits

