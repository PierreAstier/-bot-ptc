

IMAGE_REPO=/sps/lsst/SLAC/LCA-10134_Cryostat-0001/6819D/BOT_acq_recovery/v0/47896/
#	export MPLBACKEND=Agg; extremely slow, replaced by pl.ioff()

# link to the source directory
%/source :
	test -d $* || mkdir $*
	ln -s $(IMAGE_REPO) $@

% : 
	+ dir=$$(dirname $@); make $$dir/source; \
	export BFPARAMS=$$PWD/envparams.yaml; \
	export MPLBACKEND=Agg; \
	cd $$dir; make $(MAKEFAGS) -f ../Makefile.chip $$(basename $@)

tar :
	tar cvf out.tar R*/*.png R*/*.pkl 

allfits :
	group_ptcfits.py */allfits.pkl

settings :
	nohup ./collect_settings.py R??_S??


#clean: FORCE
#	echo clean
#FORCE:



.SECONDARY : 



