if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB1)
else
  INT_LAPACK_LIB=
endif

SUBDIRS := ../../libana \
		   ../../../atdyn \
		   ../../converter/crd_convert \
		   ../../trj_analysis/trj_analysis \
		   ../../trj_analysis/diffusion_analysis \
		   ../../trj_analysis/drms_analysis \
		   ../../trj_analysis/hb_analysis \
		   ../../trj_analysis/msd_analysis \
		   ../../trj_analysis/rg_analysis \
		   ../../trj_analysis/rmsd_analysis \
		   ../../clustering/kmeans_clustering \
		   ../../free_energy/mbar_analysis \
		   ../../free_energy/wham_analysis \
		   ../../free_energy/pmf_analysis \
		   ../../mode_analysis/avecrd_analysis \
		   .

lib_LTLIBRARIES = libpython_interface.la

libpython_interface_la_SOURCES = \
        conv_f_c_util.fpp \
        error_mod.fpp \
        s_molecule_c_mod.fpp \
        selection_c_mod.fpp \
        ctrl_c_mod.fpp \
        define_molecule.fpp \
        s_trajectories_c_mod.fpp \
        crd_convert_c_mod.fpp \
        trj_c_mod.fpp \
        rg_c_mod.fpp \
        rmsd_c_mod.fpp \
        drms_c_mod.fpp \
        msd_c_mod.fpp \
        hb_c_mod.fpp \
        diffusion_c_mod.fpp \
        avecrd_c_mod.fpp \
        wham_c_mod.fpp \
        mbar_c_mod.fpp \
        pmf_c_mod.fpp \
        kmeans_c_mod.fpp \
        dynamic_string_mod.fpp \
        internal_file_type_mod.fpp \
        export_molecule_c_mod.fpp \
        atdyn_c_mod.fpp

EXTRA_DIST = \
	Makefile.depends

SUFFIXES = .fpp
DEPENDS = Makefile.depends

# Ensure Position-Independent Code is generated
AM_FCFLAGS = -fPIC

libpython_interface_la_LDFLAGS = -module -avoid-version -export-dynamic -shared -no-undefined
libpython_interface_la_LDFLAGS += $(MPI_FORTRAN_LDFLAGS)
libpython_interface_la_LDFLAGS += $(PYTHON_INTERFACE_LDFLAG) -Wl,-rpath,$(libdir)
libpython_interface_la_LIBADD = \
		../../../atdyn/libatdyn.a \
		../../converter/crd_convert/libcrd_convert.a \
		../../trj_analysis/trj_analysis/libtrj_analysis.a \
		../../trj_analysis/rg_analysis/librg_analysis.a \
		../../trj_analysis/rmsd_analysis/librmsd_analysis.a \
		../../trj_analysis/drms_analysis/libdrms_analysis.a \
		../../trj_analysis/msd_analysis/libmsd_analysis.a \
		../../trj_analysis/hb_analysis/libhb_analysis.a \
		../../trj_analysis/diffusion_analysis/libdiffusion_analysis.a \
		../../mode_analysis/avecrd_analysis/libavecrd_analysis.a \
		../../clustering/kmeans_clustering/libkmeans_clustering.a \
		../../free_energy/wham_analysis/libwham_analysis.a \
		../../free_energy/mbar_analysis/libmbar_analysis.a \
		../../free_energy/pmf_analysis/libpmf_analysis.a \
		../../libana/libana.a ../../../lib/lib.a ../../../ \
		$(INT_LAPACK_LIB) \
		$(LAPACK_LIBS) \
		-lgfortran $(LIBS) \
		$(MPI_FORTRAN_LIBS)

libpython_interface_la_FCFLAGS = $(FCFLAGS) $(PYTHON_CPPFLAGS)
# The error guard (error_mod / run_guarded) hands C function pointers to
# Fortran procedures. An internal (contains) procedure there would need a
# trampoline, i.e. an executable stack, which dlopen() rejects on glibc >= 2.41.
# Make gfortran fail the build instead of silently emitting one.
if GFORTRAN
libpython_interface_la_FCFLAGS += -Werror=trampolines
endif

# Add these lines
CCLD = $(FC)
LINK = $(LIBTOOL) --tag=FC --mode=link $(FC) $(AM_CFLAGS) $(CFLAGS) \
       $(AM_FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

# FC = mpif90
# CC = mpicc

.fpp.lo:
if USEKCOMP
	cp $< $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(LIBTOOL) --tag=FC --mode=compile $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana -I../../../atdyn \
		-I../../converter/crd_convert \
		-I ../../trj_analysis/trj_analysis \
		-I ../../trj_analysis/rg_analysis \
		-I ../../trj_analysis/rmsd_analysis \
		-I ../../trj_analysis/drms_analysis \
		-I ../../trj_analysis/msd_analysis \
		-I ../../trj_analysis/hb_analysis \
		-I ../../trj_analysis/diffusion_analysis \
		-I ../../mode_analysis/avecrd_analysis \
		-I ../../clustering/kmeans_clustering \
		-I ../../free_energy/wham_analysis \
		-I ../../free_energy/mbar_analysis \
		-I ../../free_energy/pmf_analysis \
		$(libpython_interface_la_FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $< $*.f90
	$(LIBTOOL) --tag=FC --mode=compile $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana -I../../../atdyn \
		-I../../converter/crd_convert \
		-I ../../trj_analysis/trj_analysis \
		-I ../../trj_analysis/rg_analysis \
		-I ../../trj_analysis/rmsd_analysis \
		-I ../../trj_analysis/drms_analysis \
		-I ../../trj_analysis/msd_analysis \
		-I ../../trj_analysis/hb_analysis \
		-I ../../trj_analysis/diffusion_analysis \
		-I ../../mode_analysis/avecrd_analysis \
		-I ../../clustering/kmeans_clustering \
		-I ../../free_energy/wham_analysis \
		-I ../../free_energy/mbar_analysis \
		-I ../../free_energy/pmf_analysis \
		$(libpython_interface_la_FCFLAGS) -c $*.f90
endif

depend: clean_depend
	python ../../../../fortdep.py *.fpp | sed -e 's/\.o/.lo/g' > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
