AUTOMAKE_OPTIONS = dist-zip foreign subdir-objects no-dependencies -Wall -Wportability

__HDF_VER_FULL=$(shell echo @HDF5_VERSION@ | tr '.' '_')
__HDF_VER_SHORT=$(shell echo @HDF5_VERSION@ | tr '.' '_' | cut -d'_' -f1,2)

HDF5_URL = https://support.hdfgroup.org/releases/hdf5/v$(__HDF_VER_SHORT)/v$(__HDF_VER_FULL)/downloads/hdf5-@HDF5_VERSION@.tar.gz

__NPROC = $(shell nproc --ignore=1)

SUBDIRS = test

DOXYGEN_FORMATS       = @WITH_DOXYGEN_FORMATS@
DOXYGEN_CONFIGS       = doc/app.conf doc/lib.conf doc/vis.conf

bin_PROGRAMS = PiCoPiC

include_HEADERS = include/defines.hpp\
				include/SMB.hpp\
				include/collisions/collisionsTA77S.hpp\
				include/collisions/collisionsSK98.hpp\
				include/collisions/collisionsP12.hpp\
				include/outWriter.hpp\
				include/domain.hpp\
				include/outEngine.hpp\
				include/maxwellSolver/externalFields.hpp\
				include/maxwellSolver/externalFields.hpp~\
				include/maxwellSolver/maxwellSolverYee.hpp\
				include/current.hpp\
				include/defines.hpp.in\
				include/phys/rel.hpp\
				include/phys/plasma.hpp\
				include/math/vector3d.hpp\
				include/math/maxwellJuettner.hpp\
				include/math/rand.hpp\
				include/outController.hpp\
				include/specieP.hpp\
				include/current/currentVB.hpp\
				include/current/currentZigZag.hpp\
				include/msg.hpp\
				include/algo/weighter.hpp\
				include/algo/grid.hpp\
				include/algo/common.hpp\
				include/algo/grid3d.hpp\
				include/outEngine/outEngineHDF5.hpp\
				include/maxwellSolver.hpp\
				include/pusher.hpp\
				include/constant.hpp\
				include/pusher/pusherHC.hpp\
				include/pusher/pusherBoris.hpp\
				include/pusher/pusherVay.hpp\
				include/collisions.hpp\
				include/cfg.hpp\
				include/beamP.hpp\
				include/timeSim.hpp\
				include/geometry.hpp

PiCoPiC_SOURCES = src/beamP.cpp\
				src/specieP.cpp\
				src/outWriter.cpp\
				src/collisions/collisionsP12.cpp\
				src/collisions/collisionsSK98.cpp\
				src/collisions/collisionsTA77S.cpp\
				src/domain.cpp\
				src/maxwellSolver/maxwellSolverYee.cpp\
				src/maxwellSolver/externalFields.cpp\
				src/outController.cpp\
				src/SMB.cpp\
				src/PiCoPiC.cpp\
				src/phys/rel.cpp\
				src/phys/plasma.cpp\
				src/math/rand.cpp\
				src/math/maxwellJuettner.cpp\
				src/msg.cpp\
				src/current/currentVB.cpp\
				src/current/currentZigZag.cpp\
				src/algo/common.cpp\
				src/outEngine/outEngineHDF5.cpp\
				src/collisions.cpp\
				src/cfg.cpp\
				src/pusher/pusherHC.cpp\
				src/pusher/pusherVay.cpp\
				src/pusher/pusherBoris.cpp

hdf5:
	test -d $(top_srcdir)/tmp/hdf5-@HDF5_VERSION@ ||\
	( mkdir -p $(top_srcdir)/tmp;\
	cd $(top_srcdir)/tmp && wget -qO- $(HDF5_URL) | tar xzf -;\
	cd hdf5-@HDF5_VERSION@ &&\
	./autogen.sh &&\
	./configure --enable-build-mode=production --prefix=@prefix@ $(HDF5_ADDITIONAL_OPTIONS) &&\
	make -j$(__NPROC) @HDF5_ADDITIONAL_OPTIONS@ )

hdf5-install: hdf5
	$(MAKE) -C $(top_srcdir)/tmp/hdf5-@HDF5_VERSION@ install

compile_commands.json:
	$(MAKE) clean
	bear -- $(MAKE)

.install-hooks:
	cp -a tools/git_hooks/pre-commit .git/hooks

.PHONY: lib/python/picopic
lib/python/picopic:
	cd $(top_srcdir)/lib/python/picopic; \
	python setup.py install

deps: lib/python/picopic .install-hooks compile_commands.json

.PHONY: doc
doc:
	$(MAKE) -C doc
	cd $(top_srcdir)/doc/out && zip -r $(abs_top_builddir)/picopic_doc.zip html latex/refman.pdf
