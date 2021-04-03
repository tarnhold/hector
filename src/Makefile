# Makefile for the hector-1.0 programs
#
# 2/10/2012 Machiel Bos
#===============================================================================

include ../Makefile.inc

#--- Directory where programs will be installed
BINDIR=$(PREFIX)/bin

#--- Include directories
INCLUDE_DIRS = -I./ -I/usr/include/atlas

#--- convenient group names
likelihood_objects = Likelihood.o LevinsonGap.o LevinsonNoGap.o \
	AmmarGrag.o FullCov.o
noisemodel_objects = NoiseModel.o White.o Powerlaw.o PowerlawApprox.o \
	ARFIMA.o GenGaussMarkov.o
programs = estimatetrend modelspectrum estimatespectrum removeoutliers \
	date2mjd mjd2date

all : $(programs)

clean:
	rm -f *.o
	rm -f $(programs)

install:
	mkdir -p $(BINDIR)
	install $(programs) $(BINDIR)

estimatetrend : EstimateTrend.cpp Control.o Minimizer.o \
	$(noisemodel_objects) $(likelihood_objects) Observations.o \
	DesignMatrix.o HyperGeo.o
	$(CPP) $(CPPFLAGS) EstimateTrend.cpp Control.o Minimizer.o \
	$(noisemodel_objects) $(likelihood_objects) Observations.o \
	DesignMatrix.o HyperGeo.o \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas -lfftw3 -lgsl \
	-o estimatetrend

removeoutliers : RemoveOutliers.cpp RemoveOutliers.h Observations.o \
        DesignMatrix.o Control.o
	$(CPP) $(CPPFLAGS) RemoveOutliers.cpp Observations.o DesignMatrix.o \
	Control.o \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas \
	-o removeoutliers

modelspectrum : ModelSpectrum.cpp ModelSpectrum.h $(noisemodel_objects) \
	Control.o Observations.o HyperGeo.o
	$(CPP) $(CPPFLAGS) ModelSpectrum.cpp $(noisemodel_objects) \
        Control.o Observations.o HyperGeo.o \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas -lfftw3 -lgsl \
        -o modelspectrum

estimatespectrum : EstimateSpectrum.cpp Observations.o Control.o
	$(CPP) $(CPPFLAGS) EstimateSpectrum.cpp Observations.o Control.o \
	$(INCLUDE_DIRS) -lfftw3 \
        -o estimatespectrum

date2mjd : date2MJD.cpp
	$(CPP) $(CPPFLAGS) date2MJD.cpp -o date2mjd

mjd2date : MJD2date.cpp
	$(CPP) $(CPPFLAGS) MJD2date.cpp -o mjd2date


#---------- Subroutines

MLE.o : MLE.cpp MLE.h
	$(CPP) $(CPPFLAGS) -c MLE.cpp \
        $(INCLUDE_DIRS)

Observations.o : Observations.cpp Observations.h Control.o
	$(CPP) $(CPPFLAGS) -c Observations.cpp 

Control.o : Control.cpp Control.h
	$(CPP) $(CPPFLAGS) -c Control.cpp

NoiseModel.o : NoiseModel.cpp NoiseModel.h Control.o 
	$(CPP) $(CPPFLAGS) -c NoiseModel.cpp \
        $(INCLUDE_DIRS)

Minimizer.o : Minimizer.cpp Minimizer.h Control.o NoiseModel.o \
	Likelihood.o
	$(CPP) $(CPPFLAGS) -c Minimizer.cpp \
        $(INCLUDE_DIRS)

Likelihood.o : Likelihood.cpp Likelihood.h Control.o DesignMatrix.o \
	Observations.o 
	$(CPP) $(CPPFLAGS) -c Likelihood.cpp \
        $(INCLUDE_DIRS)

LevinsonGap.o : LevinsonGap.cpp LevinsonGap.h Likelihood.o NoiseModel.o
	$(CPP) $(CPPFLAGS) -c LevinsonGap.cpp \
        $(INCLUDE_DIRS)

LevinsonNoGap.o : LevinsonNoGap.cpp LevinsonNoGap.h Likelihood.o NoiseModel.o
	$(CPP) $(CPPFLAGS) -c LevinsonNoGap.cpp \
        $(INCLUDE_DIRS)

DesignMatrix.o : DesignMatrix.cpp DesignMatrix.h Control.o Observations.o
	$(CPP) $(CPPFLAGS) -c DesignMatrix.cpp \
        $(INCLUDE_DIRS)

White.o : White.cpp White.h NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c White.cpp

AmmarGrag.o : AmmarGrag.cpp AmmarGrag.h Likelihood.o NoiseModel.o
	$(CPP) $(CPPFLAGS) -c AmmarGrag.cpp \
        $(INCLUDE_DIRS)

FullCov.o : FullCov.cpp FullCov.h Likelihood.o NoiseModel.o
	$(CPP) $(CPPFLAGS) -c FullCov.cpp \
        $(INCLUDE_DIRS)

Powerlaw.o : Powerlaw.cpp Powerlaw.h NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c Powerlaw.cpp

PowerlawApprox.o : PowerlawApprox.cpp PowerlawApprox.h \
	NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c PowerlawApprox.cpp \
        $(INCLUDE_DIRS)

ARFIMA.o : ARFIMA.cpp ARFIMA.h HyperGeo.o NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c ARFIMA.cpp \
        $(INCLUDE_DIRS)

GenGaussMarkov.o : GenGaussMarkov.cpp GenGaussMarkov.h HyperGeo.o
	$(CPP) $(CPPFLAGS) -c GenGaussMarkov.cpp \
	$(INCLUDE_DIRS)

HyperGeo.o : HyperGeo.cpp HyperGeo.h 
	$(CPP) $(CPPFLAGS) -c HyperGeo.cpp 
