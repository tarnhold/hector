# Makefile for the hector-1.3 programs
#
#  2/10/2012 Machiel Bos
# 13/ 1/2013 Machiel Bos
# 11/ 6/2014 Machiel Bos
#===============================================================================

include ../Makefile.inc

VERSION=1.3

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
	simulatenoise date2mjd mjd2date

all : $(programs)

clean:
	rm -f *.o $(programs)

install:
	mkdir -p $(BINDIR)
	install $(programs) $(BINDIR)

uninstall :
	@rm -f $(addprefix $(BINDIR)/,$(programs))


estimatetrend : EstimateTrend.cpp Control.o Minimizer.o \
	$(noisemodel_objects) $(likelihood_objects) Observations.o \
	DesignMatrix.o HyperGeo.o Calendar.o Statistics.o
	$(CPP) $(CPPFLAGS) EstimateTrend.cpp Control.o Minimizer.o \
	$(noisemodel_objects) $(likelihood_objects) Observations.o \
	DesignMatrix.o HyperGeo.o Calendar.o Statistics.o \
	 -DVERSION=$(VERSION) \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas -lfftw3 -lgsl \
	-o estimatetrend

removeoutliers : RemoveOutliers.cpp RemoveOutliers.h \
	Observations.o DesignMatrix.o Control.o Calendar.o
	$(CPP) $(CPPFLAGS) RemoveOutliers.cpp Observations.o DesignMatrix.o \
	Control.o Calendar.o -DVERSION=$(VERSION) \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas \
	-o removeoutliers

modelspectrum : ModelSpectrum.cpp ModelSpectrum.h \
	$(noisemodel_objects) Control.o Observations.o HyperGeo.o Calendar.o
	$(CPP) $(CPPFLAGS) ModelSpectrum.cpp $(noisemodel_objects) \
        Control.o Observations.o HyperGeo.o Calendar.o -DVERSION=$(VERSION) \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas -lfftw3 -lgsl \
        -o modelspectrum

estimatespectrum : EstimateSpectrum.cpp Observations.o Control.o \
	Calendar.o
	$(CPP) $(CPPFLAGS) EstimateSpectrum.cpp Observations.o Control.o \
	Calendar.o -DVERSION=$(VERSION) \
	$(INCLUDE_DIRS) -lfftw3 \
        -o estimatespectrum

simulatenoise : SimulateNoise.cpp SimulateNoise.h $(noisemodel_objects) \
	Control.o HyperGeo.o
	$(CPP) $(CPPFLAGS) SimulateNoise.cpp $(noisemodel_objects) \
        Control.o HyperGeo.o -DVERSION=$(VERSION) \
        $(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas -lfftw3 -lgsl \
        -o simulatenoise

date2mjd : date2MJD.cpp Calendar.o
	$(CPP) $(CPPFLAGS) date2MJD.cpp Calendar.o -o date2mjd

mjd2date : MJD2date.cpp Calendar.o
	$(CPP) $(CPPFLAGS) MJD2date.cpp Calendar.o -o mjd2date


test_HyperGeo : test_HyperGeo.cpp HyperGeo.o
	$(CPP) $(CPPFLAGS) test_HyperGeo.cpp HyperGeo.o \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas -lgsl -lgslcblas \
	-o test_HyperGeo

test_GaussRelation : test_GaussRelation.cpp HyperGeo.o
	$(CPP) $(CPPFLAGS) test_GaussRelation.cpp HyperGeo.o \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas -lgsl -lgslcblas \
	-o test_GaussRelation

test_ARFIMA : test_ARFIMA.cpp ARFIMA.o NoiseModelBaseClass.h HyperGeo.o \
	Control.o
	$(CPP) $(CPPFLAGS) test_ARFIMA.cpp ARFIMA.o HyperGeo.o Control.o \
	$(INCLUDE_DIRS) -latlas -lcblas -llapack_atlas -lgsl -lgslcblas \
	-o test_ARFIMA

#---------- Subroutines

MLE.o : MLE.cpp MLE.h
	$(CPP) $(CPPFLAGS) -c MLE.cpp \
        $(INCLUDE_DIRS)

Observations.o : Observations.cpp Observations.h Control.o Calendar.o
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

DesignMatrix.o : DesignMatrix.cpp DesignMatrix.h Control.o Observations.o \
	Calendar.o
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

Calendar.o : Calendar.cpp Calendar.h
	$(CPP) $(CPPFLAGS) -c Calendar.cpp

Statistics.o : Statistics.cpp Statistics.h
	$(CPP) $(CPPFLAGS) -c Statistics.cpp
