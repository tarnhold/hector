# Makefile for the hector-1.5.1 programs
#
#  2/10/2012 Machiel Bos
# 13/ 1/2013 Machiel Bos
# 11/ 6/2014 Machiel Bos
# 16/ 3/2015 Machiel Bos
# 26/10/2015 Machiel Bos
#===============================================================================

include ../Makefile.inc

VERSION=1.5

#--- Directory where programs will be installed
BIN_DIR=$(PREFIX)/bin

#--- Include directories
ifeq "$(HOST)" "Machiels-MacBook"
  INC_DIR = -I/opt/local/include 
  LIB_DIR = -L/opt/local/lib 
  ATLAS   = -latlas -lcblas -llapack -lblas
else ifeq "$(HOST)" "segal"
  INC_DIR = -I./ -I/usr/include/atlas
  LIB_DIR = 
  ATLAS   = -llapack_atlas -latlas -lcblas
else 
  $(info ./src/Makefile: Unknown machine, I am using standard values)
  INC_DIR = -I./ -I/usr/include/atlas
  LIB_DIR = 
  ATLAS   = -llapack -lcblas -lblas -latlas 
endif

#--- convenient group names
likelihood_objects = Likelihood.o AmmarGrag.o FullCov.o 
STAnoisemodel_objects = White.o Powerlaw.o PowerlawApprox.o \
	ARFIMA.o GenGaussMarkov.o 
noisemodel_objects = NoiseModel.o $(STAnoisemodel_objects) 
programs = estimatetrend modelspectrum estimatespectrum simulatenoise \
	removeoutliers date2mjd mjd2date

all : $(programs)

clean:
	rm -f *.o $(programs)
	rm -rf *.dSYM

install:
	mkdir -p $(BIN_DIR)
	install $(programs) $(BIN_DIR)

uninstall :
	@rm -f $(addprefix $(BIN_DIR)/,$(programs))


#---- programs 

estimatetrend : EstimateTrend.cpp Control.o Minimizer.o \
	$(noisemodel_objects) $(likelihood_objects) Observations.o \
	DesignMatrix.o HyperGeo.o Calendar.o 
	$(CPP) $(CPPFLAGS) EstimateTrend.cpp Control.o Minimizer.o \
	$(noisemodel_objects) $(likelihood_objects) Observations.o \
	DesignMatrix.o HyperGeo.o Calendar.o  \
	 -DVERSION=$(VERSION) -DOMP=$(OMP) \
	$(INC_DIR) $(LIB_DIR) $(ATLAS) -lfftw3 -lgsl $(ATLAS) \
	-o estimatetrend

removeoutliers : RemoveOutliers.cpp RemoveOutliers.h \
	Observations.o DesignMatrix.o Control.o Calendar.o
	$(CPP) $(CPPFLAGS) RemoveOutliers.cpp Observations.o DesignMatrix.o \
	Control.o Calendar.o -DVERSION=$(VERSION) -DOMP=$(OMP) \
	$(INC_DIR) $(LIB_DIR) $(ATLAS) \
	-o removeoutliers

modelspectrum : ModelSpectrum.cpp ModelSpectrum.h \
	$(noisemodel_objects) Control.o Observations.o HyperGeo.o Calendar.o \
	$(likelihood_objects) DesignMatrix.o
	$(CPP) $(CPPFLAGS) ModelSpectrum.cpp $(noisemodel_objects) \
	$(likelihood_objects) DesignMatrix.o \
        Control.o Observations.o HyperGeo.o Calendar.o -DVERSION=$(VERSION) \
	$(INC_DIR) $(LIB_DIR) $(ATLAS) -lfftw3 -lgsl \
        -o modelspectrum

estimatespectrum : EstimateSpectrum.cpp Observations.o Control.o \
	Calendar.o
	$(CPP) $(CPPFLAGS) EstimateSpectrum.cpp Observations.o Control.o \
	Calendar.o -DVERSION=$(VERSION) \
	$(INC_DIR) $(LIB_DIR) -lfftw3 \
        -o estimatespectrum

simulatenoise : SimulateNoise.cpp SimulateNoise.h $(STAnoisemodel_objects) \
	Control.o HyperGeo.o NoiseModel.o $(likelihood_objects) \
	DesignMatrix.o Observations.o Calendar.o
	$(CPP) $(CPPFLAGS) SimulateNoise.cpp $(STAnoisemodel_objects) \
	$(likelihood_objects) DesignMatrix.o Observations.o Calendar.o \
        Control.o HyperGeo.o NoiseModel.o -DVERSION=$(VERSION) -DOMP=$(OMP) \
        $(INC_DIR) $(LIB_DIR) $(ATLAS) -lfftw3 -lgsl \
        -o simulatenoise

date2mjd : date2MJD.cpp Calendar.o
	$(CPP) $(CPPFLAGS) date2MJD.cpp Calendar.o -o date2mjd

mjd2date : MJD2date.cpp Calendar.o
	$(CPP) $(CPPFLAGS) MJD2date.cpp Calendar.o -o mjd2date

#---------- Subroutines

MLE.o : MLE.cpp MLE.h
	$(CPP) $(CPPFLAGS) -c MLE.cpp \
        $(INC_DIR)

Observations.o : Observations.cpp Observations.h Control.o Calendar.o
	$(CPP) $(CPPFLAGS) -c Observations.cpp 

Control.o : Control.cpp Control.h
	$(CPP) $(CPPFLAGS) -c Control.cpp

NoiseModel.o : NoiseModel.cpp NoiseModel.h Control.o 
	$(CPP) $(CPPFLAGS) -c NoiseModel.cpp \
        $(INC_DIR)

Minimizer.o : Minimizer.cpp Minimizer.h Control.o NoiseModel.o \
	Likelihood.o
	$(CPP) $(CPPFLAGS) -c Minimizer.cpp \
        $(INC_DIR)

Likelihood.o : Likelihood.cpp Likelihood.h Control.o DesignMatrix.o \
	Observations.o 
	$(CPP) $(CPPFLAGS) -c Likelihood.cpp \
        $(INC_DIR)

DesignMatrix.o : DesignMatrix.cpp DesignMatrix.h Control.o Observations.o \
	Calendar.o
	$(CPP) $(CPPFLAGS) -c DesignMatrix.cpp \
        $(INC_DIR)

White.o : White.cpp White.h NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c White.cpp

AmmarGrag.o : AmmarGrag.cpp AmmarGrag.h NoiseModel.o Likelihood.o
	$(CPP) $(CPPFLAGS) -c AmmarGrag.cpp -DOMP=$(OMP) \
        $(INC_DIR)

FullCov.o : FullCov.cpp FullCov.h Likelihood.o NoiseModel.o
	$(CPP) $(CPPFLAGS) -c FullCov.cpp \
        $(INC_DIR)

Powerlaw.o : Powerlaw.cpp Powerlaw.h NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c Powerlaw.cpp

PowerlawApprox.o : PowerlawApprox.cpp PowerlawApprox.h \
	NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c PowerlawApprox.cpp \
        $(INC_DIR)

ARFIMA.o : ARFIMA.cpp ARFIMA.h HyperGeo.o NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c ARFIMA.cpp \
        $(INC_DIR)

GenGaussMarkov.o : GenGaussMarkov.cpp GenGaussMarkov.h 
	$(CPP) $(CPPFLAGS) -c GenGaussMarkov.cpp 

HyperGeo.o : HyperGeo.cpp HyperGeo.h 
	$(CPP) $(CPPFLAGS) -c HyperGeo.cpp $(INC_DIR)

Calendar.o : Calendar.cpp Calendar.h
	$(CPP) $(CPPFLAGS) -c Calendar.cpp
