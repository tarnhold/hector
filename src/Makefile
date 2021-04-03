# Makefile for the hector-1.7 programs
#
#  This script is part of Hector 1.7.2
#
#  Hector is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  Hector is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Hector.  If not, see <http://www.gnu.org/licenses/>
#
#===============================================================================

VERSION=1.6

#--- Check on which machine this script is running
HOST=$(shell hostname)


#--- Include directories

#--- Hector can be compiled on Mac OS with MacPorts using Gnu C++. I am not
#    able to create static binaries.
ifeq "$(HOST)" "Machiels-MacBook"
  CPP     = /opt/local/bin/g++-mp-6   # Gnu compiler 
  CPPFLAGS= -O3 -fopenmp
  INC_DIR = -I/opt/local/include 
  LIB_DIR = -L/opt/local/lib 
  PREFIX  =/opt/local
  ATLAS   = -latlas -llapack -lcblas -lblas

#--- These settings are used on a CentOS machine
#    Note that I have compiled my own Atlas library from
#    source code which gave better computation times. 
else ifeq "$(HOST)" "carga.vguest.segal"
  CPP = g++
  PREFIX=/usr/local
  CPPFLAGS= -O3 -fopenmp -static
  INC_DIR = -I./ -I/usr/include/atlas
  LIB_DIR = -L/usr/local/lib64/
  ATLAS   = -lmy_lapack -lmy_cblas -lmy_atlas # use my own compiled Atlas lib

#--- Fedora Core 21 V.machine
else ifeq "$(HOST)" "FC21"
  CPP = g++
  $(info ./src/Makefile: PAULS Fedora Core 21 V.machine)
  PREFIX=/usr/local
  CPPFLAGS= -O3 -fopenmp
  INC_DIR = -I./ -I/usr/include/atlas-x86_64-base
   LIB_DIR = -L/usr/lib64/atlas/
  ATLAS   = -llapack -lcblas -lblas -latlas 

#--- Standard values which work in many cases
else 
  CPP = g++
  $(info ./src/Makefile: Unknown machine, I am using standard values)
  PREFIX=/usr/local
  CPPFLAGS= -O3 -fopenmp -static 
  INC_DIR = -I./ -I/usr/include/atlas
  LIB_DIR = -L/usr/lib/atlas-base/
  ATLAS   = -llapack_atlas -latlas -lcblas
endif

#--- Directory where programs will be installed
BIN_DIR=$(PREFIX)/bin

#--- convenient group names
likelihood_objects = Likelihood.o AmmarGrag.o FullCov.o MLEBase.o
noisemodel_objects = NoiseModel.o White.o Powerlaw.o PowerlawApprox.o \
	ARFIMA.o GenGaussMarkov.o 
programs = estimatetrend modelspectrum estimatespectrum simulatenoise \
	removeoutliers date2mjd mjd2date 

all : $(programs)

clean:
	rm -f *.o $(programs)

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
	 -DVERSION=$(VERSION) \
	$(INC_DIR) $(LIB_DIR) $(ATLAS) -lfftw3 -lgsl $(ATLAS) \
	-o estimatetrend

removeoutliers : RemoveOutliers.cpp RemoveOutliers.h \
	Observations.o DesignMatrix.o Control.o Calendar.o
	$(CPP) $(CPPFLAGS) RemoveOutliers.cpp Observations.o DesignMatrix.o \
	Control.o Calendar.o -DVERSION=$(VERSION) \
	$(INC_DIR) $(LIB_DIR) $(ATLAS) -lgsl $(ATLAS) \
	-o removeoutliers

modelspectrum : ModelSpectrum.cpp ModelSpectrum.h \
	$(noisemodel_objects) Control.o Observations.o HyperGeo.o Calendar.o \
	$(likelihood_objects) DesignMatrix.o
	$(CPP) $(CPPFLAGS) ModelSpectrum.cpp $(noisemodel_objects) \
	$(likelihood_objects) DesignMatrix.o \
        Control.o Observations.o HyperGeo.o Calendar.o -DVERSION=$(VERSION) \
	$(INC_DIR) $(LIB_DIR) $(ATLAS) -lfftw3 -lgsl $(ATLAS) \
        -o modelspectrum

estimatespectrum : EstimateSpectrum.cpp Observations.o Control.o \
	Calendar.o
	$(CPP) $(CPPFLAGS) EstimateSpectrum.cpp Observations.o Control.o \
	Calendar.o -DVERSION=$(VERSION) \
	$(INC_DIR) $(LIB_DIR) -lfftw3 \
        -o estimatespectrum

simulatenoise : SimulateNoise.cpp SimulateNoise.h $(noisemodel_objects) \
	Control.o HyperGeo.o $(likelihood_objects) \
	DesignMatrix.o Observations.o Calendar.o
	$(CPP) $(CPPFLAGS) SimulateNoise.cpp $(noisemodel_objects) \
	$(likelihood_objects) DesignMatrix.o Observations.o Calendar.o \
        Control.o HyperGeo.o -DVERSION=$(VERSION) \
        $(INC_DIR) $(LIB_DIR) $(ATLAS) -lfftw3 -lgsl \
        -o simulatenoise

date2mjd : date2MJD.cpp Calendar.o
	$(CPP) $(CPPFLAGS) $(LIB_DIR) date2MJD.cpp Calendar.o -o date2mjd

mjd2date : MJD2date.cpp Calendar.o
	$(CPP) $(CPPFLAGS) $(LIB_DIR) MJD2date.cpp Calendar.o -o mjd2date


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

Minimizer.o : Minimizer.cpp Minimizer.h Control.o $(noisemodel_objects) \
	$(likelihood_objects) DesignMatrix.o
	$(CPP) $(CPPFLAGS) -c Minimizer.cpp \
        $(INC_DIR)

Likelihood.o : Likelihood.cpp Likelihood.h Control.o DesignMatrix.o \
	Observations.o MLEBase.o AmmarGrag.o FullCov.o
	$(CPP) $(CPPFLAGS) -c Likelihood.cpp \
        $(INC_DIR)

MLEBase.o : MLEBase.cpp MLEBase.h Control.o DesignMatrix.o \
	Observations.o 
	$(CPP) $(CPPFLAGS) -c MLEBase.cpp \
        $(INC_DIR)

DesignMatrix.o : DesignMatrix.cpp DesignMatrix.h Control.o Observations.o \
	Calendar.o
	$(CPP) $(CPPFLAGS) -c DesignMatrix.cpp \
        $(INC_DIR)

White.o : White.cpp White.h NoiseModel.o NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c White.cpp

AmmarGrag.o : AmmarGrag.cpp AmmarGrag.h NoiseModel.o MLEBase.o
	$(CPP) $(CPPFLAGS) -c AmmarGrag.cpp \
        $(INC_DIR)

FullCov.o : FullCov.cpp FullCov.h MLEBase.o NoiseModel.o
	$(CPP) $(CPPFLAGS) -c FullCov.cpp \
        $(INC_DIR)

Powerlaw.o : Powerlaw.cpp Powerlaw.h NoiseModel.o NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c Powerlaw.cpp

PowerlawApprox.o : PowerlawApprox.cpp PowerlawApprox.h \
	NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c PowerlawApprox.cpp \
        $(INC_DIR)

ARFIMA.o : ARFIMA.cpp ARFIMA.h HyperGeo.o NoiseModel.o NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c ARFIMA.cpp \
        $(INC_DIR)

GenGaussMarkov.o : GenGaussMarkov.cpp GenGaussMarkov.h NoiseModel.o \
	NoiseModelBaseClass.h
	$(CPP) $(CPPFLAGS) -c GenGaussMarkov.cpp 

HyperGeo.o : HyperGeo.cpp HyperGeo.h 
	$(CPP) $(CPPFLAGS) -c HyperGeo.cpp $(INC_DIR)

Calendar.o : Calendar.cpp Calendar.h
	$(CPP) $(CPPFLAGS) -c Calendar.cpp
