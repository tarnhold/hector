# Makefile for the hector-1.9 programs
#
#  This script is part of Hector 1.9
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

VERSION=1.9

CPP      = g++
CPPFLAGS = -O2 -Xpreprocessor -fopenmp -std=c++1z 
INC_DIR  = -I/usr/include/openblas/ -I/usr/local/Cellar/openblas/0.3.13/include/
PREFIX   =/usr/local
OPENBLAS = -lopenblas
LIB_DIR  = -L/usr//local/Cellar/openblas/0.3.13/lib/ -lomp
FFTW     = -lfftw3

#--- Directory where programs will be installed
BIN_DIR=$(PREFIX)/bin

#--- convenient group names
programs = estimatetrend modelspectrum estimatespectrum simulatenoise \
	removeoutliers date2mjd mjd2date findoffset

all : $(programs)

clean:
	rm -f *.o $(programs)
	rm -rf *.dSYM

install:
	mkdir -p $(BIN_DIR)
	install $(programs) $(BIN_DIR)

uninstall :
	@rm -f $(addprefix $(BIN_DIR)/,$(programs))


include Makefile_common


test_Matern : test_Matern.cpp Control.o Observations.o NoiseModelBaseClass.h \
	Calendar.o Matern.o
	$(CPP) $(CPPFLAGS) test_Matern.cpp Control.o Observations.o \
        Calendar.o Matern.o \
        $(INC_DIR) $(OPENBLAS) $(FFTW) $(LIB_DIR) \
        -o test_Matern

test_2F1 : test_2F1.cpp
	$(CPP) $(CPPFLAGS) test_2F1.cpp \
        $(INC_DIR) $(OPENBLAS) $(LIB_DIR) \
        -o test_2F1
