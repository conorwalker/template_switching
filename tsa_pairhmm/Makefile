# Copyright (C) 2020 EMBL - European Bioinformatics Institute
# Contact: goldman@ebi.ac.uk, cwalker@ebi.ac.uk

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.

# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.


####################################
# Makefile for building: tsa_pairhmm
####################################

MAKEFILE      = Makefile


####### Compiler, tools and options
CC            = gcc
CXX           = g++
CFLAGS        = -m64 -pipe -O2 -Wall -std=c++11
CXXFLAGS      = -m64 -pipe -O2 -Wall -std=c++11
INCPATH       = -I. 
LINK          = g++
LFLAGS        = -m64 -Wl,-O1
LIBS          = -lboost_program_options $(SUBLIBS)
DEL_FILE      = rm -f

####### Output directory
OBJECTS_DIR   = ./


####### Files
SOURCES       = tsa_pairhmm.cpp 
OBJECTS       = tsa_pairhmm.o

DESTDIR       = 
TARGET        = tsa_pairhmm

first: all


####### Implicit rules
.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"


####### Build rules
all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

clean: compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core

compiler_clean: 


####### Compile
tsa_pairhmm.o: tsa_pairhmm.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o tsa_pairhmm.o tsa_pairhmm.cpp


####### Install
install:  FORCE

uninstall:  FORCE

FORCE:
