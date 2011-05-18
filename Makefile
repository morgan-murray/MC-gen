#######################
#
#   Makefile
# 
#   BH-DVCS Generator
#
#######################

CC = g++

WARN = -Wall

CPPFLAGS := $(shell root-config --cflags)  -g
#-D__TESTING__ -D__TESTING1__ -D__TESTING_BH__
# -D__COMPASS__
LIBS := $(shell root-config --libs)

SRC-MAIN = bhdvcs.cc calcBH1.cc

SRC-GEN = $(SRC-MAIN) gen_events.cc

HEADER = physics_types.h physics_const.h bhdvcs.h

gen_events: $(SRC-GEN)

	$(CC) $(WARN) $(CPPFLAGS) $(LIBS) -lgsl -o $@ $(SRC-GEN)

bhdvcs: $(SRC-MAIN) $(HEADER)

	$(CC) $(WARN) $(CPPFLAGS) $(LIBS) -o $@ $(SRC-MAIN)

clean: 
	rm bhdvcs bhdvs-gen
