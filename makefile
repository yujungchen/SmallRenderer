INCLUDE_DIR=include
SOURCE_DIR=src
CC=gcc
CXX=g++
INCLUDE=-I$(INCLUDE_DIR) 
CFLAGS=$(INCLUDE)

OBJECT_DIR=obj
LIB_DIR=lib

_DEPS = core.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = utility.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = glm.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = configure.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = bvh.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = light.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = radiometry.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = camera.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = sampler.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = microfacet.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = dipole.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = volume.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = direct.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = path.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = imgPPM.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_DEPS = mcRenderer.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_OBJ = main.o mcRenderer.o imgPPM.o glm.o configure.o core.o utility.o bvh.o light.o sampler.o camera.o radiometry.o volume.o dipole.o microfacet.o direct.o path.o
OBJ = $(patsubst %,$(OBJECT_DIR)/%,$(_OBJ))

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(DEPS)
	@echo "Building object $@"
	$(CXX) -O3 -c -msse3 -fopenmp -lm -o $@ $< $(CFLAGS)

all: SmallRenderer

SmallRenderer: $(OBJ)
	$(CXX) -O3 -msse3 -fopenmp -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(OBJECT_DIR)/*.o