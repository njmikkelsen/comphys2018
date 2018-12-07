# parameters
CXX      = g++ -std=c++14
CXXFLAGS = -Wall -O3
DEPS     = vectorlib.h vmc.h
OBJ      = vectorlib.o vmc.o

# make object files
%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# make executable
run: $(OBJ) $(PROG)
	$(CXX) -o $@.x $^ $(CXXFLAGS)

# remove object files
.PHONY: clean
clean:
	rm -f ./*.o *~ core $(INCDIR)/*~
