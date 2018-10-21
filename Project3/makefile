# parameters
CXX      = g++ -std=c++14
CXXFLAGS = -Wall -O3
DEPS     = vectorlib.h gravity.h
OBJ      = vectorlib.o gravity.o test.o
LIBS     = 

# make object files
%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# make executable
test: $(OBJ)
	$(CXX) -o $@.x $^ $(CXXFLAGS) $(LIBS)

# remove object files
.PHONY: clean
clean:
	rm -f ./*.o *~ core $(INCDIR)/*~
