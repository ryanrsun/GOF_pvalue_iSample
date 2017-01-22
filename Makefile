# $@ is the file name of the target
# $^ is the name of all dependencies separated by spaces, duplicate names removed
# -c means compile only (no executable, only make the .o file)
#

CXX=g++ -std=c++11
INCLUDES = -I/users/ryansun/Documents/Research/Paper2/Software/includes -I/usr/local/include
CXXFLAGS = -g -Wall $(INCLUDES)
LDFLAGS = -g -L/usr/local/lib
LDLIBS = -larmadillo

all: GOF_pvalue_iSample

# 
# It doesn't automatically know the name of the binary
#

GOF_pvalue_iSample: GOF_pvalue_iSample.o 
	$(CXX) $(LDFLAGS) GOF_pvalue_iSample.o $(LDLIBS) -o GOF_pvalue_iSample 

#
# It automatically knows to make the .o from the .c, that's why we only need $^
#

GOF_pvalue_iSample.o: GOF_pvalue_iSample.cpp 
	$(CXX) $(CXXFLAGS) -c $^ 

clean:
	rm -f *.o GOF_pvalue_iSample

