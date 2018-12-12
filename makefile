all: test
test: householderQR.cpp
	g++ householderQR.cpp -std=c++11 -g -o test
clean:
	rm -rf test
