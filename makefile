all:
	g++ -O3 -Wall -pthread -lpthread karect.cpp -o karect

clean:
	rm -f karect
