all:
	@g++ -O3 -w -pthread -lpthread karect.cpp -o karect

install:
	@mv karect /usr/local/bin

clean:
	@rm -f karect
