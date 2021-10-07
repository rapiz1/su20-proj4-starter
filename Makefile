CC = gcc
CFLAGS = -g -Wall $(shell python3-config --cflags) -I/usr/include/python3.9 
LDFLAGS = -fopenmp -lcunit $(shell python3-config --ldflags --embed)

install:
	if [ ! -f files.txt ]; then touch files.txt; fi
	rm -rf build
	xargs rm -rf < files.txt
	python3 setup.py install --record files.txt

uninstall:
	if [ ! -f files.txt ]; then touch files.txt; fi
	rm -rf build
	xargs rm -rf < files.txt

clean:
	rm -f *.o
	rm -f test
	rm -rf build
	rm -rf __pycache__

test:
	rm -f test
	$(CC) $(CFLAGS) $(LDFLAGS) mat_test.c matrix.c -o test
	./test

.PHONY: test
