LIBS := -L "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" -lmsmpi -ldl -lmsmpifec -lmsmpifmc
FILES := func.c main.c
INCLUDE := -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"

all: build\task

build\task: build
	gcc $(INCLUDE) $(FILES) -o build\task $(LIBS)

build:
	mkdir build

clean:
	rm -r build
