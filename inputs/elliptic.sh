#! /bin/bash

echo "Usage: elliptic.sh p N q a b"

P=$1
N=$2
Q=$3
A=$4
B=$5

echo $P $N $A $B

sed "s@^.*= @@" inputs/input3 > inputs/input
sed -i -e "s@a@$A@" \
	-e "s@b@$B@" inputs/input

sed -i \
	-e "s@^#define d1\t.*@#define d1\t1@" \
	-e "s@^#define d2\t.*@#define d2\t1@" \
	-e "s@^#define d3\t.*@#define d3\t1@" \
	-e "s@^#define p\t.*@#define p\t$P@" \
        -e "s@^#define r\t.*@#define r\t$N@" \
        -e "s@^#define q\t.*@#define q\t$Q@" data.h

make input_pol

time ./tester < inputs/input

exit 0
