#! /bin/bash

echo "Usage: elliptic.sh p N q a b"

P=$1
N=$2
Q=$3
A=$4
B=$5

echo $P $N $A $B

if [ -f inputs/tijdelijk ]; then
	echo "Please remove the file inputs/tijdelijk. Stop."
	exit 1
fi
grep -q "^#define LEX_ORDER" data.h
if [ ! $? = 0 ]; then
	echo "Make sure you use LEX_ORDER. Stop."
	exit 1
fi
DISC=$((27*$B*$B + 4*$A*$A*$A))
DISC=$((DISC % $P))
if [ $DISC = 0 ]; then
	echo "Discriminant is zero! Stop."
	exit 1
fi

cat >> inputs/tijdelijk << AISEJOHAN
1
0
0
0
0
a
0
-1
0
b
AISEJOHAN

sed -i -e "s@a@$A@" \
	-e "s@b@$B@" inputs/tijdelijk

sed -i \
	-e "s@^#define d1\t.*@#define d1\t1@" \
	-e "s@^#define d2\t.*@#define d2\t1@" \
	-e "s@^#define d3\t.*@#define d3\t1@" \
	-e "s@^#define d\t.*@#define d\t3@" \
	-e "s@^#define p\t.*@#define p\t$P@" \
        -e "s@^#define r\t.*@#define r\t$N@" \
        -e "s@^#define q\t.*@#define q\t$Q@" data.h

cat data.h
echo
echo "Does this look ok? (Hit return to continue,Ctrl-C to abort.)"
read ANTWOORD

make input_pol

time ./tester < inputs/tijdelijk

rm inputs/tijdelijk

exit 0
