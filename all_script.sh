#! /bin/bash


assign_inputs()
{
	c1=$1
	c2=$2
	c3=$3
	c4=$4
	c5=$5
	c6=$6
	c7=$7
}

change_input()
{
echo -e "$c1\n$c2\n$c3\n$c4\n$c5\n$c6\n$c7" > input
}

echo "A=vector(22664);"

COUNTER=0
while read LINE; do
	assign_inputs $LINE
#	change_input
#	./tester < input > outputs/9-13-14-126-p-5/$c1-$c2-$c3-$c4-$c5-$c6-$c7
	COUNTER=$(($COUNTER+1))
	tail -n9 outputs/"9-13-14-126-p-5"/"$c1-$c2-$c3-$c4-$c5-$c6-$c7" | \
		sed -e "s/^\[/A\[$COUNTER\]=[/" \
			-e "s/\]$/\];/" \
			-e "s/This.*frobenius.//"
done
