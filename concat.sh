SAMPLE=$(cat _sample.txt)

for ((i = 0; i != $SAMPLE; ++i))
do
	cat "$i.dat"
	echo
done
