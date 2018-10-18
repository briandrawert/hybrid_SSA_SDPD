#!/bin/bash

counter=1
while [ $counter -le 100 ]
do
	echo Create directory $counter
	rm -r $counter
	((counter++))
	
done
echo All done

