#!/bin/bash

./cluster points.csv -m complete > points_complete_test.txt

diff points_complete_test.txt points_complete.txt

if [ $? -ne 0 ]; then
	echo "Test failed"
	exit 1
else
	echo "Test passed"
fi
