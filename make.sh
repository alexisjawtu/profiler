#!/usr/bin/bash

current=20
echo "std $current"
c++ -std=c++$current ./source/*.cpp -I ./include  -o ./profile.o

echo "Compilation finished with return status $?"
