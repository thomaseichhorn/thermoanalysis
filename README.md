# thermoanalysis

A small ROOT script for analysis of data from the thermosetup.
Data is read from the ROOT tuples specified in a runlist.

Compile this script with (ROOT environment must be loaded):
g++ -I `root-config --incdir` -o test main.cc `root-config --libs` -Wall -std=c++0x -pedantic

Then run: ./test /path/to/runlist
Change test to a different executable name if you want.

Directly run this script in ROOT with:

root -l -b
.x main.cc

Then input the runlist when prompted.

