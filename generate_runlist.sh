#!/bin/bash

# This script generates a runlist for use in the analysis.
# Run with:
# sh generate_runlist.sh

# the directory we are in
location="$(pwd)"
echo "running in "$location

# create the output we write to
outputfile="outputlist.txt"
touch $outputfile
echo "# Runlist:" > $outputfile

# the sensorsorting
sensorsort="0432156789"

# any broken sensors
broken="4"

# sample thickness
thickness="50"

# dummy filename
oldfile="fail"

# read in the grep result, up to the first : to result, remainder into temp
while IFS=: read -r result temp
do
    # strip away the .xml
    file="$(echo $result | awk -F'.xml' '{print $1}')"
    appendix='.root'
    newfile=$file$appendix
    echo "File to open:" $newfile
    # filename is now in newfile

    # remove everyting up to the >
    descr="$(echo $temp | awk -F'>' '{print $2}')"
    # remove the </Log
    descr="$(echo $descr | awk -F'<' '{print $1}')"
    echo "File description:" $descr
    # descr now holds the file info

    # did we just write this file?

    if [ $newfile != $oldfile ]; then

	# write the vars to the output
	echo $newfile,$sensorsort,$broken,$thickness,$descr >> $outputfile
	echo "New file!"
	
    fi

    oldfile=$newfile

    echo " "

done <<< "`grep -rnw '../' -e "Log time" | grep -v "\thomas/*" | grep -v "\nils/*"`"
# grep one directory up, -e "" and pipe the output excluding some subdirs...



echo "Done!"
