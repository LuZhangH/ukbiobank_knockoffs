#!/bin/sh

if [ -z "$1" ] || [ -z "$2" ] ; then
    echo "Minimum and maximum Job Number arguments re required."
    exit 1
fi

minjobnum="$1"
maxjobnum="$2"

myself="$(id -u -n)"

for j in $(squeue --user="$myself" --noheader --format=%i) ; do
  if [ "$j" -gt "$minjobnum" ] ; then
    if [ "$j" -lt "$maxjobnum" ] ; then
      echo "scancel $j"
      scancel "$j"
    fi
  fi
done
