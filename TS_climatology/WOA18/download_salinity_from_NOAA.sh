#!/bin/bash
mkdir 1995-2004
mkdir 2005-2017
cwd=$(pwd)

# downloads annual (0), and monthly (1-12) salinity climatology fields
for i in `seq -w 0 12`
do
        echo $i
        # 1995-2004
        cd $cwd/1995-2004
        string1=https://data.nodc.noaa.gov/thredds/fileServer/ncei/woa/salinity/95A4/1.00/woa18_95A4_s$i
        string2="_01.nc" 
        string3=$string1$string2
        echo "\n"
        wget $string3 .

        cd ..

        # 2005-2017
        cd $cwd/2005-2017
        string1=https://data.nodc.noaa.gov/thredds/fileServer/ncei/woa/salinity/A5B7/1.00/woa18_A5B7_s$i
        string2="_01.nc" 
        string3=$string1$string2
        echo $string3
        echo "\n"
        wget $string3 .

        cd $cwd
        find .
done

