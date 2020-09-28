#!/bin/bash
mkdir 1995-2004
mkdir 2005-2017
cwd=$(pwd)

# downloads annual (0), and monthly (1-12) temperature climatology fields
for i in `seq -w 0 12`
do
        echo $i
        # 1995-2004
        cd $cwd/1995-2004
        string1=https://data.nodc.noaa.gov/thredds/fileServer/ncei/woa/temperature/95A4/1.00/woa18_95A4_t$i
        string2="_01.nc" 
        string3=$string1$string2
        echo "\n"
        wget -nc $string3 .

        cd ..

        # 2005-2017
        cd $cwd/2005-2017
        string1=https://data.nodc.noaa.gov/thredds/fileServer/ncei/woa/temperature/A5B7/1.00/woa18_A5B7_t$i
        string2="_01.nc" 
        string3=$string1$string2
        echo $string3
        echo "\n"
        wget -nc $string3 .

        cd $cwd
        find .
done

