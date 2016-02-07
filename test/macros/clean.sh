#!/bin/bash

rm -f att.txt

for file in `ls -trd1 jobs4mu/*.sh.* | head -500`; do
 echo $file
 echo "rm -f $file" >> att.txt
done

bash att.txt

for file in `ls -trd1 jobs4e/*.sh.* | head -500`; do
 echo $file
 echo "rm -f $file" >> att.txt
done

bash att.txt

for file in `ls -trd1 jobs2e2mu/*.sh.* | head -500`; do
 echo $file
 echo "rm -f $file" >> att.txt
done

bash att.txt
