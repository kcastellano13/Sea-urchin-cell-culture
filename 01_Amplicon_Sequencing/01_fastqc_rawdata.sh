#!/bin/bash
#fastqc of 18S sea urchin cell culture data

fileList="*.fastq.gz"

for file in ${fileList}
do

        fastqc ${file} -t 5

done
