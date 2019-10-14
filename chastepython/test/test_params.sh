#!/bin/bash
echo "Parameter test program "
echo $@
touch ${1#"output_dir="}/testout.txt