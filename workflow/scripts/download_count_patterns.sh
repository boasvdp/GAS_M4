#!/bin/bash

set -e

if [[ -s workflow/scripts/count_patterns.py ]]
then
	echo "Script has been downloaded before and size > 0 bytes"
else
	wget -O workflow/scripts/count_patterns.py https://raw.githubusercontent.com/mgalardini/pyseer/master/scripts/count_patterns.py
fi
