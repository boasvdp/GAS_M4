#!/bin/bash

set -e

if [[ -s workflow/scripts/maskrc-svg.py ]]
then
	echo "Script has been downloaded before and size > 0 bytes"
else
	wget -O workflow/scripts/maskrc-svg.py https://raw.githubusercontent.com/kwongj/maskrc-svg/master/maskrc-svg.py
fi
