#!/bin/sh

echo "======  pep8  ======"
pep8 $1
echo "======  pyflakes  ======"
pyflakes $1
echo "======  pylint  ======"
pylint --output-format=parseable $1