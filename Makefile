
.PHONY: all help try-% clean

all: help

help:
	echo "No operation. Try running 'make try-<filename>' to run a script."

try-%:
	python -m try.$(@:try-%=%)

clean:
	find . -type f -name ‘*.pyc’ -delete
