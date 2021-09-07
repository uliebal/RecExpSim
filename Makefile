# Why does a Makefile exist? Makefiles offer a system operative and language agnostic way to
# describe batch commands. They can hold python module calls, commands for git hooks and other
# utilities to help in development.

.PHONY: all help try-% clean

all: help

help:
	echo "No operation. Check the file named "Makefile" for possible operations."

# Run one of the python scripts inside the "try" folder.
try-%:
	python -m try.$(@:try-%=%)

# Clean all of the python cache files.
clean:
	find . -type f -name ‘*.pyc’ -delete
