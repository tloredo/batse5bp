MODULE = batse5bp
VER = 0.1

DIST = dist/$(MODULE)-$(VER).tar.gz

install:
	pip install $(DIST) 

# Note that --upgrade attempts to upgrade other packages, e.g., numpy,
# so we instead uninstall and reinstall the target.
reinstall:
	pip uninstall $(MODULE)
	python setup.py sdist
	pip install $(DIST)