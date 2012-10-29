install:
	pip install dist/batse5bp-0.1.tar.gz 

reinstall:
	pip uninstall batse5bp
	python setup.py sdist
	pip install dist/batse5bp-0.1.tar.gz