inplace:
	python setup.py build_ext -i

deb:
	python setup.py --command-packages=stdeb.command sdist_dsc \
	--depends "python-scipy (>= 0.10.1)" \
	 bdist_deb

clean:
	python setup.py clean
	rm deb_dist dist build -rf
	rm pyoptools*.tar.gz -f
	python clean.py