inplace:
	python setup.py build_ext -i
# DEB_BUILD_OPTIONS=nocheck make deb
deb:
	python setup.py --command-packages=stdeb.command sdist_dsc \
	 bdist_deb

clean:
	rm deb_dist dist build -rf
	rm pyoptools*.tar.gz -f
	python setup.py clean
	python clean.py
