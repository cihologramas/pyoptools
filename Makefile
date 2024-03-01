inplace:
	python3 setup.py build_ext -i
deb:
	dch -b Paquete de prueba generado
	dpkg-buildpackage -rfakeroot -uc -us -b

dsc:
	dch -b Paquete de prueba generado
	dpkg-source -b .

clean:
	rm deb_dist dist build -rf
	rm pyoptools*.tar.gz -f
	python3 setup.py clean
	python3 clean.py
