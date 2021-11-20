all: test1.png DistribucionIndividuos.html

test1.png: test1.py
	python3 $<

DistribucionIndividuos.html: DistribucionIndividuos.ipynb
	jupyter-nbconvert --exec DistribucionIndividuos.ipynb
