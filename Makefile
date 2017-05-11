.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help
define BROWSER_PYSCRIPT
import os, webbrowser, sys
try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

## remove all build, test, coverage and Python artifacts
clean: clean-build clean-pyc clean-test

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -f .coverage
	rm -fr htmlcov/

lint: ## check style with pylint
	pylint src/imfusion

test: clean-pyc ## run tests quickly with the default Python
	py.test tests

coverage: ## check code coverage quickly with the default Python
	py.test tests --cov=imfusion --cov-report=html
	$(BROWSER) htmlcov/index.html

tox: clean-pyc ## run tests in multiple pythons using tox
	rm -rf .tox
	cp tests/matplotlibrc ./
	docker run -v `pwd`:/app -t -i jrderuiter/tox-base
	rm -rf matplotlibrc .tox

docs: ## generate and serve Sphinx documentation
	sphinx-autobuild docs docs/_build

gh-pages:
	git checkout gh-pages
	find ./* -not -path '*/\.*' -prune -exec rm -r "{}" \;
	git checkout develop docs Makefile src AUTHORS.rst CONTRIBUTING.rst HISTORY.rst README.rst
	git reset HEAD
	(cd docs && make html)
	mv -fv docs/_build/html/* ./
	rm -rf docs Makefile src AUTHORS.rst CONTRIBUTING.rst HISTORY.rst README.rst
	touch .nojekyll
	git add -A
	git commit -m "Generated gh-pages for `git log develop -1 --pretty=short --abbrev-commit`" && git push origin gh-pages ; git checkout develop

dist: clean ## builds source and wheel package
	rm -rf build
	python setup.py sdist bdist_wheel

conda: clean-pyc ## build a conda release
	conda build --python 3.5 -c bioconda conda

conda-docker: clean-pyc
	docker run -v `pwd`:/imfusion -t -i condaforge/linux-anvil /bin/sh -c 'cd /imfusion && ./scripts/conda_build_docker.sh'

install: clean ## install the package to the active Python's site-packages
	python setup.py install
