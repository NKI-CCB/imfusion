GH_PAGES_SOURCES = docs Makefile

gh-pages:
	git checkout gh-pages
	rm -rf ./*
	git checkout develop $(GH_PAGES_SOURCES)
	git reset HEAD
	(cd docs && make html)
	mv -fv docs/_build/html/* ./
	rm -rf $(GH_PAGES_SOURCES)
	touch .nojekyll
	git add -A
	git commit -m "Generated gh-pages for `git log develop -1 --pretty=short --abbrev-commit`" && git push origin gh-pages ; git checkout develop
