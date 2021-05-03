## tools
pandoc
pip install sphinx --user
pip install sphinx_rtd_theme --user
asciidoc

## we go from asciidoc to markdown
### folder asciidocs
iconv -t utf-8 odgi_bin.xml | pandoc -f docbook -t markdown_strict --wrap=none | iconv -f utf-8 > odgi_bin.md
### copy to folder docs

## build man pages
### folder docs
sphinx-build -b man ./ sphinx_build_man
### take a look at it
man ./sphinx_build_man/odgi.1

## build html
### in folder docs
sphinx-build -b html ./ sphinx_build
