# how to (the German way)
## tools
pandoc
pip install sphinx --user
pip install sphinx_rtd_theme --user
pip install m2r --user 
asciidoc

https://stackoverflow.com/questions/45967058/render-output-from-markdown-file-inside-rst-file

## we go from asciidoc to markdown
### folder asciidocs
asciidoc -b docbook odgi_bin.adoc
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
### disable caching
rm -rf sphinx_build && sphinx-build -b html ./ sphinx_build



# how to (the Italian way)
Assuming you have Python already, install Sphinx:

```
pip install sphinx
```

Then, go into your documentation directory

```
cd docs
```

and edit your `*.rst` files. Finally, build them to see how they look:

```
make html
```
