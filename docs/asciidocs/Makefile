CUR_VER=$(tail -2 VERSIONS | head -1)
echo "$CUR_VER"

NEW_VER=$(tail -n 1 VERSIONS)
echo "$NEW_VER"

sed -i "s/$CUR_VER/$NEW_VER/g" *.adoc

asciidoctor -b html odgi_docs.adoc
asciidoctor -b docbook odgi_docs.adoc
asciidoctor -r asciidoctor-pdf -b pdf odgi_docs.adoc

asciidoctor -b manpage odgi.adoc
asciidoctor -b manpage odgi_build.adoc
asciidoctor -b manpage odgi_stats.adoc
asciidoctor -b manpage odgi_sort.adoc

mv *.1 man/