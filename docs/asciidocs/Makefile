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
asciidoctor -b manpage odgi_view.adoc
asciidoctor -b manpage odgi_kmers.adoc
asciidoctor -b manpage odgi_unitig.adoc
asciidoctor -b manpage odgi_viz.adoc
asciidoctor -b manpage odgi_paths.adoc
asciidoctor -b manpage odgi_prune.adoc
asciidoctor -b manpage odgi_unchop.adoc
asciidoctor -b manpage odgi_normalize.adoc
asciidoctor -b manpage odgi_subset.adoc
asciidoctor -b manpage odgi_bin.adoc
asciidoctor -b manpage odgi_matrix.adoc

mv *.1 man/