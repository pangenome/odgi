#CUR_VER=$(tail -2 VERSIONS | head -1)
CUR_VER=$(grep "{\"v" ../../src/version.cpp | tail -2 | head -1 | cut -d "," -f 1 | cut -d "\"" -f 2)
echo "$CUR_VER"

NEW_VER=$(grep "{\"v" ../../src/version.cpp | tail -n 1 | cut -d "," -f 1 | cut -d "\"" -f 2)
echo "$NEW_VER"

sed -i "s/$CUR_VER/$NEW_VER/g" *.adoc

asciidoctor -b html odgi_docs.adoc
asciidoctor -b docbook odgi_docs.adoc
asciidoctor -r asciidoctor-pdf -b pdf odgi_docs.adoc

asciidoctor -b manpage odgi.adoc
asciidoctor -b manpage odgi_build.adoc
asciidoctor -b manpage odgi_stats.adoc
asciidoctor -b manpage odgi_degree.adoc
asciidoctor -b manpage odgi_depth.adoc
asciidoctor -b manpage odgi_overlap.adoc
asciidoctor -b manpage odgi_sort.adoc
asciidoctor -b manpage odgi_view.adoc
asciidoctor -b manpage odgi_kmers.adoc
asciidoctor -b manpage odgi_unitig.adoc
asciidoctor -b manpage odgi_viz.adoc
asciidoctor -b manpage odgi_draw.adoc
asciidoctor -b manpage odgi_paths.adoc
asciidoctor -b manpage odgi_prune.adoc
asciidoctor -b manpage odgi_unchop.adoc
asciidoctor -b manpage odgi_normalize.adoc
asciidoctor -b manpage odgi_extract.adoc
asciidoctor -b manpage odgi_bin.adoc
asciidoctor -b manpage odgi_matrix.adoc
asciidoctor -b manpage odgi_chop.adoc
asciidoctor -b manpage odgi_layout.adoc
asciidoctor -b manpage odgi_flatten.adoc
asciidoctor -b manpage odgi_break.adoc
asciidoctor -b manpage odgi_pathindex.adoc
asciidoctor -b manpage odgi_panpos.adoc
asciidoctor -b manpage odgi_server.adoc
asciidoctor -b manpage odgi_test.adoc
asciidoctor -b manpage odgi_version.adoc
asciidoctor -b manpage odgi_cover.adoc
asciidoctor -b manpage odgi_explode.adoc
asciidoctor -b manpage odgi_squeeze.adoc
asciidoctor -b manpage odgi_validate.adoc
asciidoctor -b manpage odgi_position.adoc

mv *.1 man/