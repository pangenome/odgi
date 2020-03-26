asciidoctor -b html odgi_docs.adoc
asciidoctor -b docbook odgi_docs.adoc
asciidoctor -r asciidoctor-pdf -b pdf odgi_docs.adoc

cat odgi.adoc | sed "s/==/=/" | sed "s/ +//" > man/odgi_man.adoc
cat odgi_build.adoc | sed "s/==/=/" | sed "s/ +//" > man/odgi_build_man.adoc
cat odgi_stats.adoc | sed "s/==/=/" | sed "s/ +//" > man/odgi_stats_man.adoc
cat odgi_sort.adoc | sed "s/==/=/" | sed "s/ +//" > man/odgi_sort_man.adoc

asciidoctor -b manpage man/odgi_man.adoc
asciidoctor -b manpage man/odgi_build_man.adoc
asciidoctor -b manpage man/odgi_stats_man.adoc
asciidoctor -b manpage man/odgi_sort_man.adoc