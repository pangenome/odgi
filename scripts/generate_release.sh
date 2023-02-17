#!/bin/bash

cmd=$0" "$@
TEMP=`getopt -o i:o:D:a:p:n:s:l:K:F:k:x:f:B:XH:j:P:O:Me:t:T:vhASY:G:Q:d:I:R:NbrmZzV: --long input-fasta:,output-dir:,temp-dir:,input-paf:,map-pct-id:,n-mappings:,segment-length:,block-length-min:,mash-kmer:,mash-kmer-thres:,min-match-length:,sparse-map:,sparse-factor:,transclose-batch:,skip-normalization,n-haps:,path-jump-max:,subpath-min:,edge-jump-max:,threads:,poa-threads:,skip-viz,do-layout,help,no-merge-segments,do-stats,exclude-delim:,poa-length-target:,poa-params:,poa-padding:,run-abpoa,global-poa,write-maf,consensus-spec:,consensus-prefix:,pad-max-depth:,block-id-min:,block-ratio-min:,no-splits,resume,keep-temp-files,multiqc,compress,vcf-spec:,version -n 'pggb' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -i|--input-fasta) input_fasta=$2 ; shift 2 ;;
        -s|--segment-length) segment_length=$2 ; shift 2 ;;
        -l|--block-length) block_length=$2 ; shift 2 ;;
        -p|--map-pct-id) map_pct_id=$2 ; shift 2 ;;
        -n|--n-haplotypes) n_mappings=$2 ; shift 2 ;;
        -N|--no-splits) no_splits=true ; shift ;;
        -x|--sparse-map) sparse_map=$2 ; shift 2 ;;
        -K|--mash-kmer) mash_kmer=$2 ; shift 2 ;;
        -F|--mash-kmer-thres) mash_kmer_thres=$2 ; shift 2 ;;
        -Y|--exclude-delim) exclude_delim=$2 ; shift 2 ;;
        -k|--min-match-length) min_match_length=$2 ; shift 2 ;;
        -f|--sparse-factor) sparse_factor=$2 ; shift 2 ;;
        -B|--transclose-batch) transclose_batch=$2 ; shift 2 ;;
        -X|--skip-normalization) skip_normalization=true ; shift ;;
        -H|--n-haplotypes-smooth) n_haps=$2 ; shift 2 ;;
        -j|--path-jump-max) max_path_jump=$2 ; shift 2 ;;
        -e|--edge-jump-max) max_edge_jump=$2 ; shift 2 ;;
        -G|--poa-length-target) target_poa_length=$2 ; shift 2 ;;
        -P|--poa-params) poa_params=$2 ; shift 2 ;;
        -O|--poa-padding) poa_padding=$2 ; shift 2 ;;
        -d|--pad-max-depth) pad_max_depth=$2 ; shift 2 ;;
        -b|--run-abpoa) run_abpoa=true ; shift ;;
        -z|--global-poa) run_global_poa=true ; shift ;;
        -M|--write-maf) write_maf=true ; shift ;;
        #-C|--consensus-spec) consensus_spec=$2 ; shift 2 ;;
        -Q|--consensus-prefix) consensus_prefix=$2 ; shift 2 ;;
        -v|--skip-viz) do_viz=false ; do_layout=false; shift ;;
        -S|--do-stats) do_stats=true ; shift ;;
        -V|--vcf-spec) vcf_spec=$2 ; shift 2 ;;
        -m|--multiqc) multiqc=true ; shift ;;
        -o|--output-dir) output_dir=$2 ; shift 2 ;;
        -D|--temp-dir) input_temp_dir=$2 ; shift 2 ;;
        -a|--input-paf) input_paf=$2 ; shift 2 ;;
        -r|--resume) resume=true ; shift ;;
        -t|--threads) threads=$2 ; shift 2 ;;
        -T|--poa-threads) poa_threads=$2 ; shift 2 ;;
        -A|--keep-temp-files) keep_intermediate_files=true ; shift ;;
        -Z|--compress) compress=true ; shift ;;
        --version) show_version=true ; shift ;;
        -h|--help) show_help=true ; shift ;;
        --) shift ; break ;;
        *) echo "$2" "Internal error!" ; exit 1 ;;
    esac
done

SSE2=$(lscpu | grep "Flags" | grep "sse2")
SSE4_2=$(lscpu | grep "Flags" | grep "sse4_2")
AVX=$(lscpu | grep "Flags" | grep "avx")
AVX2=$(lscpu | grep "Flags" | grep "avx2")
AVX512=$(lscpu | grep "Flags" | grep "avx512")

mkdir release
sed -i '55 s/^/#/' CMakeLists.txt

if [[ ! -z "$SSE2" ]];
then
  echo "SSE2";
  cmake -H. -Bsse2 -DEXTRA_FLAGS="-Ofast -pipe -msse2" && cmake --build sse2 -- -j 15
  mv bin/odgi sse2/
  rm -r bin
fi

if [[ ! -z "$SSE4_2" ]];
then
  echo "SSE4_2";
  cmake -H. -Bsse4_2 -DEXTRA_FLAGS="-Ofast -pipe -msse4.2" && cmake --build sse4_2 -- -j 15
  mv bin/odgi sse4_2
  rm -r bin
  mv sse2/odgi release/odgi_sse2
  mv sse4_2/odgi release/odgi_sse4_2
fi

if [[ ! -z "$AVX" ]];
then
  echo "AVX";
  cmake -H. -Bavx -DEXTRA_FLAGS="-Ofast -pipe -mavx" && cmake --build avx -- -j 15
  mv bin/odgi avx
  rm -r bin
  mv avx/odgi release/odgi_avx
fi

if [[ ! -z "$AVX2" ]];
then
  echo "AVX2";
  cmake -H. -Bavx2 -DEXTRA_FLAGS="-Ofast -pipe -mavx2" && cmake --build avx2 -- -j 15
  mv bin/odgi avx2
  rm -r bin
  mv avx2/odgi release/odgi_avx2
fi

if [[ ! -z "$AVX512" ]];
then
  echo "AVX512";
  # TODO maybe only enable -mavx512f so it has the best compatibility - https://en.wikipedia.org/wiki/AVX-512#SIMD_modes - https://gcc.gnu.org/onlinedocs/gcc-6.1.0/gcc/x86-Options.html
  # AVX512_FLAGS=$(lscpu | grep "avx512" | sed 's/ /\n/g' | grep "avx512" | sed 's/avx/-mavx/g' | tr '\n' ' ')
  cmake -H. -Bavx512 -DEXTRA_FLAGS="-Ofast -pipe -mavx512f" && cmake --build avx512 -- -j 15
  mv bin/odgi avx512
  rm -r bin
  mv avx512/odgi release/odgi_avx512
fi

sed -i '55 s/^#//' CMakeLists.txt

ODGI="#!/bin/bash

# https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script#
SOURCE=\${BASH_SOURCE[0]}
while [ -L \"\$SOURCE\" ]; do # resolve \$SOURCE until the file is no longer a symlink
  DIR=\$( cd -P \"\$( dirname \"\$SOURCE\" )\" >/dev/null 2>&1 && pwd )
  SOURCE=\$(readlink \"\$SOURCE\")
  [[ \$SOURCE != /* ]] && SOURCE=\$DIR/\$SOURCE # if \$SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR=\$( cd -P \"\$( dirname \"\$SOURCE\" )\" >/dev/null 2>&1 && pwd )

if [[ -f odgi_avx512 ]];
then
  echo \"avx512\"
  \"\$DIR\"/odgi_avx512 \"\$@\"
elif [[ -f odgi_avx2 ]];
then
  echo \"avx2\"
  \"\$DIR\"/odgi_avx2 \"\$@\"
elif [[ -f odgi_avx ]];
then
  echo \"avx\"
  \"\$DIR\"/odgi_avx \"\$@\"
elif [[ -f odgi_sse4_2 ]];
then
  echo \"sse4_2\"
  \"\$DIR\"/odgi_sse4_2 \"\$@\"
else
  echo \"sse2\"
  \"\$DIR\"/odgi_sse2 \"\$@\"
fi"
echo "$ODGI" > release/odgi
chmod +x release/odgi

