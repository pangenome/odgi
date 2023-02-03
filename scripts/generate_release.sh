#!/bin/bash

SSE2=$(lscpu | grep "Flags" | grep "sse2")
SSE4_2=$(lscpu | grep "Flags" | grep "sse4_2")
AVX=$(lscpu | grep "Flags" | grep "avx")
AVX2=$(lscpu | grep "Flags" | grep "avx2")
AVX512=$(lscpu | grep "Flags" | grep "avx512")

# sed CMakeLists.txt 's//g'
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
  AVX512_FLAGS=$(lscpu | grep "avx512" | sed 's/ /\n/g' | grep "avx512" | sed 's/avx/-mavx/g' | tr '\n' ' ')
  cmake -H. -Bavx512 -DEXTRA_FLAGS="-Ofast -pipe $AVX512_FLAGS" && cmake --build avx512 -- -j 15
  mv bin/odgi avx512
  rm -r bin
  mv avx512/odgi release/odgi_avx512
fi

sed -i '55 s/^#//' CMakeLists.txt

ODGI="#!/bin/bash
if [[ -f odgi_avx512 ]];
then
  ./odgi_avx512
elif [[ -f odgi_avx2 ]];
then
# debugging
#  echo "SLJLEJLEJFLSJFLÖ"
  ./odgi_avx2
elif [[ -f odgi_avx ]];
then
  ./odgi_avx
elif [[ -f odgi_sse4_2 ]];
then
  ./odgi_sse4_2
else
  ./odgi_sse2
fi"
echo "$ODGI" > release/odgi
chmod +x release/odgi

