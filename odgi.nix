{ lib, stdenv, fetchFromGitHub, cmake, jemalloc, python3, pkg-config, gcc, git, zlib, htslib, gsl }:

stdenv.mkDerivation rec {
  pname = "odgi";
  version = "0.8.3";

  src = fetchFromGitHub {
    owner = "pangenome";
    repo = "odgi";
    rev = "86e62bac42d737808a83bb00a0e7a60069494dcf";
    sha256 = "sha256-IKHQyP3E01LZvui6ykUiWPr2zuD2iqq55ARG+O2KCxM=";
    fetchSubmodules = true;
  };

  nativeBuildInputs = [ cmake pkg-config ];
  
  buildInputs = [
    jemalloc
    gcc
    zlib
    htslib
    gsl
    python3
  ];

  postPatch = ''
    mkdir -p include
    echo "#define ODGI_GIT_VERSION \"${version}\"" > include/odgi_git_version.hpp
  '';

  makeFlags = [ "CC=${gcc}/bin/gcc" ];

  meta = with lib; {
    description = "odgi optimized dynamic sequence graph implementation";
    homepage = "https://github.com/pangenome/odgi";
    license = licenses.mit;
    platforms = platforms.linux;
    maintainers = [ maintainers.yourNameHere ];  # Replace with your name
  };
}
