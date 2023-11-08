{ pkgs ? import <nixpkgs> {} }:

pkgs.callPackage ./odgi.nix {
  inherit (pkgs) stdenv fetchFromGitHub cmake jemalloc pkgconfig python3 gcc zlib htslib gsl;
}
