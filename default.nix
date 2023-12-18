{ pkgs ? import <nixpkgs> {} }:

pkgs.callPackage ./odgi.nix {
  inherit (pkgs) stdenv fetchFromGitHub cmake jemalloc pkg-config python3 gcc zlib htslib gsl;
}
