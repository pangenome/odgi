;; To use this file to build HEAD of odgi:
;;
;;   guix build -f guix.scm
;;
;; To get a development container (emacs shell will work)
;;
;;   guix shell -C -D -f guix.scm
;;
;; For the tests you may need /usr/bin/env. In a container create it with
;;
;;   mkdir -p /usr/bin ; ln -s $GUIX_ENVIRONMENT/bin/env /usr/bin/env
;;
;; Note for python bindings you may need to run against gcc-11 with something
;; like
;;
;;   env LD_LIBRARY_PATH=/gnu/store/*gcc-11*lib/lib PYTHONPATH=lib python3 examples/explore.py
;;
;; otherwise you get ImportError:
;; /gnu/store/90lbavffg0csrf208nw0ayj1bz5knl47-gcc-10.3.0-lib/lib/libstdc++.so.6:
;; version `GLIBCXX_3.4.29' not found because it tries to pick up from gcc-10.

(use-modules
  ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix git-download)
  (guix build-system cmake)
  (guix utils)
  ;; (gnu packages algebra)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages commencement) ; gcc-toolchain
  (gnu packages curl)
  (gnu packages gdb)
  (gnu packages gcc)
  (gnu packages jemalloc)
  ;; (gnu packages llvm)
  (gnu packages python)
  (gnu packages python-xyz)
  ;; (gnu packages parallel)
  ;; (gnu packages perl)
  ;; (gnu packages perl6)
  (gnu packages pkg-config)
  (gnu packages tls)
  (gnu packages version-control)
  (srfi srfi-1)
  (ice-9 popen)
  (ice-9 rdelim))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public odgi-git
  (package
    (name "odgi-git")
    (version (git-version "0.6.3" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (inputs
     `(
       ("coreutils" ,coreutils)
       ;; ("pybind11" ,pybind11) - prebuilt fails on gcc version 11
       ("jemalloc" ,jemalloc)
       ("gcc" ,gcc-11)
       ("gcc-toolchain" ,gcc-toolchain)
       ("gdb" ,gdb)
       ("git" ,git)
       ("python" ,python)))
    (native-inputs
     `(("pkg-config" ,pkg-config)
       ))
    (arguments
      `(#:phases
        (modify-phases
         %standard-phases
         ;; This stashes our build version in the executable
         (add-after 'unpack 'set-version
           (lambda _
             (mkdir-p "include")
             (with-output-to-file "include/odgi_git_version.hpp"
               (lambda ()
                 (format #t "#define ODGI_GIT_VERSION \"~a\"~%" version)))
             #t))
         (delete 'check))
        #:make-flags (list ,(string-append "CC=" (cc-for-target)))))
     (synopsis "odgi optimized dynamic sequence graph implementation")
     (description
"odgi provides an efficient, succinct dynamic DNA sequence graph model, as well
as a host of algorithms that allow the use of such graphs in bioinformatic
analyses.")
     (home-page "https://github.com/vgteam/odgi")
     (license license:expat)))

odgi-git
