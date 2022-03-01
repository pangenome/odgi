;; To use this file to build HEAD of odgi:
;;
;;   guix build -f guix.scm
;;
;; To get a development container (emacs shell will work)
;;
;;   guix shell -C -D -f guix.scm
;;
;; and build
;;
;;   find -name CMakeCache.txt|xargs rm -v
;;   cd build
;;   cmake -DCMAKE_BUILD_TYPE=Debug ..
;;   cmake --build . --verbose
;;
;; For the tests you may need /usr/bin/env. In a container create it with
;;
;;   mkdir -p /usr/bin ; ln -s $GUIX_ENVIRONMENT/bin/env /usr/bin/env
;;
;; Note for python bindings you may need to run against gcc-11 with something
;; like
;;
;;   env PYTHONPATH=lib python3 -c 'import odgi'
;;   env PYTHONPATH=lib python3 examples/explore.py
;;

(use-modules
  (ice-9 popen)
  (ice-9 rdelim)
  ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix download)
  (guix git-download)
  (guix build-system cmake)
  (guix utils)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages commencement) ; gcc-toolchain
  (gnu packages curl)
  (gnu packages gdb)
  (gnu packages gcc)
  (gnu packages jemalloc)
  (gnu packages libffi)
  (gnu packages python)
  (gnu packages python-xyz)
  (gnu packages pkg-config)
  (gnu packages tls)
  (gnu packages version-control)
)

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
       ; ("cpp-httplib" ,cpp-httplib) later!
       ("pybind11" ,pybind11) ;; see libstd++ note in remarks above
       ; ("intervaltree" ,intervaltree) later!
       ("jemalloc" ,jemalloc)
       ("gcc" ,gcc-11)
       ("gcc-lib" ,gcc-11 "lib")
       ("gcc-toolchain" ,gcc-toolchain)
       ("gdb" ,gdb)
       ("git" ,git)
       ; ("lodepng" ,lodepng) later!
       ("python" ,python)
       ; ("sdsl-lite" ,sdsl-lite) later!
       ))
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
     (synopsis "odgi pangenome optimized dynamic sequence graph implementation")
     (description
"odgi pangenome graph tooling provides an efficient, succinct dynamic
DNA sequence graph model, as well as a host of algorithms that allow
the use of such graphs.")
     (home-page "https://github.com/vgteam/odgi")
     (license license:expat)))

odgi-git
