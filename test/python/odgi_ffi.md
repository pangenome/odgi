
Testing doctest

env LD_LIBRARY_PATH=$GUIX_ENVIRONMENT/lib:/odgi/lib:./libbf-prefix/lib/ LD_PRELOAD=libjemalloc.so.2:libsdsl.so:libbf.so PYTHONPATH=/odgi/lib/ python3 -m doctest /odgi/test/python/odgi_ffi.md -v

env LD_LIBRARY_PATH=$GUIX_ENVIRONMENT/lib ctest . --verbose -R ffi

-->

```python
>>> import odgi_ffi
>>> odgi_ffi.odgi_version()
'f937271'

>>> 1+2
3

```

```python
>>> 1+3
4

```
