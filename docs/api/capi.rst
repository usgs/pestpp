C API
=====

The C ABI exposed by the shared library (``pestpp-api.dylib`` / ``.so`` / ``.dll``). This is the
interface used to drive a PEST++ tool from another process rather than from a control file on
the command line - see :doc:`../users_guide` for what the tools themselves do.

Every entry point is prefixed ``pestpp_`` and returns a ``pestpp_status``. Nothing in the ABI
transfers ownership of memory across the boundary: no function returns an allocation the caller
is expected to free, and no ``FILE*`` or file descriptor appears in a signature. The one
``double**`` is a borrowed view of the library's own buffer and is handed back with
``pestpp_release_view``.

Header
------

The header carries the contract - status codes, buffer sizing conventions, and the ownership
notes above - so it is the place to start.

.. doxygenfile:: pestpp-api.h
   :project: pestpp

Only the header is documented here. Every entry point is declared there and defined in
pestpp_capi.cpp, so documenting both makes breathe emit each function twice and sphinx reports
every one as a duplicate declaration.
