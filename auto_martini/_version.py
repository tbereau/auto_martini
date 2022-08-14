try:
    from importlib import metadata
except ImportError:  # pragma: no cover
    # For Py37 and below, use the import_metadata backport
    import importlib_metadata as metadata

__version__ = metadata.version("auto_martini")
