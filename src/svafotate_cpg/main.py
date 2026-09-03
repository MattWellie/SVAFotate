"""Backwards-compatible entry point; the implementation lives in :mod:`svafotate_cpg.cli`."""

from .cli import get_parser, main

__all__ = ['get_parser', 'main']
