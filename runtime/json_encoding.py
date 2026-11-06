"""Encode artifact JSON with CPython's complete-value encoder.

Writers retain control of newlines, flush/fsync and atomic replacement.
"""
import json


def write_json_encoded(value, handle, **kwargs):
    """Preserve ``json.dump`` output while issuing one text write."""
    handle.write(json.dumps(value, **kwargs))
