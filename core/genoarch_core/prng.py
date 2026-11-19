"""Deterministic 32-bit PRNG shared bit-for-bit with ``web/js/engine.js``.

``fnv1a32`` hashes the canonical query string into a seed; ``Mulberry32`` is the
public-domain mulberry32 generator. Both use only 32-bit integer operations, so
the JavaScript port (``Math.imul`` / ``>>>``) yields identical sequences.
"""
from __future__ import annotations

_MASK = 0xFFFFFFFF


def _imul(a: int, b: int) -> int:
    return (a * b) & _MASK


def fnv1a32(text: str) -> int:
    h = 0x811C9DC5
    for byte in text.encode("utf-8"):
        h = _imul(h ^ byte, 0x01000193)
    return h


class Mulberry32:
    def __init__(self, seed: int) -> None:
        self.state = seed & _MASK

    def next(self) -> float:
        """Return a float in [0, 1)."""
        self.state = (self.state + 0x6D2B79F5) & _MASK
        t = self.state
        t = _imul(t ^ (t >> 15), t | 1)
        t = ((t + _imul(t ^ (t >> 7), t | 61)) & _MASK) ^ t
        return ((t ^ (t >> 14)) & _MASK) / 4294967296.0
