# Computing graded Betti tables of toric surfaces

This is a preliminary version of the software developed as part of the paper
*Computing graded Betti tables of toric surfaces*
by Wouter Castryck, Filip Cools, Jeroen Demeyer and Alexander Lemmens.

## Installation

This package requires a recent version of Sage and
the recommended way to install this package is within Sage.
Download the sources and then run

    sage --python setup.py install

## Usage

Start Sage the usual way and then execute the following commands
to compute the graded Betti table of 5 Sigma in just a few seconds:

    sage: from betti.all import *
    sage: P = polygon_normalized([(0,0), (5,0), (0,5)])
    sage: T = ToricSurfaceBettiTable(P, GF(40009))
    sage: T.betti_table()
    [     1      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0]
    [     0    165   1830  10710  41616 117300 250920 417690 548080 568854 464100 291720 134640  39780   4858    375      0      0      0]
    [     0      0      0      0      0      0      0      0      0      0      0      0      0   2002   4200   2160    595     90      6]

Useful keywords are `verbose` and `threads`:

    sage: T = ToricSurfaceBettiTable(P, GF(40009), verbose=2, threads=2)
    sage: T.betti_table()
    ...
    [     1      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0]
    [     0    165   1830  10710  41616 117300 250920 417690 548080 568854 464100 291720 134640  39780   4858    375      0      0      0]
    [     0      0      0      0      0      0      0      0      0      0      0      0      0   2002   4200   2160    595     90      6]
