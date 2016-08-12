#!/usr/bin/env python

import os
from subprocess import check_call
from sage.all import Polyhedron, GF
from betti.all import *

d = 4
Delta = Polyhedron([(0,0), (0,d), (d,0)])
out = "coimage{}sigma.tex".format(d)

Delta_picture = "\\tikz{\\fill"
for pt in Delta.integral_points():
    Delta_picture += " (%s,%s) circle[radius=3mm]" % tuple(pt)
Delta_picture += ";}"


T = ToricSurfaceBettiTable(Delta, GF(40009), verbose=1, threads=2, remove_points=[])

f = open(out, 'w')

f.write(
r"""\documentclass[a4paper,landscape,10pt]{article}
\usepackage[margin=15mm]{geometry}
\usepackage{tikz}
\pagestyle{empty}

\begin{document}
\small
""")

for l in range(1,7):
    c = T.c(l)
    s = format_table_c_coimage(T, l, format="tikz", options="xscale=0.75, yscale=0.5, baseline=(current bounding box.north)")
    f.write("\n\n")
    f.write("\\makebox[0mm]{\\Huge\\hspace{150mm}\\raisebox{-1.5em}{$c_%s(\\resizebox{1.6ex}{!}{%s}) = %s$}}" % (l, Delta_picture, c))
    f.write(s)
    f.write("\\clearpage\n")

f.write(r"""
\end{document}
""")

f.close()

check_call(["pdflatex", out], stdin=open(os.devnull))
