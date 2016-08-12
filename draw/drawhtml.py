#!/usr/bin/env python
# coding: utf-8

import os
from sage.all import Polyhedron, GF
from betti.all import *
from jinja2 import Template


html = Template(
r"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{{ title }}</title>
<style>
table, td, th {border: 0pt none; border-collapse: collapse;}
td, th {text-align: right;}
td {padding-left: 1ex; padding-right: 1ex;}
th.x {padding-left: 1ex; padding-right: 1ex; padding-top: 1ex;}
th.y {padding-right: 2ex;}
.nontriv {background-color: #cfc;}
.triv {background-color: #ffc;}
</style>
</head>
<body>
<h1>{{ typ }}<sub>{{ l }}</sub>({{ polygon }}) = {{ result }}</h1>
{{ table }}
</body>
</html>
"""
)

index_html = Template(
r"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Graded Betti tables of toric surfaces</title>
<style>
table {border-collapse: collapse;}
thead {border-bottom: 1px solid;}
col.head {border-right: 1px solid;}
th, td {text-align: center; padding-left: 1ex; padding-right: 1ex;}
a {text-decoration: none;}
</style>
</head>
<body>
<h1>Graded Betti tables of toric surfaces with bidegrees</h1>
{{ index }}
</body>
</html>
"""
)


def makedirs(dir):
    try:
        os.makedirs(dir)
    except OSError:
        if not os.path.isdir(dir):
            raise


def html_full_betti_table(Delta, htmlname, dir):
    makedirs(dir)

    N = len(Delta.integral_points())
    T = ToricSurfaceBettiTable(Delta, GF(40009), verbose=1, threads=2)

    # Betti table pages
    for l in range(1, N-2):
        for typ in ["b", "c"]:
            name = "{}{}".format(typ, l)
            title = "{}/{}".format(dir, name)
            out = os.path.join(dir, name + ".html")
            print(out)
            if typ == "b":
                result = T.betti_number(l, 1)
            else:
                result = T.betti_number(N - 2 - l, 2)
            table = format_table(typ, T, l, format="html", checksum=result)
            with open(out, 'w') as f:
                f.write(html.render(title=title, typ=typ, l=l, polygon=htmlname,
                    result=result, table=table))

    # Add to index page
    global index
    index += "\n<h2>Graded Betti table of {}</h2>\n<table>\n".format(htmlname)
    index += '<colgroup><col class="head"><col span="{}"></colgroup>\n'.format(N-2)
    index += "<thead>\n<tr><th></th>"
    for p in range(N-2):
        index += "<th>{}</th>".format(p)
    index += "</tr>\n</thead>\n<tbody>\n"
    for q in range(3):
        index += "<tr><th>{}</th>".format(q)
        for p in range(N-2):
            result = T.betti_number(p, q)
            if p > 0 and q == 1:
                result = '<a href="{}/b{}.html">{}</a>'.format(dir, p, result)
            elif p > 0 and q == 2:
                result = '<a href="{}/c{}.html">{}</a>'.format(dir, N - 2 - p, result)
            index += "<td>{}</td>".format(result)
        index += "</tr>\n"
    index += "</tbody>\n</table>\n"


def Sigma(d):
    return Polyhedron([(0,0), (d,0), (0,d)])

def Ypsilon(d):
    return Polyhedron([(0,0), (2*d, d), (d, 2*d)])

def Ypsilon_(d):
    return Polyhedron([(0,0), (d+1,1), (1,d+1)])


makedirs("betti")
os.chdir("betti")

index = ""
html_full_betti_table(Sigma(1), "&Sigma;", ".")
html_full_betti_table(Sigma(2), "2&thinsp;&Sigma;", "2sigma")
html_full_betti_table(Sigma(3), "3&thinsp;&Sigma;", "3sigma")
html_full_betti_table(Sigma(4), "4&thinsp;&Sigma;", "4sigma")
html_full_betti_table(Sigma(5), "5&thinsp;&Sigma;", "5sigma")
html_full_betti_table(Ypsilon(1), "&Upsilon;", "ypsilon")
html_full_betti_table(Ypsilon(2), "2&thinsp;&Upsilon;", "2ypsilon")
html_full_betti_table(Ypsilon_(2), "&Upsilon;<sub>2</sub>", "ypsilon2")
html_full_betti_table(Ypsilon_(3), "&Upsilon;<sub>3</sub>", "ypsilon3")
html_full_betti_table(Ypsilon_(4), "&Upsilon;<sub>4</sub>", "ypsilon4")
#html_full_betti_table(Ypsilon_(5), "&Upsilon;<sub>5</sub>", "ypsilon5")

with open("index.html", "w") as f:
    f.write(index_html.render(index=index))
