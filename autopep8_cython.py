# Small program to try to autoformat cython files

import glob
from subprocess import run

file_list = glob.glob("./**/*.pyx", recursive=True)
file_list = file_list + glob.glob("./**/*.pxd", recursive=True)
fixlist = [
    "E114",
    "E127",
    "E201",
    "E202",
    "E211",
    "E221",
    "E228",
    "E231",
    "E251",
    "E261",
    "E262",
    "E265",
    "E271",
    "E272",
    "E301",
    "E302",
    "E303",
    "E304",
    "E701",
    "E702",
    "E703",
    "E711",
    "W291",
    "W293",
    "W391",
]

select = []
for fix in fixlist:
    for filename in file_list:
        command = ["autopep8", filename, "-i", "-a", "--select", fix]
        run(command)
