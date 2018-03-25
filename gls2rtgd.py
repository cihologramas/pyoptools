#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from enthought.traits.api import Float, String
from os import access,  F_OK
"""
Utility to convert the Oslo Glass Data Files to the RayTrace Glass Data Format.
"""

from optparse import OptionParser

parser=OptionParser(usage="usage: %prog -i inputfile -o outputfile", 
                    description=\
                    "Utility to convert de Oslo Glass Data Files, to the RayTrace"\
                    "material library file")
 
parser.add_option("-i", dest="inputfile", help="GLS file to convert")
parser.add_option("-o",  dest="outputfile",  help="RTGL file to create")
(options, args) = parser.parse_args()

if options.inputfile==None:
     parser.error("\n\nError: Input filename must be given.\n\n")

print("\n\n\nProcessing file", options.inputfile, "\n\n\n")
try:
    f=open(options.inputfile, 'r')
except:
    print("File ", options.inputfile,  "not found")
    
lines=f.readlines()
f.close()


try:
    version, size, library=lines[0].split(None, 2)
    s=float(size)
except:
    exit("The input file seems not to be an OSLO glass description file")
    
if version not in ("5.00","5.11","5.31","5.32","5.41","6.02" , "6.03"):
    print("Library version migth not be supported. Please verify the output file")

# Check if the output file exists

if access(options.outputfile, F_OK):
    exit("The outputfile already exists. Use a different one.\n\n")
    
fo=open(options.outputfile, "w")

# Write the header of the library
# Library type RayTraceGlassLibrary Version 0
fo.write("RTGL-0\n")
# Description line always tab separated
# ref - Glass reference as given by the manufacturer 
# nd  - Refraction index for the d line
# vd  - Abbe's number for the d line
# B1,B2,B3,C1,C2,C3 - Dispersion constants for the glass

fo.write("ref\tnd\tvd\tB1\tB2\tB3\tC1\tC2\tC3\n")

for line in lines[1:]:
    
    line_data=line.split()
    
    # Skip blank lines
    if len(line_data)!=0:
        ndisp_coef=int(line_data[13])
        if ndisp_coef ==6:
            ref=line_data[0]
            nd=line_data[1]
            vd=line_data[2]
            B1=line_data[14]
            B2=line_data[15]
            B3=line_data[16]
            C1=line_data[17]
            C2=line_data[18]
            C3=line_data[19]
            fo.write("\"%s\"\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%\
                     (ref, nd, vd, B1, B2, B3, C1, C2, C3))
