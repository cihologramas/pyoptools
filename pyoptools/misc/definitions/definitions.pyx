# Definitions that help to speed up some processces

from numpy import inf, array
# Infinity vector is used in many places to show no intersection, to
# avoid runtime creation of many infiniti vectors, inf_vect should be used

inf_vect=array((inf, inf, inf))
