from tables import *
import numpy as np

h5file = open_file("tutorial1.h5", 'a')
table = h5file.root.detector.readout
pressure = [ x['pressure'] for x in table.iterrows()]
print(pressure)

#print(h5file)

for nodes in h5file:
    print(nodes)

for group in h5file.walk_groups("/"):
    for array in h5file.list_nodes(group, classname='Array'):
        print(array)
