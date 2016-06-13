from infermc.energy_groups import group_structures

num_groups = int(input('# groups: '))

groups = group_structures['CASMO']['{}-group'.format(num_groups)]

print('# groups: {}'.format(num_groups))

for i in range(num_groups):
    left = groups.group_edges[i]
    right = groups.group_edges[i+1]
    print('{} & {:0.4E} & {:0.4E} \\\\'.format(num_groups-i, left, right))

