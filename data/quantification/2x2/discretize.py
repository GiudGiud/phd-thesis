import openmoc


def discretize_geometry(self):
    """Discretize a 2x2 checkerboard of BEAVRS fuel assemblies into
    OpenMOC flat source regions.

     This method is intended to be attached directly to a Batchwise instance,
     as illustrated by the example below for a PerfectBatchwise object:

     > import types
     > import infermc
     > b = infermc.batchwise.PerfectBatchwise()
     > b._discretize_geometry = types.MethodType(discretize_geometry, b)

     This will allow the Batchwise instance to discretize the geometry for
     accurate OpenMOC simulations.

     """

    openmoc.log.py_printf('INFO', 'Discretizing the geometry...')

    openmc_geometry = self.mat_mgxslibs[0].openmc_geometry
    opencg_geometry = self.mat_mgxslibs[0].opencg_geometry
    all_cells = self.openmoc_geometry.getAllMaterialCells()
    all_openmoc_cells = self.openmoc_geometry.getAllMaterialCells()

    # Add angular sectors to all material-filled cells
    for cell_id in all_cells:
        all_cells[cell_id].setNumSectors(8)

    # Find the fuel clad outer radius zcylinder
    all_surfs = self.openmoc_geometry.getAllSurfaces()
    for surf_id in all_surfs:
        if all_surfs[surf_id].getName() == 'Fuel clad OR':
            fuel_or = openmoc.castSurfaceToZCylinder(all_surfs[surf_id])

    instr_tube_name = 'Instrument tube thimble radial 0: air'
    guide_tube_name = 'Empty GT above the dashpot radial 0: water'
    burn_abs1_name = 'BPRA rod active poison radial 0: air'
    burn_abs2_name = 'BPRA rod active poison radial 3: borosilicate'
    instr_guide_bp_tube_mod_name = 'Intermediate grid pincell radial 0: water'

    instr_tube = openmc_geometry.get_cells_by_name(instr_tube_name)
    guide_tube = openmc_geometry.get_cells_by_name(guide_tube_name)
    burn_abs1 = openmc_geometry.get_cells_by_name(burn_abs1_name)
    burn_abs2 = openmc_geometry.get_cells_by_name(burn_abs2_name)
    mod = openmc_geometry.get_cells_by_name(instr_guide_bp_tube_mod_name)

    for cell in instr_tube:
        print(cell.name, 'instr tube')
        all_openmoc_cells[cell.id].setNumRings(10)
    for cell in guide_tube:
        print(cell.name, 'guide tube')
        all_openmoc_cells[cell.id].setNumRings(10)
    for cell in burn_abs1:
        print(cell.name, 'burn abs1')
        all_openmoc_cells[cell.id].setNumRings(5)
    for cell in burn_abs2:
        print(cell.name, 'burn abs2')
        all_openmoc_cells[cell.id].setNumRings(5)
    for cell in mod:
        print(cell.name, 'mod')
        all_openmoc_cells[cell.id].setNumRings(5)
        all_openmoc_cells[cell.id].addSurface(surface=fuel_or, halfspace=+1)

    '''
    # Get the bounding box from the OpenCG geometry
    min_x = opencg_geometry.min_x
    min_y = opencg_geometry.min_y
    min_z = opencg_geometry.min_z
    max_x = opencg_geometry.max_x
    max_y = opencg_geometry.max_y
    max_z = opencg_geometry.max_z
    mid_x = (max_x + min_x) / 2.
    mid_y = (max_y + min_y) / 2.
    mid_z = (max_z + min_z) / 2.

    ###########################################################################
    # Discretize the instr. and guide tubes in assemblies along diagonal
    ###########################################################################

    # Get the middle of the top left fuel assembly
    assm_x = mid_x - 1.25984 * 8.5
    assm_y = mid_y + 1.25984 * 8.5

    # Discretize the guide and instrument tubes
    instr_tube = opencg_geometry.find_cell(x=assm_x, y=assm_y, z=mid_z)
    guide_tube = opencg_geometry.find_cell(x=assm_x, y=assm_y+3.77952, z=mid_z)
    burn_abs1 = opencg_geometry.find_cell(x=assm_x, y=assm_y+7.55904, z=mid_z)
    burn_abs2 = opencg_geometry.find_cell(x=assm_x, y=assm_y+7.55904+0.3, z=mid_z)
    instr_tube = all_cells[instr_tube.id]
    guide_tube = all_cells[guide_tube.id]
    burn_abs1 = all_cells[burn_abs1.id]
    burn_abs2 = all_cells[burn_abs2.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)
    burn_abs1.setNumRings(5)
    burn_abs2.setNumRings(5)

    print('instr tube: {} {}'.format(instr_tube.getName(), instr_tube.getId()))
    print('burn abs1: {} {}'.format(burn_abs1.getName(), burn_abs1.getId()))
    print('burn abs2: {} {}'.format(burn_abs2.getName(), burn_abs2.getId()))
    print('guide tube: {} {}'.format(guide_tube.getName(), guide_tube.getId()))

    # Discretize the guide and instrument tubes
    instr_tube_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+0.6, z=mid_z)
    guide_tube_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+3.77952+0.6, z=mid_z)
    burn_abs_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+7.55904+0.6, z=mid_z)
    instr_tube_moderator = all_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_cells[guide_tube_moderator.id]
    burn_abs_moderator = all_cells[burn_abs_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)
    burn_abs_moderator.setNumRings(5)

    print('instr tube mod: {} {}'.format(instr_tube_moderator.getName(), instr_tube.getId()))
    print('guide tube mod: {} {}'.format(guide_tube_moderator.getName(), guide_tube_moderator.getId()))
    print('burn abs mod: {} {}'.format(burn_abs_moderator.getName(), burn_abs_moderator.getId()))
    print('')

    ###########################################################################
    # Discretize the instr., guide tubes and BAs in assemblies off diagonal
    ###########################################################################

    # Get the middle of the top right fuel assembly
    assm_x = mid_x + 1.25984 * 8.5
    assm_y = mid_y + 1.25984 * 8.5

    # Discretize the guide and instrument tubes
    instr_tube = opencg_geometry.find_cell(x=assm_x, y=assm_y, z=mid_z)
    guide_tube = opencg_geometry.find_cell(x=assm_x, y=assm_y+3.77952, z=mid_z)
    burn_abs1 = opencg_geometry.find_cell(x=assm_x, y=assm_y+7.55904, z=mid_z)
    burn_abs2 = opencg_geometry.find_cell(x=assm_x, y=assm_y+7.55904+0.3, z=mid_z)
    instr_tube = all_cells[instr_tube.id]
    guide_tube = all_cells[guide_tube.id]
    burn_abs1 = all_cells[burn_abs1.id]
    burn_abs2 = all_cells[burn_abs2.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)
    burn_abs1.setNumRings(5)
    burn_abs2.setNumRings(5)

    print('instr tube: {} {}'.format(instr_tube.getName(), instr_tube.getId()))
    print('burn abs1: {} {}'.format(burn_abs1.getName(), burn_abs1.getId()))
    print('burn abs2: {} {}'.format(burn_abs2.getName(), burn_abs2.getId()))
    print('guide tube: {} {}'.format(guide_tube.getName(), guide_tube.getId()))

    # Discretize the guide and instrument tubes and burnable absorber moderator
    instr_tube_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+0.6, z=mid_z)
    guide_tube_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+3.77952+0.6, z=mid_z)
    burn_abs_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+7.55904+0.6, z=mid_z)
    instr_tube_moderator = all_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_cells[guide_tube_moderator.id]
    burn_abs_moderator = all_cells[burn_abs_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)
    burn_abs_moderator.setNumRings(5)

    print('instr tube mod: {} {}'.format(instr_tube_moderator.getName(), instr_tube.getId()))
    print('guide tube mod: {} {}'.format(guide_tube_moderator.getName(), guide_tube_moderator.getId()))
    print('burn abs mod: {} {}'.format(burn_abs_moderator.getName(), burn_abs_moderator.getId()))
    print('')

    ###########################################################################
    # Discretize the moderator cells around the fuel pins in both assemblies
    ###########################################################################

    # Find the fuel clad outer radius zcylinder
    all_surfs = self.openmoc_geometry.getAllSurfaces()
    for surf_id in all_surfs:
        if all_surfs[surf_id].getName() == 'Fuel clad OR':
            fuel_or = openmoc.castSurfaceToZCylinder(all_surfs[surf_id])

    # Discretize the moderator cell. First, add the fuel clad outer radius to each
    # cell since this is not done by the BEAVRS builder but is needed for ringify
    top_left = opencg_geometry.find_cell(x=min_x+0.1, y=max_y-0.1, z=mid_z)
    top_right = opencg_geometry.find_cell(x=min_x+1.25984*17, y=max_y-0.1, z=mid_z)
    top_left = all_cells[top_left.id]
    top_right = all_cells[top_right.id]
    top_left.addSurface(surface=fuel_or, halfspace=+1)
    top_right.addSurface(surface=fuel_or, halfspace=+1)
    top_left.setNumRings(10)
    top_right.setNumRings(10)

    print('top left: {} {}'.format(top_left.getName(), top_left.getId()))
    print('top right: {} {}'.format(top_left.getName(), top_left.getId()))
    print('')
    '''

    ###########################################################################
    # Discretize the fuel pin cells in both assemblies
    ###########################################################################

    # Find all pin cell universes - universes containing a cell filled with fuel
    all_univs = self.openmoc_geometry.getAllUniverses()
    pin_univs = set()
    for univ_id in all_univs:
        all_cells = all_univs[univ_id].getAllCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    pin_univs.add(all_univs[univ_id])

    # Discretize all cells within each pin cell universe into rings and sectors
    for pin_univ in pin_univs:
        all_cells = pin_univ.getCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    all_cells[cell_id].setNumRings(5)


def discretize_geometry_standalone(mat_mgxslib, openmoc_geometry):
    """Discretize a 2x2 checkerboard of BEAVRS fuel assemblies into
    OpenMOC flat source regions."""

    openmoc.log.py_printf('INFO', 'Discretizing the geometry...')

    opencg_geometry = mat_mgxslib.opencg_geometry
    all_cells = openmoc_geometry.getAllMaterialCells()

    # Add angular sectors to all material-filled cells
    for cell_id in all_cells:
        all_cells[cell_id].setNumSectors(8)

    # Get the bounding box from the OpenCG geometry
    min_x = opencg_geometry.min_x
    min_y = opencg_geometry.min_y
    min_z = opencg_geometry.min_z
    max_x = opencg_geometry.max_x
    max_y = opencg_geometry.max_y
    max_z = opencg_geometry.max_z
    mid_x = (max_x + min_x) / 2.
    mid_y = (max_y + min_y) / 2.
    mid_z = (max_z + min_z) / 2.

    ###########################################################################
    # Discretize the instr. and guide tubes in assemblies along diagonal
    ###########################################################################

    # Get the middle of the top left fuel assembly
    assm_x = mid_x - 1.25984 * 8.5
    assm_y = mid_y + 1.25984 * 8.5

    # Discretize the guide and instrument tubes
    instr_tube = opencg_geometry.find_cell(x=assm_x, y=assm_y, z=mid_z)
    guide_tube = opencg_geometry.find_cell(x=assm_x, y=assm_y+3.77952, z=mid_z)
    burn_abs1 = opencg_geometry.find_cell(x=assm_x, y=assm_y+7.55904, z=mid_z)
    burn_abs2 = opencg_geometry.find_cell(x=assm_x, y=assm_y+7.55904+0.3, z=mid_z)
    instr_tube = all_cells[instr_tube.id]
    guide_tube = all_cells[guide_tube.id]
    burn_abs1 = all_cells[burn_abs1.id]
    burn_abs2 = all_cells[burn_abs2.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)
    burn_abs1.setNumRings(5)
    burn_abs2.setNumRings(5)

    # Discretize the guide and instrument tubes
    instr_tube_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+0.6, z=mid_z)
    guide_tube_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+3.77952+0.6, z=mid_z)
    burn_abs_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+7.55904+0.6, z=mid_z)
    instr_tube_moderator = all_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_cells[guide_tube_moderator.id]
    burn_abs_moderator = all_cells[burn_abs_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)
    burn_abs_moderator.setNumRings(5)

    ###########################################################################
    # Discretize the instr., guide tubes and BAs in assemblies off diagonal
    ###########################################################################

    # Get the middle of the top right fuel assembly
    assm_x = mid_x + 1.25984 * 8.5
    assm_y = mid_y + 1.25984 * 8.5

    # Discretize the guide and instrument tubes
    instr_tube = opencg_geometry.find_cell(x=assm_x, y=assm_y, z=mid_z)
    guide_tube = opencg_geometry.find_cell(x=assm_x, y=assm_y+3.77952, z=mid_z)
    burn_abs1 = opencg_geometry.find_cell(x=assm_x, y=assm_y+7.55904, z=mid_z)
    burn_abs2 = opencg_geometry.find_cell(x=assm_x, y=assm_y+7.55904+0.3, z=mid_z)
    instr_tube = all_cells[instr_tube.id]
    guide_tube = all_cells[guide_tube.id]
    burn_abs1 = all_cells[burn_abs1.id]
    burn_abs2 = all_cells[burn_abs2.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)
    burn_abs1.setNumRings(5)
    burn_abs2.setNumRings(5)

    # Discretize the guide and instrument tubes and burnable absorber moderator
    instr_tube_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+0.6, z=mid_z)
    guide_tube_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+3.77952+0.6, z=mid_z)
    burn_abs_moderator = opencg_geometry.find_cell(x=assm_x+0.6, y=assm_y+7.55904+0.6, z=mid_z)
    instr_tube_moderator = all_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_cells[guide_tube_moderator.id]
    burn_abs_moderator = all_cells[burn_abs_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)
    burn_abs_moderator.setNumRings(5)

    ###########################################################################
    # Discretize the moderator cells around the fuel pins in both assemblies
    ###########################################################################

    # Find the fuel clad outer radius zcylinder
    all_surfs = openmoc_geometry.getAllSurfaces()
    for surf_id in all_surfs:
        if all_surfs[surf_id].getName() == 'Fuel clad OR':
            fuel_or = openmoc.castSurfaceToZCylinder(all_surfs[surf_id])

    # Discretize the moderator cell. First, add the fuel clad outer radius to each
    # cell since this is not done by the BEAVRS builder but is needed for ringify
    top_left = opencg_geometry.find_cell(x=min_x+0.1, y=max_y-0.1, z=mid_z)
    top_right = opencg_geometry.find_cell(x=min_x+1.25984*17, y=max_y-0.1, z=mid_z)
    top_left = all_cells[top_left.id]
    top_right = all_cells[top_right.id]
    top_left.addSurface(surface=fuel_or, halfspace=+1)
    top_right.addSurface(surface=fuel_or, halfspace=+1)
    top_left.setNumRings(10)
    top_right.setNumRings(10)

    ###########################################################################
    # Discretize the fuel pin cells in both assemblies
    ###########################################################################

    # Find all pin cell universes - universes containing a cell filled with fuel
    all_univs = openmoc_geometry.getAllUniverses()
    pin_univs = set()
    for univ_id in all_univs:
        all_cells = all_univs[univ_id].getAllCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    pin_univs.add(all_univs[univ_id])

    # Discretize all cells within each pin cell universe into rings and sectors
    for pin_univ in pin_univs:
        all_cells = pin_univ.getCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    all_cells[cell_id].setNumRings(5)