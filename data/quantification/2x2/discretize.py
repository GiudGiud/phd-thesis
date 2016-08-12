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

    opencg_geometry = self.mat_mgxslibs[0].opencg_geometry
    all_cells = self.openmoc_geometry.getAllMaterialCells()

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