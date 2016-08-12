import openmoc


def discretize_geometry(self):
    """Discretize the BEAVRS assembly with 1.6% enriched fuel pins and guide tubes
     into OpenMOC flat source regions.

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

    # Find the fuel clad outer radius zcylinder
    all_surfs = self.openmoc_geometry.getAllSurfaces()
    for surf_id in all_surfs:
        if all_surfs[surf_id].getName() == 'Fuel clad OR':
            fuel_or = openmoc.castSurfaceToZCylinder(all_surfs[surf_id])

    # Discretize the moderator cell. First, add the fuel clad outer radius to each
    # cell since this is not done by the BEAVRS builder but is needed for ringify
    bottom_left = opencg_geometry.find_cell(x=min_x+0.1, y=min_y+0.1, z=mid_z)
    bottom_left = all_cells[bottom_left.id]
    bottom_left.addSurface(surface=fuel_or, halfspace=+1)
    bottom_left.setNumRings(10)

    # Discretize the guide and instrument tubes
    instr_tube = opencg_geometry.find_cell(x=mid_x, y=mid_y, z=mid_z)
    guide_tube = opencg_geometry.find_cell(x=mid_x, y=mid_y+3.78, z=mid_z)
    instr_tube = all_cells[instr_tube.id]
    guide_tube = all_cells[guide_tube.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)

    # Discretize the moderator around the guide and instrument tubes
    instr_tube_moderator = opencg_geometry.find_cell(x=mid_x+0.6, y=mid_y+0.6, z=mid_z)
    guide_tube_moderator = opencg_geometry.find_cell(x=mid_x+0.6, y=mid_y+3.78+0.6, z=mid_z)
    instr_tube_moderator = all_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_cells[guide_tube_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)

    # Find all pin cell universes - universes containing a cell filled with fuel
    all_univs = self.openmoc_geometry.getAllUniverses()
    pin_univs = set()
    for univ_id in all_univs:
        all_cells = all_univs[univ_id].getAllCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    pin_univs.add(all_univs[univ_id])

    # Discretize all fuel cells within each pin cell universe into rings
    for pin_univ in pin_univs:
        all_cells = pin_univ.getCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    all_cells[cell_id].setNumRings(5)


def discretize_geometry_standalone(mat_mgxslib, openmoc_geometry):
    """Discretize the BEAVRS assembly with 1.6% enriched fuel pins and guide tubes
     into OpenMOC flat source regions."""

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

    # Find the fuel clad outer radius zcylinder
    all_surfs = openmoc_geometry.getAllSurfaces()
    for surf_id in all_surfs:
        if all_surfs[surf_id].getName() == 'Fuel clad OR':
            fuel_or = openmoc.castSurfaceToZCylinder(all_surfs[surf_id])

    # Discretize the moderator cell. First, add the fuel clad outer radius to each
    # cell since this is not done by the BEAVRS builder but is needed for ringify
    bottom_left = opencg_geometry.find_cell(x=min_x+0.1, y=min_y+0.1, z=mid_z)
    bottom_left = all_cells[bottom_left.id]
    bottom_left.addSurface(surface=fuel_or, halfspace=+1)
    bottom_left.setNumRings(10)

    # Discretize the guide and instrument tubes
    instr_tube = opencg_geometry.find_cell(x=mid_x, y=mid_y, z=mid_z)
    guide_tube = opencg_geometry.find_cell(x=mid_x, y=mid_y+3.78, z=mid_z)
    instr_tube = all_cells[instr_tube.id]
    guide_tube = all_cells[guide_tube.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)

    # Discretize the moderator around the guide and instrument tubes
    instr_tube_moderator = opencg_geometry.find_cell(x=mid_x+0.6, y=mid_y+0.6, z=mid_z)
    guide_tube_moderator = opencg_geometry.find_cell(x=mid_x+0.6, y=mid_y+3.78+0.6, z=mid_z)
    instr_tube_moderator = all_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_cells[guide_tube_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)

    # Find all pin cell universes - universes containing a cell filled with fuel
    all_univs = openmoc_geometry.getAllUniverses()
    pin_univs = set()
    for univ_id in all_univs:
        all_cells = all_univs[univ_id].getAllCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    pin_univs.add(all_univs[univ_id])

    # Discretize all fuel cells within each pin cell universe into rings
    for pin_univ in pin_univs:
        all_cells = pin_univ.getCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    all_cells[cell_id].setNumRings(5)