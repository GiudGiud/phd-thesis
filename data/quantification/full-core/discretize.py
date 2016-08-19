import openmoc


def discretize_geometry(self):
    """Discretize a 2x2 checkerboard of BEAVRS fuel assemblies surrounded by
    a water reflector into OpenMOC flat source regions.

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
    all_openmoc_cells = self.openmoc_geometry.getAllMaterialCells()

    # FIXME: Only do this for pin cells!!!
    # Add angular sectors to all material-filled cells
    for cell_id in all_openmoc_cells:
        all_openmoc_cells[cell_id].setNumSectors(8)

    ###########################################################################
    # Discretize the instrument and guide tubes and BPs
    ###########################################################################

    # Find the fuel clad outer radius zcylinder
    all_surfs = self.openmoc_geometry.getAllSurfaces()
    for surf_id in all_surfs:
        if all_surfs[surf_id].getName() == 'Fuel clad OR':
            fuel_or = openmoc.castSurfaceToZCylinder(all_surfs[surf_id])

    # Find cells by their string names in the BEAVRS benchmark
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

    # Discretize each cell into radial rings
    for cell in instr_tube:
        all_openmoc_cells[cell.id].setNumRings(10)
    for cell in guide_tube:
        all_openmoc_cells[cell.id].setNumRings(10)
    for cell in burn_abs1:
        all_openmoc_cells[cell.id].setNumRings(5)
    for cell in burn_abs2:
        all_openmoc_cells[cell.id].setNumRings(5)
    for cell in mod:
        all_openmoc_cells[cell.id].setNumRings(5)
        all_openmoc_cells[cell.id].addSurface(surface=fuel_or, halfspace=+1)

    ###########################################################################
    # Discretize the fuel pin cells in all assemblies
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

    # Discretize all fuel cells within each pin cell universe into rings
    for pin_univ in pin_univs:
        all_cells = pin_univ.getCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    all_cells[cell_id].setNumRings(5)

    ###########################################################################
    # Discretize the reflector cells using a lattice
    ###########################################################################

    # Store the names of each of the reflector cells in the BEAVRS benchmark
    cell_names = ['Water', 'North Baffle Outer Water',
                  'Baffle South radial outer: water',
                  'Baffle East radial outer: water',
                  'Baffle West radial outer: water',
                  'Baffle North East Tip radial outer: water',
                  'Baffle South East Tip radial outer: water',
                  'Tip Baffle Outer Water W',
                  'Tip Baffle Outer Water N',
                  'Baffle South West Tip radial outer: water',
                  'Baffle South West Corner radial outer: water',
                  'Baffle South East Corner radial outer: water',
                  'Baffle North East Corner radial outer: water',
                  'Corner Baffle Water Gap N',
                  'Corner Baffle Water Gap W',
                  'Corner Baffle Outer Water']

    # Find all of the water-filled cells which comprise the reflector
    refl_cells = []
    for cell_name in cell_names:
        cells = openmc_geometry.get_cells_by_name(cell_name, matching=True)
        refl_cells.extend(cells)

    # Create a water-filled reflector cell/universe to put in a lattice
    all_materials = self.openmoc_geometry.getAllMaterials()
    h2o = all_openmoc_cells[mod[0].id].getFillMaterial()
    reflector_cell = openmoc.Cell(name='Refined Reflector Cell')
    reflector_cell.setFill(all_materials[h2o.getId()])
    reflector = openmoc.Universe(name='Refined Reflector Universe')
    reflector.addCell(reflector_cell)

    # Sliced up water cells with a lattice
    mesh_per_pin = 4
    lattice = openmoc.Lattice(name='{} x {} Spaced Reflector'.format(mesh_per_pin))
    lattice.setWidth(width_x=1.26492 / mesh_per_pin, width_y=1.26492 / mesh_per_pin)
    template = [[reflector] * 17 * mesh_per_pin] * 17 * mesh_per_pin
    lattice.setUniverses([template])

    # Put the lattice in each of the water-filled reflector cells
    for refl_cell in refl_cells:
        all_openmoc_cells[refl_cell.id].setFill(lattice)
        all_openmoc_cells[refl_cell.id].setNumSectors(0)


def discretize_geometry_standalone(mat_mgxslib, openmoc_geometry):
    """Discretize the full core BEAVRS benchmark."""

    openmoc.log.py_printf('INFO', 'Discretizing the geometry...')

    openmc_geometry = mat_mgxslib.openmc_geometry
    all_openmoc_cells = openmoc_geometry.getAllMaterialCells()

    # Add angular sectors to all material-filled cells
    # FIXME: Only do this for pin cells!!!
    for cell_id in all_openmoc_cells:
        all_openmoc_cells[cell_id].setNumSectors(8)

    ###########################################################################
    # Discretize the instrument and guide tubes and BPs
    ###########################################################################

    # Find the fuel clad outer radius zcylinder
    all_surfs = openmoc_geometry.getAllSurfaces()
    for surf_id in all_surfs:
        if all_surfs[surf_id].getName() == 'Fuel clad OR':
            fuel_or = openmoc.castSurfaceToZCylinder(all_surfs[surf_id])

    # Find cells by their string names in the BEAVRS benchmark
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

    # Discretize each cell into radial rings
    for cell in instr_tube:
        all_openmoc_cells[cell.id].setNumRings(10)
    for cell in guide_tube:
        all_openmoc_cells[cell.id].setNumRings(10)
    for cell in burn_abs1:
        all_openmoc_cells[cell.id].setNumRings(5)
    for cell in burn_abs2:
        all_openmoc_cells[cell.id].setNumRings(5)
    for cell in mod:
        all_openmoc_cells[cell.id].setNumRings(5)
        all_openmoc_cells[cell.id].addSurface(surface=fuel_or, halfspace=+1)

    ###########################################################################
    # Discretize the fuel pin cells in all assemblies
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

    # Discretize all fuel cells within each pin cell universe into rings
    for pin_univ in pin_univs:
        all_cells = pin_univ.getCells()
        for cell_id in all_cells:
            if all_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_cells[cell_id].getFillMaterial().getName().lower():
                    all_cells[cell_id].setNumRings(5)

    ###########################################################################
    # Discretize the reflector cells using a lattice
    ###########################################################################

    # Store the names of each of the reflector cells in the BEAVRS benchmark
    cell_names = ['Water', 'North Baffle Outer Water',
                  'Baffle South radial outer: water',
                  'Baffle East radial outer: water',
                  'Baffle West radial outer: water',
                  'Baffle North East Tip radial outer: water',
                  'Baffle South East Tip radial outer: water',
                  'Tip Baffle Outer Water W',
                  'Tip Baffle Outer Water N',
                  'Baffle South West Tip radial outer: water',
                  'Baffle South West Corner radial outer: water',
                  'Baffle South East Corner radial outer: water',
                  'Baffle North East Corner radial outer: water',
                  'Corner Baffle Water Gap N',
                  'Corner Baffle Water Gap W',
                  'Corner Baffle Outer Water']

    # Find all of the water-filled cells which comprise the reflector
    refl_cells = []
    for cell_name in cell_names:
        cells = openmc_geometry.get_cells_by_name(cell_name, matching=True)
        refl_cells.extend(cells)

    # Create a water-filled reflector cell/universe to put in a lattice
    all_materials = openmoc_geometry.getAllMaterials()
    h2o = all_openmoc_cells[mod[0].id].getFillMaterial()
    reflector_cell = openmoc.Cell(name='Refined Reflector Cell')
    reflector_cell.setFill(all_materials[h2o.getId()])
    reflector = openmoc.Universe(name='Refined Reflector Universe')
    reflector.addCell(reflector_cell)

    # Sliced up water cells with a lattice
    mesh_per_pin = 4
    lattice = openmoc.Lattice(name='{} x {} Spaced Reflector'.format(mesh_per_pin))
    lattice.setWidth(width_x=1.26492 / mesh_per_pin, width_y=1.26492 / mesh_per_pin)
    template = [[reflector] * 17 * mesh_per_pin] * 17 * mesh_per_pin
    lattice.setUniverses([template])

    # Put the lattice in each of the water-filled reflector cells
    for refl_cell in refl_cells:
        all_openmoc_cells[refl_cell.id].setFill(lattice)
        all_openmoc_cells[refl_cell.id].setNumSectors(0)