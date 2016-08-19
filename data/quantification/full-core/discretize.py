import opencg
import openmoc


def find_assembly(assembly_name, opencg_geometry, wrap_geometry=True):
    """Find a fuel assembly with some string name in the BEAVRS OpenCG model."""

    # Get all OpenCG Universes
    all_univ = opencg_geometry.get_all_universes()

    # Iterate over all Universes
    fuel_assembly = None
    for univ_id, univ in all_univ.items():
        if univ.name == assembly_name:
            fuel_assembly = univ

    # Wrap lattice in a Geometry if requested by the user
    if wrap_geometry:

        # Create a root Cell
        root_cell = opencg.Cell(name='root cell')
        root_cell.fill = fuel_assembly

        # Make mixed reflective / vacuum boundaries
        min_x = opencg.XPlane(x0=root_cell.fill.min_x, boundary='reflective')
        max_x = opencg.XPlane(x0=root_cell.fill.max_x, boundary='reflective')
        min_y = opencg.YPlane(y0=root_cell.fill.min_y, boundary='reflective')
        max_y = opencg.YPlane(y0=root_cell.fill.max_y, boundary='reflective')
        max_z = opencg.ZPlane(z0=208., boundary='reflective')
        min_z = opencg.ZPlane(z0=203., boundary='reflective')

        # Add boundaries to the root Cell
        root_cell.add_surface(surface=min_x, halfspace=+1)
        root_cell.add_surface(surface=max_x, halfspace=-1)
        root_cell.add_surface(surface=min_y, halfspace=+1)
        root_cell.add_surface(surface=max_y, halfspace=-1)
        root_cell.add_surface(surface=min_z, halfspace=+1)
        root_cell.add_surface(surface=max_z, halfspace=-1)

        # Create a root Universe
        root_univ = opencg.Universe(universe_id=0, name='root universe')
        root_univ.add_cell(root_cell)

        # Create a Geometry
        fuel_assembly = opencg.Geometry()
        fuel_assembly.root_universe = root_univ

    return fuel_assembly


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
    opencg_geometry = self.mat_mgxslibs[0].opencg_geometry
    all_opencg_cells = opencg_geometry.get_all_material_cells()
    all_openmoc_cells = self.openmoc_geometry.getAllMaterialCells()

    # Add angular sectors to all material-filled cells
    # FIXME: Only do this for pin cells!!!
    # FIXME: Or extract the water and do it directly to it
    # Add angular sectors to all material-filled cells
    for cell_id in all_openmoc_cells:
        all_openmoc_cells[cell_id].setNumSectors(8)

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
                    all_cells[cell_id].setNumSectors(8)


                    
    ###########################################################################
    # Discretize the reflector cells using same methodology as for C5G7
    ###########################################################################

#    core_lattice = openmc_geometry.get_lattices_by_name('Core Lattice', case_sensitive=True)
    water = openmc_geometry.get_cells_by_name('Water', matching=True)[0]
#    water_n = openmc_geometry.get_cells_by_name('Baffle North radial outer: water', matching=True)[0]
    water_n = openmc_geometry.get_cells_by_name('North Baffle Outer Water', matching=True)[0]
    water_s = openmc_geometry.get_cells_by_name('Baffle South radial outer: water', matching=True)[0]
    water_e = openmc_geometry.get_cells_by_name('Baffle East radial outer: water', matching=True)[0]
    water_w = openmc_geometry.get_cells_by_name('Baffle West radial outer: water', matching=True)[0]
    
    water_ne = openmc_geometry.get_cells_by_name('Baffle North East Tip radial outer: water', matching=True)[0]
    water_se = openmc_geometry.get_cells_by_name('Baffle South East Tip radial outer: water', matching=True)[0]
#    water_nw = openmc_geometry.get_cells_by_name('Baffle North West Tip radial outer: water', matching=True)[0]
    water_nw = openmc_geometry.get_cells_by_name('Tip Baffle Outer Water W', matching=True)[0]
    water_nw2 = openmc_geometry.get_cells_by_name('Tip Baffle Outer Water N', matching=True)[0]
    water_sw = openmc_geometry.get_cells_by_name('Baffle South West Tip radial outer: water', matching=True)[0]
    
    water_sw_corn = openmc_geometry.get_cells_by_name('Baffle South West Corner radial outer: water', matching=True)[0]
    water_se_corn = openmc_geometry.get_cells_by_name('Baffle South East Corner radial outer: water', matching=True)[0]
#    water_nw_corn = openmc_geometry.get_cells_by_name('Baffle North West Corner radial outer: water', matching=True)[0]
    water_ne_corn = openmc_geometry.get_cells_by_name('Baffle North East Corner radial outer: water', matching=True)[0]

    watera = openmc_geometry.get_cells_by_name('Corner Baffle Water Gap N', matching=True)[0]
    waterb = openmc_geometry.get_cells_by_name('Corner Baffle Outer Water', matching=True)[0]
    waterc = openmc_geometry.get_cells_by_name('Corner Baffle Water Gap W', matching=True)[0]

#    Baffle South West Corner radial outer: water
    
    all_materials = self.openmoc_geometry.getAllMaterials()

    # Reflector
    h2o = all_openmoc_cells[mod[0].id].getFillMaterial()
    reflector_cell = openmoc.Cell(name='moderator')
    reflector_cell.setFill(all_materials[h2o.getId()])

    reflector = openmoc.Universe(name='Reflector')
    reflector.addCell(reflector_cell)

    refined_reflector_cell = openmoc.Cell(name='Semi-Finely Spaced Reflector')
    refined_reflector = openmoc.Universe(name='Semi-Finely Spaced Moderator')
    refined_reflector.addCell(refined_reflector_cell)
    
    # Sliced up water cells - semi finely spaced
    lattice = openmoc.Lattice(name='Semi-Finely Spaced Reflector')
    lattice.setWidth(width_x=1.26492, width_y=1.26492)
    template = [[reflector] * 17] * 17
    lattice.setUniverses([template])
    refined_reflector_cell.setFill(lattice)

    all_openmoc_cells[water.id].setFill(lattice)
    all_openmoc_cells[water_n.id].setFill(lattice)
    all_openmoc_cells[water_s.id].setFill(lattice)
    all_openmoc_cells[water_e.id].setFill(lattice)
    all_openmoc_cells[water_w.id].setFill(lattice)
    all_openmoc_cells[water_ne.id].setFill(lattice)
    all_openmoc_cells[water_se.id].setFill(lattice)
    all_openmoc_cells[water_nw.id].setFill(lattice)
    all_openmoc_cells[water_nw2.id].setFill(lattice)
    all_openmoc_cells[water_sw.id].setFill(lattice)
    all_openmoc_cells[water_sw_corn.id].setFill(lattice)
    all_openmoc_cells[water_se_corn.id].setFill(lattice)
#    all_openmoc_cells[water_nw_corn.id].setFill(lattice)
    all_openmoc_cells[water_ne_corn.id].setFill(lattice)
    all_openmoc_cells[waterc.id].setFill(lattice)
    all_openmoc_cells[waterb.id].setFill(lattice)
    all_openmoc_cells[watera.id].setFill(lattice)

    all_openmoc_cells[water.id].setNumSectors(0)
    all_openmoc_cells[water_n.id].setNumSectors(0)
    all_openmoc_cells[water_s.id].setNumSectors(0)
    all_openmoc_cells[water_e.id].setNumSectors(0)
    all_openmoc_cells[water_w.id].setNumSectors(0)
    all_openmoc_cells[water_ne.id].setNumSectors(0)
    all_openmoc_cells[water_se.id].setNumSectors(0)
    all_openmoc_cells[water_nw.id].setNumSectors(0)
    all_openmoc_cells[water_nw2.id].setNumSectors(0)
    all_openmoc_cells[water_sw_corn.id].setNumSectors(0)
    all_openmoc_cells[water_se_corn.id].setNumSectors(0)
#    all_openmoc_cells[water_nw_corn.id].setNumSectors(0)
    all_openmoc_cells[water_ne_corn.id].setNumSectors(0)
    all_openmoc_cells[waterc.id].setNumSectors(0)
    all_openmoc_cells[waterb.id].setNumSectors(0)
    all_openmoc_cells[watera.id].setNumSectors(0)

    # FIXME
    return


    '''

    ###########################################################################
    # Discretize the reflector cells using same methodology as for C5G7
    ###########################################################################

    all_materials = self.openmoc_geometry.getAllMaterials()

    # Reflector
    water = top_left.getFillMaterial()
    reflector_cell = openmoc.Cell(name='moderator')
    reflector_cell.setFill(all_materials[water.getId()])

    reflector = openmoc.Universe(name='Reflector')
    reflector.addCell(reflector_cell)

    # Cells
    refined_reflector_cell = openmoc.Cell(name='Semi-Finely Spaced Reflector')
    right_reflector_cell = openmoc.Cell(name='Right Reflector')
    corner_reflector_cell = openmoc.Cell(name='Bottom Corner Reflector')
    bottom_reflector_cell = openmoc.Cell(name='Bottom Reflector')

    refined_reflector = openmoc.Universe(name='Semi-Finely Spaced Moderator')
    right_reflector = openmoc.Universe(name='Right Reflector')
    corner_reflector = openmoc.Universe(name='Bottom Corner Reflector')
    bottom_reflector = openmoc.Universe(name='Bottom Reflector')

    refined_reflector.addCell(refined_reflector_cell)
    right_reflector.addCell(right_reflector_cell)
    corner_reflector.addCell(corner_reflector_cell)
    bottom_reflector.addCell(bottom_reflector_cell)

    # Initialize a list of reflector lattices
    lattices = list()

    # Sliced up water cells - semi finely spaced
    lattices.append(openmoc.Lattice(name='Semi-Finely Spaced Reflector'))
    lattices[-1].setWidth(width_x=1.25984/10., width_y=1.25984/10.)
    template = [[reflector] * 10] * 10
    lattices[-1].setUniverses([template])
    refined_reflector_cell.setFill(lattices[-1])

    # Sliced up water cells - right side of geometry
    lattices.append(openmoc.Lattice(name='Right Reflector'))
    lattices[-1].setWidth(width_x=1.25984, width_y=1.25984)
    template = [[refined_reflector] * 11 + [reflector] * 6] * 17
    lattices[-1].setUniverses([template])
    right_reflector_cell.setFill(lattices[-1])

    # Sliced up water cells for bottom corner of geometry
    lattices.append(openmoc.Lattice(name='Bottom Corner Reflector'))
    lattices[-1].setWidth(width_x=1.25984, width_y=1.25984)
    template = [[refined_reflector] * 11 + [reflector] * 6] * 11
    template += [[reflector] * 17] * 6
    lattices[-1].setUniverses([template])
    corner_reflector_cell.setFill(lattices[-1])

    # Sliced up water cells for bottom of geometry
    lattices.append(openmoc.Lattice(name='Bottom Reflector'))
    lattices[-1].setWidth(width_x=1.25984, width_y=1.25984)
    template = [[refined_reflector] * 17] * 11
    template += [[reflector] * 17] * 6
    lattices[-1].setUniverses([template])
    bottom_reflector_cell.setFill(lattices[-1])

    # Find the root cell
    all_cells = self.openmoc_geometry.getAllCells()
    for cell_id in all_cells:
        if all_cells[cell_id].getName() == 'root cell':
            root_cell = all_cells[cell_id]

    # Extract core lattice
    core_lattice = root_cell.getFillUniverse()
    core_lattice = openmoc.castUniverseToLattice(core_lattice)

    # Assign reflector universes to core lattice
    core_lattice.updateUniverse(0, 0, 0, bottom_reflector)
    core_lattice.updateUniverse(1, 0, 0, bottom_reflector)
    core_lattice.updateUniverse(2, 0, 0, corner_reflector)
    core_lattice.updateUniverse(2, 1, 0, right_reflector)
    core_lattice.updateUniverse(2, 2, 0, right_reflector)
'''

def discretize_geometry_standalone(mat_mgxslib, openmoc_geometry):
    """Discretize a 2x2 checkerboard of BEAVRS fuel assemblies surrounded by
    a water reflector into OpenMOC flat source regions."""

    openmoc.log.py_printf('INFO', 'Discretizing the geometry...')

    opencg_geometry = mat_mgxslib.opencg_geometry
    all_opencg_cells = opencg_geometry.get_all_material_cells()
    all_openmoc_cells = openmoc_geometry.getAllMaterialCells()

    # Add angular sectors to all material-filled cells
    for cell_id in all_opencg_cells:
        all_openmoc_cells[cell_id].setNumSectors(8)

    # FIXME: Discretize a 1.6% enriched assembly
    fuel_assembly = find_assembly('Fuel 1.6% enr instr no BAs', opencg_geometry)

    # Get the bounding box from the OpenCG geometry
    min_x = fuel_assembly.min_x
    min_y = fuel_assembly.min_y
    min_z = fuel_assembly.min_z
    max_x = fuel_assembly.max_x
    max_y = fuel_assembly.max_y
    max_z = fuel_assembly.max_z
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
    bottom_left = fuel_assembly.find_cell(x=min_x+0.1, y=min_y+0.1, z=mid_z)
    bottom_left = all_openmoc_cells[bottom_left.id]
    bottom_left.addSurface(surface=fuel_or, halfspace=+1)
    bottom_left.setNumRings(10)

    # Discretize the guide and instrument tubes
    instr_tube = fuel_assembly.find_cell(x=mid_x, y=mid_y, z=mid_z)
    guide_tube = fuel_assembly.find_cell(x=mid_x, y=mid_y+3.78, z=mid_z)
    instr_tube = all_openmoc_cells[instr_tube.id]
    guide_tube = all_openmoc_cells[guide_tube.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)


    # FIXME: Try to get the guide tube
    # Discretize the guide and instrument tubes
    print(min_x, min_y, min_z)
    print(max_x, max_y, max_z)
    print(mid_x, mid_y, mid_z)
    hmm_tube = fuel_assembly.find_cell(x=mid_x+1.25984*3, y=mid_y+1.25984*3, z=mid_z)
    hmm_tube = all_openmoc_cells[hmm_tube.id]
    hmm_tube.setNumRings(10)


    # Discretize the moderator around the guide and instrument tubes
    instr_tube_moderator = fuel_assembly.find_cell(x=mid_x+0.6, y=mid_y+0.6, z=mid_z)
    guide_tube_moderator = fuel_assembly.find_cell(x=mid_x+0.6, y=mid_y+3.78+0.6, z=mid_z)
    instr_tube_moderator = all_openmoc_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_openmoc_cells[guide_tube_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)



    # FIXME: Discretize a 3.1% enriched assembly w 20 BPs
    fuel_assembly = find_assembly('Fuel 3.1% enr instr 20', opencg_geometry)

    # Get the bounding box from the OpenCG geometry
    min_x = fuel_assembly.min_x
    min_y = fuel_assembly.min_y
    min_z = fuel_assembly.min_z
    max_x = fuel_assembly.max_x
    max_y = fuel_assembly.max_y
    max_z = fuel_assembly.max_z
    mid_x = (max_x + min_x) / 2.
    mid_y = (max_y + min_y) / 2.
    mid_z = (max_z + min_z) / 2.

    # Discretize the moderator cell. First, add the fuel clad outer radius to each
    # cell since this is not done by the BEAVRS builder but is needed for ringify
    bottom_left = fuel_assembly.find_cell(x=min_x+0.1, y=min_y+0.1, z=mid_z)
    bottom_left = all_openmoc_cells[bottom_left.id]
    bottom_left.addSurface(surface=fuel_or, halfspace=+1)
    bottom_left.setNumRings(10)

    # Discretize the guide and instrument tubes and burnable absorbers
    instr_tube = fuel_assembly.find_cell(x=mid_x, y=mid_y, z=mid_z)
    guide_tube = fuel_assembly.find_cell(x=mid_x, y=mid_y+3.78, z=mid_z)
    burn_abs1 = fuel_assembly.find_cell(x=mid_x+3.78, y=mid_y+7.56, z=mid_z)
    burn_abs2 = fuel_assembly.find_cell(x=mid_x+3.78, y=mid_y+7.56+0.3, z=mid_z)
    instr_tube = all_openmoc_cells[instr_tube.id]
    guide_tube = all_openmoc_cells[guide_tube.id]
    burn_abs1 = all_openmoc_cells[burn_abs1.id]
    burn_abs2 = all_openmoc_cells[burn_abs2.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)
    burn_abs1.setNumRings(5)
    burn_abs2.setNumRings(5)

    # Discretize the guide and instrument tubes and burnable absorber moderator
    instr_tube_moderator = fuel_assembly.find_cell(x=mid_x+0.6, y=mid_y+0.6, z=mid_z)
    guide_tube_moderator = fuel_assembly.find_cell(x=mid_x+0.6, y=mid_y+3.78+0.6, z=mid_z)
    burn_abs_moderator = fuel_assembly.find_cell(x=mid_x+3.78+0.6, y=mid_y+7.56+0.6, z=mid_z)
    instr_tube_moderator = all_openmoc_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_openmoc_cells[guide_tube_moderator.id]
    burn_abs_moderator = all_openmoc_cells[burn_abs_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)
    burn_abs_moderator.setNumRings(5)



    # FIXME: Discretize a 2.4% enriched assembly w 12 BPs
    fuel_assembly = find_assembly('Fuel 2.4% enr instr 12', opencg_geometry)

    # Get the bounding box from the OpenCG geometry
    min_x = fuel_assembly.min_x
    min_y = fuel_assembly.min_y
    min_z = fuel_assembly.min_z
    max_x = fuel_assembly.max_x
    max_y = fuel_assembly.max_y
    max_z = fuel_assembly.max_z
    mid_x = (max_x + min_x) / 2.
    mid_y = (max_y + min_y) / 2.
    mid_z = (max_z + min_z) / 2.

    # Discretize the moderator cell. First, add the fuel clad outer radius to each
    # cell since this is not done by the BEAVRS builder but is needed for ringify
    bottom_left = fuel_assembly.find_cell(x=min_x+0.1, y=min_y+0.1, z=mid_z)
    bottom_left = all_openmoc_cells[bottom_left.id]
    bottom_left.addSurface(surface=fuel_or, halfspace=+1)
    bottom_left.setNumRings(10)

    # Discretize the guide and instrument tubes and burnable absorbers
    instr_tube = fuel_assembly.find_cell(x=mid_x, y=mid_y, z=mid_z)
    guide_tube = fuel_assembly.find_cell(x=mid_x, y=mid_y+3.78, z=mid_z)
    burn_abs1 = fuel_assembly.find_cell(x=mid_x+3.78, y=mid_y+7.56, z=mid_z)
    burn_abs2 = fuel_assembly.find_cell(x=mid_x+3.78, y=mid_y+7.56+0.3, z=mid_z)
    instr_tube = all_openmoc_cells[instr_tube.id]
    guide_tube = all_openmoc_cells[guide_tube.id]
    burn_abs1 = all_openmoc_cells[burn_abs1.id]
    burn_abs2 = all_openmoc_cells[burn_abs2.id]
    instr_tube.setNumRings(10)
    guide_tube.setNumRings(10)
    burn_abs1.setNumRings(5)
    burn_abs2.setNumRings(5)

    # Discretize the guide and instrument tubes and burnable absorber moderator
    instr_tube_moderator = opencg_geometry.find_cell(x=mid_x+0.6, y=mid_y+0.6, z=mid_z)
    guide_tube_moderator = opencg_geometry.find_cell(x=mid_x+0.6, y=mid_y+3.78+0.6, z=mid_z)
    burn_abs_moderator = opencg_geometry.find_cell(x=mid_x+3.78+0.6, y=mid_y+7.56+0.6, z=mid_z)
    instr_tube_moderator = all_openmoc_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_openmoc_cells[guide_tube_moderator.id]
    burn_abs_moderator = all_openmoc_cells[burn_abs_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)
    burn_abs_moderator.setNumRings(5)



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
