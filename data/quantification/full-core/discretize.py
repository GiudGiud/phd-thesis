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


def discretize_geometry(mat_mgxslib, openmoc_geometry):
    """Discretize a 2x2 checkerboard of BEAVRS fuel assemblies surrounded by
    a water reflector into OpenMOC flat source regions."""

    openmoc.log.py_printf('INFO', 'Discretizing the geometry...')

    opencg_geometry = mat_mgxslib.opencg_geometry
    all_openmoc_cells = openmoc_geometry.getAllMaterialCells()

    # FIXME
#    return

    # FIXME: Discretize a 1.6% enriched assembly
    fuel_assembly = find_assembly('Fuel 1.6% enr instr no BAs', opencg_geometry)
    all_opencg_cells = fuel_assembly.get_all_material_cells()

    # Add angular sectors to all material-filled cells
    for cell_id in all_opencg_cells:
        all_openmoc_cells[cell_id].setNumSectors(8)

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

    # Discretize the moderator around the guide and instrument tubes
    instr_tube_moderator = fuel_assembly.find_cell(x=mid_x+0.6, y=mid_y+0.6, z=mid_z)
    guide_tube_moderator = fuel_assembly.find_cell(x=mid_x+0.6, y=mid_y+3.78+0.6, z=mid_z)
    instr_tube_moderator = all_openmoc_cells[instr_tube_moderator.id]
    guide_tube_moderator = all_openmoc_cells[guide_tube_moderator.id]
    instr_tube_moderator.setNumRings(5)
    guide_tube_moderator.setNumRings(5)

    '''
    # Find all pin cell universes - universes containing a cell filled with fuel
    all_univs = openmoc_geometry.getAllUniverses()
    pin_univs = set()
    for univ_id in all_univs:
        all_cells = all_univs[univ_id].getAllCells()
        for cell_id in all_openmoc_cells:
            if all_openmoc_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_openmoc_cells[cell_id].getFillMaterial().getName().lower():
                    pin_univs.add(all_univs[univ_id])

    # Discretize all fuel cells within each pin cell universe into rings
    for pin_univ in pin_univs:
        all_cells = pin_univ.getCells()
        for cell_id in all_openmoc_cells:
            if all_openmoc_cells[cell_id].getType() == openmoc.MATERIAL:
                if 'fuel' in all_openmoc_cells[cell_id].getFillMaterial().getName().lower():
                    all_openmoc_cells[cell_id].setNumRings(5)
    '''