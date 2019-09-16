#!/usr/bin/env python3.5

"""
vtugeo aims to convert  2D (x-y) vtu unstructured grids to geodataframe 
structures and vice versa. Furthermore, output as polygon shapefile is provided
All node and element fields from vtu are transfered to the geodataframe 
and can be written out to shapefile as attributes. Node based fields
start with an "N_{i}" where i is the number of corresponding node
LIMITATION: UNTIL NOW ONLY PURE TRIANGULAR VTU GRIDS ARE ALLOWED

V02: exception of field data is vector data catched
    improved readability
"""

__author__ = "Erik Nixdorf"
__propertyof__ = "Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ. "
__email__ = "erik.nixdorf@ufz.de"
__version__ = "0.2"

# Import functions

import pandas as pd
import geopandas as gpd
import sys
import meshio
import numpy as np
from operator import itemgetter
import shapely.geometry


def unique_nd(a):
    """
    A function which allows to find unique sets in n dimension
    inspired by https://stackoverflow.com/questions/8560440/
    removing-duplicate-columns-and-rows-from-a-numpy-2d-array
    """
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1)
    return a[ui]


def vtu_to_geo(fctn_inpt='testmesh.vtu', inpt_flfrmt='vtu-ascii',
               srs_crs='EPSG:31468', elmt_type='triangle', Output=None):
    """
    vtu unstructured grid to geodataframe
    """
    # check whether the input is string or instance of vtk objects
    if isinstance(fctn_inpt, str):
        filename = fctn_inpt
        # read the mesh
        mesh = meshio.read(filename, file_format=inpt_flfrmt)
    else:
        sys.exit('Only explicite VTK and VTU files allowed as input!')

    # get triangles nodes, quads could exist as well, not covered yet
    if len(mesh.cells.keys()) == 1:
        try:
            mesh_nds = mesh.cells[elmt_type]
        except:
            sys.exit('Only triangular 2D Meshes supported ')
    else:
        sys.exit('More than one type of elements not supported')

    # Get the Nodes corresponding to the elements
    # https://www.geeksforgeeks.org/python-accessing-all-elements-at-given-list-of-indexes/
    nds_cord = list(itemgetter(*mesh_nds)(mesh.points))
    tri_geom = [shapely.geometry.Polygon(zip(nds_cord[i][:, 0],
                                             nds_cord[i][:, 1],
                                             nds_cord[i][:, 2]))
                for i in range(0, len(nds_cord))]
    # next we create the dataframe for the element fields
    df_elmnt_flds = pd.DataFrame()
    for key in mesh.cell_data[elmt_type]:
        # check whether scalar or vector input
        if len(mesh.cell_data[elmt_type][key].shape) == 1:
            df_elmnt_flds[key] = mesh.cell_data[elmt_type][key]
        else:
            print(key, 'could not used, no support of vector data input')

    # finally for the node_based informationn get the NodalValues as Element
    # based attributes with particular column for each node
    meshpnt_fldnms = list(mesh.point_data.keys())
    # create the dataframe
    df_pnt_flds = pd.DataFrame()
    # we have to define the order of the nodes in a polygon , in vtk it is
    # counter clockwise, however in shapefile it is clockwise
    if isinstance(Output, str):
        NodeOrder = [0, 2, 1]
    else:
        NodeOrder = [0, 1, 2]
    # for each of the fields we convert and append
    for meshpnt_fldnm in meshpnt_fldnms:
        fld_data = np.stack(itemgetter(*mesh_nds)
                            (mesh.point_data[meshpnt_fldnm]))
        try:
            df_pnt_flds[['N_'+str(i)+'_'+meshpnt_fldnm
                         for i in NodeOrder]] = pd.DataFrame(fld_data)
        except:
            print(meshpnt_fldnm, ' is vector data input, so support yet')
    print('... done writing fields')
    # Finally we merge all the dataframes to one geodataframe    
    # create a geopandas dataframe
    gdf_mesh = gpd.GeoDataFrame(pd.concat([df_elmnt_flds, df_pnt_flds], axis=1),
                                geometry=tri_geom, crs=srs_crs)
    if isinstance(Output, str):
        gdf_mesh.to_file(Output)
    return gdf_mesh


def geo_to_vtu(fctn_inpt='testmesh.shp',
               elmt_fldnms_long=['MaterialIDs', 'bulk_element_ids', 'bulk_face_ids'],
               pnt_fldnms_long=['bulk_node_ids'],
               Output='outputmsh.vtu',
               Output_flfrmt='vtu-binary'):
    """
    This Function takes a geodataframe/shapefile and converts to a VTU mesh
    conventions: only triangles supported
    All fields starting with a number from 0 (to 2) belong to nodal values
    BULK NODE IDS ARE PROBABLY WRONGLY ASSIGNED
    elmt_fldnms_long and pnt_fldnms_long define some (more than 10 digit)
    OGS field names which we need later to replace
    CHECK BY INDEX
    @author: Erik
    """
    if isinstance(fctn_inpt, str):
        gdf_inpt = gpd.read_file(fctn_inpt)
    else:
        gdf_inpt = fctn_inpt

    # extract cell_data arrays, which neither start with "N_" or geometry
    cell_fld_nms = [col_nm for col_nm in gdf_inpt.columns
                    if col_nm.startswith('N_') is False
                    and col_nm != 'geometry']
    print('...cell fieldnames ', cell_fld_nms, ' have been identified')
    # extract the nodal_array_nms and reshape datasets
    # We check for all point_datafields and append them to existing array
    pnt_fld_uniqs = set()
    for col_nm in gdf_inpt.columns:
        if col_nm.startswith('N_'):
            node_array_nms = col_nm
            pnt_fld_uniqs.add(col_nm[4:])
    print('...point fieldnames ', list(pnt_fld_uniqs), ' have been identified')
    # we check whether we have 2D or 3D Data
    if gdf_inpt['geometry'][0].has_z:
        dimension = 3
        dim_clm_name = ['X', 'Y', 'Z']
        clm_nms = ['Element_ID', 'X', 'Y', 'Z']
    else:
        dimension = 2
        dim_clm_name = ['X', 'Y']
        clm_nms = ['Element_ID', 'X', 'Y']
    # %% Extract all information from the shape which are relevant for VTU Mesh
    # first coordiantes for all nodes forming a polygon (mesh element)
    # https://stackoverflow.com/questions/20474549/
    # extract-points-coordinates-from-a-polygon-in-shapely
    pnts_per_elmnt = len(gdf_inpt['geometry'][0].exterior.coords[:-1])
    df_pnt_cell_flds_lst = [np.hstack((np.ones((pnts_per_elmnt, 1))*i,
                            np.array(gdf_inpt['geometry'][i].exterior.coords)[:-1]))
                            for i in range(0, gdf_inpt['geometry'].shape[0])]
    # create pandas array which contains all elements and nodes
    df_pnt_cell_flds = pd.DataFrame(np.vstack(df_pnt_cell_flds_lst),
                                    columns=clm_nms)
    # create a new df with all fields belonging to elements are repeated two t
    df_cell_flds = pd.concat([gdf_inpt[cell_fld_nms]] *
                             pnts_per_elmnt).sort_index()
    # append the cell fld to the pnt_cell_flds
    df_pnt_cell_flds = pd.concat([df_pnt_cell_flds,
                                  df_cell_flds.reset_index(drop=True)], axis=1)
    # We check for all point_datafields and append to existing numpy array
    for pnt_fld_uniq in pnt_fld_uniqs:
        node_array_nms = ['N_'+str(i)+'_'+pnt_fld_uniq
                          for i in range(0, dimension)]
        df_pnt_cell_flds[pnt_fld_uniq] = gdf_inpt[node_array_nms].to_numpy().flatten().reshape(-1, 1)
    print('... done processing point fields')
    # set the last column to new node_id
    # finding the unique points using our new unique2d function
    df_uniq_pts = pd.DataFrame(unique_nd(df_pnt_cell_flds[dim_clm_name].to_numpy()),
                               columns=dim_clm_name)
    df_uniq_pts['Node_ID'] = np.linspace(0, len(df_uniq_pts)-1, len(df_uniq_pts))
    # merge the uniquepoints with the element based dataframe
    # in order to get the NodeID for each point in each element
    df_pnt_cell_flds = pd.merge(df_uniq_pts, df_pnt_cell_flds,
                                how="inner", on=dim_clm_name)
    print('... done merging point and cell flds')
    # %%Create the dataframe with the points
    df_msh_pnts = df_pnt_cell_flds.drop_duplicates(subset=['Node_ID']).sort_values(by=['Node_ID'])
    # delete element relevant fields
    df_msh_pnts.drop(columns=cell_fld_nms, inplace=True)
    # delete the element ID
    df_msh_pnts.drop(columns=['Element_ID', 'Node_ID'], inplace=True)
    # create the numpy array for point coordinates
    msh_pnt_arr = df_msh_pnts[dim_clm_name].to_numpy()
    # create the dictionary with pnt fields
    msh_pnt_flds = {}
    for pnt_fld_uniq in pnt_fld_uniqs:
        msh_pnt_fld = df_msh_pnts[pnt_fld_uniq].to_numpy()
        # check if the cell_array name to be replaced (10 digit bug in shp)
        try:
            pnt_fld_uniq = [s for s in pnt_fldnms_long if pnt_fld_uniq in s][0]
        except:
            pass
        # add to dictionary
        msh_pnt_flds[pnt_fld_uniq] = msh_pnt_fld
    # Next we create the msh cell dictionary and
    # the array which points element to node

    # sort by element id
    df_msh_cells = df_pnt_cell_flds.sort_values(by=['Element_ID'])
    # get node index per cell
    msh_cell_pntids_arr = df_msh_cells['Node_ID'].to_numpy().reshape(
            int(np.max(df_pnt_cell_flds['Element_ID']))+1, 3).astype(int)
    msh_cell_pntids = {"triangle": msh_cell_pntids_arr}
    # reduce the number of columns
    df_msh_cell_flds = df_msh_cells.drop_duplicates(subset=['Element_ID'])
    df_msh_cell_flds = df_msh_cell_flds[cell_fld_nms]
    # write the dictionary
    msh_cell_fld_dict = {}
    for cell_fld_nm in cell_fld_nms:
        msh_cell_fld = df_msh_cell_flds[cell_fld_nm].to_numpy()
        # check if the cell_array name to be replaced (10 digit bug in shp)
        try:
            cell_fld_nm = [s for s in elmt_fldnms_long if cell_fld_nm in s][0]
        except:
            pass
        # add to dictionary
        msh_cell_fld_dict[cell_fld_nm] = msh_cell_fld
    # %% create object and write
    mesh = meshio.Mesh(msh_pnt_arr, msh_cell_pntids)
    # add cellfield data
    mesh.cell_data['triangle'] = msh_cell_fld_dict
    # add point field data
    mesh.point_data = msh_pnt_flds
    # write out
    if isinstance(Output, str):
        meshio.write(Output, mesh, file_format=Output_flfrmt)
