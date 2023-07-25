import sys

import numpy as np

from typing import Tuple
from mayavi import mlab

n_0         = 5
thickness_0 = 775
zetas_0     = np.linspace(thickness_0, thickness_0/2, n_0 + 2)[1:-1]
zetas_1     = [2750.0, 2450.0, 2175.0, 1900.0, 1625.0]


def make_layers (n, thickness = thickness_0, zetas = zetas_0):
    correction = 0
    answer     = []
    for x in zetas:
        new         = quadr(x, thickness) + correction
        answer     += [abs(new)]
        correction += answer[-1]

    return list(zip(zetas, thickness - np.array(zetas), answer))


def quadr (z, t):
    c = 1/t
    return c * (z - t) * z


def draw_isolated_points(vertices_file: str, folder: str = "."):
    vertices_file     = folder + "/" + vertices_file
    vertices: np.array = np.loadtxt(vertices_file, delimiter=",")
    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]

    fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1))


    mlab.points3d(
            x, y, z, 
            scale_factor=.061, 
            color=(1, 0, 0)
    )

    mlab.show()


def plot_all_tetrahedra(
        vertices_file: str,      # nodes.dat
        connectivity_file: str,  # elements.dat
        folder: str = ".",
        isolated_points: str = None,
        number_of_isolated_points_to_draw: int = 3,
        vert_delim: str = ",",
        elem_delim: str = ",",
        colors: Tuple = (.2,.3,.4)
    ) -> None:

    vertices_file     = folder + "/" + vertices_file
    connectivity_file = folder + "/" + connectivity_file

    vertices: np.array = np.loadtxt(vertices_file, delimiter=vert_delim)

    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]

    del vertices
    
    all_elements: np.array = np.loadtxt(connectivity_file, delimiter=elem_delim)

    # Here we have all tetrahedra, which have 6 edges.
    n_con = 0
    for e in all_elements:
        n_con = n_con + 6

    connections = np.zeros((n_con,2))

    links: np.array = np.array([
                            [0,1,2,3,0,1],
                            [1,2,3,0,2,3]
                        ]).T

    last: int = 0

    try:
        for i in range(len(all_elements)):
            row = np.array(all_elements[i])
            connections[last: last + 6, :] = row[links]
            last = last + 6

    except IndexError:
        print(i, row)
        input()

    except KeyError:
        print(i, row)
        input()

    # plot:
    fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1))

    src = mlab.pipeline.scalar_scatter(x, y, z)
    
    src.mlab_source.dataset.lines = connections 

    mlab.pipeline.surface(src, color=colors)

    src.update()

    # Put some balls to visualize the axes
    # mlab.points3d(
    #     [0, 1, 0, 0],
    #     [0, 0, 1, 0],
    #     [0, 0, 0, 1],
    #     scale_factor=.06,
    #     color=(.5, 1, .5)
    # )

    # # These two are GenMesh::Point pA and GenMesh::Point pD and Point B
    # mlab.points3d(
    #     [.15, -.15, .05],
    #     [14, 14, 14],
    #     [0, 0, 0],
    #     scale_factor=.04,
    #     color=(.7,.5,.7)
    # )
    # # Inner layer scatter
    # mlab.points3d(
    #     list(x)[:122],
    #     list(y)[:122],
    #     [0] * 122,
    #     scale_factor=.01,
    #     color=(.4,.5,.4)
    # )
    # # bdry_adj_true
    # mlab.points3d(
    #     [0.294000, -0.294000],
    #     [14.000000, 14.000000],
    #     [0.077500, 0.077500],
    #     scale_factor=.02,
    #     color=(.2,.3,.4)
    # )

    if isolated_points:

        isolated_points = folder + "/" + isolated_points
        isolated_points = np.loadtxt(isolated_points, delimiter=vert_delim)  # [0:number_of_isolated_points_to_draw,:]

        number_of_isolated_points_to_draw = len(isolated_points)
        
        mlab.points3d(
            isolated_points[:,0],
            isolated_points[:,1],
            np.zeros(number_of_isolated_points_to_draw), 
            scale_factor=.061, 
            color=(1, 0, 0)
        )

    mlab.show()


def plot_separate(
        vertices_file_inner: str,
        elements_file_inner: str,
        vertices_file_surrounding: str = None,
        elements_file_surrounding: str = None,
        isolated_points: str = None,
        vert_delim: str = None,
        nodes_delim: str = None,
        colors_inner: Tuple = (.2,.7,.2),
        colors_surrounding: Tuple = (.7,.2,.2)
    ) -> None:
    
    """
    plot_separate(
        "experiments/wafer24ago22/nodes_inner.dat",
        "experiments/wafer24ago22/elements_inner.dat",
        "experiments/wafer24ago22/physical_vertices_16613627703281198.dat",
        "experiments/wafer24ago22/wafer_profile_elements_16613627703281198.dat"
    ) 
    """

    vertices_inner = np.loadtxt(vertices_file_inner, delimiter=",")

    x = vertices_inner[:, 0]
    y = vertices_inner[:, 1]
    z = vertices_inner[:, 2]

    del vertices_inner
    
    elements_inner: np.array = np.loadtxt(elements_file_inner, delimiter=",")

    # all_elements: np.array = np.vstack((elements_inner, elements))

    # Here we have all tetrahedra, which have 6 edges.
    cant_edges = {4: 6}
    n_con = 0
    for e in elements_inner:
        n_con = n_con + 6

    connections: np.array = np.zeros((n_con, 2))

    connections_tetra: np.array = np.array([
        [0,1,2,3,0,1],
        [1,2,3,0,2,3]
    ]).T

    new_connections = {4: connections_tetra}

    last = 0

    for i in range(len(elements_inner)):
        row = np.array(elements_inner[i])
        #add connections depending on the type of elelement

        # TODO: this [[]+1]-1 at the end of the line should be depending of lang = C/Octave
        connections[last: last + cant_edges[row[0]], :] = row[new_connections[row[0]] + 1] - 1
        last = last + cant_edges[row[0]]

    #plot:
    fig = mlab.figure(1, size=(400, 400), bgcolor=(.92, .92, .92))

    src = mlab.pipeline.scalar_scatter(x, y, z)
    
    src.mlab_source.dataset.lines = connections 

    mlab.pipeline.surface(src, color=colors_inner)

    src.update()

    # Put some balls to visualize the axes
    mlab.points3d(
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        scale_factor=.06,
        color=(.5, 1, .5)
    )

    if vertices_file_surrounding:
        vertices_surrounding: np.array = np.loadtxt(vertices_file_surrounding, delimiter=",")
        elements_surrounding: np.array = np.loadtxt(elements_file_surrounding)
        # prepend a legacy 4
        rows: int = elements_surrounding.shape[0]
        elements_surrounding = np.hstack((4 * np.ones((rows, 1)), elements_surrounding))

        x_surr = vertices_surrounding[:, 0]
        y_surr = vertices_surrounding[:, 1]
        z_surr = vertices_surrounding[:, 2]

        del vertices_surrounding
        
        n_con_surr: int = 6 * elements_surrounding.shape[0]

        connections_surr = np.zeros((n_con_surr, 2))

        new_connections_surr = {4: connections_tetra}

        last_surr: int = 0

        for i in range(len(elements_surrounding)):
            row = np.array(elements_surrounding[i])

            # TODO: this [[]+1]-1 at the end of the line should be depending of lang = C/Octave

            try:
                connections_surr[last_surr: last_surr + 6, :] = row[connections_tetra + 1] - 1
                last_surr = last_surr + 6
            except:
                print(i, row, row.shape)
                exit()

        #plot:
        fig = mlab.figure(1, size=(400, 400), bgcolor=(.92, .92, .92))

        src = mlab.pipeline.scalar_scatter(x_surr, y_surr, z_surr)
        
        src.mlab_source.dataset.lines = connections_surr 

        mlab.pipeline.surface(src, color=colors_surrounding)

        src.update()

    if isolated_points:
        isolated_points = np.loadtxt(isolated_points, delimiter=vert_delim)
        
        mlab.points3d(
            isolated_points[:,0],
            isolated_points[:,1],
            isolated_points[:,2], 
            np.linspace(.001, .3, len(isolated_points)),
            scale_factor=1, 
            color=(1, 0, 0)
        )

    mlab.show()
