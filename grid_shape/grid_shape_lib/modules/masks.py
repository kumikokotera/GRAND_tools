import numpy as np
import os
import matplotlib.pyplot as plt

from grid_shape_lib.modules import grids as grids
from grid_shape_lib.modules import hexy as hx


def plot_mask(mask, nring, plot_suffix):
    pos, _ = grids.create_grid_univ("trihex", 250, input_n_ring=nring)
    plt.figure()
    plt.scatter(
        pos[0],
        pos[1],
        c = 1-mask[:,0],
        s = 3
    )
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('{} antennas'.format(mask.sum()))
    plt.axis('equal')
    plt.savefig(
        os.path.join("layout_{}.png".format(plot_suffix))
    )

def get_spiral(a, b, theta):
    x = a*b**theta*np.cos(theta)
    y = a*b**theta*np.sin(theta)
    return np.array([x, y])


def get_closest_antenna(x, y, pos):
    d = np.sqrt((pos[0, :]-x)**2 + (pos[1, :]-y)**2)
    return np.argmin(d)


def make_spiral_mask_with_minipose(a, b, ntheta, shape="hexhex", minipos_nring=0, minipos_shape='trihex'):
    theta = np.arange(0, ntheta) / ntheta * 4 * np.pi
    spi1 = get_spiral(a, b, theta)
    spi1[0] -= a
    th = 2*np.pi/3
    r = np.array(
        [
            [np.cos(th), -np.sin(th)],
            [np.sin(th), np.cos(th)]
        ]
    )

    spi2 = np.dot(r, spi1)
    spi3 = np.dot(r, spi2)

    pos, _ = grids.create_grid_univ(
        "trihex",
        375,
        angle=0,
        input_n_ring=25    
    )

    posfine, _  = grids.create_grid_univ(
        "trihex",
        250,
        angle=0,
        input_n_ring = 20    
    )

    plt.scatter(
        pos[0],
        pos[1],
        c = "c",
        s = 3
    )

    p1 = pos.transpose()
    mask_new = np.zeros(p1.shape[0], dtype=bool) 

    for spi in [spi1, spi2, spi3]:

        for i in range(spi.shape[1]):
            xy = spi[:,i]
            id = get_closest_antenna(xy[0], xy[1], pos)
            mask_new[id] = 1


    minipos, _ = grids.create_grid_univ(
        minipos_shape, 250, do_prune=False, input_n_ring=minipos_nring
    )

    minipos2 = minipos.copy()

    idd = np.where(mask_new)[0]

    for id in idd:
        minipos2 = np.hstack([
            minipos2,
            minipos + pos[:,id].reshape(3,1)]
        )


    mask = prune_pos(
        posfine,
        minipos2
    )

    mask_hex1000 = make_coarse_out_of_fine_trihex(
    1000,
    shape, 
    250,
    20
    )

    new_pos = np.hstack(
        [
            posfine[:, mask[:,0]],
            posfine[:,mask_hex1000[:,0]]
        ]
    )
    mask2 = prune_pos(posfine, new_pos)
    return mask2


def make_spiral_mask(a, b):

    theta = np.arange(0, 50)/50 *4*np.pi
    spi1 = get_spiral(a, b, theta)
    spi1[0] -= a
    th = 2*np.pi/3
    r = np.array(
        [
            [np.cos(th), -np.sin(th)],
            [np.sin(th), np.cos(th)]
        ]
    )

    spi2 = np.dot(r, spi1)
    spi3 = np.dot(r, spi2)

    plt.figure()
    plt.plot(spi1[0], spi1[1], "r.")
    plt.plot(spi2[0], spi2[1], "g.")
    plt.plot(spi3[0], spi3[1], "b.")
    
    pos, _  = grids.create_grid_univ(
        "trihex",
        250,
        angle=0,
        input_n_ring = 25    
    )

    posfine, _  = grids.create_grid_univ(
        "trihex",
        250,
        angle=0,
        input_n_ring = 20    
    )

    plt.scatter(
        pos[0],
        pos[1],
        c = "c",
        s = 3
    )

    p1 = pos.transpose()
    mask_new = np.zeros(p1.shape[0], dtype=bool) 

    for spi in [spi1, spi2, spi3]:

        for i in range(spi.shape[1]):
            xy = spi[:,i]
            id = get_closest_antenna(xy[0], xy[1], pos)
            mask_new[id] = 1



    plt.figure()
    plt.scatter(
        pos[0],
        pos[1],
        c = 1-mask_new,
        s = 3
    )

    mask = prune_pos(
        posfine,
        pos[:,mask_new]
    )

    plt.figure()
    plt.scatter(
        posfine[0],
        posfine[1],
        c = 1-mask[:,0],
        s = 3
    )

    print(sum(mask))

    mask_hex1000 = make_coarse_out_of_fine_trihex(
    1000,
    "hexhex", 
    250,
    20
    )

    new_pos = np.hstack(
        [
            posfine[:, mask[:,0]],
            posfine[:,mask_hex1000[:,0]]
        ]
    )
    mask2 = prune_pos(posfine, new_pos)


    plt.figure()
    plt.scatter(
        posfine[0],
        posfine[1],
        c = 1-mask2[:,0],
        s = 3
    )
    return mask2


def make_all_mask(input_n_ring=5):
    pos, _  = grids.create_grid_univ(
        "trihex",
        1000,
        angle=0,
        input_n_ring = input_n_ring      
    )
    n_antennas = pos.shape[1]

    mask3 = np.ones(n_antennas, dtype = bool)
    return np.expand_dims(mask3, 1) 


def make_mask_random(input_n_ring=5, n_keep_ratio=0.1):
    """keep n_keep_ratio antennas from a trihex grid with size given by
    input_n_ring
    """

    pos, _  = grids.create_grid_univ(
        "trihex",
        1000,
        angle=0,
        input_n_ring = input_n_ring      
    )
    n_antennas = pos.shape[1]

    a = np.arange(n_antennas)
    np.random.shuffle(a)
    mask3 = np.ones(n_antennas, dtype = bool) * False
    mask3[a[0:np.int32(n_antennas * (n_keep_ratio))]] = True
    
    return np.expand_dims(mask3, 1) * True


def prune_pos(pos1, pos2):
        """
        create the mask of the commun positions bewtween 
        pos1 and pos2, mask to be applied to pos1.
        """

        p1 = pos1.transpose()
        p2 = pos2.transpose()
            
        grid_x = p1[:, 0]
        grid_y = p1[:, 1]

        x_pos_flat_fl = grid_x
        y_pos_flat_fl = grid_y

        grid_x2 = p2[:, 0]
        grid_y2 = p2[:, 1]

        x_pos_flat_fl2 = grid_x2
        y_pos_flat_fl2 = grid_y2

        xmin = np.min([x_pos_flat_fl.min(), x_pos_flat_fl2.min()])
        ymin = np.min([y_pos_flat_fl.min(), y_pos_flat_fl2.min()])
        xmax = np.max([x_pos_flat_fl.max(), x_pos_flat_fl2.max()])

        scal = (x_pos_flat_fl - xmin) + (xmax - xmin) * (y_pos_flat_fl - ymin)
        scal = np.floor(scal)

        scal2 = (x_pos_flat_fl2 - xmin) + (xmax - xmin) * (y_pos_flat_fl2 - ymin)
        scal2 = np.floor(scal2)
        
        mask_new = np.zeros(p1.shape[0], dtype=bool) 
        
        for i in range(scal.shape[0]):
            mask_new[i] = np.any(np.isclose(scal[i], scal2, rtol=1e-8))
        
        mask_new = np.expand_dims(mask_new, 1)
        return mask_new


def make_coarse_out_of_fine_trihex(
    step_coarse,
    shape_coarse,
    step_fine,
    input_n_ring
):

    pos_coarse, _ = grids.create_grid_univ(
        shape_coarse,
        step_coarse,
        do_prune=False,
        input_n_ring=input_n_ring
    )

    pos_fine, _ =  grids.create_grid_univ(
        "trihex",
        step_fine,
        do_prune=False,
        input_n_ring=input_n_ring
    )

    mask_new = prune_pos(pos_fine, pos_coarse)
    return mask_new 


def make_coarse_out_of_fine_trihex_2(
    step_coarse,
    shape_coarse,
    n_ring_coarse,
    step_fine,
    input_n_ring
):

    pos_coarse, _ = grids.create_grid_univ(
        shape_coarse,
        step_coarse,
        do_prune=False,
        input_n_ring=n_ring_coarse
    )

    pos_fine, _ =  grids.create_grid_univ(
        "trihex",
        step_fine,
        do_prune=False,
        input_n_ring=input_n_ring
    )

    mask_new = prune_pos(pos_fine, pos_coarse)
    return mask_new 


def make_trihex_new_out_of_125(pos_125, step, input_n_ring):
    pstep, offsetstep = grids.create_grid_univ(
        "trihex",
        step,
        do_prune=False,
        input_n_ring=input_n_ring
    )

    p125 = pos_125.transpose()
    pstep = pstep.transpose()

    p3 = np.vstack([pstep, p125])

    grid_x = p3[:,0]
    grid_y = p3[:,1]

    x_pos_flat_fl = grid_x
    y_pos_flat_fl = grid_y
    scal = (x_pos_flat_fl - x_pos_flat_fl.min()) + (x_pos_flat_fl.max() - x_pos_flat_fl.min()) * (y_pos_flat_fl - y_pos_flat_fl.min())
    scal = np.floor(scal)
    unique, index, inverse, count = np.unique(scal, return_index=True, return_inverse=True, return_counts=True)

    mask_dum = np.ones(p3.shape[0], dtype = bool)
    mask_dum[index] = False

    ind = np.where(mask_dum ==True)

    ind2 = np.array(ind) - pstep.shape[0]

    mask_new = np.zeros(p125.shape[0], dtype=bool)
    mask_new[ind2] = True

    mask_new = np.expand_dims(mask_new, 1)

    return mask_new


def make_mini_island():

    minipos, _ = grids.create_grid_univ(
        "trihex", 125, do_prune=False, input_n_ring=0
    )
    minipos2 = minipos.copy()
    pos, offset = grids.create_grid_univ(
        "trihex", 125, do_prune=False, input_n_ring=10
    )

    mask_tri1000 = make_trihex_new_out_of_125(pos, 1000, 10) 

    id1000 = np.where(mask_tri1000)[0]

    for id in id1000:
        minipos2 = np.hstack([minipos2,  minipos + pos[:,id].reshape(3,1)])

    mini_island_mask = prune_pos(pos, minipos2)
    return mini_island_mask


def make_mini_islands(fine_step, fine_n_ring, with_center=False):

    if with_center:
        minipos, _ = grids.create_grid_univ(
            "trihex", fine_step, do_prune=False, input_n_ring=0
        )
    else:
        minipos, _ = grids.create_grid_univ(
            "hexhex", fine_step, do_prune=False, input_n_ring=0
        )
    minipos2 = minipos.copy()
    pos, offset = grids.create_grid_univ(
        "trihex", fine_step, do_prune=False, input_n_ring=fine_n_ring
    )

    mask_tri1000 = make_coarse_out_of_fine_trihex(
    1000,
    "trihex", 
    fine_step,
    fine_n_ring
    )

    id1000 = np.where(mask_tri1000)[0]

    for id in id1000:
        minipos2 = np.hstack([minipos2,  minipos + pos[:,id].reshape(3,1)])

    mini_island_mask = prune_pos(pos, minipos2)
    return mini_island_mask


def make_limited_mini_island(
    island_shape,
    island_step,
    coarse_step,
    coarse_shape,
    limiting_shape,
    limiting_nring,
    limiting_step,
    fine_step, 
    fine_nring
):
    """
    create a masks of minislands of [tri,hex]hex (n_ring = 0), with given step
    at coarse positions within a given n_ring, prune on a fine grid
    """
    island_pos, _ = grids.create_grid_univ(
        island_shape, island_step, do_prune=False, input_n_ring=0
    )
    island_pos2 = island_pos.copy()

    fine_pos, _ = grids.create_grid_univ(
        "trihex", fine_step, do_prune=False, input_n_ring=fine_nring
    )

    limiting_pos, _ = grids.create_grid_univ(
        limiting_shape,
        limiting_step,
        do_prune=False,
        input_n_ring=limiting_nring
    )
    coarse_pos, _ = grids.create_grid_univ(
        coarse_shape, coarse_step, do_prune=False, input_n_ring=fine_nring
    )

    mask_coarse = make_coarse_out_of_fine_trihex(coarse_step, coarse_shape, fine_step, fine_nring)
    id_coarse = np.where(mask_coarse)[0]
    for id in id_coarse:
        island_pos2 = np.hstack([island_pos2,  island_pos + fine_pos[:,id].reshape(3,1)])


    m_temp = prune_pos(limiting_pos, island_pos2)

    pos_temp = np.hstack([limiting_pos[:, m_temp[:,0]], coarse_pos])

    mask_limited_mini_island = prune_pos(
        fine_pos,
        pos_temp
    )
    return mask_limited_mini_island

def make_mini_island2():

    minipos, _ = grids.create_grid_univ(
        "hexhex", 125, do_prune=False, input_n_ring=0
    )
    minipos2 = minipos.copy()
    pos, offset = grids.create_grid_univ(
        "trihex", 125, do_prune=False, input_n_ring=10
    )

    mask_tri1000 = make_trihex_new_out_of_125(pos, 1000, 10) 

    id1000 = np.where(mask_tri1000)[0]

    for id in id1000:
        minipos2 = np.hstack([minipos2,  minipos + pos[:,id].reshape(3,1)])

    mini_island_mask = prune_pos(pos, minipos2)
    return mini_island_mask


def make_central_island_out_of_fine(
    step_coarse,
    shape_coarse,
    step_island,
    shape_island,
    n_ring_island,
    step_fine,
    n_ring_fine
    ):
    """
    make a mask consisting of a coarse [tri,hex]hex of a given step
    with a central island [tri,hex]hex of n_ring of a given step,  out of a 
    fine trihex grid of a given n_ring
    """

    pos_coarse, _ = grids.create_grid_univ(
        shape_coarse,
        step_coarse,
        do_prune=False,
        input_n_ring=n_ring_fine
    )

    pos_fine, _ =  grids.create_grid_univ(
        "trihex",
        step_fine,
        do_prune=False,
        input_n_ring=n_ring_fine
    )
    pos_island, _ = grids.create_grid_univ(
        shape_island,
        step_island,
        do_prune=False,
        input_n_ring=n_ring_island
    )
    pos_to_trim = np.hstack([pos_coarse, pos_island])


    mask_new = prune_pos(pos_fine, pos_to_trim)
    return mask_new 


def make_centralisland_out_of_125(pos_125, input_n_ring):

    pstep, offsetstep = grids.create_grid_univ(
        "trihex",
        1000,
        do_prune=False,
        input_n_ring=10
    )
    pstep2, offsetstep2 = grids.create_grid_univ(
        "trihex",
        125,
        do_prune=False,
        input_n_ring=2
    )

    pstep = np.hstack([pstep, pstep2])
    
    mask_new = prune_pos(pos_125, pstep)

    return mask_new




def make_centralisland_out_of_125_v2(pos_125, input_n_ring):

    pstep, offsetstep = grids.create_grid_univ(
        "trihex",
        1000,
        do_prune=False,
        input_n_ring=10
    )
    pstep2, offsetstep2 = grids.create_grid_univ(
        "trihex",
        125,
        do_prune=False,
        input_n_ring=3
    )

    pstep = np.hstack([pstep2])


    mask_new = prune_pos(pos_125, pstep)

    return mask_new


def make_mask_gp13(dense_step, radius, input_n_ring=10):

    pos, _ = grids.create_grid_univ(
        "trihex",
        dense_step,
        do_prune=False,
        input_n_ring=input_n_ring
    )
    area = hx.get_area(radius) * 3
    cube0 = np.array(
        [
            [0, 0, 0],
            [1, 0, -1],
            [0, 1, -1],
        ]
    )

    pos_13 = hx.cube_to_pixel(cube0, radius)
    corners = hx.get_corners(pos_13, radius)

    sh = np.array(corners).shape
    corners = corners.transpose(0, 2, 1)
    corners = np.array(corners).reshape((sh[0]*sh[2], sh[1]))

    corners_x, corners_y = grids.remove_redundant_point(corners[:, 0], corners[:, 1])

    pos_gp13 = np.array([corners_x, corners_y]).transpose()
    pos_gp13 = pos_gp13.T
    print(pos_gp13.shape)
    print(pos.shape)
    mask_gp13 = prune_pos(pos, pos_gp13)
    return mask_gp13


def make_mask_gp13_2_with_cubeinput(dense_step, radius, cube0, input_n_ring=10):

    pos, _ = grids.create_grid_univ(
        "trihex",
        dense_step,
        do_prune=False,
        input_n_ring=input_n_ring
    )
    area = hx.get_area(radius) * 3
    
    pos_13 = hx.cube_to_pixel(cube0, radius)
    corners = hx.get_corners(pos_13, radius)

    sh = np.array(corners).shape
    corners = corners.transpose(0, 2, 1)
    corners = np.array(corners).reshape((sh[0]*sh[2], sh[1]))

    corners_x, corners_y = grids.remove_redundant_point(corners[:, 0], corners[:, 1])

    pos_gp13 = np.array([corners_x, corners_y]).transpose()
    pos_gp13 = pos_gp13.T
    print(pos_gp13.shape)
    print(pos.shape)
    mask_gp13 = prune_pos(pos, pos_gp13)
    return mask_gp13


def make_mask_gaa10_2_with_cubeinput(dense_step, radius, cube0, input_n_ring=10):

    pos, _ = grids.create_grid_univ(
        "trihex",
        dense_step,
        do_prune=False,
        input_n_ring=input_n_ring
    )
    area = hx.get_area(radius) * 3
    
    pos_13 = hx.cube_to_pixel(cube0, radius)
    corners = hx.get_corners(pos_13, radius)

    sh = np.array(corners).shape
    corners = corners.transpose(0, 2, 1)
    corners = np.array(corners).reshape((sh[0]*sh[2], sh[1]))

    corners_x, corners_y = grids.remove_redundant_point(corners[:, 0], corners[:, 1])

    pos_gaa10 = np.array([corners_x, corners_y]).transpose()
    pos_gaa10 = pos_gaa10.T
    print(pos_gaa10.shape)
    print(pos.shape)
    mask_gaa10 = prune_pos(pos, pos_gaa10)
    return mask_gaa10






def make_mask_gp13_2(dense_step, radius, input_n_ring=10):

    pos, _ = grids.create_grid_univ(
        "trihex",
        dense_step,
        do_prune=False,
        input_n_ring=input_n_ring
    )
    area = hx.get_area(radius) * 3
    cube0 = np.array(
        [
            [0, 0, 1],
            [1, 0, 1],
            [1, 0, 0],
        ]
    )

    pos_13 = hx.cube_to_pixel(cube0, radius)
    corners = hx.get_corners(pos_13, radius)

    sh = np.array(corners).shape
    corners = corners.transpose(0, 2, 1)
    corners = np.array(corners).reshape((sh[0]*sh[2], sh[1]))

    corners_x, corners_y = grids.remove_redundant_point(corners[:, 0], corners[:, 1])

    pos_gp13 = np.array([corners_x, corners_y]).transpose()
    pos_gp13 = pos_gp13.T
    print(pos_gp13.shape)
    print(pos.shape)
    mask_gp13 = prune_pos(pos, pos_gp13)
    return mask_gp13


def make_mask_gaA10(dense_step, radius, input_n_ring=10):

    pos, _ = grids.create_grid_univ(
        "trihex",
        dense_step,
        do_prune=False,
        input_n_ring=input_n_ring
    )
    #radius = 1000  # [m]

    #radius = 1291.472  # [m]

    area = hx.get_area(radius) * 3
    cube0 = np.array(
        [
            #[0, 0, 0],
            [1, 0, -1],
            [0, 1, -1],
        ]
    )

    pos_13 = hx.cube_to_pixel(cube0, radius)
 
    corners = hx.get_corners(pos_13, radius)


    sh = np.array(corners).shape
    corners = corners.transpose(0, 2, 1)
    corners = np.array(corners).reshape((sh[0]*sh[2], sh[1]))

    corners_x, corners_y = grids.remove_redundant_point(corners[:, 0], corners[:, 1])

    pos_gp13 = np.array([corners_x, corners_y]).transpose()
    pos_gp13 = pos_gp13.T
    print(pos_gp13.shape)
    print(pos.shape)
    mask_gp13 = prune_pos(pos, pos_gp13)
    return mask_gp13



def make_mask_gaa10_entangled_with_cubeinput(dense_step, radius, cube0, input_n_ring=10):

    pos, _ = grids.create_grid_univ(
        "trihex",
        dense_step,
        do_prune=False,
        input_n_ring=input_n_ring
    )

    pix_ = hx.cube_to_pixel(cube0, radius)
    corners_ = hx.get_corners(pix_, radius)

    sh = np.array(corners_).shape
    corners_ = corners_.transpose(0, 2, 1)
    corners_ = np.array(corners_).reshape((sh[0]*sh[2], sh[1]))

    corners_x, corners_y = grids.remove_redundant_point(corners_[:, 0], corners_[:, 1])

    pos_ = np.array([corners_x, corners_y]).transpose()
    pos_ = pos_.T
    print(pos_.shape)
    print(pix_)
    id = np.where(pos_[1,:]==pos_[1,:].max())[0]
    print(id)
    print(pos_[:,id])
    
    offset = pos_[:,id].T - pix_
    pos1 = np.hstack([pos_, pos_ +offset.reshape(2,1)])
    mask_gaa10 = prune_pos(pos, pos1)

    
    return mask_gaa10
