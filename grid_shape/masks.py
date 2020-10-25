import numpy as np
from grid_shape import grids as grids

def make_all_mask(input_n_ring = 5):
    pos, _  = grids.create_grid_univ(
        "trihex",
        1000,
        angle=0,
        input_n_ring = input_n_ring      
    )
    n_antennas = pos.shape[1]

    mask3 = np.ones(n_antennas, dtype = bool)
    return np.expand_dims(mask3, 1) 


def make_mask_random(input_n_ring = 5, n_keep_ratio=0.1):
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
