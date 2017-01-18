
def conIntTheta(con, X, Y, Z, dm, src, con_cri):
    import numpy as np
    """Concentration intergrate along constant R"""

    ### Main Body ###
    ## Resolution
    L = np.array([dm.Lx, dm.Ly])
    delta = np.array([dm.Dx, dm.Dy])
    src_norm = src / L
    rxy_max = np.abs((1-np.floor(src_norm/0.5)) * (L-delta) - src)
    r_max = np.sqrt(np.sum(rxy_max ** 2))
    r_min = np.min(L-delta-rxy_max)
    n_r = int(np.floor(r_max / np.amin(delta)))
    ang_res_rad = np.amin(delta) / r_max
    n_ang = int(np.ceil(2*np.pi/ang_res_rad))
    ang_res_rad = 2*np.pi / n_ang

    ## Variable used to store
    int_con = np.zeros((n_r, dm.Nz))
    r = np.arange(np.amin(delta),
            (np.floor(r_max/np.amin(delta))+1)*np.amin(delta), np.amin(delta))
    r_in = r_max*np.ones((dm.Nz, 1)) #after r_in, avg_con cannot be trusted
    x_cmax = np.zeros((dm.Nz, 1))
    y_cmax = np.zeros((dm.Nz, 1))

    ## Loop
    for iz in range(dm.Nz):
        # Determine r_in
        #- ix=0
        r_bdy = np.sqrt((X[iz,:,0]-src[0])**2 + (Y[iz,:,0]-src[1])**2)
        r_cri = r_bdy[con[iz,:,0] >= con_cri]
        if (r_cri.size != 0):
            r_in[iz] = np.amin(r_cri)
        #- ix=nx-1
        r_bdy = np.sqrt((X[iz,:,dm.Nx-1]-src[0])**2 +
                (Y[iz,:,dm.Nx-1]-src[1])**2)
        r_cri = r_bdy[con[iz,:,dm.Nx-1] >= con_cri]
        if (r_cri.size != 0):
            r_in[iz] = np.amin((np.amin(r_cri), r_in[iz]))
        #- iy=0
        r_bdy = np.sqrt((X[iz,0,:]-src[0])**2 + (Y[iz,0,:]-src[1])**2)
        r_cri = r_bdy[con[iz,0,:] >= con_cri]
        if (r_cri.size != 0):
            r_in[iz] = np.amin((np.amin(r_cri), r_in[iz]))
        #- iy=Ny-1
        r_bdy = np.sqrt((X[iz,dm.Ny-1,:]-src[0])**2 +
                (Y[iz,dm.Ny-1,:]-src[1])**2)
        r_cri = r_bdy[con[iz,dm.Ny-1,:] >= con_cri]
        if (r_cri.size != 0):
            r_in[iz] = np.amin((np.amin(r_cri), r_in[iz]))
        if (r_in[iz] == r_max):
            r_in[iz] = r_min
        # integrate the concentration
        for i_r in range(n_r):
            for i_ang in range(n_ang):
                # coordinate
                ang = ang_res_rad*i_ang;
                xc = np.cos(ang)*r[i_r] + src[0]
                yc = np.sin(ang)*r[i_r] + src[1]
                ix = xc / dm.Dx
                iy = yc / dm.Dy
                # closest share point
                ixc = np.floor(ix)
                iyc = np.floor(iy)
                if ((ixc >= 0 and ixc < dm.Nx-1) and (iyc >= 0 and iyc < dm.Ny-1)):
                    cx = ix - ixc
                    cy = iy - iyc
                    con_xy = (1-cy)*((1-cx)*con[iz,iyc,ixc] + cx*con[iz,iyc,ixc+1])\
                            + cy*((1-cx)*con[iz,iyc+1,ixc] + cx*con[iz,iyc+1,ixc+1])
                    int_con[i_r, iz] += con_xy
            int_con[i_r, iz] = int_con[i_r, iz] / n_ang * 2*np.pi*r[i_r]
    return int_con,r
    ### END Main Body ###
