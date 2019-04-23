# Some tools to help with budget calculations in mom5 (I tried to keep the naming consistent with the MOM5 manual [link?])
import xarray as xr
from xgcm import Grid
import numpy as np
import os

# from xarrayutils.cm26_codebucket import drop_all_coords

# I think this should go to the EUC-shape processing, since it is quite specialized
def add_split_tendencies(ds):
    """Reconstructs various fluxes and tendencies (x-y_x) from the monthly averaged output. Hardcoded to o2 atm."""
    ds = ds.copy()
    rho = 1035  # reference density
    grid = Grid(ds)

    # This should be calculated upstream (see: add_vertical_spacing), it is possibly totally wrong
    if "dzwt" not in ds.coords:
        print(
            "Vertical spacing for vertical vel cell is approximated!!! Use with caution"
        )
        #         ds.coords['dzwt'] = grid.interp(ds['dzt'], 'Z')
        # This avoids computation. Somehow things are triggered when these fields are coordinates.
        ds["dzwt"] = grid.interp(ds["dzt"], "Z")
    # These should be less critical, but should be given nontheless
    if "dxte" not in ds.coords:
        print("Spacing for `dxte` is approximated!!! Use with caution")
        ds.coords["dxte"] = grid.interp(ds["dxt"], "X")
    if "dytn" not in ds.coords:
        print("Spacing for `dytn` is approximated!!! Use with caution")
        ds.coords["dytn"] = grid.interp(ds["dyt"], "Y")

    # Calculate thickness weighted mass transports from velocities according to MOM5 formulation
    uhrho_et, vhrho_nt, wrhot = reconstruct_hrho_trans(
        ds.u, ds.v, ds.wt, ds.dzt * rho, ds.dzu * rho, grid, rho
    )

    # Reconstruct the flux terms
    ds["o2_xflux_adv_recon"], ds["o2_yflux_adv_recon"], ds[
        "o2_zflux_adv_recon"
    ] = approximate_transport_op(
        uhrho_et, vhrho_nt, wrhot, ds["o2"], grid, boundary="extend"
    )

    # Calculate tendencies from all advective fluxes
    for suffix in ["", "_recon"]:
        ds["o2_advection_x" + suffix], ds["o2_advection_y" + suffix], ds[
            "o2_advection_z" + suffix
        ] = tend_from_fluxes(
            ds["o2_xflux_adv" + suffix],
            ds["o2_yflux_adv" + suffix],
            ds["o2_zflux_adv" + suffix],
            grid,
        )

    # Reconstruct tendency terms seperately for changes in tracer and changes in transport
    udiv, vdiv, wdiv = tend_from_fluxes(
        uhrho_et * grid._ds.dyte,
        vhrho_nt * grid._ds.dxtn,
        wrhot * grid._ds.area_t,
        grid,
    )
    ds["o2_advection_x_recon_vel"] = ds["o2"] * udiv
    ds["o2_advection_y_recon_vel"] = ds["o2"] * vdiv
    ds["o2_advection_z_recon_vel"] = ds["o2"] * wdiv

    ds["o2_advection_x_recon_tracer"] = -grid.interp(
        grid.diff(ds["o2"], "X") / grid._ds.dxte * uhrho_et, "X"
    )
    ds["o2_advection_y_recon_tracer"] = -grid.interp(
        grid.diff(ds["o2"], "Y") / grid._ds.dytn * vhrho_nt, "Y"
    )
    ds["o2_advection_z_recon_tracer"] = -grid.interp(
        grid.diff(ds["o2"], "Z") / grid._ds.dzwt * wrhot, "Z"
    )

    #     # some interpolations... might want to remap these in the future, to get more accurate estimates
    #     u = grid.interp(ds['u'], 'Y')
    #     v = grid.interp(ds['v'], 'X')
    #     wt = ds['wt']
    # #     dxte = grid.interp(ds.dxu, 'Y')
    # #     dytn = grid.interp(ds.dyu, 'X')

    #     # o2_flux reconstructed from tracer and velocity field (vel*tracer*dyt*dzt*rho)
    #     # Are the values actually interpolated or do they take the center tracer value? Read up in MOM5 manual and correct if needed.
    #     # This will need some more advanced testing...in the end we cannot really reproduce the complex advection scheme, but it is worth trying
    # #     # to get as close as possible.
    # #     ds['o2_xflux_adv_recon'] = grid.interp(ds['o2'], 'X') * u * ds.dzt * ds.dyte * rho
    # #     ds['o2_yflux_adv_recon'] = grid.interp(ds['o2'], 'Y') * v * ds.dzt * ds.dxtn * rho
    # #     ds['o2_zflux_adv_recon'] = grid.interp(ds['o2'], 'Z') * wt * ds.dxt * ds.dyt * rho

    #     # Reconstruct the advective tendencies (as (tracer* s^-1) * dzt * rho)
    #     # also not sure about the numerics here...this implements finite difference approach...which mom used for some variables but not all...

    #     ds['o2_advection_x_recon_full'] = - (grid.diff(ds['o2_xflux_adv_recon'], 'X') / ds.area_t)
    #     ds['o2_advection_x_recon_du'] = - (ds['o2'] * grid.diff(u, 'X') / ds.dxt * ds.dzt * rho)
    #     ds['o2_advection_x_recon_do2'] = - (u * grid.diff(ds['o2'], 'X') / dxte * ds.dzt * rho)

    #     ds['o2_advection_y_recon_full'] = - (grid.diff(ds['o2_yflux_adv_recon'], 'Y') / ds.area_t)
    #     ds['o2_advection_y_recon_dv'] = - (ds['o2'] * grid.diff(v, 'Y') / ds.dyt * ds.dzt * rho)
    #     ds['o2_advection_y_recon_do2'] = - (v * grid.diff(ds['o2'], 'Y') / dytn * ds.dzt * rho)

    #     ds['o2_advection_z_recon_full'] =  (grid.diff(ds['o2_zflux_adv_recon'], 'Z') / ds.area_t)
    #     ds['o2_advection_z_recon_dwt'] =  (ds['o2'] * grid.diff(wt, 'Z') * rho)
    #     ds['o2_advection_z_recon_do2'] =  (wt * grid.diff(ds['o2'], 'Z') * rho)
    return ds


################# more high level convenience functions


def add_vertical_spacing(ds):
    grid = Grid(ds)
    ds.coords["dst"] = calculate_ds(ds, dim="st")
    ds.coords["dswt"] = calculate_ds(ds, dim="sw")
    ds.coords["dzt"] = calculate_dz(ds["eta_t"], ds["ht"], ds["dst"])
    ds.coords["dzu"] = grid.min(grid.min(ds["dzt"], "X"), "Y")
    #     # Avoids suspected computation midway during the dask graph creation...Should probably raise an xarray issue about this.
    #     ds['dzt'] = calculate_dz(ds['eta_t'], ds['ht'], ds['dst'])
    #     ds['dzu'] = grid.min(grid.min(ds['dzt'], 'X'),'Y')
    #     # Lets try this as a workaround...
    #     for vv in ['dzt', 'dzu']:
    #         ds.coords[vv] = ds[vv]
    # the dzwt value is dependent on the model version (finite vol vs. engergetically consistent?; See MOM5 manual section 10.4.2)
    return ds


def split_adv_budget(ds):
    print("I think this is outdated....")
    ds = ds.copy()
    if "o2_xflux_adv" in list(ds.data_vars):
        grid = Grid(ds)
        area = ds.area_t
        div_x = -grid.diff(ds.o2_xflux_adv, "X", boundary="fill") / area
        div_y = -grid.diff(ds.o2_yflux_adv, "Y", boundary="fill") / area
        div_z = grid.diff(ds.o2_zflux_adv, "Z", boundary="fill") / area

        for data, name in zip(
            [div_x, div_y, div_z],
            ["o2_advection_%s" % a for a in ["x", "y", "z"]],
        ):
            ds[name] = data
    return ds


def budget_prep(ds, tracer="o2"):
    """Does some simplifications of the budget to compare high res and low res better
    This works for o2 but namingconventions might be different for other tracers...use with care"""
    ds = ds.copy()
    # check if budget terms are in ds
    if "o2_tendency" in list(ds.data_vars):

        # combine vertical terms
        ds["%s_diff_ver" % tracer] = (
            ds["%s_nonlocal_KPP" % tracer] + ds["%s_vdiffuse_impl" % tracer]
        )

        # combine fields from the two resolutions into comparable fields
        if "neutral_diffusion_%s" % tracer in ds.data_vars:
            ds["%s_diff_hor" % tracer] = (
                ds["neutral_diffusion_%s" % tracer]
                + ds["neutral_gm_%s" % tracer]
            )
            ds["%s_diff" % tracer] = (
                ds["%s_diff_ver" % tracer]
                + ds["%s_diff_hor" % tracer]
                + ds["%s_rivermix" % tracer]
            )
        else:
            ds["%s_diff" % tracer] = ds["%s_diff_ver" % tracer]

        ds["%s_residual" % tracer] = ds["%s_tendency" % tracer] - (
            ds["%s_advection" % tracer]
            + ds["j%s" % tracer]
            + ds["%s_diff" % tracer]
            + ds["%s_submeso" % tracer]
        )
    return ds


##################################
# Bare bones budget calculations #
##################################

##### Vertical coordinates #######
def calculate_dzstar(ds, dim="st"):
    print("Outdated function. Use `calculate_ds` instead")
    return calculate_ds(ds, dim=dim)


def calculate_ds(ds, dim="st", partial_bottom=True):
    """Creates static thickness (ds) for MOM5 tracer cells ()."""
    z = ds["%s_ocean" % dim]
    diff_dim = "%s_edges_ocean" % dim
    edges = ds[diff_dim]
    if partial_bottom:
        edges = calculate_partial_bottom(edges, ds["kmt"], ds["ht"])
        z = z * xr.ones_like(ds["kmt"])

    dz_raw = edges.diff(diff_dim).data
    # I am basing this on the edge values, assuming these are the dz values at time =0 (see MOM5_manual; 10.4.2)
    return xr.DataArray(dz_raw, coords=z.coords)


def calculate_partial_bottom(edges, bot_idx, bot_dep):
    """Creates array of cell edges, accounting for partial cells at the bottom

    Parameters
    ----------
    edges : xr.DataArray
        1D array of cell edges (usual names: `st_edges_ocean`, 'swt_edges_ocean' or similar).
    bot_idx : xr.DataArray
        2D array of number of ocean cells in the vertical (usual names: `kmt` or `kmu`)
    bot_idx : xr.DataArray
        2D array of static bottom depth (bathymetry) (usual names: `ht` or `hu`)

    Returns
    -------
    edges_partial_bottom: xr.DataArray
        3D array of cell edges, accounting for the bathymetry.
        True if successful, False otherwise.

    .. _PEP 484:
        https://www.python.org/dev/peps/pep-0484/
    """
    full_edges = edges * (xr.ones_like(bot_dep))

    edge_count = xr.ones_like(edges) * range(len(edges))
    full_count = edge_count * xr.ones_like(
        full_edges
    )  # full array of depth indicies

    # mask values larger than bottom depth
    rough_mask = full_count <= bot_idx
    full_edges_masked = full_edges.where(rough_mask)

    # replace values at bottom of last cell with actual depth
    bottom_mask = full_count == bot_idx
    bottom_correction = (
        full_edges_masked.where(bottom_mask) - bot_dep
    ).fillna(0)
    edges_partial_bottom = full_edges_masked - bottom_correction
    return edges_partial_bottom


def calculate_dz(eta, h, ds):
    """reconstructs cell thickness dz(x,y,z,t) from sea surface elevation ('eta') and the full ocean depth ('h')
    and fixed thickness levels (ds)."""
    dz = ds * (1 + (eta / h))
    return dz


# # Legacy version (deprecate soon)
# def reconstruct_thickness(eta, h, mask, dz_star, ref_array):
#     """reconstructs cell thickness dz(x,y,z,t) from sea surface elevation ('eta') and the full ocean depth ('h')
#     and fixed thickness levels (dz_star). ref has to be a full (x,y,z,t) array to carry 3d masking"""
#     print('Old version with masking implemented. Change to `calculate_dzt` and mask afterwards')
#     ones  = (ref_array*0+1).where(mask)
#     dz = ones*dz_star*(1+(eta/h))
#     return dz

######### Transports #########


def reconstruct_hrho_trans(u, v, wt, rho_dzt, rho_dzu, grid, rho):
    """Reconstruct thickness weighted mass transport in x,y,z direction.
    Units: uhrho_et/vhrho_nt [kg * s^-1 * m^-1]
           wrhot [kg * s^-1 * m^-2]"""
    uhrho_et = remap_u_2_et((u * rho_dzu), grid)
    vhrho_nt = remap_v_2_nt((v * rho_dzu), grid)
    wrhot = wt * rho
    return uhrho_et, vhrho_nt, wrhot


# I need to write a function to get the wrhot from contiunity! Otherwise the reconstructed fluxes might actually be divergent.


def approximate_transport_op(
    uflux, vflux, wflux, tracer, grid, boundary="extend"
):
    """Approximate the MOM5 transport operator. Fluxes are given as total tracer mass flux across x,y,z face
    Units: [tracer_units * kg * s^-1]"""
    tracer_xflux_adv = (
        grid.interp(tracer, "X", boundary=boundary) * uflux * grid._ds.dyte
    )
    tracer_yflux_adv = (
        grid.interp(tracer, "Y", boundary=boundary) * vflux * grid._ds.dxtn
    )
    tracer_zflux_adv = (
        grid.interp(tracer, "Z", boundary=boundary) * wflux * grid._ds.area_t
    )
    return tracer_xflux_adv, tracer_yflux_adv, tracer_zflux_adv


def horizontal_t_cell_div(u, v, w, grid, boundary="extend"):
    """Computes a simple divergence operator and total change in mass/tracer
    Units:  div_u/div_v/div_w [tracer_units * kg * s^-1 * m^-1]
    """
    div_u = -grid.diff(u, "X", boundary=boundary)
    div_v = -grid.diff(v, "Y", boundary=boundary)
    div_w = grid.diff(w, "Z", boundary=boundary)
    return div_u, div_v, div_w


def tend_from_fluxes(xflux, yflux, zflux, grid):
    """depth weighted tracer tendencies (x-y-z) from tracer fluxes (fully integrated over cell face) at the t-cell boundary"""
    x, y, z = horizontal_t_cell_div(
        xflux, yflux, zflux, grid, boundary="extend"
    )
    return x / grid._ds.area_t, y / grid._ds.area_t, z / grid._ds.area_t


# These are not modular enough...keep em for convenience but phase out at some point
def t_cell_tendency_split(uflux, vflux, wflux, tracer, grid):
    """converts thicknesswighted (horizontal) mass transport into tracer tendencies (x,y,z) using the """

    tracer_xflux_adv, tracer_yflux_adv, tracer_zflux_adv = approximate_transport_op(
        uflux, vflux, wflux, tracer, grid
    )

    uf_div, vf_div, wf_div = horizontal_t_cell_div(
        tracer_xflux_adv, tracer_yflux_adv, tracer_zflux_adv, grid
    )

    return (
        (uf_div / grid._ds.area_t),
        (vf_div / grid._ds.area_t),
        (wf_div / grid._ds.area_t),
    )


def t_cell_tendency(uflux, vflux, wflux, tracer, grid):
    """converts thicknesswighted (horizontal) mass trans"""

    u_tendency, v_tendency, w_tendency = t_cell_tendency_split(
        uflux, vflux, wflux, tracer, grid
    )
    return u_tendency + v_tendency, w_tendency


# Not needed anymore since grid.min is implemented in xgcm
# def grid_shift(da, axis, grid, boundary='extend'):
#     ref = grid.interp(da, axis, boundary=boundary)
#     return xr.DataArray(da.data, coords=ref.coords, dims=ref.coords)

# # Reconstruct dzu (minimum of sourrounding grid cells)
# def min_at_u(t_array, xdim='xt_ocean', ydim='yt_ocean'):
#     ll = t_array
#     lr = t_array.roll(**{xdim:-1})
#     ul = t_array.roll(**{ydim:-1})
#     ur = t_array.roll(**{xdim:-1, ydim:-1})
#     u_min = xr.ufuncs.minimum(ul,ur)
#     l_min = xr.ufuncs.minimum(ll,lr)
#     return xr.ufuncs.minimum(u_min, l_min)


def remap_u_2_et(u, grid, boundary="extend"):
    dyu = grid._ds.dyu
    dyte = grid._ds.dyte
    u_on_et = grid.interp((u * dyu), "Y", boundary=boundary) / dyte
    return u_on_et


def remap_v_2_nt(v, grid, boundary="extend"):
    dxu = grid._ds.dxu
    dxtn = grid._ds.dxtn
    v_on_nt = grid.interp((v * dxu), "X", boundary=boundary) / dxtn
    return v_on_nt
