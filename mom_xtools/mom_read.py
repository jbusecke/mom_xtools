### Tools to read in MOM5 data cleanly
# Adapted from gist of rabernat

from glob import glob
import os
import xarray as xr

# How to deal with multiple filenames vs only one?

# time inconsistency ( that has to be dealt with after this step)
def time_add_refyear(ds, timedim="time", refyear=2000):
    # I am still not partiularly happy with this...I thought xarray can now deal with years from 0...
    ds_clean = ds.copy()
    ds = ds.copy()
    # Fix the time axis (I added 1900 years, since otherwise the stupid panda
    # indexing does not work)
    new_time_units = "days since %s-01-01 00:00:00" % str(refyear + 1)
    ds.time.attrs["units"] = new_time_units
    ds = xr.decode_cf(ds)
    ds.attrs["refyear_shift"] = refyear
    # I could just compute the time and then past it into a fresh (chunked copy of the ds)
    ds_clean["time"] = ds["time"]
    return ds_clean


def process_coords(
    ds,
    concat_dim="time",
    drop=True,
    hard=True,
    extra_coord_vars=[
        "time_bounds",
        "average_T1",
        "average_T2",
        "average_DT",
        "nv",
    ],
):
    """Preprocessor function to drop all non-dim coords (all coords if hard is true)
    to speed concatenation."""

    # Remove unneeded attrs from ds (Some of these might be usefull to keep,
    # but I am worried the ones that differ will slow down reading)
    remove_attrs = ["filename", "title", "grid_type", "grid_tile"]
    for ra in remove_attrs:
        try:
            del ds.attrs[ra]
        except ValueError:
            pass

    # This detects coords that are declared as variables (usually not a problem at GFDL...)
    coord_vars = [v for v in ds.data_vars if concat_dim not in ds[v].dims]
    non_dim_coords = [v for v in ds.coords if v not in ds.dims]

    target_coords = coord_vars + non_dim_coords

    for ecv in extra_coord_vars:
        if ecv in ds:
            target_coords += extra_coord_vars

    if hard:
        target_coords += [d for d in ds.dims if d != concat_dim]

    # make sure all target_coords are in the dataset
    target_coords = [a for a in target_coords if a in list(ds.keys())]

    if drop:
        return ds.drop(target_coords)
    else:
        return ds.set_coords(target_coords)


def open_coord_ds(fname):
    """open a 'reference' dataset that is used to overwrite coords on main dataset."""
    dsc = xr.open_dataset(fname, decode_cf=False, chunks={"time": 1})
    varnames = [v for v in dsc.data_vars if "time" in dsc[v].dims]
    dsc = dsc.drop(varnames)
    dsc = dsc.drop("time")
    try:
        dsc = dsc.drop("time_bounds")
    except ValueError:
        pass
    dsc = dsc.set_coords(dsc.data_vars)
    dsc = xr.decode_cf(dsc, decode_times=False)
    return dsc


def clean(ds):
    # get rid of coords attr, which causes problems with serialization (not sure if all need to be removed
    # but coordinates, definitely caused problems)
    remove_attrs = [
        "time_avg_info",
        "cell_methods",
        "valid_range",
        "coordinates",
    ]
    for ra in remove_attrs:
        for var in ds.variables:
            try:
                del ds[var].attrs[ra]
            except KeyError:
                pass
            # don't need original encoding cluttering things up
            ds[var].encoding.clear()
    return ds


def merge_and_clean(ds_vars, ds_coords):
    ds = xr.merge([ds_vars, ds_coords])
    ds.attrs.update(ds_coords.attrs)
    ds = clean(ds)
    return ds


def parse_filenames_from_dir(ddir, debug=False):
    """detects filenames for different variables in ddir"""
    flist = sorted(os.listdir(ddir))
    # remove files or dirs that do not end in '.nc'
    flist = [f for f in flist if ".nc" in f]
    # Get a unique sorted list of years
    step_list = sorted(list(set([a.split(".")[0] for a in flist])))
    # This assumes that all variables are given for the last time step.
    var_list = [b.split(".")[-2] for b in flist if step_list[-1] in b]
    return var_list


def data_vars_consistency(path1, path2, **kwargs):
    """returns list of variables that are found in the file at path1 but not path2
    and vice versa"""
    ds1 = xr.open_mfdataset(
        path1, decode_cf=False, chunks={"time": 1}, **kwargs
    )
    ds2 = xr.open_mfdataset(
        path2, decode_cf=False, chunks={"time": 1}, **kwargs
    )
    return list(set(ds1.variables) ^ set(ds2.variables))


def open_mom5_single_var(
    ddir,
    var,
    years=None,
    yearfmt="%040101",
    filenamefmt=None,
    drop_inconsistent=True,
    add_coords=True,
    add_reftime=True,
    switch_year_name=False,
    gridfile=None,
    parallel=False,
    debug=False,
    **kwargs
):
    """Open files for a single filename(`var`) in the directory `ddir`. Inconsistent variables (not written in the last and first file) are automaticlly dropped. Adding coordinates from a reference file (first file in the list) can be supressed by setting `add_coords` to false.
    Document the filefmt stuff...
    This is an example for files found in /archive/Julius.Busecke/CM2.6/CM2.6_A_Control-1860_V03/pp/ocean/annual_1yr...

    fmt = "{var}.{year}.ann.nc"

    Eventually replace all other logic with that...
    """
    if years is None:
        all_fnames = sorted(glob(os.path.join(ddir, "*%s*.nc" % var)))
    else:
        if filenamefmt:
            # this shit just got fancy!
            filled_years = [str(year).zfill(4) for year in years]
            pathspecs = [
                os.path.join(ddir, filenamefmt.format(var=var, year=year))
                for year in filled_years
            ]
        else:
            if switch_year_name:
                pathspecs = [
                    os.path.join(ddir, "%s*%04i*.nc" % (var, year))
                    for year in years
                ]
            else:
                pathspecs = [
                    os.path.join(ddir, "%04i*%s*.nc" % (year, var))
                    for year in years
                ]

        flists = [glob(pathspec) for pathspec in pathspecs]
        if debug:
            print(flists)
            print(pathspecs[0])
        all_fnames = sorted([j for i in flists for j in i])
        # if some errors arise, ill have to tinker with the yearfmt
    if not all_fnames:
        raise RuntimeError("No files found for Variable(%s)" % (var))
    drop_vars = []
    print
    if drop_inconsistent:
        incon_var = data_vars_consistency(
            all_fnames[0], all_fnames[-1], **kwargs
        )
        if incon_var:
            print(
                "%s - Variable inconsistencies between files found for (%s). \
            Variables are dropped. Use `drop_inconsistent=False` to disable."
                % (var, incon_var)
            )
        drop_vars += incon_var
    else:
        print(
            "This might take a long time to read the files if incosistencies are present."
        )

    # There might be a more pythonic way to do this.
    preprocess = kwargs.pop("preprocess", None)
    if preprocess:

        def pp_func(ds):
            return preprocess(process_coords(ds))

    else:

        def pp_func(ds):
            return process_coords(ds)

    # The reftime is a legacy addition, not most of the things are working nicely with the ncdatetime in xarray
    dsc = open_coord_ds(all_fnames[0])
    if add_reftime:
        decode_times = False
    else:
        decode_times = True

    ds = xr.open_mfdataset(
        all_fnames,
        drop_variables=drop_vars,
        decode_times=decode_times,
        #                            autoclose=True,
        concat_dim="time",
        preprocess=pp_func,
        chunks={"time": 1},
        parallel=parallel,
        **kwargs
    )

    if add_reftime:
        print("adding offset of 2000 to time")
        print("decode_times: %s" % decode_times)
        ds = time_add_refyear(ds, timedim="time", refyear=2000)

    if add_coords:
        ds_out = merge_and_clean(ds, dsc)
    else:
        ds_out = clean(ds)

    if gridfile:
        ds_grid = xr.open_mfdataset(gridfile)
        ds_out.update(ds_grid)

    return ds_out


def open_mom5_all_vars(ddir, varlist=None, gridfile=None, **kwargs):
    """auto open all files in `ddir`, if not otherwise specified in `varlist` other kwargs are passed to `open_mom5_single_var` """
    if varlist is None:
        vars = parse_filenames_from_dir(ddir)
    else:
        vars = varlist
    datasets = []
    for var in vars:
        ds = open_mom5_single_var(ddir, var, **kwargs)
        datasets.append(ds)

    ds_out = xr.merge(datasets)

    if gridfile:
        ds_grid = xr.open_mfdataset(gridfile)
        ds_out.update(ds_grid)

    return ds_out


def open_mom5_CM_ESM(
    basedir,
    timespec="monthly_1yr",
    subfolder="av",
    varfolderlist=None,
    force_varlist=False,
    reftime=None,
    years=None,
    gridfile=None,
    **kwargs
):

    """Wrap and merge function to read all GFDL MOM5 based CM/ESM files stored in GFDL filestructure.


    Parameters
    ----------

    basedir : path
        typically a ./pp/ directory

    timespec: str
        specifies the subdirectory to read (e.g. 'monthly_1yr', 'annual_1yr')

    subfolder: str
        specifies the subfolder for each subfolder (typically 'av' or 'ts')

    varfolderlist: {None, list}
        specifies the folders to be read (default is reading all existing ones)

    force_varlist: bool
        `force_varlist` is currently necessary for ESM2.6 since the files are saved in different naming convention.
<<<<<<< HEAD
        This is not well explained...

=======

>>>>>>> 291f420a89e3071b23e39907fe7576cb18b515a8
    reftime: {None, str}
        Physical and Biogeochemical output can have different timestamps (leading to problems when merging).
        `reftime` uses the time dimension of the specified fields and overwrites all others. If folder
        corresponding to `reftime is not found, defaults to first dataset.

    years: {range, np.array}

    gridfile: {None, path}
        `gridfile' defines a static file which overwrites the existing coordinates, if not specified fields are returned without coordinates...
        #TODO return them with coordinates if gridfile is not specified...this could lead to merge problems since the values are not exactly the same
        (similar as the time).


    Returns
    -------
    Dataset
    """

    if varfolderlist is None:
        varfolderlist = os.listdir(basedir)

    datasets = []
    for va in varfolderlist:
        if force_varlist:
            varlist = [va]
        else:
            varlist = None

        path = os.path.join(basedir, va, subfolder, timespec)
        if os.path.isdir(path):
            print("reading %s" % va)
            print(varlist)
            ds = open_mom5_all_vars(
                path, varlist, add_coords=False, years=years, **kwargs
            )
            print("timesteps: %i" % len(ds.time))
            datasets.append(ds)
        else:
            if os.path.isdir(os.path.join(basedir, va, subfolder)):
                print('No time resolution "%s" found for %s' % (timespec, va))
            elif os.path.isdir(os.path.join(basedir, va)):
                print("No `subfolder`(%a) found for %s" % (subfolder, va))
            else:
                print(
                    "something is wrong with the path for %s (%s)" % (va, path)
                )

    # Replace time from ref dataset. This is necessary since the BGC and Physics have different time
    # stamps, leading to broadcasting of enormous merged datsets
    if reftime is not None:
        if not all(len(x.time) == len(datasets[0].time) for x in datasets):
            raise RuntimeError(
                "Reference Time replacement can only be done when all datasets\
            have the same length. Current lengths: %s"
                % [len(d.time) for d in datasets]
            )
        if reftime not in varfolderlist:
            ref = 0
        else:
            ref = varfolderlist.index(reftime)
        for dd in range(len(datasets)):
            datasets[dd].time.data = datasets[ref].time.data

    ds_out = xr.merge(datasets)
    if gridfile:
        ds_grid = xr.open_mfdataset(gridfile)
        ds_out.update(ds_grid)
    return ds_out


def parse_xgcm_attributes(
    ds,
    xc="xt_ocean",
    xg="xu_ocean",
    yc="yt_ocean",
    yg="yu_ocean",
    zc="st_ocean",
    zg="sw_ocean",
):
    """ Adds axis attributes needed for xgcm to recognize the grid"""
    if (xc is not None) and (xc in ds.dims):
        ds[xc] = ds[xc].assign_attrs(axis="X")

    if (xg is not None) and (xg in ds.dims):
        ds[xg] = ds[xg].assign_attrs(axis="X")
        ds[xg] = ds[xg].assign_attrs(c_grid_axis_shift=0.5)

    if (yc is not None) and (yc in ds.dims):
        ds[yc] = ds[yc].assign_attrs(axis="Y")

    if (yg is not None) and (yg in ds.dims):
        ds[yg] = ds[yg].assign_attrs(axis="Y")
        ds[yg] = ds[yg].assign_attrs(c_grid_axis_shift=0.5)

    if (zc is not None) and (zc in ds.dims):
        ds[zc] = ds[zc].assign_attrs(axis="Z")

    if (zg is not None) and (zg in ds.dims):
        ds[zg] = ds[zg].assign_attrs(axis="Z")
        ds[zg] = ds[zg].assign_attrs(c_grid_axis_shift=0.5)
    return ds
