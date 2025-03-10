from zero_point import zpt as _zpt


def zero_point_correction(sample = None, *zargs):
    """
    Computes the parallax zero point correction for a given sample of stars,
    following the "recipe" given in the dited paper
    > `Lindergren, et al., A&A 649, A4 (2021)`

    The user must provide either a `grasp sample` object, containing the parameters:
    - `phot_g_mean_mag`
    - `nu_eff_used_in_astrometry`
    - `pseudocolour`
    - `ecl_lat`
    - `astrometric_params_solved`
    
    or, alternatively, these quantities separately, as singular arguments.

    NOTE:
    -----
    for 5-p solutions (ra-dec-parallax-pmra-pmdec), the field 
    `astrometric_params_solved` equals 31 and the `pseudocolour` variable can take
    any arbitrary values (even None). On the other hand, for 6-p solutions (ra-dec
    -parallax-pmra-pmdec-pseudocolour), the field `astrometric_params_solved` equals
    95 and the `nu_eff_used_in_astrometry` variable can take any arbitrary values
    (even None). For 2-p solutions (ra-dec), the field `astrometric_params_solved`
    equals 3 and the zero point correction cannot be computed.
    
    Parameters
    ----------
    sample : grasp._utility.sample.Sample, optional
        The sample of stars to compute the zero point correction for, containing
        the needed parameters.
    *zargs : list, optional
        The parameters needed to compute the zero point correction, in the following
        order:
        - phot_g_mean_mag
        - nu_eff_used_in_astrometry
        - pseudocolour
        - ecl_lat
        - astrometric_params_solved

    Returns
    -------
    numpy.ndarray or Sample object
        The zero point correction for the given sample of stars. If the sample
        object was provided for the computation, this will automatically add the
        `zero_point_correction` column in the sample. If not, then the array of 
        the computed values are returned.
    """
    _zpt.load_tables()
    print("""
Using the "Zero Point Correction" tool from Pau Ramos
(Lindergren, et al., A&A 649, A4 (2021)"""
    )
    if sample is None:
        zpc = _zpt.zero_point_correction(*zargs)
        return zpc
    else:
        from pandas import DataFrame
        from grasp._utility.sample import Sample as _Sample
        if isinstance(sample, _Sample):
            zpc = _zpt.zpt_wrapper(sample.to_pandas())
            sample.sample.add_column("zero_point_correction", zpc)
        elif isinstance(sample, DataFrame):
            zpc = _zpt.zpt_wrapper(sample)
            sample["zero_point_correction"] = zpc
        else:
            raise ValueError(
                "The sample provided is not a valid Sample object."
            )
        return sample
