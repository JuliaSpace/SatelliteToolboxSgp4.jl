## Description #############################################################################
#
# Precompilation.
#
############################################################################################

import PrecompileTools

PrecompileTools.@compile_workload begin

    # == Orbit Propagation =================================================================

    # We select multiple TLEs to make sure all functions are precompiled.
    tles = tles"""
    #                       # TEME example
    1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753
    2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667
    #   MOLNIYA 2-14        # 12h resonant ecc in 0.65 to 0.7 range
    1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813
    2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656
    #   MINOTAUR R/B        # Sub-orbital case - Decayed 2005-11-29
    #                       #(perigee = -51km), lost in 50 minutes
    1 28872U 05037B   05333.02012661  .25992681  00000-0  24476-3 0  1534
    2 28872  96.4736 157.9986 0303955 244.0492 110.6523 16.46015938 10708
    """

    for tle in tles
        sgp4(10.0, tle; sgp4c = SGP4C_WGS84)
        sgp4(10.0f0, tle; sgp4c = Sgp4Constants{Float32}(SGP4C_WGS84))
    end

    # Exercise the compact and the rich representations of the propagator.
    sgp4d = sgp4_init(first(tles))
    show(IOBuffer(), sgp4d)
    show(IOBuffer(), MIME("text/plain"), sgp4d)

    # == Orbit Propagation Using OMM =======================================================

    # Exercise the initialization from an Orbit Mean-Elements Message.
    omm = parse_omm(
        """
        <?xml version="1.0" encoding="utf-8"?>
        <ndm><omm id="CCSDS_OMM_VERS" version="3.0">
        <header><CREATION_DATE>2025-12-30T23:36:37</CREATION_DATE><ORIGINATOR>18 SPCS</ORIGINATOR></header>
        <body><segment>
        <metadata><OBJECT_NAME>AMAZONIA 1</OBJECT_NAME><OBJECT_ID>2021-015A</OBJECT_ID><CENTER_NAME>EARTH</CENTER_NAME><REF_FRAME>TEME</REF_FRAME><TIME_SYSTEM>UTC</TIME_SYSTEM><MEAN_ELEMENT_THEORY>SGP4</MEAN_ELEMENT_THEORY></metadata>
        <data>
        <meanElements><EPOCH>2025-12-30T18:12:04.533984</EPOCH><MEAN_MOTION>14.40772474</MEAN_MOTION><ECCENTRICITY>0.00011240</ECCENTRICITY><INCLINATION>98.3721</INCLINATION><RA_OF_ASC_NODE>75.0877</RA_OF_ASC_NODE><ARG_OF_PERICENTER>97.3772</ARG_OF_PERICENTER><MEAN_ANOMALY>262.7545</MEAN_ANOMALY></meanElements>
        <tleParameters><EPHEMERIS_TYPE>0</EPHEMERIS_TYPE><CLASSIFICATION_TYPE>U</CLASSIFICATION_TYPE><NORAD_CAT_ID>47699</NORAD_CAT_ID><ELEMENT_SET_NO>999</ELEMENT_SET_NO><REV_AT_EPOCH>25439</REV_AT_EPOCH><BSTAR>0.00015330000000</BSTAR><MEAN_MOTION_DOT>0.00000447</MEAN_MOTION_DOT><MEAN_MOTION_DDOT>0.0000000000000</MEAN_MOTION_DDOT></tleParameters>
        </data>
        </segment></body>
        </omm></ndm>
        """,
    )

    sgp4(10.0, omm; sgp4c = SGP4C_WGS84)
    sgp4(10.0f0, omm; sgp4c = Sgp4Constants{Float32}(SGP4C_WGS84))

    # == TLE Fitting =======================================================================

    tle_input = tle"""
        AMAZONIA 1
        1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
        2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652"""

    vjd     = [2.46002818657856e6]
    vr_teme = [@SVector [-6792.402703741442, 2192.6458461287293, 0.18851758695295118]]
    vv_teme = [@SVector [0.3445760107690598, 1.0395135806993514, 7.393686131436984]]

    redirect_stdout(devnull) do
        # Fit a TLE without and with an initial guess.
        fit_sgp4_mean_elements(
            TLE, vjd, vr_teme, vv_teme; estimate_bstar = false, max_iterations = 1
        )

        fit_sgp4_mean_elements(
            TLE,
            vjd,
            vr_teme,
            vv_teme;
            estimate_bstar = false,
            initial_guess  = tle_input,
            max_iterations = 1,
        )

        # Fit an OMM with the covariance section.
        fit_sgp4_mean_elements(
            OrbitMeanElementsMessage,
            vjd,
            vr_teme,
            vv_teme;
            estimate_bstar = false,
            max_iterations = 1,
        )
    end

    # == Mean Elements Epoch Update ========================================================

    tle = tle"""
        AMAZONIA 1
        1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
        2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652"""

    redirect_stdout(devnull) do
        update_sgp4_mean_elements_epoch(tle, 2.46002818657856e6 + 1; max_iterations = 1)
        update_sgp4_mean_elements_epoch(omm, 2.46104025838581e6 + 1; max_iterations = 1)
    end

    # == Structure Copying =================================================================

    sgp4d = sgp4_init(
        tle"""
        AMAZONIA 1
        1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
        2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652""";
        sgp4c = SGP4C_WGS84,
    )
    copy(sgp4d)

    sgp4d_f32 = sgp4_init(
        tle"""
        AMAZONIA 1
        1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
        2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652""";
        sgp4c = Sgp4Constants{Float32}(SGP4C_WGS84),
    )
    copy(sgp4d_f32)
end
