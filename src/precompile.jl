# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Precompilation.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import SnoopPrecompile

SnoopPrecompile.@precompile_all_calls begin

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
        sgp4(10.0,   tle; sgp4c = sgp4c_wgs84)
        sgp4(10.0f0, tle; sgp4c = sgp4c_wgs84_f32)
    end
end
