# -bot-ptc
Analysis code of the flat collected at SLAC for the LSST camera tests.

This code measures variances and cvariances of flats in order to constrain
the electrostatics inside the sensor, because the later determines the
details of the brighter-fatter effect, something that affects pixel sensors
and in particular CCDs.

This code depends on bfptc, which was used to produce the analysis published in 
 https://arxiv.org/abs/1905.08677
