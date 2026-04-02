MATLAB scripts specifics

ts_mat: it generates series resistance, input resistance and capacitance values from test pulse recordings

autoinh_vc_ana: it analyses two auto-evoked inhibition current recordings by averaging sweeps to produce mean current traces and calculating the response charge (pC) for each.
It then subtracts the second from the first to generate a difference trace and the corresponding difference charge (area_pC_diff).

autoinh_vc_ana_modif: it analyses an auto-evoked inhibition current recording by averaging sweeps to produce a mean current trace and calculating the response charge (pC). It outputs the mean trace and the integrated charge over a defined post-stimulus time window.

AP_10ms_mat: This script takes recordings from short current injections, finds a spike, and measures some basic properties of that spike (like when it starts, how big it is, etc.). But for this paper, we mostly used it to turn each spike into a phase plane plot (dV/dt vs voltage) and look at its shape. Following the approach used in the eLife paper you linked, the key thing was that axon-bearing neurons show a biphasic phase plane (a clear “kink”/two-stage rise), while anaxonic neurons show a smooth, monophasic curve