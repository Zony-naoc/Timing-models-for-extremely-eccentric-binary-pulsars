# Timing-models-for-extremely-eccentric-binary-pulsars
We present new timing models in tempo2 format, adapted from DD and DDGR models for highly eccentric binary pulsars (\sqrt(1-ecc*ecc)<<1). The new timing models avoid high correlations between orbital parameters x, ecc and om, and seperate the periastron advance into relativistic and non-relativistic parts to measure them separately.

# model describtions
The new models have these key features:
1. In these timing models, we have replaced orbital parameters x, ecc and om with zeta1=xsin(om)/ecc, zeta2=xcos(om)\sqrt(1-ecc*ecc) and zeta3=\sqrt(1-ecc*ecc), which provide much better describition of the Keplarian motion for exteremly eccentric binary pulsars.
2. For exteremly eccentric binary pulsars, their periastron advance rates have lage variations in one orbital period. Fitting for such variations and the long-term changes of om allow the independent measurement of the relativistic periastron advance rate and the total periastron advance rate, providing constraints on the non-relativistic contributions on the periastron advance.

# usage
1. In DDe/DDe+ model, if XOMODT parameter is included, the timing model will fit for both the relativistic periastron advance rate and the total periastron advance rate. Otherwise, the model assumes that the periastron advance is fully relativistic the same as the DD model.
2. In DDeGR+ model, the timing model will always fit for both the total mass (which corresponds to the relativistic periastron advance rate) and the total periastron advance rate. Besides that, sini can be used in this model to replace mc (optional).





