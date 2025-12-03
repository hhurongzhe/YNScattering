# YNScattering

Hyperon-Nucleon (YN) Coupled-Channel Scattering Calculator, written in C++.

YN scattering observables which can be computed:
phase shift, colomb-matched phase shift, differential cross section, total cross section.

Chiral hyperon-nucleon interaction at leading-order is used, see YN.ini.
Other types of interaction can also be easily implemented.

Phase shifts are calculated by solving coupled-channel LS equation with matrix inversion method in momentum space.
Two-body states are defined by SYM-LSJ convention, under which matrix elements are calculated.

The evaluation of interaction matrix elements are by aPWD method extened to YN interaction case.

The code is compiled with xmake, see contents in xmake.lua.
Of course a naive Makefile would be possible, too.

If you have any futher needs or questions, just contact me: rongzhe_hu@pku.edu.cn
