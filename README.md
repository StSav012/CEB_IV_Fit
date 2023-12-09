# CEB_IV_Fit
A complicated program for automatic fitting of CEB IV-curves.
It needs an IV-curve *.txt file and a file with start parameters.

#### Parameters in `startparams.txt`:
0. incident power for all structure [pW]
1. returned power ratio, < 1
2. exponent for Te-ph
3. (unused) gap smearing
4. volume of absorber [um³]
5. (unused) volume of superconductor [um³]
6. heat exchange in normal metal [nW/(K⁵×um³)]
7. (unused) sigma of superconductor
8. critical temperature of the superconductor [K]
9. total normal resistance of all bolometers [Ohm]
10. total leakage resistance of all bolometers [Ohm]
11. transparency of the barrier
12. depairing energy
13. coefficient for Andreev current
14. normal resistance of 1 absorber [Ohm]
15. number of bolometers in series
16. number of bolometers in parallel
17. phonon temperature [K]
18. upper voltage limit relative to the gap
19. lower voltage limit relative to the gap
20. voltage step relative to the gap
