#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <cmath>

/*
 * Ensure the math constants are defined one way or another
 */
#ifndef M_PI
const double M_PI = std::atan2(1.0, 1.0) * 4.0;
#endif

#ifndef M_SQRT2
const double M_SQRT2 = std::sqrt(2.0);
#endif

#ifndef M_SQRT1_2
const double M_SQRT1_2 = 1.0 / M_SQRT2;
#endif

/*
 * 2018 CODATA recommended values
 *
 * see https://physics.nist.gov/cuu/Constants/index.html
 */
const double ME = 9.1093837015e-31;         // ± 0.0000000028e-31, [kg]
const double EL = 1.602176634e-19;          // exact,              [C]
const double HPLANCK = 6.62607015e-34;      // exact,              [J / Hz]
const double HBAR = HPLANCK / (2.0 * M_PI); //                     [J / Hz]
const double KBOL = 1.380649e-23 / EL;      // exact,              [eV / K]

#define THREADS 56

/*
 * The scales to use when calculating integrals
 *
 * The variable to integrate over makes the steps of this size
 * relative to a characteristic value.
 */
// the step of an integral calculation
const double INTEGRATION_SCALE = 0.2e-3;
// the large enough value
const double ESSENTIALLY_INFINITY_SCALE = 20.0;

/*
 * Voltage and current noise for some amplifiers
 */
const double vn = (3.2e-9) * M_SQRT2;    //V/sqrt(Hz), for 2 amp * sqrt(2) AD745 amplifier
const double in = (6.9e-15) * M_SQRT1_2; //A/sqrt(Hz), for 2 amp / sqrt(2)

//		vn = (8.0e-9) * M_SQRT2;					//V/sqrt(Hz), for 2 amp * sqrt(2) OPA111 amplifier
//		in = (0.8e-15) * M_SQRT1_2;					//A/sqrt(Hz), for 2 amp / sqrt(2)

//		vn = (0.9e-9) * M_SQRT2;					//V/sqrt(Hz), for 2 amp * sqrt(2) AD797 amplifier
//		in = (2.0e-12) * M_SQRT1_2;					//A/sqrt(Hz), for 2 amp / sqrt(2)

//		vn = (1.1e-9) * M_SQRT2;					//V/sqrt(Hz), for 2 amp * sqrt(2) IFN146 amplifier
//		in = (0.3e-15) * M_SQRT1_2;					//A/sqrt(Hz), for 2 amp / sqrt(2)

//		vn = (5.1e-9) * M_SQRT2;					//V/sqrt(Hz), for 2 amp * sqrt(2) OPA1641 amplifier
//		in = (0.8e-15) * M_SQRT1_2;					//A/sqrt(Hz), for 2 amp / sqrt(2)

#endif
