/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CURVE_H_
#define CURVE_H_

typedef struct /* fuer Kurven (Stuetzstellen) */
{
	double punkt;
	double wert;
} StuetzStellen;

typedef struct /* fuer Kurven (Kurven) */
{
	long anz_stuetzstellen;
	/* Anzahl der Stuetzstellen */
	StuetzStellen* stuetzstellen;
	/* Feld mit den eingelesenen Stuetzstellen */
} Kurven;

extern double GetCurveValue(int kurve, int methode, double punkt, int* gueltig);
extern double GetCurveValueInverse(int kurve,
                                   int methode,
                                   double wert,
                                   int* gueltig);
extern double GetCurveDerivative(int kurve,
                                 int methode,
                                 double punkt,
                                 int* gueltig);
extern double GetCurveInverseDerivative(int kurve,
                                        int methode,
                                        double wert,
                                        int* gueltig);
extern Kurven* kurven; /* Feld mit Kurven */
extern int anz_kurven; /* Anzahl der Kurven */


#endif // CURVE_H_
