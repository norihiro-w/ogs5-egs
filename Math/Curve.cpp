/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Curve.h"

#include <cmath>

#include "makros.h"
#include "display.h"

#include "MathTools.h"

Kurven* kurven = NULL;
int anz_kurven = 0;


/**************************************************************************
   ROCKFLOW - Funktion: GetCurveValue

   Aufgabe:
   Liefert Wert aus einer Kurve fuer angegebenen Punkt.
   Liegt der Punkt ausserhalb des durch die Kurve beschriebenen
   Bereichs, wird der letzte bzw. erste Wert zurueckgeliefert und
   der Flag fuer die Gueltigkeit auf 0 gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int kurve: Kurvennummer, >= 0
   E int methode  : Interpolationsmethode
   E double punkt: Punkt
   R int gueltig: Flag fuer die Gueltigkeit des zurueckgelieferten Wertes

   Ergebnis:
   s.o.

   Programmaenderungen:
   Basiert auf "GetRBZValue"

   04/1999     C.Thorenz     Gueltigkeit und Methode eingefuehrt
   09/2000     C.Thorenz     Fehlermeldung bei falscher Kurvennummer

**************************************************************************/
double GetCurveValue(int kurve, int methode, double punkt, int* gueltig)
{
	static long anz;
	register long i;
	static StuetzStellen* s;

	if (kurve == 0)
	{
		*gueltig = 1;
		return 1.0;
	}

#ifdef ERROR_CONTROL
	if ((kurve < 0) || (kurve >= anz_kurven))
	{
		DisplayMsgLn("");
		DisplayMsg("PANIC! Curve ");
		DisplayLong(kurve);
		DisplayMsgLn(" is requested but not defined!");
		abort();
	}
#endif

	anz = kurven[kurve].anz_stuetzstellen;
	s = kurven[kurve].stuetzstellen;
	*gueltig = 1;

	i = 1l;
	while (punkt > s[i].punkt)
		i++;
	//
	// Check curve bounds
	if (punkt < s[0].punkt)
	{
		*gueltig = 0;
		return s[0].wert;
	}
	if (punkt > s[anz - 1l].punkt)
	{
		*gueltig = 0;
		return s[anz - 1l].wert;
	}
	//
	// Otherwise, get interpolated value
	switch (methode)
	{
		default:
			ScreenMessage("ERROR: GetCurveValue() --> Invalid curve.\n");
			return 0.0;
		//
		case 0:  // Linear Interpolation
			return s[i - 1].wert +
			       (s[i].wert - s[i - 1l].wert) /
			           (s[i].punkt - s[i - 1l].punkt) *
			           (punkt - s[i - 1l].punkt);
		//
		case 1:                    // Piece wise constant
			return s[i - 1].wert;  // BG changed from i to i-1, 2011
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: GetCurveValueInverse

   Aufgabe:
   Liefert zu einem Wert aus einer Kurve den zugehoerigen Punkt.
   Liegt der Punkt ausserhalb des durch die Kurve beschriebenen
   Bereichs, wird der letzte bzw. erste Wert zurueckgeliefert und
   der Flag fuer die Gueltigkeit auf 0 gesetzt.

   Kurven muessen streng monoton fallend oder steigend sein.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int kurve: Kurvennummer, >= 0
   E int methode  : Interpolationsmethode
   E double wert: Wert
   R int gueltig: Flag fuer die Gueltigkeit des zurueckgelieferten Wertes

   Ergebnis:
   s.o.

   Programmaenderungen:
   Basiert auf "GetRBZValue"

   12/1999     C. Thorenz Erste Version

**************************************************************************/
double GetCurveValueInverse(int kurve, int methode, double wert, int* gueltig)
{
	static long anz;
	register long i;
	static StuetzStellen* s;

#ifdef ERROR_CONTROL
	if ((kurve < 0) || (kurve >= anz_kurven))
	{
		DisplayMsgLn("");
		DisplayMsg("PANIC! Curve ");
		DisplayLong(kurve);
		DisplayMsgLn(" is requested but not defined!");
		abort();
	}
#endif

	anz = kurven[kurve].anz_stuetzstellen;
	s = kurven[kurve].stuetzstellen;
	*gueltig = 1;
	i = 1l;

	if (s[0].wert < s[anz - 1l].wert)
	{
		/* Monoton steigend */
		if (wert < s[0].wert)
		{
			*gueltig = 0;
			return s[0].punkt;
		}
		if (wert > s[anz - 1].wert)
		{
			*gueltig = 0;
			return s[anz - 1].punkt;
		}
		/* Suchen der Stuetzstelle. Vorraussetzung: Zeitpunkte aufsteigend
		 * geordnet */
		while (wert > s[i].wert)
			i++;
	}
	else
	{
		/* Monoton fallend */
		if (wert > s[0].wert)
		{
			*gueltig = 0;
			return s[0].punkt;
		}
		if (wert < s[anz - 1].wert)
		{
			*gueltig = 0;
			return s[anz - 1].punkt;
		}
		/* Suchen der Stuetzstelle. Vorraussetzung: Zeitpunkte aufsteigend
		 * geordnet */
		while (wert < s[i].wert)
			i++;
	}

	switch (methode)
	{
		default:
			ScreenMessage("ERROR: GetCurveValue() --> Invalid curve.\n");
			return 0.0;
		//
		case 0:  // Lineare Interpolation
			return s[i - 1].punkt +
			       (s[i].punkt - s[i - 1l].punkt) /
			           (s[i].wert - s[i - 1l].wert) * (wert - s[i - 1l].wert);
		//
		case 1:  // Piece wise constant
			return s[i].punkt;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: GetCurveDerivative

   Aufgabe:
   Liefert die Ableitung zu einem Punkt auf einer Kurve.
   Liegt der Punkt ausserhalb des durch die Kurve beschriebenen
   Bereichs, wird die letzte bzw. erste moegliche Ableitung
   zurueckgeliefert und der Flag fuer die Gueltigkeit auf 0 gesetzt.

   Die Ableitung kann ueber mehrere Methoden bestimmt werden:
     0: Stueckweise konstant
   1: Gleitend

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int kurve: Kurvennummer, >= 0
   E int methode  : Ableitungsmethode
   E double punkt: Punkt
   R int gueltig: Flag fuer die Gueltigkeit des zurueckgelieferten Wertes

   Ergebnis:
   s.o.

   Programmaenderungen:

   3/2002   C. Thorenz Erste Version
**************************************************************************/
double GetCurveDerivative(int kurve, int methode, double punkt, int* gueltig)
{
	static long anz;
	register long i;
	static StuetzStellen* s;
	static double w, s1, s2;

	if (kurve == 0)
	{
		*gueltig = 1;
		return 1.0;
	}

#ifdef ERROR_CONTROL
	if ((kurve < 0) || (kurve >= anz_kurven))
	{
		DisplayMsgLn("");
		DisplayMsg("PANIC! Curve ");
		DisplayLong(kurve);
		DisplayMsgLn(" is requested but not defined!");
		abort();
	}
#endif

	anz = kurven[kurve].anz_stuetzstellen;
	s = kurven[kurve].stuetzstellen;
	*gueltig = 1;
	i = 1l;

	if (punkt < s[0].punkt)
	{
		*gueltig = 0;
		i = 1;
		punkt = s[0].punkt;
	}
	else if (punkt > s[anz - 1l].punkt)
	{
		*gueltig = 0;
		i = anz - 1;
		punkt = s[anz - 1].punkt;
	}
	else
		/* Suchen der Stuetzstelle. Vorraussetzung: Zeitpunkte aufsteigend
		 * geordnet */
		while (punkt > s[i].punkt)
			i++;

	switch (methode)
	{
		default:
		case 0:
			/* Stueckweise konstant */
			if (fabs(s[i].punkt - s[i - 1].punkt) > DBL_MIN)
				return (s[i].wert - s[i - 1].wert) /
				       (s[i].punkt - s[i - 1].punkt);
			else
				return Signum(s[i + 1].wert - s[i].wert) / DBL_EPSILON;
		case 1:
			/* Gleitend */
			if ((i > 1) && (i < anz - 2))
			{
				s1 = (0.5 * s[i].wert - 0.5 * s[i - 2].wert) /
				     (0.5 * s[i].punkt - 0.5 * s[i - 2].punkt);

				s2 = (0.5 * s[i + 1].wert - 0.5 * s[i - 1].wert) /
				     (0.5 * s[i + 1].punkt - 0.5 * s[i - 1].punkt);

				w = (punkt - s[i - 1].punkt) / (s[i].punkt - s[i - 1].punkt);

				return (1. - w) * s1 + w * s2;
			}
			else
			{
				/* Stueckweise konstant */
				if (fabs(s[i].punkt - s[i - 1].punkt) > DBL_MIN)
					return (s[i].wert - s[i - 1].wert) /
					       (s[i].punkt - s[i - 1].punkt);
				else
					return Signum(s[i + 1].wert - s[i].wert) / DBL_EPSILON;
			}
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: GetCurveValueInverse

   Aufgabe:
   Liefert zu einem Wert aus einer Kurve den zugehoerigen Punkt.
   Liegt der Punkt ausserhalb des durch die Kurve beschriebenen
   Bereichs, wird der letzte bzw. erste Wert zurueckgeliefert und
   der Flag fuer die Gueltigkeit auf 0 gesetzt.

   Kurven muessen streng monoton fallend oder steigend sein.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int kurve: Kurvennummer, >= 0
   E int methode  : Interpolationsmethode
   E double wert: Wert
   R int gueltig: Flag fuer die Gueltigkeit des zurueckgelieferten Wertes

   Ergebnis:
   s.o.

   Programmaenderungen:
   Basiert auf "GetRBZValue"

   03/2012 JT: Implementation. Based on other curve methods existing

**************************************************************************/
double GetCurveInverseDerivative(int kurve,
                                 int methode,
                                 double wert,
                                 int* gueltig)
{
	static long anz;
	register long i;
	static StuetzStellen* s;
	static double w, s1, s2;

	if (kurve == 0)
	{
		*gueltig = 1;
		return 1.0;
	}

#ifdef ERROR_CONTROL
	if ((kurve < 0) || (kurve >= anz_kurven))
	{
		DisplayMsgLn("");
		DisplayMsg("PANIC! Curve ");
		DisplayLong(kurve);
		DisplayMsgLn(" is requested but not defined!");
		abort();
	}
#endif

	anz = kurven[kurve].anz_stuetzstellen;
	s = kurven[kurve].stuetzstellen;
	*gueltig = 1;
	i = 1l;

	if (s[0].wert < s[anz - 1l].wert)
	{
		/* Monoton steigend */
		if (wert < s[0].wert)
		{
			*gueltig = 0;
			i = 1;
			wert = s[0].wert;
		}
		else if (wert > s[anz - 1].wert)
		{
			*gueltig = 0;
			i = anz - 1;
			wert = s[anz - 1].wert;
		}
		else
		{ /* Suchen der Stuetzstelle. Vorraussetzung: Zeitpunkte aufsteigend
		     geordnet */
			while (wert > s[i].wert)
				i++;
		}
	}
	else
	{
		/* Monoton fallend */
		if (wert > s[0].wert)
		{
			*gueltig = 0;
			i = 1;
			wert = s[0].wert;
		}
		else if (wert < s[anz - 1].wert)
		{
			*gueltig = 0;
			i = anz - 1;
			wert = s[anz - 1].wert;
		}
		else
		{ /* Suchen der Stuetzstelle. Vorraussetzung: Zeitpunkte aufsteigend
		     geordnet */
			while (wert < s[i].wert)
				i++;
		}
	}

	switch (methode)
	{
		default:
		case 0:
			/* Stueckweise konstant */
			if (fabs(s[i].wert - s[i - 1].wert) > DBL_MIN)
				return (s[i].punkt - s[i - 1].punkt) /
				       (s[i].wert - s[i - 1].wert);
			else
				return Signum(s[i + 1].punkt - s[i].punkt) / DBL_EPSILON;
		case 1:
			/* Gleitend */
			if ((i > 1) && (i < anz - 2))
			{
				s1 = (0.5 * s[i].punkt - 0.5 * s[i - 2].punkt) /
				     (0.5 * s[i].wert - 0.5 * s[i - 2].wert);

				s2 = (0.5 * s[i + 1].punkt - 0.5 * s[i - 1].punkt) /
				     (0.5 * s[i + 1].wert - 0.5 * s[i - 1].wert);

				w = (wert - s[i - 1].wert) / (s[i].wert - s[i - 1].wert);

				return (1. - w) * s1 + w * s2;
			}
			else
			{
				/* Stueckweise konstant */
				if (fabs(s[i].wert - s[i - 1].wert) > DBL_MIN)
					return (s[i].punkt - s[i - 1].punkt) /
					       (s[i].wert - s[i - 1].wert);
				else
					return Signum(s[i + 1].punkt - s[i].punkt) / DBL_EPSILON;
			}
	}
}

