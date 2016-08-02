/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: files0.c
 */
/* Aufgabe:
   Enthaelt die uebergeordneten Datei- Ein- und Ausgaberoutinen, sowie
   das Speichern der Durchbruchskurven.
 */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   02/1999     CT         Bugfix: Anpassung fuer mit 0 beginnende Elemente.
                          Kennung fuer Version im Kopf. Versionsabhaengiges
                          lesen.
   02/1999     CT         Korrekturen fuer herausgeloeste Materialgruppen
   03/1999     AH         Korrekturen nicht noetig, da die Materialzuweisung
                          nicht mehr auf dieser Ebene erfolgt. Die Abfrage ist
                          jetzt auskommentiert. Sie war hier nur um
   Kompatibilitaeten
                          zum Rockflow aufzubewahren.
   03/1999     CT         anz_matxx, start_mat entfernt
   02/2000     CT         Restart wieder hergestellt
   07/2000     AH         Vorbereitungen zum HGM.
   9/2000     CT         Neu: RefreshNodeOutputData, Warnungen beseitigt
   10/2001     AH         Inverse Modellierung
   Trennung und Anpassung (CreateFileData)
   In DestroyFileData Datenfeld auskommentiert.
   Neues Konzept fuer Datenbank-Verwaltung ah inv
   01/2002     MK         DisplayMsgX-Umleitung in *.msg-Datei: OpenMsgFile
   08/2002     MK         GetPathRFDFile
   ConfigFileData aus CreateFileData herausgeloest
   03/2003     RK         Quellcode bereinigt, Globalvariablen entfernt
 */
/**************************************************************************/

#include "FileToolsRF.h"

#include <cmath>
#include <cstdlib>

#include "makros.h"
#include "Configure.h"
#include "MemWatch.h"
#include "display.h"
#include "memory.h"
#include "FileTools.h"


#define KEYWORD '#'
#define SUBKEYWORD '$'


using namespace std;


/**************************************************************************/
/* ROCKFLOW - Funktion: GetLineFromFile1
 */
/* Aufgabe:
   Liest aus dem Eingabefile *ein die n�chste Zeile
   F�ngt die Zeile mit ";" an oder ist sie leer, wird sie ausgelassen
   R�ckgabe ist ist ein string mit dem Zeileninhalt ab dem ersten
   Nicht-Leerzeichen
   bis zum ersten Auftreten des Kommentartzeichens ";"
 */
/* Programmaenderungen:
    05/2004     SB  First Version
 */
/*  09/2005     CC move from fem to geo
 **************************************************************************/
string GetLineFromFile1(ifstream* ein)
{
	//	return readNonBlankLineFromInputStream(*ein);
	string line, zeile = "";
	int fertig = 0, i = 0, j = 0;
	char zeile1[MAX_ZEILE];
	line = "";  // WW
	//----------------------------------------------------------------------
	while (fertig < 1)
	{
		if (ein->getline(zeile1, MAX_ZEILE))  // Zeile lesen
		{
			line = zeile1;  // character in string umwandeln
			i = (int)line.find_first_not_of(
			    " ", 0);  // Anf�ngliche Leerzeichen �berlesen, i=Position des
			              // ersten Nichtleerzeichens im string
			j = (int)line.find(";", i);  // Nach Kommentarzeichen ; suchen. j =
			                             // Position des Kommentarzeichens, j=-1
			                             // wenn es keines gibt.
			if (j != i) fertig = 1;      // Wenn das erste nicht-leerzeichen ein
			// Kommentarzeichen ist, zeile �berlesen. Sonst ist
			// das eine Datenzeile
			if ((i != -1))
				zeile = line.substr(i, j - i);  // Ab erstem nicht-Leerzeichen
			                                    // bis Kommentarzeichen
			                                    // rauskopieren in neuen
			                                    // substring, falls Zeile nicht
			                                    // leer ist
			i = (int)zeile.find_last_not_of(" ");  // Suche nach dem letzten
			                                       // Zeichen, dass kein
			                                       // Leerzeichen ist
			if (i >= 0)
			{
				//		  line.clear(); // = "";
				line = zeile.substr(
				    0, i + 1);  // Leerzeichen am Ende rausschneiden
				//		  zeile.clear(); // = "";
				zeile = line;
			}
		}
		else  // end of file found

			fertig = 1;
	}  // end while(...)
	//----------------------------------------------------------------------
	return zeile;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: FilePrintString
 */
/* Aufgabe:
   Schreibt Zeichenkette ohne Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
   E char *s: Zeichenkette
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int FilePrintString(FILE* f, const char* s)
{
	if ((int)fprintf(f, "%s", s) != (int)strlen(s)) return 0;
	return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: FilePrintInt
 */
/* Aufgabe:
   Schreibt Integer-Wert ohne Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
   E int x: Integer-Wert
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   05/1994     MSR        Erste Version
 */
/**************************************************************************/
int FilePrintInt(FILE* f, int x)
{
	if (fprintf(f, " %i ", x) < 0) return 0;
	return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: FilePrintLong
 */
/* Aufgabe:
   Schreibt Long-Wert ohne Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
   E long x: Long-Wert
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   05/1994     MSR        Erste Version
 */
/**************************************************************************/
int FilePrintLong(FILE* f, long x)
{
	if (fprintf(f, " %ld ", x) < 0) return 0;
	return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: FilePrintDouble
 */
/* Aufgabe:
   Schreibt Double-Wert ohne Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
   E double x: Double-Wert
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   05/1994     MSR        Erste Version
   12/1995     cb         E-Format
 */
/**************************************************************************/
int FilePrintDouble(FILE* f, double x)
{
#ifdef FORMAT_DOUBLE
	if (fprintf(f, " % #*.*g ", FPD_GESAMT, FPD_NACHKOMMA, x) < 0) return 0;
#else
	if (fprintf(f, " % #g ", x) < 0) return 0;
#endif
	return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrReadDouble
 */
/* Aufgabe:
   Liest Double-Wert aus String und schreibt Protokoll in Datei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R double *x: gelesener Double-Wert
   E char *s: Zeichenkette, aus der gelesen werden soll
   E FILE *f: Dateizeiger fuer Protokolldatei
   E FctTestDouble func: Funktionszeiger auf die Funktion, die den
                         eingelesenen Wert auf Gueltigkeit testet
   R int *pos: Anzahl der bisher gelesenen Zeichen
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int StrReadDouble(double* x, char* s, FILE* f, int* pos)
{
	*x = 0.0;
	if (sscanf(s, " %lf%n", x, pos) <= 0)
	{
		*pos = 0; /* nichts sinnvolles gelesen */
		fprintf(f,
		        "\n %f      *** Fehler: Kein Wert eingelesen (double) !!!\n",
		        *x);
		return 0;
	}
	else
	{
		/* CT: Protokolformat geaendert */
		if ((fabs(*x) < 100000.) && (fabs(*x) >= 0.1))
			fprintf(f, " %f ", *x);
		else
			fprintf(f, " %e ", *x);

		return 1;
	}
}



/**************************************************************************
   STRLib-Method: SubKeyword
   Task:
   Programing:
   09/2004 OK Implementation
   last modification:
**************************************************************************/
bool SubKeyword(const std::string& line)
{
	if (line.find(SUBKEYWORD) != std::string::npos)
		return true;
	else
		return false;
}

/**************************************************************************
   STRLib-Method: SubKeyword
   Task:
   Programing:
   09/2004 OK Implementation
   last modification:
**************************************************************************/
bool Keyword(const std::string& line)
{
	if (line.find(KEYWORD) != std::string::npos)
		return true;
	else
		return false;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrUp
 */
/* Aufgabe:
   wandelt Zeichenkette in Grossbuchstaben um
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   X char *s: umzuwandelnde Zeichenkette
 */
/* Ergebnis:
   umgewandelte Zeichenkette
 */
/* Programmaenderungen:
   03/1994   MSR   Erste Version
 */
/**************************************************************************/
char* StrUp(const char* s)
{
	int i;
	int l = (int)strlen(s);
	char* tmp = new char[l];
	strcpy(tmp, s);
	for (i = 0; i < l; i++)
		if (islower((int)s[i])) tmp[i] = (char)toupper((int)s[i]);
	return tmp;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StringReadStr
 */
/* Aufgabe:
   Liest Zeichenkette aus String
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R char **x: Zeiger auf Adresse der gelesenen Zeichenkette; eine vorher
               vorhandene wird geloescht, es sollten nur mit malloc erzeugte
               Zeichenketten verwendet werden.
   E char *s: Zeichenkette, aus der gelesen werden soll
   R int *pos: Anzahl der bisher gelesenen Zeichen
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   06/1999   OK   aus StrReadString
 */
/**************************************************************************/
int StringReadStr(char** x, char* s, int* pos)
{
	*x = NULL;
	//  *x = (char *) Malloc(256);
	*x = new char[256];  // CC
	*x[0] = '\0';
	if (sscanf(s, " %s%n", *x, pos) <= 0)
	{
		int a = (int)strlen(*x);  // CC
		// delete[] *x;//CC
		*x = new char[a + 1];  // CC
		//*x = (char *) Realloc(*x,((int)strlen(*x)+1));
		*pos = 0; /* nichts sinnvolles gelesen */
		return 0;
	}
	else
		return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: LineFeed
 */
/* Aufgabe:
   Schreibt Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int LineFeed(FILE* f)
{
	if (fprintf(f, "\n") < 0) return 0;
	return 1;
}

// int TFDouble ( double *x, FILE *f )
//{
//   return 1;
//}

// int TFString ( char* x, FILE* f )
//{
//	return 1;
//}

/**************************************************************************
   STRLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   last modification:
**************************************************************************/
void remove_white_space(std::string* buffer)
{
	int pos = 0;
	while (pos >= 0)
	{
		pos = (int)buffer->find_first_of(" ");
		if (pos < 0) break;
		buffer->erase(pos, 1);
	}
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrReadStr
 */
/* Aufgabe:
   Liest Zeichenkette aus String und schreibt Protokoll in Datei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R char **x: Adresse der gelesenen Zeichenkett, es muss vorher Speicher
               allokiert werden
   E char *s: Zeichenkette, aus der gelesen werden soll
   E FILE *f: Dateizeiger fuer Protokolldatei
   E FctTestString func: Funktionszeiger auf die Funktion, die den
                         eingelesenen Wert auf Gueltigkeit testet
   R int *pos: Anzahl der bisher gelesenen Zeichen
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   08/2000     CT        Erste Version
 */
/**************************************************************************/
int StrReadStr(char* x, char* s, FILE* f, /*FctTestString func,*/ int* pos)
{
	//   int test;
	x[0] = '\0';
	if (sscanf(s, " %s%n", x, pos) <= 0)
	{
		*pos = 0; /* nichts sinnvolles gelesen */
		fprintf(
		    f, "\n %s      *** Fehler: Kein Wert eingelesen (string) !!!\n", x);
		return 0;
	}
	else
	{
		//      test = func(x,f);
		//      fprintf(f,"%s ",x);
		//      return test;
		fprintf(f, "%s ", x);
		return 1;
	}
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrTestDouble
 */
/* Aufgabe:
   Testet, ob in s noch ein Double kommt;
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
 */
/* Ergebnis:
   0: nein; 1: ja
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int StrTestDouble(char* s)
{
	double i;
	if (sscanf(s, " %lf", &i) <= 0)
		return 0;
	else
		return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrTestHash
 */
/* Aufgabe:
   Testet, ob in s ein # folgt
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
   R int *pos: gelesene Zeichen bis nach dem # (wenn gefunden)
 */
/* Ergebnis:
   0: nein; 1: ja
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int StrTestHash(char* s, int* pos)
{
	int p;
	char h[256];
	if (sscanf(s, " %s%n", h, &p) <= 0)
		return 0;
	else
	{
		if (strcmp(h, "#") == 0)
		{
			*pos = p;
			return 1;
		}
		else
			return 0;
	}
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrOnlyReadStr
 */
/* Aufgabe:
   Liest Zeichenkette aus String aber schreibt Protokoll in Datei nicht
   (for Phreeqc read function)
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R char **x: Adresse der gelesenen Zeichenkett, es muss vorher Speicher
               allokiert werden
   E char *s: Zeichenkette, aus der gelesen werden soll
   E FILE *f: Dateizeiger fuer Protokolldatei
   E FctTestString func: Funktionszeiger auf die Funktion, die den
                         eingelesenen Wert auf Gueltigkeit testet
   R int *pos: Anzahl der bisher gelesenen Zeichen
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   06/2003     MX        Erste Version
 */
/**************************************************************************/
/*MX*/
int StrOnlyReadStr(char* x,
                   char* s,
                   FILE* /*f*/,
                   /*FctTestString func,*/ int* pos)
{
	//   int test;

	x[0] = '\0';
	if (sscanf(s, " %s%n", x, pos) <= 0)
	{
		*pos = 0; /* nichts sinnvolles gelesen */
		return 0;
	}
	else
		//      test = func(x,f);
		//      return test;
		return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrReadSubKeyword
 */
/* Aufgabe:
   Liest ein Sub-Keyword (eingeleitet mit "$") aus Keyword-String. Nur bis zum
   naechsten Hash (naechstes Keyword) oder Stringende
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R char *sub :  Zeichenkette, in die das Sub-Keywords geschrieben wird
   E char *s   : Zeichenkette, aus der gelesen werden soll
   E int begin : Erstes Zeichen, ab dem gesucht werden soll
   R int *found: Erstes Zeichen des Sub-Keywords im String s
   R int *ende : Letztes Zeichen des Sub-Keywords im String s, oder
                 Beginn des naechsten Keywords, oder Ende des Strings

 */
/* Ergebnis:
   1 bei gefundenem Sub-Keyword, sonst 0
 */
/* Programmaenderungen:
   08/2000 C.Thorenz  Erste Version
 */
/**************************************************************************/
int StrReadSubKeyword(char* sub, char* s, int beginn, int* found, int* ende)
{
	int i, xi = 0;

	*found = -1;
	*ende = (int)strlen(s);

	for (i = beginn; i < (int)strlen(s); i++)
	{
		if (s[i] == '$')
		{
			if (*found < 1) /* Anfang des Sub-Keywords merken */
				*found = i;
			else
			{
				/* Ende des Sub-Keywords merken (neues Sub-Keyword folgt) */
				*ende = i;
				break;
			}
		}

		if (s[i] == '#')
		{
			/* Ende des Sub-Keywords merken (neues Keyword folgt) */
			*ende = i;
			break;
		}

		if (*found >= 0)
		{
			sub[xi] = s[i];
			xi++;
		}
	}

	if (*found >= 0) sub[xi] = '\0';

	return *found >= 0;
}

/**************************************************************************
   STRLib-Method: get_sub_string
   Task: sub_string between pos1 and delimiter
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
std::string get_sub_string(const std::string& buffer,
                           const std::string& delimiter,
                           int pos1,
                           int* pos2)
{
	int pos = 0;
	std::string empty_string("");
	// string sub_string_this;
	*pos2 = (int)buffer.find(delimiter, pos1);
	if (*pos2 < 0) return empty_string;
	while (*pos2 <= pos1)
	{
		pos1++;
		*pos2 = (int)buffer.find(delimiter, pos1);
		if (*pos2 < 0)
		{
			*pos2 = (int)buffer.size();
			break;
		}
		if (pos1 >= (int)buffer.size()) break;
	}
	string sub_string_this = buffer.substr(pos1, *pos2);
	while (pos >= 0)
	{
		pos = (int)sub_string_this.find_first_of(" ");
		if (pos < 0) break;
		sub_string_this.erase(pos, 1);
	}
	return sub_string_this;
}

/**************************************************************************
   STRLib-Method: get_sub_string2
   Task:
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
std::string get_sub_string2(const std::string& buffer,
                            const std::string& delimiter,
                            std::string* tmp)
{
	int pos2 = (int)buffer.find_first_of(delimiter);
	std::string sub_string = buffer.substr(0, pos2);
	*tmp = buffer.substr(pos2 + delimiter.size());
	return sub_string;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetUncommentedLine
 */
/* Aufgabe:
   R�ckgabe ist ist ein string mit dem Zeileninhalt ab dem ersten
   Nicht-Leerzeichen
   bis zum ersten Auftreten des Kommentartzeichens ";"
   Abgeleitet aus GetLineFromFile1() */
/* Programmaenderungen:
    06/2009     SB  First Version
 **************************************************************************/
std::string GetUncommentedLine(std::string line)
{
	std::string zeile = "";
	int i = 0, j = 0;
	//----------------------------------------------------------------------
	i = (int)line.find_first_not_of(" ", 0);  // Anf�ngliche Leerzeichen
	                                          // �berlesen, i=Position des
	                                          // ersten Nichtleerzeichens im
	                                          // string
	j = (int)line.find(";", i);  // Nach Kommentarzeichen ; suchen. j = Position
	                             // des Kommentarzeichens, j=-1 wenn es keines
	                             // gibt.
	if ((i != -1))
		zeile = line.substr(i, j - i);  // Ab erstem nicht-Leerzeichen bis
	                                    // Kommentarzeichen rauskopieren in
	                                    // neuen substring, falls Zeile nicht
	                                    // leer ist
	i = (int)zeile.find_last_not_of(
	    " ");  // Suche nach dem letzten Zeichen, dass kein Leerzeichen ist
	if (i >= 0)
	{
		line = zeile.substr(0, i + 1);  // Leerzeichen am Ende rausschneiden
		zeile = line;
	}

	return zeile;
}
