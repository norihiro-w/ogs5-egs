/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   rf_kinreact.cpp

   KINETIC REACTIONS

   FEMLib-Object KinReact

   Programming:
   01/2004    Dirk Schaefer      Original Implementation
   02/2006    Sebastian Bauer    Adaption to C++ Class structure, new FEM
 concept
   05/2007    Dirk Schaefer      Addition of NAPL-dissolution
 ***************************************************************************/

#include "rf_kinreact.h"

#include <cfloat>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

#include "display.h"
#include "FileToolsRF.h"
#include "makros.h"
#include "StringTools.h"

#include "msh_lib.h"

#include "rfmat_cp.h"
#include "rf_mfp_new.h"
#include "rf_mmp_new.h"
#include "rf_msp_new.h"
#include "rf_tim_new.h"
#include "tools.h"

using namespace std;
using SolidProp::CSolidProperties;
using namespace Math_Group;

vector<CKinReact*> KinReact_vector;  // declare instance CKinReact_vector
vector<CKinReactData*>
    KinReactData_vector;           // declare instance CKinReact_vector
vector<CKinBlob*> KinBlob_vector;  // declare extern instance of class Blob

static double dmaxarg1, dmaxarg2;
#define DMAX(a, b)                   \
	(dmaxarg1 = (a), dmaxarg2 = (b), \
	 (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
#define DMIN(a, b)                   \
	(dmaxarg1 = (a), dmaxarg2 = (b), \
	 (dmaxarg1) < (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

/* Constructor for MonodSubstruct */
MonodSubstruct::MonodSubstruct(void)
{
	species = "blubb";
	speciesnumber = -1;
	concentration = -1.0e9;
	order = -99.0;
	isotopecouplenumber = -1;  // CB isotope fractionation
	threshhold = false;
	threshConc = -1.0e9;
	threshOrder = -99;
}

/* Destructor for MonodSubstruct */
MonodSubstruct::~MonodSubstruct(void)
{
}

/***************************************************************************
   FEMLib-Method:
   Task: CKinReact constructor
   Programing:
   02/2006 SB Implementation
   06/2007 DS NAPL-dissolution added
***************************************************************************/
CKinReact::CKinReact(void)
{
	name = "NULL";
	type = "NULL";
	number_reactionpartner = 0;
	reactionpartner.clear();
	stochmet.clear();
	rateconstant = 0.0;
	rateorder = 0.0;
	number_monod = 0;
	number_inhibit = 0;
	number_production = 0;
	number_isotope_couples = 0;  // CB isotope fractionation
	monod.clear();
	inhibit.clear();
	production.clear();
	bacteria_name = "NULL";
	bacteria_number = -1;
	ProductionStoch.clear();
	grow = -1;
	isoenfac = 0;      // CB Isotope fractionation
	degType = "NULL";  // CB Isotope fractionation

	// CB Not this particular reaction on specified GEO-Objects
	switched_off_node.clear();

	ProdStochhelp.clear();
	ex_param.clear();
	ex_species.clear();
	ex_species_names.clear();
	exSurfaceID = -1;
	exType = "NULL";
	// NAPL-dissolution
	blob_name = "NULL";
	blob_ID = -1;
	Csat_pure = 0.;
	current_Csat = 0.;
	Density_NAPL = 0.;
	//
	typeflag_monod = 0;
	typeflag_exchange = 0;
	typeflag_exchange_linear = 0;
	typeflag_exchange_langmuir = 0;
	typeflag_exchange_freundlich = 0;
	typeflag_napldissolution = 0;
	typeflag_iso_fract = 0;  // CB isotope fractionation
}

/***************************************************************************
   FEMLib-Method:
   Task: CKinReact destructor
   Programing:
   02/2006 SB Implementation
***************************************************************************/
CKinReact::~CKinReact(void)
{
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   02/2004 SB Implementation
**************************************************************************/
bool KRRead(const std::string& file_base_name,
            const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
	// File handling
	std::string krc_file_name(file_base_name + KRC_FILE_EXTENSION);
	std::ifstream krc_file(krc_file_name.data(), ios::in);
	if (!krc_file.good()) return false;

	CKinReact* m_kr = NULL, * m_kr1 = NULL;
	CKinReactData* m_krd = NULL;
	CKinBlob* m_bp = NULL;
	char line[MAX_ZEILE];
	string sub_line;
	string line_string;
	string database_file_name;
	int found;
	size_t length;
	ios::pos_type position;
	string m_name, sp_name;

	KRCDelete();
	//========================================================================
	// Keyword loop
	cout << "KinReact Read" << endl;
	while (!krc_file.eof())
	{
		krc_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos) break;
		//----------------------------------------------------------------------
		// keyword found  // Read kinetic Reactions
		if (line_string.find("#REACTION") != string::npos)
		{
			m_kr = new CKinReact();
			m_kr->Read(&krc_file, geo_obj, unique_name);
			m_kr->number = (int)KinReact_vector.size();
			KinReact_vector.push_back(m_kr);
		}  // keyword found
		   // keyword found  // Read KinReactData
		if (line_string.find("#KINREACTIONDATA") != string::npos)
		{
			m_krd = new CKinReactData();
			m_krd->Read(&krc_file, geo_obj, unique_name);
			KinReactData_vector.push_back(m_krd);
		}  // keyword found
		   // keyword found  // Read BlobProperties
		if (line_string.find("#BLOB_PROPERTIES") != string::npos)
		{
			m_bp = new CKinBlob();
			m_bp->Read(&krc_file, geo_obj, unique_name);
			KinBlob_vector.push_back(m_bp);
		}  // keyword found
	}      // eof
	// Close input file
	krc_file.close();

	if (!m_krd)
	{
		cout << " No keyword #KINREACTIONDATA specified - setting defaults"
		     << endl;
		m_krd = new CKinReactData();
		KinReactData_vector.push_back(m_krd);
	}

	//========================================================================

	/* Check, if data has to be completed from database file */
	//  cout << "checking database" << endl;
	if (!KinReact_vector.empty())
	{
		// File handling, open reaction database file
		database_file_name = "reactions.dbf";
		ifstream dbf_file(database_file_name.data(), ios::in);
		//	if (!dbf_file.good()) cout << " No database file found " << endl;

		// go through all reactions in reaction vector
		length = KinReact_vector.size();
		for (size_t i = 0; i < length; i++)
		{
			m_kr = KinReact_vector[i];
			// if no type is given in input file, only name, then read from
			// database file
			if (m_kr->type == "NULL")
			{
				found = 0;
				dbf_file.seekg(0L, ios::beg);  // rewind database file
				while (!dbf_file.eof() && found == 0)
				{
					if (!GetLineFromFile(line, &dbf_file)) break;
					line_string = line;
					// keyword found
					if (line_string.find("#REACTION") != string::npos)
					{
						m_kr1 = new CKinReact();
						position = m_kr1->Read(&dbf_file, geo_obj, unique_name);
						// check if this is the right one
						if (m_kr->name == m_kr1->name)
						{
							// Insert in Reaction vector and remove old reaction
							// (replacement)
							KinReact_vector.insert(
							    KinReact_vector.begin() + i + 1, m_kr1);
							KinReact_vector.erase(KinReact_vector.begin() + i);
							found = 1;
						}
						else
							// not the right equation:
							delete m_kr1;
					}  // end if(line_str..)  keyword found
				}      // eof dbf-file

				if (found == 0)
				{
					// reaction not complete in input file and also not
					// specified in dbf_file
					cout << " ERROR! Reaction " << m_kr->name
					     << "  not found in database file" << endl;
					cout << " Reaction removed from set of reaction " << endl;
					// remove reaction from reaction vector
					KinReact_vector.erase(KinReact_vector.begin() + i);
				}
			}  // end of if(m_kr=NULL
		}      // end for
		// close database file
		dbf_file.close();
	}  // end if kinreaction_vec.size() >
	/* end data base file read for reactions  */

	/* check reaction data consistency */
	std::cout << " Checking reaction data consistency for "
	          << KinReact_vector.size() << " specified reactions "
	          << "\n";
	length = KinReact_vector.size();
	for (size_t i = 0; i < length; i++)
	{
		m_kr = KinReact_vector[i];
		if (!m_kr->CheckReactionDataConsistency())
		{
			cout << " ERROR! Reaction " << m_kr->name
			     << "  removed from set of reactions due to data inconsistency"
			     << endl;
			KinReact_vector.erase(KinReact_vector.begin() + i);
		}
	}  // end for(i=0;.. consistency check
	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ configure function
   Programing:
   02/2004 SB Implementation
   05/2007 DS NAPL dissolution added
   07/2009 CB Isotope fractionation
   08/2009 CB Reaction deactivation
**************************************************************************/
void KRConfig(const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
	CKinReact* m_kr = NULL;  //, *m_kr1=NULL;
	CKinReactData* m_krd = NULL;
	CompProperties* m_cp = NULL;
	CRFProcess* m_pcs = NULL;
	vector<long> nodes_vector;
	CMediumProperties* m_mat_mp = NULL;
	bool ok = true;

	// CB reaction deactivation
	MeshLib::CElem* m_dnele = NULL;
	MeshLib::CElem* m_snele = NULL;
	MeshLib::CNode* m_dnnod = NULL;
	vector<int> ReactNeighborNodes;
	vec<long> secnnodesindices(8);
	vec<long> primnnodesindices(8);

	if (!KinReactData_vector.empty())
	{
		m_krd = KinReactData_vector[0];
		if (m_krd == NULL) cout << " Error - no m_krd data " << endl;
		// Set up additional data structures for calculations

		// Set number of reactions
		m_krd->NumberReactions = (int)KinReact_vector.size();
		size_t length(cp_vec.size());

		// Check if all reaction partners are specified as processes
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			// Check presence of bacteria species
			if (m_kr->bacteria_name.compare("NULL") != 0)
			{
				m_pcs = PCSGet("MASS_TRANSPORT", m_kr->bacteria_name);
				if (m_pcs == NULL)
				{
					cout << " Warning: Component " << m_kr->bacteria_name
					     << " specified in KinReact as biomass but not given "
					        "as transport process " << endl;
					ok = false;
				}
			}
			// Check if monod species are specified as components
			for (int i = 0; i < m_kr->number_monod; i++)
			{
				m_pcs = PCSGet("MASS_TRANSPORT", m_kr->monod[i]->species);
				if (m_pcs == NULL)
				{
					cout << " Warning: Component " << m_kr->monod[i]->species
					     << " specified in KinReact as monod species but not "
					        "given as transport process " << endl;
					ok = false;
				}
			}
			// Check if inhibition species are specified as components
			for (int i = 0; i < m_kr->number_inhibit; i++)
			{
				m_pcs = PCSGet("MASS_TRANSPORT", m_kr->inhibit[i]->species);
				if (m_pcs == NULL)
				{
					cout << " Warning: Component " << m_kr->inhibit[i]->species
					     << " specified in KinReact as inhibition species but "
					        "not given as transport process " << endl;
					ok = false;
				}
			}
			// Check productionstoch
			for (size_t i = 0; i < m_kr->ProdStochhelp.size(); i++)
			{
				m_pcs =
				    PCSGet("MASS_TRANSPORT", m_kr->ProdStochhelp[i]->species);
				if (m_pcs == NULL)
				{
					cout << " Warning: Component "
					     << m_kr->ProdStochhelp[i]->species
					     << " specified in KinReact as produced species but "
					        "not given as transport process " << endl;
					ok = false;
				}
			}
			if (m_kr->type.compare("exchange") == 0)
				for (int i = 0; i < m_kr->number_reactionpartner; i++)
				{
					m_pcs = PCSGet("MASS_TRANSPORT", m_kr->reactionpartner[i]);
					if (m_pcs == NULL)
					{
						cout << " Warning: Component "
						     << m_kr->reactionpartner[i]
						     << " specified in KinReact as reaction partner "
						        "but not given as transport process " << endl;
						ok = false;
					}
				}
			// check isotope couples
			// todo CB isotope fract
			for (int i = 0; i < m_kr->number_isotope_couples; i++)
			{
				m_pcs = PCSGet("MASS_TRANSPORT", m_kr->Isotope_light);
				if (m_pcs == NULL)
				{
					cout << " Warning: Component " << m_kr->Isotope_light
					     << " specified in KinReact in isotope couple but not "
					        "given as transport process " << endl;
					ok = false;
				}
			}
		}  // loop over m_krd->NumberReactions

		//
		if (ok == false)
		{
			cout << " Components missing, Stopping" << endl;
			cout.flush();
			exit(1);
		}

		// Set vector is_a_bacterium
		//   cout << " Length of cp_vec: " << length << endl;
		for (size_t j = 0; j < length; j++)
			m_krd->is_a_bacterium.push_back(0);  // initialize
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->type.compare("monod") == 0)
				for (size_t i = 0; i < length; i++)
				{
					m_cp = cp_vec[i];
					if (m_kr->bacteria_name.compare(m_cp->compname) == 0)
						m_krd->is_a_bacterium[i] = 1;
				}
		}
		// Set Vector ProductionStoch for each reaction
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->type.compare("monod") == 0)
			{
				for (size_t i = 0; i < length; i++)
					// initialize
					m_kr->ProductionStoch.push_back(0.0);

				// Get Stochiometry
				/*
				   for(k=0;k<(int)m_kr->reactionpartner.size();k++){
				   m_name = m_kr->reactionpartner[k];
				   for(i=0;i<length;i++)
				   if(m_name.compare(cp_vec[i]->compname) == 0)
				   m_kr->ProductionStoch[i] = m_kr->stochmet[k];
				   }
				 */
				for (size_t k = 0; k < m_kr->ProdStochhelp.size(); k++)
				{
					std::string const& m_name(m_kr->ProdStochhelp[k]->species);
					for (size_t i = 0; i < length; i++)
						if (m_name.compare(cp_vec[i]->compname) == 0)
							m_kr->ProductionStoch[i] =
							    m_kr->ProdStochhelp[k]->concentration;
				}
			}  // if type == monod
		}      // vector ProductionStoch

		// Set numbers for monod species for each reaction
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->type.compare("monod") == 0)
				for (size_t i = 0; i < m_kr->monod.size(); i++)
				{
					std::string const& m_name(m_kr->monod[i]->species);
					for (size_t k = 0; k < length; k++)
						if (m_name.compare(cp_vec[k]->compname) == 0)
							m_kr->monod[i]->speciesnumber = k;
				}
		}  // monod substructure numbers

		// CB Isotope fractionation
		// Set number for isotope couple in monod reaction
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->degType.compare("isotope_fractionation") == 0)
				for (size_t i = 0; i < m_kr->monod.size(); i++)
				{
					std::string m_name(m_kr->monod[i]->species);
					if (m_name.compare(m_kr->Isotope_light) == 0)
					{
						m_name = m_kr->Isotope_heavy;
						for (size_t k = 0; k < length; k++)
							if (m_name.compare(cp_vec[k]->compname) == 0)
							{
								m_kr->monod[i]->isotopecouplenumber = k;
								break;
							}
					}
					else if (m_name.compare(m_kr->Isotope_heavy) == 0)
					{
						m_name = m_kr->Isotope_light;
						for (size_t k = 0; k < length; k++)
							if (m_name.compare(cp_vec[k]->compname) == 0)
							{
								m_kr->monod[i]->isotopecouplenumber = k;
								break;
							}
					}
				}
		}  // monod isotope substructure numbers

		// Set numbers for inhibition species for each reaction
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->type.compare("monod") == 0)
				for (size_t i = 0; i < m_kr->inhibit.size(); i++)
				{
					std::string const& m_name(m_kr->inhibit[i]->species);
					for (size_t k = 0; k < length; k++)
						if (m_name.compare(cp_vec[k]->compname) == 0)
							m_kr->inhibit[i]->speciesnumber = k;
				}
		}  // inhibition substructure numbers

		// Set bacteria number for each reaction
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->type.compare("monod") == 0)
			{
				std::string const& m_name(m_kr->bacteria_name);
				for (size_t k = 0; k < length; k++)
				{
					m_cp = cp_vec[k];
					if (m_name.compare(m_cp->compname) == 0)
						m_kr->bacteria_number = k;
				}
			}
		}  // bacteria numbers

		// Set flags type_monod and type_exchange
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			//			dummy = m_kr->type;
			if (m_kr->type.compare("monod") == 0)  // CB Isotope fractionation
			{
				m_kr->typeflag_monod = 1;
				if (m_kr->degType.compare("isotope_fractionation") == 0)
					m_kr->typeflag_iso_fract = 1;
			}
			if (m_kr->type.compare("exchange") == 0)
			{
				m_kr->typeflag_exchange = 1;
				if (m_kr->exType.compare("linear") == 0)
					m_kr->typeflag_exchange_linear = 1;
				if (m_kr->exType.compare("langmuir") == 0)
					m_kr->typeflag_exchange_langmuir = 1;
				if (m_kr->exType.compare("freundlich") == 0)
					m_kr->typeflag_exchange_freundlich = 1;
			}
			if (m_kr->type.compare("NAPLdissolution") == 0)
				m_kr->typeflag_napldissolution = 1;
		}  // typeflags

		// exchange reactions
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->type.compare("exchange") == 0)
				// move names from equation to vector with names
				for (size_t i = 0; i < m_kr->reactionpartner.size(); i++)
				{
					std::string const& m_name(m_kr->reactionpartner[i]);
					m_kr->ex_species_names.push_back(m_name);
					//	find species numbers for soecies names
					for (size_t k = 0; k < length; k++)
						if (m_name.compare(cp_vec[k]->compname) == 0)
							m_kr->ex_species.push_back(k);
				}
		}  // exchange species numbers

		//#ds NAPLdissolution reaction
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->type.compare("NAPLdissolution") == 0)
				// move names from equation to vector with names
				for (size_t i = 0; i < m_kr->reactionpartner.size(); i++)
				{
					std::string const& m_name(m_kr->reactionpartner[i]);
					m_kr->ex_species_names.push_back(m_name);
					//	find species numbers for species names
					for (size_t k = 0; k < length; k++)
						if (m_name.compare(cp_vec[k]->compname) == 0)
							m_kr->ex_species.push_back(k);
				}

			// set index numbers for blob_id
			CKinBlob* m_kb = NULL;
			for (size_t i = 0; i < KinBlob_vector.size(); i++)
			{
				m_kb = KinBlob_vector[i];
				if (m_kr->blob_name.compare(m_kb->name) == 0) m_kr->blob_ID = i;
			}
		}  // NAPLdissolution species numbers

		// exchange reactions numbers for m_krd
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if (m_kr->exType.compare("linear") == 0) m_krd->NumberLinear++;
			if (m_kr->exType.compare("langmuir") == 0) m_krd->NumberLangmuir++;
			if (m_kr->exType.compare("freundlich") == 0)
				m_krd->NumberFreundlich++;
			if (m_kr->type.compare("monod") == 0) m_krd->NumberMonod++;
			if (m_kr->type.compare("NAPLdissolution") == 0)
				m_krd->NumberNAPLdissolution++;
		}
		/*
		   // set up vector sp_index
		   length = (int)cp_vec.size();
		   for(j=0; j < length; j++)
		   m_krd->sp_index.push_back(-10); //initialize
		   // Go through all reactions and set sp_index[x] active, if species x
		   is found
		   for(j=0;j<m_krd->NumberReactions;j++){
		   m_kr = KinReact_vector[j];
		   if(m_kr->type.compare("monod")==0){
		   if(m_kr->bacteria_number > 0) m_krd->sp_index[m_kr->bacteria_number]
		   = 1;
		   for(k=0;k<m_kr->number_monod;k++)
		   if(m_kr->monod[k]->speciesnumber > 0)
		   m_krd->sp_index[m_kr->monod[k]->speciesnumber] = 1;
		   for(k=0;k<m_kr->number_inhibit;k++)
		   if(m_kr->inhibit[k]->speciesnumber > 0)
		   m_krd->sp_index[m_kr->inhibit[k]->speciesnumber] = 1;
		   for(k=0;k<m_kr->number_production;k++)
		   if(m_kr->production[k]->speciesnumber > 0)
		   m_krd->sp_index[m_kr->production[k]->speciesnumber] = 1;
		   for(k=0;k<(int)m_kr->ProdStochhelp.size();k++)
		   if(m_kr->ProdStochhelp[k]->speciesnumber > 0)
		   m_krd->sp_index[m_kr->ProdStochhelp[k]->speciesnumber] = 1;
		   for(k=0;k<length;k++)
		   if(m_kr->ProductionStoch[k] !=0) m_krd->sp_index[k] = 1;
		   }
		   if(m_kr->type.compare("exchange")==0){
		   for(k=0;k<(int)m_kr->ex_species.size();k++)
		   if(m_kr->ex_species[k] > 0) m_krd->sp_index[m_kr->ex_species[k]] = 1;

		   }
		   }

		   // renumber in sp_index
		   count =0;
		   for(j=0;j<length;j++)
		   if(m_krd->sp_index[j] > 0){
		   m_krd->sp_index[j] = count;
		   count++;
		   }
		   // output for checking
		   m_krd->kr_active_species = count;
		   cout << "m_krd->sp_index : " << endl;
		   for(j=0;j<length;j++) cout << j << ", " << m_krd->sp_index[j] <<
		   endl;
		   cout << " total number of active species: " <<
		   m_krd->kr_active_species << endl;
		 */
		// set up vectors sp_pcs and sp_varind
		for (size_t j = 0; j < length; j++)
		{
			std::string const& sp_name(cp_vec[j]->compname);
			// HS, PCSGet not needed any more.
			// m_pcs = PCSGet("MASS_TRANSPORT", sp_name);
			// new style:
			m_pcs = cp_vec[cp_name_2_idx[sp_name]]->getProcess();
			if (m_pcs != NULL)
				// HS, not needed any more
				// idummy = PCSGetPCSIndex("MASS_TRANSPORT", sp_name);
				// m_krd->sp_pcsind.push_back(idummy);
				// new timelevel
				m_krd->sp_varind.push_back(m_pcs->GetNodeValueIndex(sp_name) +
				                           1);
			//	cout << " PCS: " << j << ", " << sp_name << ", " << idummy << ",
			//" << m_pcs->GetNodeValueIndex(sp_name) + 1 << endl;
		}
		/*
		   // set up short versions of vectors is_a_bacterium2 and
		   ProductionStoch for each reaction
		   for(j=0;j<m_krd->kr_active_species;j++)
		   m_krd->is_a_bacterium2.push_back(0);
		   for(j=0;j<length;j++)
		   if(m_krd->sp_index[j] > -1)
		   m_krd->is_a_bacterium2[m_krd->sp_index[j]] =
		   m_krd->is_a_bacterium[j];
		   cout << "is_a_bacterium2:" << endl;
		   for(j=0;j<m_krd->kr_active_species;j++)
		   cout << m_krd->is_a_bacterium2[j] << endl;
		   // Set Vector ProductionStoch2 for each reaction
		   for(j=0;j<m_krd->NumberReactions;j++){
		   m_kr = KinReact_vector[j];
		   if(m_kr->type.compare("monod")==0){
		   // initialize
		   for(i=0;i<m_krd->kr_active_species;i++)
		   m_kr->ProductionStoch2.push_back(0.0); //initialize
		   // Get Stochiometry from ProductionStoch
		   for(i=0;i<length;i++)
		   if(m_krd->sp_index[i] > -1)
		   m_kr->ProductionStoch2[m_krd->sp_index[i]] =
		   m_kr->ProductionStoch[i];
		   cout << " ProductionStoch2 for reaction " << m_kr->name << endl;
		   for(i=0;i<m_krd->kr_active_species;i++)
		   cout << m_kr->ProductionStoch2[i] << endl;
		   cout << endl << " ProductionStoch for reaction " << m_kr->name <<
		   endl;
		   for(i=0;i<length;i++)
		   cout << m_kr->ProductionStoch[i] << endl;
		   } //if type == monod
		   } // vector ProductionStoch2
		 */

		// CB isotope fractionation
		// modify rateconstant for Monod-type reactions and isotope
		// fractionation
		for (int j = 0; j < m_krd->NumberReactions; j++)
		{
			m_kr = KinReact_vector[j];
			if ((m_kr->typeflag_monod == 1) && (m_kr->typeflag_iso_fract == 1))
				m_kr->rateconstant *= (1 + m_kr->isoenfac / 1000);
		}
		/********************************************************/
		// check global requirements for kinetic reactions:
		// check if porosities for phases are set
		if (m_krd->NumberMonod > 0)
			for (size_t j = 0; j < mmp_vector.size(); j++)
				if (mmp_vector[j]->vol_bio < MKleinsteZahl)
					std::cout << "Warning: Porosity of bio-phase is 0.0 ! "
					             "Change Settings in *.mmp file "
					          << "\n";
		//#ds  k= m_krd->NumberLinear + m_krd->NumberFreundlich +
		// m_krd->NumberLangmuir;
		int k = m_krd->NumberLinear + m_krd->NumberFreundlich +
		        m_krd->NumberLangmuir + m_krd->NumberNAPLdissolution;
		if (k > 0)
			for (size_t j = 0; j < mmp_vector.size(); j++)
				if (mmp_vector[j]->vol_mat < MKleinsteZahl)
					std::cout << "Warning: Porosity of solid phase is 0.0 ! "
					             "Change Settings in *.mmp file "
					          << "\n";

		/********************************************************/
		// Set up vector is_a_CCBC
		CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
		if (m_msh == NULL)
		{
			cout << "No mesh in KRConfig" << endl;
			exit(1);
		}

		const size_t mesh_node_vector_size(m_msh->nod_vector.size());
		// Initialize vector is_a_CCBC
		for (size_t l = 0; l < mesh_node_vector_size; l++)
			m_krd->is_a_CCBC.push_back(false);
		// Go through specified geometry elements
		std::string s_geo_name, s_geo_type;
		size_t s_geo_id;
		for (size_t j = 0; j < m_krd->NoReactGeoName.size(); j++)
		{
			s_geo_name = m_krd->NoReactGeoName[j];
			s_geo_type = m_krd->NoReactGeoType[j];
			s_geo_id = m_krd->NoReactGeoID[j];
			//------------------------------------------------------------------
			if (s_geo_type.compare("POINT") == 0)
			{
				const std::vector<GEOLIB::Point*>* pnt_vec(
				    geo_obj.getPointVec(unique_name));
				if (pnt_vec)
					m_krd->is_a_CCBC[m_msh->GetNODOnPNT(
					    (*pnt_vec)[m_krd->NoReactGeoID[j]])] = true;
			}
			//------------------------------------------------------------------
			if (s_geo_type.compare("POLYLINE") == 0)
			{
				//			CGLPolyline *ply = NULL;
				//			ply = GEOGetPLYByName(s_geo_name);// get Polyline by
				// name

				CGLPolyline* ply(polyline_vector[s_geo_id]);
				//            if (ply) {
				//               m_msh->GetNODOnPLY(ply, nodes_vector);
				//               for (size_t i = 0; i < nodes_vector.size();
				//               i++) {
				//                  ll = nodes_vector[i];
				//                  l = ll; //+ShiftInNodeVector;
				//                  m_krd->is_a_CCBC[l] = true;
				//               }
				//            }

				std::vector<GEOLIB::Polyline*> const* const ply_vec(
				    geo_obj.getPolylineVec(unique_name));
				GEOLIB::Polyline const* const polyline((*ply_vec)[s_geo_id]);
				double msh_min_edge_length(m_msh->getMinEdgeLength());
				m_msh->setMinEdgeLength(ply->epsilon);
				m_msh->GetNODOnPLY(polyline, nodes_vector);
				m_msh->setMinEdgeLength(msh_min_edge_length);
				for (size_t i = 0; i < nodes_vector.size(); i++)
					//					m_krd->is_a_CCBC[nodes_vector[i]+ShiftInNodeVector]
					//= true;
					m_krd->is_a_CCBC[nodes_vector[i]] = true;
			}  // if(POLYLINE)
			//------------------------------------------------------------------
			if (s_geo_type.compare("SURFACE") == 0)
			{
				Surface* m_surface = NULL;
				m_surface = GEOGetSFCByName(s_geo_name);
				if (m_surface)
				{
					m_msh->GetNODOnSFC(m_surface, nodes_vector);
					for (size_t i = 0; i < nodes_vector.size(); i++)
					{
						long ll = nodes_vector[i];  //+ShiftInNodeVector;
						m_krd->is_a_CCBC[ll] = true;
					}
				}
			}
		}  // end of for(j=0;j<m_krd->NoReactGeoName.size()...

		nodes_vector.clear();

		/********************************************************/
		// Set up vectors switched_off_node for individual reactions
		m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
		if (m_msh == NULL)
		{
			std::cout << "No mesh in KRConfig"
			          << "\n";
			exit(1);
		}
		// for all reactions
		for (int k = 0; k < m_krd->NumberReactions; k++)
		{
			m_kr = KinReact_vector[k];
			if (!m_kr->NotThisReactGeoName.empty())
			{
				// Initialize vector is_a_CCBC
				for (size_t l = 0; l < mesh_node_vector_size; l++)
					m_kr->switched_off_node.push_back(false);
				// Go through specified geometry elements
				for (size_t j = 0; j < m_kr->NotThisReactGeoName.size(); j++)
				{
					s_geo_name = m_kr->NotThisReactGeoName[j];
					s_geo_type = m_kr->NotThisReactGeoType[j];
					s_geo_id = m_kr->NotThisReactGeoID[j];
					//------------------------------------------------------------------
					if (s_geo_type.compare("POINT") == 0)
					{
						const std::vector<GEOLIB::Point*>* pnt_vec(
						    geo_obj.getPointVec(unique_name));
						m_kr->switched_off_node[m_msh->GetNODOnPNT(
						    (*pnt_vec)[m_kr->NotThisReactGeoID[j]])] = true;
					}
					//------------------------------------------------------------------
					if (s_geo_type.compare("POLYLINE") == 0)
					{
						//					CGLPolyline
						//*ply(polyline_vector[s_geo_id]);
						//					if (ply) {
						//						if (ply->getType() == 100) //WW
						//							m_msh->GetNodesOnArc(ply,
						// nodes_vector);
						//						else m_msh->GetNODOnPLY(ply,
						// nodes_vector);
						//						for (i = 0; i < (long)
						// nodes_vector.size(); i++) {
						//							ll = nodes_vector[i];
						//							l = ll;
						////+ShiftInNodeVector;
						//							m_kr->switched_off_node[l] =
						// true;
						//						}
						//					}

						std::vector<GEOLIB::Polyline*> const* const ply_vec(
						    geo_obj.getPolylineVec(unique_name));
						GEOLIB::Polyline const* const polyline(
						    (*ply_vec)[s_geo_id]);
						m_msh->GetNODOnPLY(polyline, nodes_vector);
						for (size_t i = 0; i < nodes_vector.size(); i++)
						{
							long ll = nodes_vector[i];  //+ShiftInNodeVector;
							m_krd->is_a_CCBC[ll] = true;
						}
					}

					if (s_geo_type.compare("SURFACE") == 0)
					{
						Surface* m_surface = NULL;
						m_surface = GEOGetSFCByName(s_geo_name);
						if (m_surface)
						{
							m_msh->GetNODOnSFC(m_surface, nodes_vector);
							for (size_t i = 0; i < nodes_vector.size(); i++)
							{
								long ll =
								    nodes_vector[i];  //+ShiftInNodeVector;
								m_kr->switched_off_node[ll] = true;
							}
						}
					}
					//------------------------------------------------------------------
				}  // loop over j NotThisReactGeoName elements
				nodes_vector.clear();
			}  // if NotThisReactGeoName.size > 0
		}      // loop over k nreactions

		/********************************************************/
		// Get foc average of connected elements
		MeshLib::CNode* m_nod = NULL;
		MeshLib::CElem* m_ele = NULL;
		for (size_t l = 0; l < mesh_node_vector_size; l++)
			m_krd->node_foc.push_back(1.0);
		for (size_t l = 0; l < mesh_node_vector_size; l++)
		{
			m_nod = m_msh->nod_vector[l];
			double ww(0.0);
			double w(0.0);
			for (size_t i = 0; i < m_nod->getConnectedElementIDs().size(); i++)
			{
				m_ele = m_msh->ele_vector[m_nod->getConnectedElementIDs()[i]];
				m_mat_mp = mmp_vector[m_ele->GetPatchIndex()];
				ww += m_mat_mp->foc * m_ele->GetVolume();
				w += m_ele->GetVolume();
			}
			double foc(ww / w);
			//// CB dirty fix Brand project
			// if(foc>0.1)
			//  foc = 1;
			// else
			//  foc = 1e-100;
			m_krd->node_foc[l] = foc;
		}

		/********************************************************/
		// CB Set neighborhood relations for reaction deactivation switch
		if ((m_krd->ReactDeactEpsilon) < 0)
			m_krd->ReactDeactFlag = false;
		else
			m_krd->ReactDeactFlag = true;
		if (m_krd->ReactDeactFlag)
		{
			for (size_t l = 0; l < mesh_node_vector_size; l++)
			{
				// initialize all nodes as active
				m_krd->ReactDeact.push_back(false);
				// initialize all reaction terms as zero
				m_krd->React_dCdT.push_back(0);
				m_nod = m_msh->nod_vector[l];  // get current node object

				//  // loop over direct neighbour nodes
				// //this routine neglects some of the secondary neighbor
				// element nodes
				//  for(ll=0;ll<(long)m_nod->connected_nodes.size();ll++){
				//    dnnode = m_nod->connected_nodes[ll];          // get index
				//    of direct neighbour
				//    /* // check vs. all previously identified nodes to prevent
				//    duplicate entries , not necessary, as node is contained in
				//    connected node list of neighbours
				//    duplicate = 0; // flag
				//    for(nn=0;nn<(long)ReactNeighborNodes.size();nn++){
				//      if(dnnode == ReactNeighborNodes[nn]){ // node has
				//      already been found, not fresh
				//        duplicate = 1;
				//        break;
				//      }
				//    }
				//    if(duplicate==0) // it is fresh node
				//      ReactNeighborNodes.push_back(dnnode); // push back
				//      secondary neighbour index
				//    else
				//      duplicate = 0;  // no fresh direct neighbour node,
				//      reinitialize
				//    m_dnnod = m_msh->nod_vector[dnnode];          // in any
				//    case, get direct neighbour node object
				//    // loop over secondary neighbour nodes linked to a direct
				//    neighbour node
				//    for(lll=0;lll<(long)m_dnnod->connected_nodes.size();lll++){
				//      snnode = m_dnnod->connected_nodes[lll];     // get index
				//      of secondary neighbour
				//      // check vs. all previously identified nodes to prevent
				//      duplicate entries
				//      duplicate2 = 0; // flag, initialize
				//      for(nn=0;nn<(long)ReactNeighborNodes.size();nn++){
				//        if(snnode == ReactNeighborNodes[nn]){ // node has
				//        already been found
				//          duplicate2 = 1;
				//          break;
				//        }
				//      }
				//      if(duplicate2==0) // fresh node
				//        ReactNeighborNodes.push_back(snnode); // push back
				//        secondary neighbour index
				//      else
				//        duplicate2 = 0; // no fresh secondary neighbour node,
				//        reinitialize
				//    }
				//  }

				// loop over PRIMARY NEIGHBOR elements
				// this routine takes longer, but gets all the nodes
				for (size_t ll = 0; ll < m_nod->getConnectedElementIDs().size();
				     ll++)
				{
					// get index of direct neighbour element
					size_t dnele_idx = m_nod->getConnectedElementIDs()[ll];
					// get direct neighbour element object
					m_dnele = m_msh->ele_vector[dnele_idx];
					// get the neighbor element node indices
					m_dnele->GetNodeIndeces(primnnodesindices);
					// get the neighbor element number of nodes
					int nnodpneigh = m_dnele->GetNodesNumber(false);
					// loop over primary neighbour element number of nodes
					for (int lll = 0; lll < nnodpneigh; lll++)
					{
						// get the current node index
						long pneighnod_idx = primnnodesindices[lll];
						// get the node object
						m_dnnod = m_msh->nod_vector[pneighnod_idx];
						// loop over the connected elements, this now includes
						// the SECONDARY NEIGHBORS
						for (size_t llll = 0;
						     llll < m_dnnod->getConnectedElementIDs().size();
						     llll++)
						{
							// get the current connected element indices
							size_t snele_idx =
							    m_dnnod->getConnectedElementIDs()[llll];
							// get secondary neighbour element object
							m_snele = m_msh->ele_vector[snele_idx];
							// get the neighbor element node indices
							m_snele->GetNodeIndeces(secnnodesindices);
							// get the neighbor element number of nodes
							int nnodsneigh = m_snele->GetNodesNumber(false);
							// loop over secondary neighbour element number of
							// nodes
							for (int lllll = 0; lllll < nnodsneigh; lllll++)
							{
								bool duplicate(false);  // flag, initialize
								int annode_idx = secnnodesindices[lllll];
								// check vs. all previously identified nodes to
								// prevent duplicate entries
								for (size_t nn = 0;
								     nn < ReactNeighborNodes.size();
								     nn++)
									// node has already been found
									if (annode_idx == ReactNeighborNodes[nn])
									{
										duplicate =
										    true;  // set flag to "not fresh"
										break;     // skip rest of comparisons
									}
								if (!duplicate)  // fresh node
									// push back node index
									ReactNeighborNodes.push_back(annode_idx);
								else
									duplicate = false;  // was no fresh node,
								                        // reinitialize flag
							}
						}
					}
				}
				// Add local node neighborhood indices vector to global vector
				// This is the same for both above neighborhood search routines
				m_krd->ReactNeighborhood.push_back(ReactNeighborNodes);
				ReactNeighborNodes.clear();
			}
			// debug
			// for (nn=0;nn<m_krd->ReactNeighborhood.size();nn++){
			//  cout << nn << " " << m_krd->ReactNeighborhood[nn].size() << " "
			//  ;
			//  for (ll=0;ll<m_krd->ReactNeighborhood[nn].size();ll++)
			//    cout << m_krd->ReactNeighborhood[nn][ll] << " " ;
			//  cout << endl;
			//}
			m_krd->concentrationmatrix = (double**)malloc(
			    ((long)mesh_node_vector_size) * sizeof(double*));
			for (size_t ix = 0; ix < mesh_node_vector_size; ix++)
				m_krd->concentrationmatrix[ix] =
				    (double*)malloc((length) * sizeof(double));
			for (size_t ix = 0; ix < mesh_node_vector_size; ix++)
				for (size_t ixx = 0; ixx < length; ixx++)
					m_krd->concentrationmatrix[ix][ixx] = 0;
		}

		if (m_krd->debugoutflag)
			m_krd->debugoutfilename = FileName + "_KR_Debug_out.txt";
	}  // if KinReactData_vector.size() > 0
}

/**************************************************************************
   FEMLib-Method:
   Task: Clear KinReact_vector
   Programing:
   02/2006 SB Implementation
   last modified:
**************************************************************************/
void KRCDelete(void)
{
	long i;
	int no_krc = (int)KinReact_vector.size();
	for (i = 0; i < no_krc; i++)
		delete KinReact_vector[i];
	KinReact_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task: Reaction class read function
   Programing:
   05/2004 SB Implementation - adapted from OK rf_bc_new
   02/2006 SB New Class structure, IO, C++, FEM
   06/2010 TF formated, removed unused variables, changed signature for new
GEOLIB data structures
**************************************************************************/
bool CKinReact::Read(std::ifstream* rfd_file,
                     const GEOLIB::GEOObjects& geo_obj,
                     const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string, line_str1;
	bool new_keyword = false, new_subkeyword = false;
	std::stringstream in;
	MonodSubstruct* m_monod = NULL, * m_inhibit = NULL, * m_production = NULL;
	size_t index, index1;
	double dval;
	// CB 10/09
	std::string s_geo_type, s_geo_name;
	string thresh_species;
	double thresh_conc, thresh_ord;
	bool found = false;
	int i;
	// clear vectors
	monod.clear();
	inhibit.clear();
	production.clear();
	reactionpartner.clear();
	stochmet.clear();

	while (!new_keyword)
	{
		index = rfd_file->tellg();
		if (!GetLineFromFile(line, rfd_file)) break;
		line_string = line;
		if (line_string.find("#") != string::npos)
		{
			new_keyword = true;
			rfd_file->seekg(index);  // Dateipointer zur�cksetzen, sonst ist das
			                         // n�chste keyword weg
			break;
		}
		/* Keywords nacheinander durchsuchen */
		//....................................................................
		// subkeyword found
		if (line_string.find("$NAME") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> name;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$TYPE") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> type;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$BACTERIANAME") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> bacteria_name;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$EQUATION") != string::npos)
		{
			line_str1 = GetLineFromFile1(rfd_file);
			ReadReactionEquation(line_str1);
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$RATECONSTANT") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> rateconstant >> rateorder;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$GROWTH") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> grow;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$MONODTERMS") != string::npos)
			while (!new_subkeyword)
			{
				index1 = rfd_file->tellg();
				line_str1 = GetLineFromFile1(rfd_file);
				// Check for end of data block
				if ((line_str1.find("#") != string::npos) ||
				    (line_str1.find("$") != string::npos))
				{
					if (line_str1.find("#") != string::npos) new_keyword = true;
					rfd_file->seekg(index1);  // Dateipointer zur�cksetzen,
					                          // sonst ist das n�chste
					                          // subkeyword weg
					break;
				}
				m_monod = new MonodSubstruct();
				number_monod++;
				in.str(line_str1);
				in >> m_monod->species >> m_monod->concentration >>
				    m_monod->order;
				if ((m_monod->order != -99.0) &&
				    (m_monod->concentration != -1.0e9))  // check for read in
					monod.push_back(m_monod);
				else
				{
					DisplayMsgLn(" ERROR reading Monod Terms  - skipping");
					number_monod--;
					delete m_monod;
				}
				in.clear();
			}
		//....................................................................
		// subkeyword found
		if (line_string.find("$THRESHHOLDTERMS") != string::npos)
			while (!new_subkeyword)
			{
				index1 = rfd_file->tellg();
				line_str1 = GetLineFromFile1(rfd_file);
				// Check for end of data block
				if ((line_str1.find("#") != string::npos) ||
				    (line_str1.find("$") != string::npos))
				{
					if (line_str1.find("#") != string::npos) new_keyword = true;
					rfd_file->seekg(index1);  // Dateipointer zur�cksetzen,
					                          // sonst ist das n�chste
					                          // subkeyword weg
					break;
				}
				in.str(line_str1);
				in >> thresh_species >> thresh_conc >> thresh_ord;
				// check for read in
				if ((thresh_ord != -99.0) && (thresh_conc != -1.0e9))
				{
					found = false;
					for (i = 0; i < number_monod; i++)
						// find the respective monod species
						if (thresh_species.compare(monod[i]->species) == 0)
						{
							monod[i]->threshhold = true;  // and store the data
							monod[i]->threshConc = thresh_conc;
							monod[i]->threshOrder = thresh_ord;
							found = true;
							break;
						}
					if (found == false)
					{
						cout << " WARNING: no matching MONOD SPECIES found in "
						        "reaction ";
						cout << name << " for THRESHHOLD TERM SPECIES "
						     << thresh_species << endl;
					}
				}
				else
					DisplayMsgLn(" ERROR reading Threshhold Terms  - skipping");
				in.clear();
			}
		//....................................................................
		// subkeyword found
		if (line_string.find("$INHIBITIONTERMS") != string::npos)
			while (!new_subkeyword)
			{
				index1 = rfd_file->tellg();
				line_str1 = GetLineFromFile1(rfd_file);
				// Check for end of data block
				if ((line_str1.find("#") != string::npos) ||
				    (line_str1.find("$") != string::npos))
				{
					if (line_str1.find("#") != string::npos) new_keyword = true;
					rfd_file->seekg(index1);  // Dateipointer zur�cksetzen,
					                          // sonst ist das n�chste
					                          // subkeyword weg
					break;
				}
				in.str(line_str1);
				m_inhibit = new MonodSubstruct();
				number_inhibit++;
				in >> m_inhibit->species;
				in >> m_inhibit->concentration;
				in >> m_inhibit->order;
				if ((m_inhibit->order != -99.0) &&
				    (m_inhibit->concentration != -1.0e9))  // check for read in
					inhibit.push_back(m_inhibit);
				else
				{
					delete m_inhibit;
					DisplayMsgLn(" ERROR reading Inhibition Terms  - skipping");
					number_inhibit--;
				}

				in.clear();
			}
		//....................................................................
		// subkeyword found
		if (line_string.find("$PRODUCTIONTERMS") != string::npos)
			while (!new_subkeyword)
			{
				index1 = rfd_file->tellg();
				line_str1 = GetLineFromFile1(rfd_file);
				// Check for end of data block
				if ((line_str1.find("#") != string::npos) ||
				    (line_str1.find("$") != string::npos))
				{
					if (line_str1.find("#") != string::npos) new_keyword = true;
					rfd_file->seekg(index1);  // Dateipointer zur�cksetzen,
					                          // sonst ist das n�chste
					                          // subkeyword weg
					break;
				}
				in.str(line_str1);
				m_production = new MonodSubstruct();
				number_production++;
				in >> m_production->species >> m_production->concentration >>
				    m_production->order;
				if ((m_production->order != -99.0)
				    // check for read in
				    &&
				    (m_production->concentration != -1.0e9))
					production.push_back(m_production);
				else
				{
					DisplayMsgLn(" ERROR reading Production Terms  - skipping");
					number_production--;
				}
				in.clear();
			}

		//....................................................................
		// subkeyword found
		if (line_string.find("$PRODUCTIONSTOCH") != string::npos)
			while (!new_subkeyword)
			{
				index1 = rfd_file->tellg();
				line_str1 = GetLineFromFile1(rfd_file);
				// Check for end of data block
				if ((line_str1.find("#") != string::npos) ||
				    (line_str1.find("$") != string::npos))
				{
					if (line_str1.find("#") != string::npos)
						new_subkeyword = true;
					rfd_file->seekg(index1);  // Dateipointer zur�cksetzen,
					                          // sonst ist das n�chste
					                          // subkeyword weg
					break;
				}
				in.str(line_str1);
				m_production = new MonodSubstruct();
				in >> m_production->species >> m_production->concentration;
				ProdStochhelp.push_back(m_production);
				in.clear();
			}
		//....................................................................
		// CB Iso routine for reading in the isotope couples
		// subkeyword found
		if (line_string.find("$ISOTOPE_FRACTIONATION") != string::npos)
		{
			number_isotope_couples++;
			in.str(GetLineFromFile1(rfd_file));
			in >> Isotope_light >> Isotope_heavy >> isoenfac;
			in.clear();
			degType = "isotope_fractionation";  // CB besser in KRConfig??
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$EXCHANGE_PARAMETERS") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> dval;  // velocity
			ex_param.push_back(dval);
			in >> dval;  // kd
			ex_param.push_back(dval);
			if (exType.compare("langmuir") == 0) in >> exSurfaceID;
			if (exType.compare("freundlich") == 0)
			{
				in >> dval;
				ex_param.push_back(dval);
			}

			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$SORPTION_TYPE") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> exType;
			in.clear();
		}

		//....................................................................
		// subkeyword found
		if (line_string.find("$NAPL_PROPERTIES") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> blob_name >> Csat_pure >> Density_NAPL;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$SWITCH_OFF_GEOMETRY") != string::npos)
			while (!new_subkeyword)
			{
				index1 = rfd_file->tellg();
				line_str1 = GetLineFromFile1(rfd_file);
				// Check for end of data block
				if ((line_str1.find("#") != string::npos) ||
				    (line_str1.find("$") != string::npos))
				{
					if (line_str1.find("#") != string::npos) new_keyword = true;
					rfd_file->seekg(index1);  // Dateipointer zurueksetzen,
					                          // sonst ist das naehste
					                          // subkeyword weg
					break;
				}
				in.str(line_str1);
				in >> s_geo_type >> s_geo_name;
				NotThisReactGeoType.push_back(s_geo_type);
				NotThisReactGeoName.push_back(s_geo_name);
				size_t geo_obj_idx(std::numeric_limits<size_t>::max());

				if (s_geo_type.find("POINT") != std::string::npos)
					// TF 06/2010 - get the point vector and set the geo_obj_idx
					if (!((geo_obj.getPointVecObj(unique_name))
					          ->getElementIDByName(s_geo_name, geo_obj_idx)))
					{
						std::cerr
						    << "error in CKinReact::Read: point name not found!"
						    << "\n";
						exit(1);
					}
				if (s_geo_type.find("POLYLINE") != std::string::npos)
					// TF 07/2010 - get the polyline vector and set the
					// geo_obj_idx
					if (!((geo_obj.getPolylineVecObj(unique_name))
					          ->getElementIDByName(s_geo_name, geo_obj_idx)))
					{
						std::cerr << "error in CKinReact::Read: polyline name "
						             "not found!"
						          << "\n";
						exit(1);
					}
				NotThisReactGeoID.push_back(geo_obj_idx);

				in.clear();
			}
	}
	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Check consistency of reactions read in from input or database file
   Programing:
   02/2006 SB Implementation
   last modified:
**************************************************************************/
int CKinReact::CheckReactionDataConsistency()
{
	int i, ok = 1, length;
	bool found = false;
	string name, name1;
	//    CRFProcess *m_pcs = NULL;
	//    int no_processes =(int)pcs_vector.size();

	// Check number of reaction partners
	if (number_reactionpartner < 2)
	{
		// SB todo		ok=0;
		// SB todo		cout << " Less than two reaction partners in equation
		// found"
		// << endl;
	}
	// Ratekonstant zero?
	if ((fabs(rateconstant) < MKleinsteZahl) && (type.compare("monod") == 0))
	{
		//#ds erstetzt statt:	if((fabs(rateconstant) < MKleinsteZahl) &&
		//(type.compare("exchange")!=0) ){
		ok = 0;
		cout << " Rateconstant is zero" << endl;
	}
	// Number of partners higher than in arrays?
	length = (int)reactionpartner.size();  // Muss so gemacht werden, weil sonst
	                                       // signed/unsigned warning kommt
	if (number_reactionpartner != length)
	{
		ok = 0;
		cout << " Not enough reaction partners" << endl;
		for (i = 0; i < length; i++)
			cout << reactionpartner[i] << "  ";
		cout << endl;
	}
	length = (int)stochmet.size();
	if (stochmet.size() != reactionpartner.size())
	{
		ok = 0;
		cout << " Not enough stochieometric coefficients for reaction partners"
		     << endl;
		for (i = 0; i < length; i++)
			cout << stochmet[i] << "  ";
		cout << endl;
	}
	// check type
	if (type.compare("monod") != 0)
	{
		if (type.compare("exchange") != 0)
			if (type.compare("NAPLdissolution") != 0)
			{
				ok = 0;
				cout << "Unknown reaction type" << endl;
			}
	}

	/* Check Monod-, Inhibition and Production terms */
	length = (int)monod.size();
	if (length < number_monod) number_monod = (int)monod.size();
	length = (int)inhibit.size();
	if (length < number_inhibit) number_inhibit = (int)inhibit.size();
	length = (int)production.size();
	if (length < number_production) number_production = (int)production.size();

	if ((number_monod + number_inhibit + number_production == 0) &&
	    (type == "monod"))
	{
		//		cout << "Warning: no monod terms specified for monod equation
		//"<<
		// endl;
	}

	// Test concentrations and constants for > 0.0
	length = (int)monod.size();
	for (i = 0; i < length; i++)
		if (monod[i]->concentration < MKleinsteZahl)
		{
			cout << " Monod Concentration of reaction " << name
			     << " smaller than 0.0  " << endl;
			ok = 0;
		}
	length = (int)inhibit.size();
	for (i = 0; i < length; i++)
		if (inhibit[i]->concentration < MKleinsteZahl)
		{
			cout << " Inhibition Concentration of reaction " << name
			     << " smaller than 0.0  " << endl;
			ok = 0;
		}
	//#ds check NAPL dissolution reactions
	if (type.compare("NAPLdissolution") == 0)
	{
		// does blob-name exist?
		CKinBlob* m_kb = NULL;
		for (i = 0; i < (int)KinBlob_vector.size(); i++)
		{
			m_kb = KinBlob_vector[i];
			if (blob_name.compare(m_kb->name) == 0) found = true;
		}
		if (found == false)
		{
			ok = 0;
			cout << " Blob_Name " << blob_name << " defined in Reaction "
			     << name << " does not exist " << endl;
		}
		// Csat_pure >= 0
		if (Csat_pure < -MKleinsteZahl)
		{
			ok = 0;
			cout << " Invalid maximum solubility Csat_pure: " << Csat_pure
			     << " in Reaction " << name << endl;
		}
		// Density_NAPL >0
		if (Density_NAPL < MKleinsteZahl)
		{
			ok = 0;
			cout << " Invalid NAPL density Density_NAPL: " << Density_NAPL
			     << " in Reaction " << name << endl;
		}
	}  // end type=NAPLdissolution

	return ok;
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction class read function
   Programing:
   05/2004 SB Implementation - adapted from OK rf_bc_new
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReact::Write(ofstream* rfe_file)
{
	int i, flag = 0, length;

	// Write Keyword
	*rfe_file << "#REACTION" << endl;
	// Name of reaction
	*rfe_file << "$NAME" << endl
	          << name << endl;
	// Type of reaction
	*rfe_file << "$TYPE" << endl
	          << type << endl;
	// bacteria name
	if (type == "monod")
		*rfe_file << "$BACTERIANAME" << endl
		          << bacteria_name << endl;
	if (type == "exchange")
		*rfe_file << "$SORPTION_TYPE" << endl
		          << exType << endl;
	// ReactionEquation
	*rfe_file << "$EQUATION" << endl;
	for (i = 0; i < number_reactionpartner; i++)
	{
		if (stochmet[i] < 0.0)  // left side of equation
		{
			if (i == 0)
				*rfe_file << " " << fabs(stochmet[i]) << " "
				          << reactionpartner[i];
			else
				*rfe_file << " + " << fabs(stochmet[i]) << " "
				          << reactionpartner[i];
		}
		if (stochmet[i] > 0 && (flag > 0))  // remaining right hand side
			*rfe_file << " + " << fabs(stochmet[i]) << " "
			          << reactionpartner[i];
		if (stochmet[i] > 0 &&
		    (flag == 0))  // " = " Sign and first term on right hand side
		{
			*rfe_file << " = " << fabs(stochmet[i]) << " "
			          << reactionpartner[i];
			flag = 1;
		}
	}
	*rfe_file << endl;
	// Rateconstant and order
	if (type == "monod")
	{
		*rfe_file << "$RATECONSTANT" << endl
		          << rateconstant << "   " << rateorder << endl;
		*rfe_file << "$GROWTH" << endl
		          << grow << endl;
		// Monod terms
		*rfe_file << "$MONODTERMS" << endl;  //<< number_monod << endl;
		for (i = 0; i < number_monod; i++)
			*rfe_file << monod[i]->species << "  " << monod[i]->concentration
			          << "  " << monod[i]->order << endl;
		// Inhibition terms
		*rfe_file << "$INHIBITIONTERMS" << endl;  // << number_inhibit << endl;
		for (i = 0; i < number_inhibit; i++)
			*rfe_file << inhibit[i]->species << "  "
			          << inhibit[i]->concentration << "  " << inhibit[i]->order
			          << endl;
		// Production Terms
		//*rfe_file << "$PRODUCTIONTERMS" << endl << number_production << endl;
		// for(i=0;i<number_production;i++)
		//	*rfe_file << production[i]->species << "  " <<
		// production[i]->concentration << "  " << production[i]->order << endl;
		// Production Terms
		length = (int)ProdStochhelp.size();
		*rfe_file << "$PRODUCTIONSTOCH" << endl;  // << length << endl;
		for (i = 0; i < length; i++)
			*rfe_file << ProdStochhelp[i]->species << "  "
			          << ProdStochhelp[i]->concentration << "  " << endl;
		if (degType == "isotope_fractionation")
		{
			*rfe_file << "$ISOTOPE_FRACTIONATION" << endl;
			*rfe_file << Isotope_light << "  " << Isotope_heavy << "  "
			          << isoenfac << endl;
		}
		//#ds output f�r NAPL-dissolution
		//	*rfe_file << "$NAPL_PROPERTIES" << endl;
		//	*rfe_file << "blob_name " << blob_name << " Csat_pure " << Csat_pure
		//<< " Density_NAPL " << Density_NAPL << endl;
	}
	if (type == "exchange")
	{
		*rfe_file << "EXCHANGE_PARAMETERS" << endl;
		length = (int)ex_param.size();
		for (i = 0; i < length; i++)
			*rfe_file << ex_param[i] << "  ";
		*rfe_file << endl;
	}
	*rfe_file << endl;
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction class read function
   Programing:
   05/2004 SB Implementation
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReact::ReadReactionEquation(string line_string_all)
{
	string line_string, name, helpstring, c, calt;
	string substring, subsubstring;
	int indexlow, indexhigh, help1, strl, linelength, ih1, ih2, ih3,
	    indexgleich, i, partners, p;
	double zahl, vz = -1.0;

	stochmet.clear();
	reactionpartner.clear();

	// Anf�ngliche Leerzeichen, Kommentar am Ende raus
	ih1 = (int)line_string_all.find_first_not_of(" ", 0);
	ih2 = (int)line_string_all.find_first_of(";", ih1);
	ih3 = (int)line_string_all.length();
	line_string = line_string_all.substr(ih1, ih2 - ih1);

	// Auf mehrfache Leerzeichen und Tabs �berpr�fen - ersetzt durch je ein
	// Leerzeichen
	ih3 = (int)line_string.length();
	for (i = 0; i < ih3; i++)
	{
		c = line_string[i];
		if (c == "	") c = " ";  // replace tab by space
		if ((c == " ") && (calt == " "))
			ih1++;  // nothing
		else
		{
			helpstring.append(c);
			calt = c;
		}
	}

	// Leerzeichen am Ende raus
	ih1 = (int)helpstring.find_first_not_of(" ", 0);
	ih2 = (int)helpstring.find_last_of(" ");
	ih3 = (int)helpstring.length();
	if (ih2 == ih3 - 1)  // Leerzeichen am Ende
		line_string = helpstring.substr(ih1, ih2 - ih1);
	else
		line_string = helpstring.substr(ih1, ih3 - ih1);

	// Count reaction partners by looking for " + "
	linelength = (int)line_string.length();
	indexhigh = 0;
	partners = 0;
	ih1 = (int)line_string.find(" = ");
	if (ih1 < 0) DisplayMsgLn(" Error in keyword REACTION");
	// Exception handling
	while (indexhigh < linelength)
	{
		ih1 = (int)line_string.find(" + ", indexhigh);
		if (ih1 > 0)
		{
			partners++;
			indexhigh = ih1 + 3;
		}
		else
			// no " + " found
			indexhigh = linelength;
	}
	// Number of partners is 2 for " = " and one additional for each " + "
	number_reactionpartner = partners + 2;

	/* Zerlegen der Gleichung in Bl�che mit "Zahl Namen"*/

	indexhigh = 0;
	p = 0;  // Zaehlt partner hoch
	indexlow = indexhigh;
	linelength = (int)line_string.length();

	indexgleich = (int)line_string.find(" = ");

	while (indexhigh < linelength)
	{
		/* einen Block holen  - find next occurence of +, -, = */
		ih1 = (int)line_string.find(" + ", indexlow);
		if (ih1 < 0) ih1 = linelength;
		ih2 = (int)line_string.find(" - ", indexlow);
		if (ih2 < 0) ih2 = linelength;
		ih3 = (int)line_string.find(" = ", indexlow);
		if (ih3 < 0) ih3 = linelength;
		indexhigh = min(ih1, min(ih2, ih3));
		if (indexhigh > indexgleich)
			vz = 1.0;  // rechte Seite der Gleichung: Vorzeichenwechsel
		substring = line_string.substr(indexlow, indexhigh - indexlow);

		/* Leerzeichen drin ? */
		help1 = (int)substring.find(" ", 0);
		strl = (int)substring.length();
		if (help1 > 0)  // Leerzeichen gefunden
		{
			subsubstring = substring.substr(0, help1);
			//	if(subsubstring.find(".",0)>0) //double value given
			zahl = atof(subsubstring.data()) * vz;
			stochmet.push_back(zahl);
			subsubstring = substring.substr(help1 + 1, strl - help1 - 1);
			// reactionpartner[p] = subsubstring.data();
			//	sprintf(reactionpartner[p],"%s", subsubstring.data());
			reactionpartner.push_back(subsubstring);
		}  // kein Leerzeichen gefunden, substring enth�lt nur den Namen => Zahl
		   // wird 1 gesetzt
		else
		{
			zahl = 1.0 * vz;
			stochmet.push_back(zahl);
			//	reactionpartner[p] = substring.data();
			//	sprintf(reactionpartner[p],"%s", subsubstring.data());
			reactionpartner.push_back(substring);
		}
		p++;  // Reactionpartner hochzaehlen
		// cout << zahl <<" "<< name << " ";

		indexlow = indexhigh + 3;  // add length of " + "
	}                              // end while

	if (p != number_reactionpartner)
		cout << "Error: Parser for kinetic equations found varying number of "
		        "reaction partners" << endl;
}

/***************************************************************************
   FEMLib-Method:
   Task: CKinBlob constructor
   Programing:
   02/2007 DS Implementation
***************************************************************************/
CKinBlob::CKinBlob(void)
{
	// default= nonsense values for input, will be checked in CKinCheck()
	d50 = -1.;
	Sh_factor = -1.;
	Re_expo = -1.;
	Sc_expo = -1.;
	Geometry_expo = -1.;
	Mass = 0.;
	Volume = 0.;
	Masstransfer_k = 0.;
	current_Interfacial_area = 0.;
	BlobGeoType.clear();
	BlobGeoName.clear();
	Area_Value.clear();
	Interfacial_area.clear();
}

/***************************************************************************
   FEMLib-Method:
   Task: CKinBlob destructor
   Programing:
   02/2007 DS Implementation
***************************************************************************/
CKinBlob::~CKinBlob(void)
{
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function for CKinBlob-Structure
   Programing:
   02/2007 DS Implementation
**************************************************************************/
bool CKinBlob::Read(ifstream* rfd_file,
                    const GEOLIB::GEOObjects& geo_obj,
                    const std::string& unique_name)
{
	char line[MAX_ZEILE];
	string line_string, line_str1, s_geo_type, s_geo_name;
	string hash("#"), dollar("$");
	bool new_keyword = false, OK = true;
	size_t index, index1;
	double d_inivalue;
	std::stringstream in;

	//========================================================================
	while (!new_keyword)
	{
		index = rfd_file->tellg();
		//    if(!rfd_file->getline(line,MAX_ZEILE)) break;
		if (!GetLineFromFile(line, rfd_file)) break;
		line_string = line;
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			rfd_file->seekg(index);  // Dateipointer zur�cksetzen, sonst ist das
			                         // n�chste keyword weg
			break;
		}
		/* Keywords nacheinander durchsuchen */
		//....................................................................
		// subkeyword found
		if (line_string.find("$NAME") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> name;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$D50") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> d50;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$CALC_SHERWOOD") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> Sh_factor >> Re_expo >> Sc_expo;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$GEOMETRY") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> Geometry_expo;
			in.clear();
		}
		// subkeyword found
		if (line_string.find("$INTERFACIAL_AREA") != string::npos)
			while (OK)
			{
				index1 = rfd_file->tellg();
				if (!GetLineFromFile(line, rfd_file)) break;
				line_str1 = line;
				if ((line_str1.find(hash) != string::npos) ||
				    (line_str1.find(dollar) != string::npos))
				{
					OK = false;
					rfd_file->seekg(index1);  // Dateipointer zur�cksetzen,
					                          // sonst ist das n�chste keyword
					                          // weg
					break;
				}
				in.str(line_str1);
				//		in >> s_geo_type >> s_geo_name >> d_inivalue ;
				in >> s_geo_type;

				BlobGeoType.push_back(s_geo_type);

				size_t geo_obj_idx(std::numeric_limits<size_t>::max());
				if (s_geo_type.compare("DOMAIN") == 0)
				{
					BlobGeoName.push_back("domain");
					in >> d_inivalue;
					Area_Value.push_back(d_inivalue);
				}
				else
				{
					in >> s_geo_name >> d_inivalue;
					BlobGeoName.push_back(s_geo_name);

					if (s_geo_type.find("POINT") != std::string::npos)
						// TF 06/2010 - get the point vector and set the
						// geo_obj_idx
						if (!((geo_obj.getPointVecObj(unique_name))
						          ->getElementIDByName(s_geo_name,
						                               geo_obj_idx)))
						{
							std::cerr << "error in CKinBlob::Read: point name "
							             "not found!"
							          << "\n";
							exit(1);
						}

					if (s_geo_type.find("POLYLINE") != std::string::npos)
						// TF 06/2010 - get the polyline vector and set the
						// geo_obj_idx
						if (!((geo_obj.getPolylineVecObj(unique_name))
						          ->getElementIDByName(s_geo_name,
						                               geo_obj_idx)))
						{
							std::cerr << "error in CKinBlob::Read: polyline "
							             "name not found!"
							          << "\n";
							exit(1);
						}

					Area_Value.push_back(d_inivalue);
				}
				BlobGeoID.push_back(geo_obj_idx);
				in.clear();
			}
	}  // end while keyword

	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function for CKinBlob-Structure
   Programing:
   02/2007 DS Implementation
   06/2010 TF changed signature
    (in order to use the Write method for general output streams, for example
std::cout)
**************************************************************************/
void CKinBlob::Write(std::ostream& rfe_file) const
{
	rfe_file << endl;
	rfe_file << "#BLOB_PROPERTIES" << endl;
	rfe_file << "$NAME" << endl
	         << name << endl;
	rfe_file << "$D50" << endl
	         << d50 << endl;
	rfe_file << "$CALC_SHERWOOD" << endl;
	rfe_file << " Sh_factor : " << Sh_factor << endl;
	rfe_file << " Re_expo   : " << Re_expo << endl;
	rfe_file << " Sc_expo   : " << Sc_expo << endl;
	rfe_file << "$GEOMETRY" << endl
	         << Geometry_expo << endl;
	rfe_file << "$INTERFACIAL_AREA" << endl;
	for (size_t j = 0; j < BlobGeoName.size(); j++)
	{
		rfe_file << " Geotype : " << BlobGeoType[j];
		rfe_file << "   Geoname : " << BlobGeoName[j];
		rfe_file << "   Value   : " << Area_Value[j] << endl;
	}
	rfe_file << endl;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ configure function
   Programing:
   05/2007 DS Implementation
**************************************************************************/
void KBlobConfig(const GEOLIB::GEOObjects& geo_obj,
                 const std::string& unique_name)
{
	std::string s_geo_name, s_geo_type;
	size_t s_geo_id;
	std::vector<long> nodes_vector;

	// create vector for Interfacial_area and initialization with input value
	// (Interfacial_area[0]=input value)
	CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
	if (m_msh == NULL)
	{
		cout << "No mesh in KRConfig" << endl;
		exit(1);
	}

	for (size_t i = 0; i < KinBlob_vector.size(); i++)
	{
		CKinBlob* m_kb(KinBlob_vector[i]);
		for (size_t l = 0; l < m_msh->nod_vector.size(); l++)
			m_kb->Interfacial_area.push_back(-1.);  // Vorbelegung mit Area=-1

		for (size_t j = 0; j < m_kb->BlobGeoName.size(); j++)
		{
			s_geo_name = m_kb->BlobGeoName[j];
			s_geo_type = m_kb->BlobGeoType[j];
			s_geo_id = (m_kb->getBlobGeoID())[j];
			if (s_geo_type.compare("POINT") == 0)
			{
				// 06/2010 TF - switch to new GEOLIB - REMOVE CANDIDATE
				//				CGLPoint* m_geo_point = NULL; // make new GEO
				// point
				//				m_geo_point = GEOGetPointByName(s_geo_name);//Get
				//GEO
				// point by name
				//				if (m_geo_point)
				//					l = m_msh->GetNODOnPNT(m_geo_point); // +
				// ShiftInNodeVector; // find MSH point number stored in l

				const std::vector<GEOLIB::Point*>* pnt_vec(
				    geo_obj.getPointVec(unique_name));
				size_t msh_node_id =
				    m_msh->GetNODOnPNT((*pnt_vec)[(m_kb->getBlobGeoID())[j]]);

				m_kb->Interfacial_area[msh_node_id] = m_kb->Area_Value[j];
			}  // end if POINT

			if (s_geo_type.compare("POLYLINE") == 0)
			{
				//            CGLPolyline *ply (polyline_vector[s_geo_id]);
				//            if (ply)
				//            {
				//               if (ply->getType() == 100)         //WW
				//                  m_msh->GetNodesOnArc(ply, nodes_vector);
				//               else
				//                  m_msh->GetNODOnPLY(ply, nodes_vector);
				//               for (size_t k = 0; k < nodes_vector.size();
				//               k++) {
				//					size_t msh_node_id = nodes_vector[k];
				//					m_kb->Interfacial_area[msh_node_id] =
				// m_kb->Area_Value[j];
				//				}
				//            }

				std::vector<GEOLIB::Polyline*> const* const ply_vec(
				    geo_obj.getPolylineVec(unique_name));
				GEOLIB::Polyline const* const polyline((*ply_vec)[s_geo_id]);
				m_msh->GetNODOnPLY(polyline, nodes_vector);
				for (size_t k = 0; k < nodes_vector.size(); k++)
					m_kb->Interfacial_area[nodes_vector[k]] =
					    m_kb->Area_Value[j];
			}  // end if POLYLINE

			if (s_geo_type.compare("SURFACE") == 0)
			{
				Surface* m_surface = NULL;
				m_surface = GEOGetSFCByName(s_geo_name);
				if (m_surface)
				{
					m_msh->GetNODOnSFC(m_surface, nodes_vector);
					for (size_t k = 0; k < nodes_vector.size(); k++)
					{
						//+ShiftInNodeVector;
						size_t msh_node_id = nodes_vector[k];
						m_kb->Interfacial_area[msh_node_id] =
						    m_kb->Area_Value[j];
					}
				}
			}  // end if SURFACE

			if (s_geo_type.compare("DOMAIN") == 0)
				for (size_t k = 0; k < m_msh->nod_vector.size(); k++)
					m_kb->Interfacial_area[k] = m_kb->Area_Value[j];
			// end if SURFACE
		}
	}  // end loop over i = Number of blobs

	nodes_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task: Check OBJ configuration for input errors
   Programing:
   05/2007 DS Implementation
**************************************************************************/
void KBlobCheck(void)
{
	int i;
	long k, l;

	CKinBlob* m_kb = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt

	cout << endl
	     << "Checking defined blob classes" << endl;
	for (i = 0; i < (int)KinBlob_vector.size(); i++)
	{
		m_kb = KinBlob_vector[i];
		cout << m_kb->name << endl;

		if (m_kb->d50 <= 0.)
		{
			cout << "  Warning: D50 <= 0., D50 is set to 1.E-4" << endl;
			m_kb->d50 = 1.E-4;
		}
		if (m_kb->Sh_factor < 0.)
		{
			cout << "  Warning: Sh_factor < 0., Sh_factor is set to 1.15"
			     << endl;
			m_kb->Sh_factor = 1.15;
		}
		if (m_kb->Re_expo < 0.)
		{
			cout << "  Warning: Re_expo < 0., Re_expo is set to 0.654" << endl;
			m_kb->Re_expo = 0.654;
		}
		if (m_kb->Sc_expo < 0.)
		{
			cout << "  Warning: Sc_expo < 0., Sc_expo is set to 0.486" << endl;
			m_kb->Sc_expo = 0.486;
		}
		if (m_kb->Geometry_expo < 0.)
		{
			cout
			    << "  Warning: Geometry_expo < 0., Geometry_expo is set to 0.66"
			    << endl;
			m_kb->Geometry_expo = 0.66;
		}

		k = 0;
		for (l = 0; l < (long)m_msh->nod_vector.size(); l++)
		{
			if (m_kb->Interfacial_area[l] < 0.)
			{
				cout << "  Warning: no value for node " << l
				     << " , Interfacial area is set to zero." << endl;
				m_kb->Interfacial_area[l] = 0.;
			}
			else
				k += 1;
		}
		cout << "  Values for interfacial Surfaces have been defined for " << k
		     << " nodes by user." << endl;
	}  // end loop i over blob classes
}

/**************************************************************************
   Reaction-Method:
   Task: Class constructor
   Programing:
   02/2006 SB Implementation
**************************************************************************/
CKinReactData::CKinReactData(void)
{
	SolverType = 1;
	relErrorTolerance = 1.0e-10;
	minTimestep = 1.0e-10;
	initialTimestep = 1.0e-10;
	NumberReactions = 0;
	NumberFreundlich = 0;
	NumberLangmuir = 0;
	NumberLinear = 0;
	NumberMonod = 0;
	NumberNAPLdissolution = 0;
	usedt = -1.0;
	maxBacteriaCapacity = -1.0;
	is_a_bacterium.clear();
	testoutput = false;  // true;
	maxSurfaces = 3;
	exSurface.clear();
	// sp_index.clear();
	// kr_active_species = -1;
	sp_varind.clear();
	// sp_pcsind.clear(); //HS
	// das Surface-Array k�nnte auch dynamisch allokiert werden
	for (int i = 0; i < maxSurfaces; i++)
		exSurface.push_back(-1.0);
	NoReactGeoName.clear();
	NoReactGeoType.clear();
	is_a_CCBC.clear();
	node_foc.clear();

	ReactDeactMode = -1;
	ReactDeactEpsilon = -1.0;  // CB ReactDeact
	ReactDeactFlag = false;
	ReactDeactPlotFlag = 0;
	ReactDeact.clear();  // flags for individual nodes
	React_dCdT.clear();  // Sum of reaction rates for individual nodes
	ReactNeighborhood
	    .clear();  // node indices of local neighborhood around individual nodes

	debugoutflag = false;
}

/**************************************************************************
   Reaction-Method:
   Task: Class destructor
   Programing:
   02/2006 SB Implementation
**************************************************************************/
CKinReactData::~CKinReactData(void)
{
}

/**************************************************************************
   FEMLib-Method:
   Task: Reaction class data read function
   Programing:
   02/2006 SB New Class structure, IO, C++, FEM
**************************************************************************/
bool CKinReactData::Read(ifstream* rfd_file,
                         const GEOLIB::GEOObjects& geo_obj,
                         const std::string& unique_name)
{
	char line[MAX_ZEILE];
	string sub_line;
	string line_string, line_str1, help("1");
	string delimiter(" ");
	bool new_keyword = false, OK = true;
	string hash("#"), dollar("$");
	std::stringstream in;
	size_t index, index1;
	int /* count_surf, */ surf_id;
	string s_geo_type, s_geo_name;

	//========================================================================
	while (!new_keyword)
	{
		index = rfd_file->tellg();
		if (!GetLineFromFile(line, rfd_file)) break;
		line_string = line;
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			rfd_file->seekg(index);  // Dateipointer zur�cksetzen, sonst ist das
			                         // n�chste keyword weg
			break;
		}
		/* Keywords nacheinander durchsuchen */
		//....................................................................
		// subkeyword found
		if (line_string.find("$SOLVER_TYPE") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> SolverType;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$RELATIVE_ERROR") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> relErrorTolerance;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$MIN_TIMESTEP") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> minTimestep;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$INITIAL_TIMESTEP") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> initialTimestep;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$BACTERIACAPACITY") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> maxBacteriaCapacity;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$SURFACES") != string::npos)
			while (OK)
			{
				index1 = rfd_file->tellg();
				if (!GetLineFromFile(line, rfd_file)) break;
				line_str1 = line;
				if ((line_str1.find(hash) != string::npos) ||
				    (line_str1.find(dollar) != string::npos))
				{
					OK = false;
					rfd_file->seekg(index1);  // Dateipointer zur�cksetzen,
					                          // sonst ist das n�chste keyword
					                          // weg
					break;
				}
				in.str(line_str1);
				in >> surf_id;
				in >> exSurface[surf_id];
				in.clear();
			}
		// in.str(GetLineFromFile1(rfd_file));
		//    in >> count_surf;
		// in.clear();
		// for(int i=0;i<count_surf;i++){
		// in.str(GetLineFromFile1(rfd_file));
		// in >> surf_id;
		// in >> exSurface[surf_id];
		// in.clear();
		// }
		//....................................................................
		// subkeyword found
		if (line_string.find("$NO_REACTIONS") != string::npos)
			while (OK)
			{
				index1 = rfd_file->tellg();
				if (!GetLineFromFile(line, rfd_file)) break;
				line_str1 = line;
				if ((line_str1.find(hash) != string::npos) ||
				    (line_str1.find(dollar) != string::npos))
				{
					OK = false;
					rfd_file->seekg(index1);  // Dateipointer zur�cksetzen,
					                          // sonst ist das n�chste keyword
					                          // weg
					break;
				}
				in.str(line_str1);
				in >> s_geo_type >> s_geo_name;
				NoReactGeoType.push_back(s_geo_type);
				NoReactGeoName.push_back(s_geo_name);

				// TF 06/2010 - for the change string-ID to size_t-ID
				size_t geo_obj_idx(std::numeric_limits<size_t>::max());

				if (s_geo_type.find("POINT") != std::string::npos)
					// get the point vector and set the geo_obj_idx
					if (!((geo_obj.getPointVecObj(unique_name))
					          ->getElementIDByName(s_geo_name, geo_obj_idx)))
					{
						std::cerr << "error in CKinReactData::Read: (type="
						          << s_geo_type << "): " << s_geo_name
						          << " point name not found!"
						          << "\n";
						exit(1);
					}
				if (s_geo_type.find("POLYLINE") != std::string::npos)
					// get the point vector and set the geo_obj_idx
					if (!((geo_obj.getPolylineVecObj(unique_name))
					          ->getElementIDByName(s_geo_name, geo_obj_idx)))
					{
						std::cerr
						    << "error in CKinReactData::Read: polyline name "
						    << s_geo_name << " not found!"
						    << "\n";
						exit(1);
					}
				NoReactGeoID.push_back(geo_obj_idx);

				in.clear();
			}
		//....................................................................
		// subkeyword found
		if (line_string.find("$REACTION_DEACTIVATION") != string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> ReactDeactMode >> ReactDeactEpsilon >> ReactDeactPlotFlag;
			in.clear();
		}
		// subkeyword found
		if (line_string.find("$DEBUG_OUTPUT") != string::npos)
			debugoutflag = true;
		//....................................................................
	}
	usedt = DMAX(initialTimestep, minTimestep);
	return true;
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction data class write function
   Programing:
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReactData::Write(ofstream* rfe_file)
{
	int i;
	// Write Keyword
	*rfe_file << endl;
	*rfe_file << "#KINREACTIONDATA" << endl;
	*rfe_file << "$SOLVER_TYPE" << endl
	          << SolverType << endl;
	*rfe_file << "$RELATIVE_ERROR" << endl
	          << relErrorTolerance << endl;
	*rfe_file << "$MIN_TIMESTEP" << endl
	          << minTimestep << endl;
	*rfe_file << "$INITIAL_TIMESTEP" << endl
	          << initialTimestep << endl;
	*rfe_file << "$BACTERIACAPACITY" << endl
	          << maxBacteriaCapacity << endl;
	//*rfe_file << " max Surfaces : " << maxSurfaces << endl;
	*rfe_file << "$SURFACES" << endl;  //<< (int)exSurface.size() << endl;
	for (i = 0; i < (int)exSurface.size(); i++)
		*rfe_file << i + 1 << "  " << exSurface[i] << endl;
	*rfe_file << "NO_REACTIONS"
	          << endl;  // << (int)NoReactGeoName.size() << endl;
	for (i = 0; i < (int)NoReactGeoName.size(); i++)
		*rfe_file << NoReactGeoType[i] << "  " << NoReactGeoName[i] << endl;
	*rfe_file << "$REACTION_DEACTIVATION	" << endl
	          << ReactDeactMode << " " << ReactDeactEpsilon << " "
	          << ReactDeactPlotFlag << endl;
	*rfe_file << "$DEBUG_OUTPUT	" << endl
	          << debugoutflag << endl;
	//*rfe_file << " Number of reactions: " << NumberReactions << endl;
	//*rfe_file << " Number of linear exchange reactions: " << NumberLinear <<
	// endl;
	//*rfe_file << " Number of freundlich exchange reactions: " <<
	// NumberFreundlich << endl;
	//*rfe_file << " Number of langmuir exchange reactions: " << NumberLangmuir
	//<< endl;
	//*rfe_file << " is_a_bacterium: "  << endl;
	// for(i=0;i<is_a_bacterium.size();i++) *rfe_file << is_a_bacterium[i] << "
	// ";
	//*rfe_file << endl;
	*rfe_file << endl;
	/*
	   *rfe_file << " usedt "<< usedt << endl;
	   *rfe_file << " Number Reactions "<< NumberReactions << endl;
	   *rfe_file << " is_a_bacterium " << endl;
	   for(int i=0; i < (int)is_a_bacterium.size();i++)
	   *rfe_file <<  is_a_bacterium[i];
	   *rfe_file << endl;
	 */
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction data class write function
   Programing:
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReactData::TestWrite(void)
{
	int i;  // CB
	// Write Keyword
	cout << "#KINREACTIONDATA" << endl;
	cout << "$SOLVER_TYPE" << endl
	     << SolverType << endl;
	cout << "$RELATIVE_ERROR" << endl
	     << relErrorTolerance << endl;
	cout << "$MIN_TIMESTEP" << endl
	     << minTimestep << endl;
	cout << "$INITIAL_TIMESTEP" << endl
	     << initialTimestep << endl;
	cout << "$BACTERIACAPACITY" << endl
	     << maxBacteriaCapacity << endl;
	cout << "$REACTION_DEACTIVATION	" << endl
	     << ReactDeactMode << " " << ReactDeactEpsilon << " "
	     << ReactDeactPlotFlag << endl;
	cout << "$DEBUG_OUTPUT	" << endl
	     << debugoutflag << endl;

	cout << endl;
	cout << " usedt " << usedt << endl;
	cout << " Number Reactions " << NumberReactions << endl;
	cout << " is_a_bacterium " << (int)is_a_bacterium.size() << endl;
	for (i = 0; i < (int)is_a_bacterium.size(); i++)
		cout << is_a_bacterium[i] << " ";
	cout << endl;
	cout << " Exchange reactions : " << endl;
	// ACHTUNG: die folgende Zeile wird ausgegeben BEVOR NumberXXX berechnet
	// wurde !
	cout << " Linear exchange : " << NumberLinear << endl
	     << " Freundlich exchange : " << NumberFreundlich << endl
	     << " Langmuir exchange : " << NumberLangmuir << endl
	     << " NAPL dissolution : " << NumberNAPLdissolution << endl;
	for (i = 0; i < (int)exSurface.size(); i++)
		cout << " " << exSurface[i];
	cout << endl;
}

/**************************************************************************
   Reaction-Method:
   Task: KinReaction write function - echo of input values to rfe - file
   Programing:
   05/2004 SB Implementation
   02/2006 SB C++, IO
**************************************************************************/
bool KRWrite(std::string const& prot_name)
{
	CKinReact* m_kr = NULL;
	CKinReactData* m_krd = NULL;
	CKinBlob* m_kp = NULL;
	string rfe_file_name;

	if (KinReact_vector.empty() && KinBlob_vector.empty() &&
	    KinReactData_vector.empty())
		return true;  // NW skip writing a file if no KinReaction data are
	                  // loaded

	//========================================================================
	// File handling
	rfe_file_name = prot_name + "_echo";
	ofstream rfe_file(rfe_file_name.data(), ios::app);
	if (!rfe_file.good()) return false;
	rfe_file.seekp(0L, ios::end);  // go to end
	//========================================================================
	rfe_file << endl
	         << "; Reactions "
	            "--------------------------------------------------------------"
	            "--- " << endl;
	// Output all Reactions
	size_t length(KinReact_vector.size());
	for (size_t i = 0; i < length; i++)
	{
		m_kr = KinReact_vector[i];
		m_kr->Write(&rfe_file);
		if (KinReactData_vector[0]->testoutput) m_kr->TestWrite();
	}
	// Output all BlobProperties
	length = KinBlob_vector.size();
	for (size_t i = 0; i < length; i++)
	{
		m_kp = KinBlob_vector[i];
		m_kp->Write(rfe_file);
		m_kp->Write(std::cout);
	}
	// Output all kinetic reaction data
	if (!KinReactData_vector.empty()) m_krd = KinReactData_vector[0];
	if (m_krd != NULL)
	{
		m_krd->Write(&rfe_file);
		if (KinReactData_vector[0]->testoutput) m_krd->TestWrite();
	}
	rfe_file.close();

	return true;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ExecuteKinReact()                                 */
/* DS-TBC                                                                 */
/*                                                                        */
/* Aufgabe:                                                               */
/* Haupt-Subroutine zur Berechnung der mikrobiellen kinetischen Reaktionen*/
/* Aufruf in void ExecuteReactions(void)  (rf_react.cpp)                  */
/*                                                                        */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 05/2003     DS         Erste Version                                   */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/*                                                                        */
/**************************************************************************/
void CKinReactData::ExecuteKinReact(void)
{
	std::cout << " ExecuteKineticReactions"
	          << "\n";

	if (debugoutflag) debugoutstr.open(debugoutfilename.c_str());

	if (time_vector.size() > 0) dt = time_vector[0]->CalcTimeStep();

	if ((dt > 1.E-20) && (aktueller_zeitschritt > 0))
	{
		const size_t nnodes(
		    fem_msh_vector[0]->nod_vector.size());  // SB: ToDo hart gesetzt
		size_t save_node(0);
		int save_nok = 0, save_nbad = 0;
		double usedttmp = 1.E+30;
		/* Einstellungen Gleichungsl�ser f�r alle Knoten gleich */
		/* relative Genauigkeit des Gleichungsl�sers (eps< fehler/c0) */
		double eps = relErrorTolerance;
		/* min zulaessiger Zeitschritt*/
		double hmin = minTimestep;

		//	cout << " NumberReactions: " << NumberReactions << endl;
		if (NumberReactions > 0)
		{
			// CB Reaction deactivation for this time step
			if (ReactDeactFlag)
				if (aktueller_zeitschritt > 2) ReactionDeactivation(nnodes);
			// Check if nodes should be deactivated or activated for this time
			// step

			unsigned count = 0;
			for (size_t node = 0; node < nnodes; node++)
			{
				// cout << node << endl;
				// no reactions at Concentration BCs
				if (is_a_CCBC[node] == true)
				{
				}
				// CB no reactions at deactivated nodes
				else if ((ReactDeactFlag) && (ReactDeact[node] == true))
				{
				}
				else
				{
					double usedtneu(0.0);
					int nok(0);
					int nbad(0);
					Biodegradation(node, eps, hmin, &usedtneu, &nok, &nbad);
					if (usedtneu < usedttmp)
					{
						usedttmp = usedtneu;
						save_nok = nok;
						save_nbad = nbad;
						save_node = node;
					}
					count++;
				}
			}  // end for(node...

			// CB Reaction deactivation for next time step
			if (ReactDeactFlag)
			{
				cout << "    Kinetic reactions executed at " << count << " of "
				     << nnodes << " nodes." << endl;
				if (ReactDeactMode == 2)  // For mode 2 the C_new must be
					                      // updated by C_old (C after eraction
					                      // of last time step)
					ReactDeactSetOldReactionTerms(nnodes);
				Aromaticum(nnodes);
			}
		}  // end if(NumberRactions>0)

		if (usedttmp < usedt)
		{
			cout << endl
			     << "Kinetics in node " << save_node
			     << " limit integration step - nok: ";
			cout << save_nok << " nbad: " << save_nbad << endl;
		}

		// update des zul�ssigen Integrationsschritts, verwendet beim Aufruf von
		// odeint
		// kleinster Wert, der in einem der Knoten zu einer zuverl�ssigen
		// Integration gef�hrt hat
		// konservative, aber stabile Annahme
		usedttmp = DMAX(usedttmp, hmin);
		usedt = DMIN(usedttmp, dt);
		// cout << endl << " Next suggested integration step " << usedt << endl;
	}  // end if((dt>1.E-20)&&(aktueller_zeitschritt>0)){

	if (debugoutflag) debugoutstr.close();
	if (ReactDeactFlag)
		if (ReactDeactPlotFlag == 1) ReactDeactPlotFlagsToTec();
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Biodegradation(node, )                            */
/* DS-TBC                                                                 */
/*                                                                        */
/* Aufgabe:                                                               */
/* Berechnet Bioabbau im Knoten node in dt                                */
/* Tr�gt berechnete Senke in ?? ein                                       */
/*                                                                        */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)                 */
/* E: node (Knotennummer), tstart (Zeitpunkt Rechenbeginn)                */
/*    tend (Zeitpunkt Rechenende)                                         */
/*                                                                        */
/* Ergebnis:                                                              */
/* - void -                                                               */
/* Ergebnisse werden in Datenstruktur zur�ckgeschrieben                   */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 05/2003     DS         Erste Version                                   */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/* 05/2007     DS         NAPL-dissolution added                          */
/*                                                                        */
/**************************************************************************/

void CKinReactData::Biodegradation(long node, double eps, double hmin,
                                   double* usedtneu, int* nok, int* nbad)
{
	double* Concentration;
	double* newVolume;
	double nexth = 0.;
	long sp;
	//  int nok=0, nbad=0, Number_of_Components;
	int Number_of_Components, nreactions, r, Sp1, blob, Number_of_blobs;
	double Csat_max, DensityNAPL, DensityAQ, DiffusionAQ, ViscosityAQ,
	    // OK411
	    PoreVelocity = 0.0, d50, Reynolds, Schmidt, Sherwood;
	double tstart, tend, dt = 0.0;
	double baditerations;
	//  CRFProcess* m_pcs = NULL;
	CTimeDiscretization* m_tim = NULL;
	//  CompProperties *m_cp = NULL;
	string speciesname = " dummy";

	if (debugoutflag) debugoutstr << "Biodegradation node timestep " << flush;
	// debugoutstr << " --> Biodegradation" << endl << flush;

	CKinReact* m_kr = NULL;
	CKinBlob* m_kb = NULL;
	CKinReactData* m_krd = NULL;
	m_krd = KinReactData_vector[0];

	if (time_vector.size() > 0)
	{
		m_tim = time_vector[0];
		dt = m_tim->CalcTimeStep();
	}
	// Number_of_Components = kr_active_species; //
	Number_of_Components = (int)cp_vec.size();
	// Get storage for vector of concentrations
	Concentration = dvector(1, Number_of_Components);

	/* Konzentrationen aller Substanzen aus Datenstruktur auslesen und in neuem
	 * Array speichern */
	for (sp = 0; sp < Number_of_Components; sp++)
	{
		// HS, old:
		// Concentration[sp + 1] = pcs_vector[sp_pcsind[sp]]->GetNodeValue(node,
		// sp_varind[sp]);
		// new:
		Concentration[sp + 1] =
		    cp_vec[sp]->getProcess()->GetNodeValue(node, sp_varind[sp]);
		if (fabs(Concentration[sp + 1]) < 1.e-19) Concentration[sp + 1] = 0.0;
		// SB todo - ist abewr gerade eh noch ein dummy
		// ExchangeTerm[sp]=TBCGetExchange(node,sp)/(dt);
	}

	if (debugoutflag)
	{
		// debugoutstr << " Concentrations before odeint: " << endl << " " <<
		// flush;
		for (sp = 0; sp < Number_of_Components; sp++)
			// debugoutstr << pcs_vector[sp_pcsind[sp]]->nod_val_name_vector[0]
			debugoutstr << cp_vec[sp]->getProcess()->nod_val_name_vector[0]
			            << " " << flush;
		debugoutstr << "baditerations" << endl;

		debugoutstr << " " << node << " " << aktueller_zeitschritt << " "
		            << flush;
		for (sp = 0; sp < Number_of_Components; sp++)
			debugoutstr << Concentration[sp + 1] << " " << flush;
		debugoutstr << "-" << endl
		            << flush;
	}
	//#ds
	/* PREPARE PARAMETERS FOR NAPL-DISSOLUION*/
	/* calculate Mass and Volume of blobs for this node */

	// set all mass(blob)=0, volume(blob)=0
	Number_of_blobs = (int)KinBlob_vector.size();
	for (r = 0; r < Number_of_blobs; r++)
	{
		m_kb = KinBlob_vector[r];
		m_kb->Mass = 0.;
		m_kb->Volume = 0.;
	}

	nreactions = m_krd->NumberReactions;
	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;
		if (m_kr->typeflag_napldissolution)  // dissolution reaction identified
		{
			Sp1 = m_kr->ex_species[0] + 1;  // Sp1 = NAPL-species
			blob = m_kr->blob_ID;
			DensityNAPL = cp_vec[Sp1 - 1]->molar_density;
			// DensityNAPL = m_kr->Density_NAPL; // CB: this should be obtained
			// from comp properties
			m_kb = KinBlob_vector[blob];  // pointer to blob-properties set in
			                              // the reaction r
			if (Concentration[Sp1] > 0.)
			{
				m_kb->Mass += DMAX(Concentration[Sp1], 0.);
				m_kb->Volume += DMAX(Concentration[Sp1], 0.) / DensityNAPL;
				// Sb todo Achtung - das wird ja gar nicht zur�ckgespeichert...
			}
		}
	}  // end for nreactions

	/* calculate current Csat depending on Raoults law for this node */
	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;
		if (m_kr->typeflag_napldissolution)  // dissolution reaction identified
		{
			Sp1 = m_kr->ex_species[0] + 1;  // Sp1 = NAPL-species
			blob = m_kr->blob_ID;
			Csat_max = cp_vec[Sp1 - 1]->max_solubility;
			// Csat_max    = m_kr->Csat_pure;    // CB: this should be obtained
			// from comp properties
			m_kb = KinBlob_vector[blob];  // pointer to blob-properties set in
			                              // the reaction r
			if (m_kb->Mass > 0.)
				m_kr->current_Csat =
				    Csat_max * DMAX(Concentration[Sp1], 0.) / m_kb->Mass;
			else
				m_kr->current_Csat = Csat_max;  // keine NAPL-Masse vorhanden,
			                                    // NAPL-Bildung m�glich wenn
			                                    // c(singleSubstance) > Csat
		}
	}  // end for nreactions
	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;
		if (m_kr->typeflag_napldissolution)
			PoreVelocity = m_kr->GetNodePoreVelocity(node);
	}

	/* calculate current Masstransfer-coefficient k for this node */
	for (r = 0; r < Number_of_blobs; r++)
	{
		m_kb = KinBlob_vector[r];
		d50 = m_kb->d50;
		DiffusionAQ = mfp_vector[0]->getDiffusion();  // CB Todo: this should be
		                                              // a component property =>
		                                              // Sherwood is component
		                                              // dependent
		DensityAQ = mfp_vector[0]->Density();
		ViscosityAQ = mfp_vector[0]->Viscosity();
		Reynolds = DensityAQ * PoreVelocity * d50 / ViscosityAQ;
		Schmidt = ViscosityAQ / DiffusionAQ / DensityAQ;
		Sherwood = m_kb->Sh_factor * pow(Reynolds, m_kb->Re_expo) *
		           pow(Schmidt, m_kb->Sc_expo);
		// k in m/s
		m_kb->Masstransfer_k = Sherwood * DiffusionAQ / d50;
	}

	/* save current Interfacial areas for this node */
	for (r = 0; r < Number_of_blobs; r++)
	{
		m_kb = KinBlob_vector[r];
		m_kb->current_Interfacial_area = m_kb->Interfacial_area[node];
	}

	tstart = DMAX(aktuelle_zeit - dt, 0.);
	tend = aktuelle_zeit;
//  tstart=tstart/86400.0; tend = tend/86400.0 ;  // alte Version: hier wurde
//  nur im kinetischen Teil mit Tagen gerechnet
//  cout << " times: " << tstart << ", " << tend << endl;
/* Aufruf Gleichungsl�ser */
/* eigentliche Rechenroutinen sind "derivs" und "jacobn" (namen fest
   vorgegeben),
   die vom Gleichungsl�ser aufgerufen werden */
#if 0  // TODO fix compililaion error
	odeint(Concentration, Number_of_Components, tstart, tend, eps, usedt, hmin,
	       &nexth, nok, nbad, derivs, stifbs, node);
#endif

	baditerations = double(*nbad) / double(*nok + *nbad);

	if (baditerations < 0.001)
	{
		// fehlerfreie Integration, zeitschritt kann vergr��ert werden
		if (nexth > usedt)
			*usedtneu = DMAX(nexth, usedt * 2.);
		else
			*usedtneu = usedt * 1.5;
	}
	else
	{
		// Integrationsfehler, zeitschritt beibehalten oder verkleinern
		if (*nbad == 1)
			*usedtneu = DMAX(nexth, usedt * 1.10);
		else if (*nok > *nbad * 2)
			*usedtneu = DMAX(nexth, usedt * 1.01);
		else
			*usedtneu = DMAX(nexth, usedt / 5.);
	}

	// update results
	for (sp = 0; sp < Number_of_Components; sp++)
	{
		// Notl�sung gegen das vollst�ndige Absterben der Bakterien
		if ((is_a_bacterium[sp]) && (Concentration[sp + 1] < 1.E-30))
			Concentration[sp + 1] = 1.E-30;
		// Konzentrationen aller Substanzen in Datenstruktur zur�ckschreiben
		cp_vec[sp]->getProcess()->SetNodeValue(node, sp_varind[sp],
		                                       Concentration[sp + 1]);
		// save exchange term SB todo
	}

	if (debugoutflag)
	{
		// debugoutstr << " Concentrations after odeint: " << endl << " " <<
		// flush;
		debugoutstr << " " << node << " " << aktueller_zeitschritt << " "
		            << flush;
		for (sp = 0; sp < Number_of_Components; sp++)
			debugoutstr << Concentration[sp + 1] << " " << flush;
		debugoutstr << baditerations << endl
		            << flush;
	}

	/* #ds calculate Interfacial areas for this node after dissolution for next
	 * time step */
	newVolume = dvector(0, Number_of_blobs);
	for (r = 0; r < Number_of_blobs; r++)
		newVolume[r] = 0.;
	nreactions = m_krd->NumberReactions;
	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;
		if (m_kr->typeflag_napldissolution)  // dissolution reaction identified
		{
			Sp1 = m_kr->ex_species[0] + 1;  // Sp1 = NAPL-species
			blob = m_kr->blob_ID;
			DensityNAPL = cp_vec[Sp1 - 1]->molar_density;
			// DensityNAPL = m_kr->Density_NAPL; // CB: this should be obtained
			// from comp properties
			newVolume[blob] += DMAX(Concentration[Sp1], 0.) / DensityNAPL;
		}
	}  // end for nreactions
	//  double dummy;
	for (r = 0; r < Number_of_blobs; r++)
	{
		m_kb = KinBlob_vector[r];
		if ((newVolume[r] > 0.) && (m_kb->Volume > 0.))
			m_kb->Interfacial_area[node] =
			    m_kb->current_Interfacial_area *
			    pow((newVolume[r] / m_kb->Volume), m_kb->Geometry_expo);
		else
			m_kb->Interfacial_area[node] = 1.E-20;  // residual interfacial area
		                                            // to allow re-building of
		                                            // phase
	}

	free_dvector(newVolume, 0, Number_of_blobs);
	free_dvector(Concentration, 1, Number_of_Components);
}

/*************************************************************************************/
/* Routine f�r Bulirsch-Stoer Gleichungsl�ser (steife ODEs) */
/* DS-TBC */
/*                                                                                   */
/* Berechnung der ersten Ableitungen �ber die Zeit dC/dt=... */
/*                                                                                   */
/* Input: */
/* t = aktuelle Zeit */
/* c[1..Number_of_Components] = aktuelle Konzentration der Substanzen in Zelle
 */
/* n =! Number_of_Components */
/*                                                                                   */
/* Output: */
/* dcdt[1..Number_of_Components] = aktuelle Ableitung der Konzentrationen nach
 * der   */
/*            Zeit, = Wachstum der Bakterien und Verbrauch der Substanzen */
/*                                                                                   */
/* Programmaenderungen: */
/* 05/2003     DS         Erste Version */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/* 05/2007     DS         NAPL-dissolution added */
/*                                                                                   */
/*************************************************************************************/

void derivs(double t, double c[], double dcdt[], int n, long node)
{
	t = t;  // OK411
	int i, r, nreactions, BacteriaNumber;
	int Sp1, Sp2, surfaceID = -1, blob;
	double BacteriaMass, BacGrowth, Yield, sumX = 0., maxkap;
	double porosity1, porosity2, exchange = 0.0, exch, kd, density1,
	                             saturation2, kadsorb, kdesorb, totalSurface,
	                             exponent, parameter,
	                             chochexp;  // OK411
	double dt;
	double foc;
	//#ds
	//	int blob;
	double Csat;
	//	double occupiedSurface[m_krd->maxSurfaces+1];
	vector<double> occupiedSurface;

	CKinReact* m_kr = NULL;
	// OK411 CKinReact *m_kr1 = NULL;
	CKinBlob* m_kb = NULL;
	CKinReactData* m_krd = NULL;
	m_krd = KinReactData_vector[0];

	// if(m_krd->debugoutflag)
	//  m_krd->debugoutstr << " derivs" << endl << flush;

	CTimeDiscretization* m_tim = NULL;
	m_tim = time_vector[0];
	dt = m_tim->CalcTimeStep();

	/* reset array with derivatives */
	/* ACHTUNG, unterschiedliche Indizierung der Arrays, c[1..n]
	 * BioDegradation[0..n-1] */
	for (i = 0; i < n; i++)
		// SBtodo		dcdt[i+1]=ExchangeTerm[i];
		dcdt[i + 1] = 0.0;

	/* calculate present bacteria capacity */
	sumX = 0.0;  // SB added
	maxkap = m_krd->maxBacteriaCapacity;
	if (maxkap > 1.E-30)
		for (i = 0; i < n; i++)
			if (m_krd->is_a_bacterium[i]) sumX = sumX + c[i + 1];

	/**********************************************************************************************/
	/* Anzahl der mikrobiellen Reaktionen aus Datenstruktur auslesen */
	nreactions = m_krd->NumberReactions;  // BioDegradation.NumberReactions;
	/* loop over reactions dX/dt= nymax * X * monodterms * inhibitionterms */
	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;
		if (m_kr->typeflag_monod)
		{
			BacteriaNumber = m_kr->bacteria_number + 1;
			BacteriaMass = c[BacteriaNumber];
			porosity1 = m_kr->GetReferenceVolume(BacteriaNumber - 1, node);
			if (BacteriaMass > 1.E-40)
			{
				m_kr->currentnode = node;  // CB 19/10/09 This is eclusively for
				                           // Brand model to allow porosity in
				                           // Inhibition constant calculation

				// This is where growth rate is computed
				BacGrowth = m_kr->BacteriaGrowth(r, c, sumX, -1);
				if (m_kr->grow) dcdt[BacteriaNumber] += BacGrowth;
				/* microbial consumption of substances */
				for (i = 0; i < n; i++)
				{
					Yield = m_kr->ProductionStoch[i];
					if (fabs(Yield) > 1.E-30)
					{
						porosity2 = m_kr->GetReferenceVolume(i, node);
						dcdt[i + 1] +=
						    BacGrowth * Yield * porosity1 / porosity2;
					}
				}
			}
		}  // type == monod
	}      // nreactions

	/**********************************************************************************************/
	/* Berechnung der Austauschprozesse */

	// calculate already occupied surface
	if (m_krd->NumberLangmuir > 0)
	{
		// Initialise Surfaces for langmuir isotherms
		for (i = 0; i < m_krd->maxSurfaces; i++)
			occupiedSurface.push_back(0.0);

		for (r = 0; r < nreactions; r++)
		{
			m_kr = KinReact_vector[r];
			// CB new reaction switch for individual reactions
			if (m_kr->switched_off_node.size() > 0)
				if (m_kr->switched_off_node[node] == true) continue;
			if ((m_kr->typeflag_exchange) && (m_kr->typeflag_exchange_langmuir))
			{
				Sp1 = m_kr->ex_species[0] + 1;
				surfaceID = m_kr->exSurfaceID;
				occupiedSurface[surfaceID] += c[Sp1];
			}
		}
	}  // if NumberLangmuir > 0

	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;
		if (m_kr->typeflag_exchange)
		{
			/* linearer Austausch ggf. mit kd */
			if (m_kr->typeflag_exchange_linear)
			{
				// Matrix
				Sp1 = m_kr->ex_species[0] + 1;
				porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
				density1 = m_kr->GetDensity(Sp1 - 1, node);
				// geloest
				Sp2 = m_kr->ex_species[1] + 1;
				porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);

				exch = m_kr->ex_param[0];
				kd = m_kr->ex_param[1];

				if (fabs(kd) < MKleinsteZahl)
				{
					// no kd, exchange between two species in solution
					exchange = exch * (c[Sp2] - c[Sp1]);
					dcdt[Sp1] += exchange / porosity1;
					dcdt[Sp2] += -exchange / porosity2;
				}
				else
				{
					// with kd, exchange between matrix (mol/kg) and solution
					// (mol/l)
					foc = m_krd->node_foc[node];

					if (foc > MKleinsteZahl) kd = kd * foc;
					// else
					//  kd = 0;
					exchange = exch * (c[Sp2] * kd - c[Sp1]);
					/* Die Abfrage verringert die Desorptionsgeschwindigkeit,
					 * wenn absehbar ist, dass Csorbiert im Negativen landet */
					if (-exchange * dt > c[Sp1]) exchange = -c[Sp1] / dt;

					dcdt[Sp1] += exchange;
					dcdt[Sp2] += -exchange * porosity1 / porosity2 * density1;
				}
			}  // ende if exType == linear

			/* Freundlich Kinetik */
			if (m_kr->typeflag_exchange_freundlich)
			{
				// Matrix
				Sp1 = m_kr->ex_species[0] + 1;
				porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
				density1 = m_kr->GetDensity(Sp1 - 1, node);
				// geloest
				Sp2 = m_kr->ex_species[1] + 1;
				porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);

				exponent = m_kr->ex_param[2];
				parameter = m_kr->ex_param[1];
				exch = m_kr->ex_param[0];

				if (c[Sp2] > residual)
					// no linearisation required
					chochexp = pow(c[Sp2], exponent);
				else
					// linearisation required due to instability of c^x if
					// c<residual
					chochexp =
					    (1. - exponent) * pow(residual, exponent) +
					    exponent * pow(residual, (exponent - 1)) * c[Sp2];

				exchange = exch * (parameter * chochexp - c[Sp1]);
				/* Die Abfrage verringert die Desorptionsgeschwindigkeit, wenn
				 * absehbar ist, dass Csorbiert im Negativen landet */
				if (-exchange * dt > c[Sp1]) exchange = -c[Sp1] / dt;

				dcdt[Sp1] += exchange;
				dcdt[Sp2] += -exchange * porosity1 / porosity2 * density1;
			}  // if freundlich

			/* Langmuir Kinetik */
			if (m_kr->typeflag_exchange_langmuir)
			{
				// Surfaces were initialized above
				//	for (i=0; i<nexchange; i++)	{
				Sp1 = m_kr->ex_species[0] + 1;
				porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
				density1 = m_kr->GetDensity(Sp1 - 1, node);
				Sp2 = m_kr->ex_species[1] + 1;
				porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);

				kadsorb = m_kr->ex_param[0];
				kdesorb = m_kr->ex_param[1];
				// SB_langmuir		surfaceID	= m_kr->exSurfaceID;
				// //Exchange.Langmuir[i].SurfaceID;
				totalSurface = m_krd->exSurface[m_kr->exSurfaceID];

				//      occupiedSurface was calculated above
				//		double occsurf = occupiedSurface[surfaceID];

				//#ds ACHTUNG hier muss sicher gestellt sein, dass Sp1 die
				// adsorbierte und Sp2 die gel�ste Species ist !
				exchange = kadsorb *
				               (totalSurface - occupiedSurface[surfaceID]) *
				               c[Sp2] -
				           kdesorb * c[Sp1];

				/* Die Abfrage verringert die Desorptionsgeschwindigkeit, wenn
				 * absehbar ist, dass Csorbiert im Negativen landet */
				if (-exchange * dt > c[Sp1]) exchange = -c[Sp1] / dt;

				dcdt[Sp1] += exchange;
				dcdt[Sp2] += -exchange * porosity1 / porosity2 * density1;
				//	}
			}  // ende if exType == langmuir
		}      // if type == exchange
	}          // for r

	//#ds
	/**********************************************************************************************/
	/* Berechnung NAPL-L�sung */
	/**********************************************************************************************/

	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;

		if (m_kr->typeflag_napldissolution)
		{
			Sp1 = m_kr->ex_species[0] +
			      1;  // Exchange.Linear[i].Species1;    Sp1 muss NAPL sein
			//		porosity1	= m_kr->GetPorosity(Sp1-1,node);
			//		density1	= m_kr->GetDensity(Sp1-1,node);
			////GetDensity(phase);
			blob = m_kr->blob_ID;
			m_kb = KinBlob_vector[blob];  // pointer to blob-properties set in
			                              // the reaction r

			Sp2 = m_kr->ex_species[1] +
			      1;  // Exchange.Linear[i].Species2;    Sp2 = mobile Phase
			// CB this includes the saturation
			porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
			//#ds TODO	    saturation2 = ??
			saturation2 = 1.;

			/**************************/
			// CB 040808 Saturation for aqueous phase
			// Is this necessary? Saturation is included in GetReferenceVolume!!
			// I think it's wrong ,so I deactivate 09/2009
			// CRFProcess *m_pcs = NULL;
			// m_pcs = PCSGet("TWO_PHASE_FLOW");
			// if(m_pcs->pcs_type_number==0)
			//  m_pcs = pcs_vector[m_pcs->pcs_number+1]; // this is the
			//  saturation equation
			// int idxs1 = m_pcs->GetNodeValueIndex("SATURATION1"); // Sat of
			// water phase
			// saturation2 = m_pcs->GetNodeValue(node, idxs1);
			/**************************/
			Csat = m_kr->current_Csat;  // Csat externally calculated in
			                            // Function Biodegradation
			// k * A externally calculated in Function Biodegradation
			exch = m_kb->Masstransfer_k * m_kb->current_Interfacial_area;

			if (exch > 0) exchange = exch * (Csat - c[Sp2]);

			/* Die Abfrage verringert die L�sungsgeschwindigkeit, wenn absehbar
			   ist, dass CNAPL im Negativen landet
			   Verhindert negative CNAPL-Konzentrationen */
			if (exchange * dt > c[Sp1]) exchange = c[Sp1] / dt;

			// no dissolution or NAPL mass present
			if ((exchange < 0.) || (c[Sp1] > MKleinsteZahl))
			{
				dcdt[Sp1] += -exchange;  // NAPL
				// concentration in solution, refers to mass / volume of water
				dcdt[Sp2] += exchange / porosity2 / saturation2;
			}
		}  // ifNAPL dissolution
	}      // loop ofer reactions r
	//#ds

	occupiedSurface.clear();
}

/**************************************************************************/
/* Phase bestimmen, in der sich Substanz i befindet                       */
/* DS-TBC                                                                 */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/**************************************************************************/
int CKinReact::GetPhase(int species)
{
	int phase = -1;
	CompProperties* cp_m = NULL;

	cp_m = cp_vec[species];
	if (cp_m != NULL)
		phase = cp_m->transport_phase;
	else
		cout << " Error: component does not exist !" << endl;
	return phase;
}

/**************************************************************************/
/* Porosit�t einer Phase bestimmen                                        */
/* DS-TBC                                                                 */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/* 05/2009     CB         Replaced by Function GetReferenceVolume         */
/**************************************************************************/
//
// double   CKinReact::GetPorosity( int comp, long index ){
//
// double poro = 0.0, theta = 1.0;
// long group, phase;
// CMediumProperties *m_mat_mp = NULL;
// CRFProcess *m_pcs=NULL;
// CKinReactData *m_krd = NULL;
//
// m_krd = KinReactData_vector[0];
//
//// Get process
// m_pcs = PCSGet("MASS_TRANSPORT",cp_vec[comp]->compname); //SB todo check
// m_pcs = pcs_vector[m_krd->sp_pcsind[comp]];
// theta = m_pcs->m_num->ls_theta;
// phase = cp_vec[comp]->transport_phase;
//
//// Get material properties of element
// group = 0; //SB todo m_pcs->m_msh->ele_vector[index]->GetPatchIndex();
// m_mat_mp = mmp_vector[group];
// if(phase ==0){ // water phase
//	poro = m_mat_mp->Porosity(index, theta);
//}
// else if(phase == 1){ // solid phase
//	poro = m_mat_mp->vol_mat;
//}
// else if (phase == 2){ // bio phase
//	poro = m_mat_mp->vol_bio;
//}
// else if (phase == 3){ // NAPL phase (refers to REV)
//	poro = 1.;
//}
// else
//	cout << " Error: No porosity found for phase " << phase << endl;
//
//// cout << " Get Porosity returns: vol_water: " << m_mat_mp->Porosity(index,
/// NULL, theta) << ", vol_bio: " << m_mat_mp->vol_bio <<", vol_mat: " <<
/// m_mat_mp->vol_mat << ", phase: " << phase << ", poro: " << poro << endl;
// return poro;
//}

/*****************************************************************************************/
/* Calculate the reference volume of a phase at a node */
/* DS-TBC */
/* 02/2006     SB         Introduced new C++ concept, Data structures */
/* 08/2008     DS         Consider saturation of water phase in case of
 * multiphase flow  */
/* 09/2009     CB         Heterogeneous Porosities update */
/*****************************************************************************************/
double CKinReact::GetReferenceVolume(int comp, long index)
{
	// Get process
	CRFProcess const* const pcs(cp_vec[comp]->getProcess());
	double theta(pcs->m_num->ls_theta), refvol;
	long phase = cp_vec[comp]->transport_phase;

	if (phase == 0)
	{
		int idx;
		refvol = GetPhaseVolumeAtNode(index, theta, phase);
		// water phase, reference volume might be less than pore space in case
		// of multiphase or richards flow
		// --> Get node saturation of mobile (water) phase, required for all
		// exchange processes
		double saturation = 1.0;  // default
		CRFProcess* pcs_flow = PCSGetFlow();
		if (pcs_flow->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
		{
			if (pcs_flow->pcs_type_number == 0)
				// this is the saturation equation
				pcs_flow = pcs_vector[pcs_flow->pcs_number + 1];
			// Sat of water phase
			idx = pcs_flow->GetNodeValueIndex("SATURATION1");
			saturation = pcs_flow->GetNodeValue(index, idx);
		}
		else if (pcs_flow->getProcessType() == FiniteElement::RICHARDS_FLOW)
		{
			// Sat of water phase
			idx = pcs_flow->GetNodeValueIndex("SATURATION1");
			saturation = pcs_flow->GetNodeValue(index, idx);
		}
		refvol *= saturation;
	}
	else if (phase == 3)  // NAPL phase (refers to REV)
		refvol = 1.0;
	else
		// solid or bio phase, 1 and 2
		refvol = GetPhaseVolumeAtNode(index, theta, phase);

	return refvol;
}

/**************************************************************************/
/* Return the volume fraction of a particular phase at a node             */
/* 0 pore space, 1 solid phase, 2 bio phase                               */
/* DS-TBC                                                                 */
/* 09/2009     CB         Introduced new C++ concept, Data structures     */
/* 09/2011	TF hanged access to coordinates of mesh node,
    - substituted access to mesh_element from pointer to direct access into the
   vector
    - made the mesh node a const pointer
    - made the pointer to the mesh const, made the mesh itself const
    - substituted pow(x,2) by x*x
    - reduced scope of loop variable i
    - moved declaration of gravity center to appropriate place			  */
/**************************************************************************/
double CKinReact::GetPhaseVolumeAtNode(long node_number, double theta,
                                       int phase)
{
	CFEMesh const* const mesh(fem_msh_vector[0]);  // SB: ToDo hart gesetzt

	size_t idx = 0;
	long elem;  // OK411
	double distance, weight, sum_w = 0;
	double vol = 0, poro = 0;

	// get Indices for phase 1 or 2, only if heterogeneous porosity model = 11,
	// i.e. vol_mat_model = vol_bio_model = 2
	long group = 0;  // SB todo group = m_ele->GetPatchIndex(); Todo CB
	CMediumProperties* m_mat_mp(mmp_vector[group]);
	if (m_mat_mp->vol_bio_model == 2 && m_mat_mp->vol_mat_model == 2)
	{
		switch (phase)
		{
			case 1:  // solid phase
				// Get VOL_MAT index
				for (idx = 0; idx < mesh->mat_names_vector.size(); idx++)
					if (mesh->mat_names_vector[idx].compare("VOL_MAT") == 0)
						break;
				break;
			case 2:  // bio phase
				// Get VOL_BIO index
				for (idx = 0; idx < mesh->mat_names_vector.size(); idx++)
					if (mesh->mat_names_vector[idx].compare("VOL_BIO") == 0)
						break;
				break;
			default:
				break;
		}
	}

	// Get node coordinates
	MeshLib::CNode const* const node(mesh->nod_vector[node_number]);
	double const* const coord = node->getData();  // Coordinates(coord);

	for (size_t el = 0; el < node->getConnectedElementIDs().size(); el++)
	{
		// initialize for each connected element
		distance = weight = poro = 0;
		// Get the connected element
		elem = node->getConnectedElementIDs()[el];  // element index
		MeshLib::CElem const* const m_ele(mesh->ele_vector[elem]);
		// get the phase volume of current element elem
		group = 0;  // group = m_ele->GetPatchIndex(); Todo CB
		m_mat_mp = mmp_vector[group];
		switch (phase)
		{
			case 0:  // pore space
				// CB Now provides also heterogeneous porosity, model 11
				poro = m_mat_mp->Porosity(elem, theta);
				break;
			case 1:                                // solid phase
				if (m_mat_mp->vol_mat_model == 1)  // homogeneous
					poro = m_mat_mp->vol_mat;
				else if (m_mat_mp->vol_mat_model == 2)  // CB heterogeneous
					poro = m_ele->mat_vector(idx);
				else
					cout << "Warning! No valid VOL_MAT model in "
					        "CKinReact::GetPhaseVolumeAtNode, vol_mat_model ="
					     << m_mat_mp->vol_mat_model << endl;
				break;
			case 2:                                // bio phase
				if (m_mat_mp->vol_bio_model == 1)  // homogeneous
					poro = m_mat_mp->vol_bio;
				else if (m_mat_mp->vol_bio_model == 2)  // CB heterogeneous
					poro = m_ele->mat_vector(idx);
				else
					cout << "Warning! No valid VOL_BIO model in "
					        "CKinReact::GetPhaseVolumeAtNode, vol_bio_model ="
					     << m_mat_mp->vol_bio_model << endl;
				break;
			case 3:  // NAPL phase (refers to REV)
				poro = 1.0;
				break;
			default:
				cout << "Error in CKinReact::GetPhaseVolumeAtNode: no valid "
				        "phase" << endl;
				break;
		}
		// calculate distance node <-> element center of gravity
		double const* grav_c(m_ele->GetGravityCenter());
		for (size_t i = 0; i < 3; i++)
			distance +=
			    (coord[i] - grav_c[i]) *
			    (coord[i] - grav_c[i]);  // pow((coord[i] - grav_c[i]), 2);
		// linear inverse distance weight = 1/(distance)
		distance =
		    sqrt(distance);  // for quadratic interpolation uncomment this line
		weight = (1 / distance);
		sum_w += weight;
		// add the weighted phase volume
		vol += poro * weight;
	}  // loop over connected elements

	// normalize weighted sum by sum_of_weights sum_w
	vol *= 1 / sum_w;

	return vol;
}

/**************************************************************************/
/* Dichte einer Phase bestimmen                                        */
/* DS-TBC                                                                 */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/**************************************************************************/

double CKinReact::GetDensity(int comp, long index)
{
	double dens = 0.;
	long group, phase;
	// WW CRFProcess *m_pcs = NULL;

	group = index;  // avoid warning
	                // Get process
	// WW m_pcs = PCSGet("MASS_TRANSPORT", cp_vec[comp]->compname);
	// theta = m_pcs->m_num->ls_theta;
	phase = cp_vec[comp]->transport_phase;

	// Get material properties of element
	group = 0;  // m_pcs->m_msh->ele_vector[index]->GetPatchIndex(); //SB todo

	if (phase == 0)  // component in water phase, return liquid density

		dens = mfp_vector[0]->Density();
	// component in solid phase, return solid phase density
	else if (phase == 1)
		dens = msp_vector[group]->Density();

	return dens;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: BacteriaGrowth(reaction)                          */
/* DS-TBC                                                                 */
/*                                                                        */
/* Aufgabe:                                                               */
/* Berechnet mikrobielles Wachstum aufgrund Reaktion reaction             */
/*                                                                        */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)                 */
/* E: reaction - Nummer der Wachstumsreaktion                             */
/*    Concentration - Array aktuelle Konzentrationen aus Gleichungsl�ser  */
/*    sumX - Bakteriendichte f�r Ber�cksichtigung MaxKapazit�t            */
/*                                                                        */
/* Ergebnis:                                                              */
/* R: Growth - Wachstum der zur Reaktion gehoerigen Bakteriengruppe       */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 05/2003     DS         Erste Version                                   */
/*                                                                        */
/**************************************************************************/

double CKinReact::BacteriaGrowth(int r, double* c, double sumX, int exclude)
{
	r = r;  // OK411
	int i, BacteriaNumber, MonodSpecies, InhibitionSpecies, Isotopespecies,
	    NumberMonod, NumberInhibition;
	double Growth, BacteriaMass, maxVelocity, maxkap;
	double MonodConcentration, InhibitionConcentration, C;
	double Ctot = 0;  // CB Isotope fractionation
	double MonodOrder;
	CKinReactData* m_krd = NULL;

	m_krd = KinReactData_vector[0];

	BacteriaNumber = bacteria_number + 1;
	maxVelocity = rateconstant;
	BacteriaMass = c[BacteriaNumber];

	Growth = maxVelocity * BacteriaMass;

	/* Hemmung durch Bakteriendichte sumX nur bei Wachstum */
	if ((sumX > 1.E-30) && (maxVelocity > 0.))
	{
		/* Max. Bakterienkapazit�t aus Datenstruktur auslesen */
		maxkap = m_krd->maxBacteriaCapacity;
		Growth = Growth * maxkap / (sumX + maxkap);
	}

	/*  FOR-Schleife �ber vorhandene Monodterme */

	NumberMonod = number_monod;
	for (i = 0; i < NumberMonod; i++)
		/* M�glichkeit zum Weglassen Monodterm f�r partielle Ableitungen
		 * erforderlich */
		if (i != exclude)
		{
			MonodSpecies = monod[i]->speciesnumber;
			MonodConcentration = monod[i]->concentration;
			MonodOrder = monod[i]->order;  // CB higher order Monod terms
			C = c[MonodSpecies + 1];
			// Growth= Growth * Monod(MonodConcentration,C);  // old formulation
			// without isotope fractionation
			// CB Isotope fractionation: hier muss Ctot mit �bergeben werden
			if ((typeflag_iso_fract == 1) &&
			    (monod[i]->isotopecouplenumber >= 0))
			{
				Isotopespecies = monod[i]->isotopecouplenumber;
				Ctot = C + c[Isotopespecies + 1];
				Growth =
				    Growth * Monod(MonodConcentration, C, Ctot, MonodOrder);
			}  // this is the standard case without fractionation
			else
			{
				Growth = Growth * Monod(MonodConcentration, C, C, MonodOrder);
				if (monod[i]->threshhold == true)
				{
					// now multiply by additional Threshhold Term, technically
					// this is the same as a Monod term
					MonodConcentration = monod[i]->threshConc;
					MonodOrder =
					    monod[i]->threshOrder;  // CB higher order Monod terms
					Growth =
					    Growth * Monod(MonodConcentration, C, C, MonodOrder);
				}
			}
		}
	//  for NumberMonod

	/*	FOR-Schleife �ber vorhandene Inhibitionsterme */
	NumberInhibition = number_inhibit;
	for (i = 0; i < NumberInhibition; i++)
	{
		InhibitionSpecies = inhibit[i]->speciesnumber;
		InhibitionConcentration = inhibit[i]->concentration;
		// ATTENTION!!!
		// CB 16/10/09 Hardcode fix of Fe3 inhibition concentration for Brand
		// model,
		// this parameter depends on porosity as in Min3P Fe3 is in solid phase
		// and Inhibition concentrations are expressed in terms of Volume
		// fraction
		// Vol_Fe3/Vol_BulkAquifer [m?m�]
		// while in Geosys, Fe3 is considered an immobile aqueous species
		// and concentrations were converted to [mol/L_water] by
		// C_gs = C_min3p * rho_Fe3 / molweight / porosity
		// for inhibition concentrations, division by porosity is still required
		if (inhibit[i]->species.compare("Fe3") == 0)
			InhibitionConcentration *=
			    1 / GetPhaseVolumeAtNode(currentnode, 1, 0);
		// CB  further changes in derivs (1), jacbn (2), class CKinReact{}
		C = c[InhibitionSpecies + 1];
		Growth = Growth * Inhibition(InhibitionConcentration, C);
	}

	return Growth;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Monod(MC, C)                                      */
/* DS-TBC                                                                 */
/*                                                                        */
/* Aufgabe:                                                               */
/* Berechnet modifizierten Monod-Term                                     */
/*                                                                        */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)                 */
/* E: MC - MonodConcentration                                             */
/*    C  - Concentration of substance                                     */
/*                                                                        */
/* Ergebnis:                                                              */
/* A: Monod - Wert des Monodterms                                         */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 05/2003     DS         Erste Version                                   */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/*                                                                        */
/**************************************************************************/
double CKinReact::Monod(double MC, double C, double Ctot, double order)
// double CKinReact::Monod ( double MC, double C )
{
	double Monod;
	// CKinReactData *m_krd = KinReactData_vector[0];
	if (C > 0.)
		// Monod=C/(MC+C);        // normaler Monodterm
		// Monod=C/(MC+Ctot);  // CB Isotope fractionation : Ctot = C in case of
		// no fractionation
		Monod = pow((C / (MC + Ctot)), order);  // CB higher order Monod terms
	                                            // --> factor order for partial
	                                            // derivatives
	else
		/* linearisierter Term fuer c<=0, with very small slope due to high
		 * concentrations */
		// Monod=(C/1000.)/MC;    //Monod=(C)/MC;   CB changed on behalf of SB,
		// 10.07.07
		Monod = pow(((C / 1000.) / MC), order);  // CB higher order Monod terms

	return Monod;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Inhibition(IC, C)
   DS-TBC                                                                 */
/* Aufgabe:
   Berechnet modifizierten Inhibitions-Term
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: IC - InhibitionConcentration
   C  - Concentration of substance
 */
/* Ergebnis:
   A: Inhibition - Wert des Inhibitionsterms
 */
/* Programmaenderungen:
   05/2003     DS         Erste Version									  */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/*                                                                          */
/**************************************************************************/

double CKinReact::Inhibition(double IC, double C)
{
	double Inhibition;
	//  CKinReactData *m_krd = KinReactData_vector[0];

	if (C > 0.) /* normaler Inhibitionsterm */
		Inhibition = IC / (IC + C);
	else
		/* linearisierter Term fuer C<=0 */
		Inhibition = 1.0;  // Inhibition=1.-C/IC;   CB changed due to stimulance
	                       // of growth for neg conc.

	return Inhibition;
}

/*************************************************************************************/
/* Routine f�r Bulirsch-Stoer Gleichungsl�ser (steife ODEs) */
/* DS-TBC */
/*                                                                                   */
/* Berechnung der Jacobi-Matrix, partielle Ableitungen von dC/dt aus Funktion */
/* derive nach */
/* a) der Zeit d2C/d2t = dfdt[] = 0. f�r mikrobiellen Abbau */
/* b) den Konzentrationen d2C/dt*dC = dfdc[][] */
/*                                                                                   */
/* Input: */
/* t = aktuelle Zeit */
/* c[1..Number_of_Components] = aktuelle Konzentration der Substanzen in Zelle
 */
/* n =! Number_of_Components */
/*                                                                                   */
/* Output: */
/* dfdt[1..Number_of_Components] = aktuelle Ableitung der Funktion nach der Zeit
 */
/* dfdc[1..Number_of_Components][1..Number_of_Components] = aktuelle Ableitung
 */
/*             der Funktion nach den Konzentrationen */
/*                                                                                   */
/* Programmaenderungen: */
/* 05/2003     DS         Erste Version */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/*                                                                                   */
/*************************************************************************************/

void jacobn(double t, double c[], double dfdt[], double** dfdc, int n,
            long node)
{
	t = t;  // OK411
	int i, j, r, nreactions, BacteriaNumber, NumberMonod, MonodSpecies,
	    NumberInhibition, InhibitionSpecies;
	int Sp1, Sp2, SpX, surfaceID = -1, surfaceID2, blob;
	double maxkap, BacteriaMass, sumX = 0., BacGrowth, maxVelocity, *d2X_dtdS;
	double CMonodSpecies, MonodConcentration, CInhibitionSpecies,
	    InhibitionConcentration, Yield;
	double porosity1, porosity2, exch, kd, density1, saturation2, kadsorb;
	double kdesorb, totalSurface, adsorb, exponent, parameter;
	// SBtodo	double occupiedSurface[maxSurfaces+1];
	vector<double> occupiedSurface;
	int IsotopeSpecies;
	double Ciso;
	double MonodOrder;
	double ThreshOrder;
	double ThreshConc;

	CKinReact* m_kr = NULL, * m_kr1 = NULL;  // OK411 *m_kr2=NULL;
	CKinBlob* m_kb = NULL;
	CKinReactData* m_krd = NULL;
	//	CMediumProperties *m_mat_mp = NULL;
	double foc;

	m_krd = KinReactData_vector[0];

	// if(m_krd->debugoutflag)
	//  m_krd->debugoutstr << " jacobn" << endl << flush;

	/* Hilfsvektor f�r partielle Ableitung des Bakterienwachstums nach Species S
	 */
	d2X_dtdS = dvector(1, n);

	/* weitere Ableitungen nach t dfdt[] alle null */
	/* Ableitungen nach c dfdc[][] werden inkrementiv berechnet, also erst alles
	 * null setzen */
	/* ACHTUNG, unterschiedliche Indizierung der Arrays, c[1..n]
	 * BioDegradation[0..n-1] */
	for (i = 0; i < n; i++)
	{
		dfdt[i + 1] = 0.;
		for (j = 0; j < n; j++)
			dfdc[i + 1][j + 1] = 0.;
	}

	/* calculate present bacteria capacity */
	maxkap = m_krd->maxBacteriaCapacity;
	// F�r Berechnung der Ableitungen f�r den Fall dass eine maximale Kapazit�t
	// ber�cksichtigt werden muss
	// Muss sein, weil Ableitungen h�here Ordnung haben (Bakterienmasse steckt
	// auch in Kapazit�tsgleichung)
	sumX = 0.;  // added CB
	if (maxkap > 1.E-30)
	{
		for (i = 0; i < n; i++)
			if (m_krd->is_a_bacterium[i])
			{
				BacteriaMass = c[i + 1];
				sumX += BacteriaMass;
			}
	}

	/* Anzahl der mikrobiellen Reaktionen aus Datenstruktur auslesen */
	nreactions = m_krd->NumberReactions;

	/* loop over reactions dX/dt= nymax * X * monodterms * inhibitionterms */
	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;

		if (m_kr->typeflag_monod)
		{
			BacteriaNumber = m_kr->bacteria_number + 1;
			BacteriaMass = c[BacteriaNumber];

			if (BacteriaMass > 1.E-40)
			{
				/* Ableitungen werden aus dX/dt = BacGrowth berechnet */
				// sumX is different for case with (>0) or without (==0) maxkap

				m_kr->currentnode = node;  // CB This is eclusively for Brand
				                           // model to allow porosity in
				                           // Inhibition constant calculation

				BacGrowth = m_kr->BacteriaGrowth(r, c, sumX, -1);
				for (i = 0; i < n; i++)
					d2X_dtdS[i + 1] = 0.;

				// Berechnung der Bakterien Ableitungen f�r den Fall dass eine
				// maximale Kapazit�t ber�cksichtigt werden muss
				// Muss sein, weil Ableitungen h�here Ordnung haben
				// (Bakterienmasse steckt auch in Kapazit�tsgleichung)
				if (maxkap > 1.E-30)
				{
					maxVelocity = m_kr->rateconstant;
					/* Wachstumsterm, ber�cksichtige Kapazit�tsterm */
					if (maxVelocity > 1.E-30)
					{
						// Erst Berechnen der Ableitungen des Bakterienwachstums
						// nach allen anderen Substanzen
						//   Ableitung nach der Bakterienmasse (mit
						//   Ber�cksichtigung Kapazit�tsterm)
						//   d2Xi / dt*dXi = BacGrowth *(sumx+maxkap-Xi) /
						//   (Xi*(sumx+maxkap))
						d2X_dtdS[BacteriaNumber] =
						    BacGrowth * (sumX + maxkap - BacteriaMass) /
						    (BacteriaMass * (sumX + maxkap));
						for (i = 0; i < n; i++)
							if (m_krd->is_a_bacterium[i] &&
							    (i + 1 != BacteriaNumber))
								// Ableitung nach den anderen Bakterienmassen im
								// Kapazit�tsterm
								//   d2Xi / dt*dXj = BacGrowth / -(sumx+maxkap)
								d2X_dtdS[i + 1] = BacGrowth / -(sumX + maxkap);
					} /* Sterbeterm, grunds�tzlich keine Ber�cksichtigung des
					     Kapazit�tsterms */
					else
						/* d2Xi / dt*dXi = BacGrowth / Xi */
						d2X_dtdS[BacteriaNumber] = BacGrowth / BacteriaMass;
					/* d2Xi / dt*dXj = 0 for decay */
				}
				// Berechnung der Bakterien Ableitungen f�r den Fall dass KEINE
				// maximale Kapazit�t ber�cksichtigt werden muss
				else  // maxkap = 0

					/* d2Xi / dt*dXi = BacGrowth / Xi */
					d2X_dtdS[BacteriaNumber] = BacGrowth / BacteriaMass;

				/* Schleife f�r Ableitungen nach Substanzen in Monodtermen;
				 * unabh�ngig von maxkap */
				NumberMonod = m_kr->number_monod;
				for (i = 0; i < NumberMonod; i++)
				{
					// d2X / dt*dS_j =      S_j = monod-species
					//   S_j may be ZERO or below !
					MonodSpecies = m_kr->monod[i]->speciesnumber + 1;
					MonodConcentration = m_kr->monod[i]->concentration;
					MonodOrder = m_kr->monod[i]->order;
					CMonodSpecies = c[MonodSpecies];
					if (CMonodSpecies > 1.E-20)
					{
						// S_j > 0, normal Monod Term used
						//   divide BacGrowth through Monod-Term of spec. j
						//   and multiplicate with partial derivate of
						//   Monod-Term

						// In case of isotope fractionation of substrate Si,Sj
						// (i,j=l,h ; light,heavy)
						// - the partial derivative d2X/dtdSi is different:
						//      d2X_dtdS[Si] = BacGrowth *
						//      (MonodConcentration+CisotopePartner)
						//           / CMonodSpecies /
						//           (MonodConcentration+CMonodSpecies+CisotopePartner)
						// - an additional partial derivative d2X/dtdSj with
						// respect to the isotope partner Sj appears
						//   and must be accounted for and added in the vector
						//   d2X_dtdS[...] for the current reaction:
						//      d2X_dtdS[Sj] = BacGrowth * (-1) /
						//      (MonodConcentration+CMonodSpecies+CisotopePartner)...

						// if((m_kr->typeflag_iso_fract==1) &&
						// (m_kr->monod[i]->isotopecouplenumber>=0))
						if ((m_kr->monod[i]->species == m_kr->Isotope_heavy) ||
						    (m_kr->monod[i]->species == m_kr->Isotope_light))
						{
							IsotopeSpecies =
							    m_kr->monod[i]->isotopecouplenumber + 1;
							Ciso = c[IsotopeSpecies];
							// this is the term for the Monod-Species
							d2X_dtdS[MonodSpecies] =
							    BacGrowth * MonodOrder *
							    (MonodConcentration + Ciso) / CMonodSpecies /
							    (MonodConcentration + CMonodSpecies +
							     Ciso);  // no isofrac
							// now get the partial derivative d2x/dtdSj with
							// respect to isotope partner
							d2X_dtdS[IsotopeSpecies] =
							    BacGrowth * MonodOrder * (-1) /
							    (MonodConcentration + CMonodSpecies + Ciso);
						}  // no isofrac standard case: d2X_dtdS[MonodSpecies] =
						   // BacGrowth * MonodConcentration / CMonodSpecies /
						   // (MonodConcentration+CMonodSpecies); // no isofrac
						else
						{
							d2X_dtdS[MonodSpecies] =
							    BacGrowth * MonodOrder * MonodConcentration /
							    CMonodSpecies
							    // no isofrac
							    /
							    (MonodConcentration + CMonodSpecies);
							// If a threshhold term exists for the Monod
							// species, the partial derivative is modified for
							// this species
							//   - the ODE with Monod and threshold term for
							//   Monod species C is:
							//     dX/dt = my*R * [C/(C*K)]^n * [C/(C*T)]^m
							//      - with n and K the Order and
							//      Monod-concentration of the Monod term
							//      - with m and T the Order and
							//      Threshhold-concentration of the threshhold
							//      term
							//   - The partial derivative with respect to C
							//   (taking into account Division by the Monod and
							//   Threshhold term) is
							//     d2X/dtdC = my*R * [n*K/C/(C+K) - p*T/C/(C+T)]
							//   - The latter term hence must be substracted
							//   from the previously calculated first term
							if (m_kr->monod[i]->threshhold == true)
							{
								ThreshConc = m_kr->monod[i]->threshConc;
								ThreshOrder = m_kr->monod[i]->threshOrder;
								d2X_dtdS[MonodSpecies] +=
								    BacGrowth * ThreshOrder * ThreshConc /
								    CMonodSpecies /
								    (ThreshConc + CMonodSpecies);  // no isofrac
							}
						}
					}  // Todo CB isofrac special case necessary?
					else if (CMonodSpecies < -1.E-20)
					{
						/* S_j << 0, linear Monod Term used */
						// d2X_dtdS[MonodSpecies] = BacGrowth / CMonodSpecies ;
						d2X_dtdS[MonodSpecies] =
						    BacGrowth * MonodOrder / CMonodSpecies /
						    1000;  // Changed monod term with smaller slope CB
						if (m_kr->monod[i]->threshhold == true)
						{
							ThreshConc = m_kr->monod[i]->threshConc;
							ThreshOrder = m_kr->monod[i]->threshOrder;
							d2X_dtdS[MonodSpecies] +=
							    BacGrowth * ThreshOrder * ThreshConc /
							    CMonodSpecies /
							    (ThreshConc + CMonodSpecies);  // no isofrac
						}
					}  // Todo CB isofrac special case necessary? Threshhold
					   // terms??
					else
					{
						// S_j near 0 numerically instable
						//   recompute BacGrowth without S_j
						//   (hope, that will only sometimes occur)

						m_kr->currentnode =
						    node;  // CB 19/10/09 This is eclusively for Brand
						           // model to allow porosity in Inhibition
						           // constant calculation

						d2X_dtdS[MonodSpecies] =
						    m_kr->BacteriaGrowth(r, c, sumX, MonodSpecies) /
						    MonodConcentration;
						// if(m_kr->monod[i]->threshhold==true){
						//  ThreshConc = m_kr->monod[i]->threshConc;
						//  ThreshOrder = m_kr->monod[i]->threshOrder;
						//  d2X_dtdS[MonodSpecies] += BacGrowth * ThreshOrder *
						//  ThreshConc / CMonodSpecies /
						//  (ThreshConc+CMonodSpecies); // no isofrac
						//}
					}
				}  // for NumberMonod

				/* Schleife f�r Ableitungen nach Substanzen in
				 * Inhibitionstermen, unabh�ngig von maxkap */
				NumberInhibition = m_kr->number_inhibit;
				for (i = 0; i < NumberInhibition; i++)
				{
					// d2X / dt*dS_j =      S_j = inhibition-species
					//   S_j may be Zero without any problem
					InhibitionSpecies = m_kr->inhibit[i]->speciesnumber + 1;
					InhibitionConcentration = m_kr->inhibit[i]->concentration;

					// ATTENTION!!!
					// CB 16/10/09 Hardcode fix of Fe3 inhibition concentration
					// for Brand model,
					// this parameter depends on porosity as in Min3P Fe3 is in
					// solid phase
					// and Inhibition concentrations are expressed in terms of
					// Volume fraction
					// Vol_Fe3/Vol_BulkAquifer [m?m�]
					// while in Geosys, Fe3 is considered an immobile aqueous
					// species
					// and concentrations were converted to [mol/L_water] by
					// C_gs = C_min3p * rho_Fe3 / molweight / porosity
					// for inhibition concentrations, division by porosity is
					// still required
					if (m_kr->inhibit[i]->species.compare("Fe3") == 0)
						InhibitionConcentration *=
						    1 / m_kr->GetPhaseVolumeAtNode(node, 1, 0);
					// CB  further changes in derivs (1), jacbn (2), class
					// CKinReact{}

					CInhibitionSpecies = c[InhibitionSpecies];
					if (CInhibitionSpecies > 0.)
						// S_j > 0, normal Inhibition Term used
						//   divide BacGrowth through Inhibition-Term of spec. j
						//   and multiplicate with partial derivate of Inhi-Term
						d2X_dtdS[InhibitionSpecies] =
						    -BacGrowth /
						    (InhibitionConcentration + CInhibitionSpecies);
					else
						/* S_j <= 0, linear Inhibition Term used */
						// d2X_dtdS[InhibitionSpecies] = - BacGrowth /
						// (InhibitionConcentration-CInhibitionSpecies);// CB
						// changed as in next line
						d2X_dtdS[InhibitionSpecies] =
						    -BacGrowth /
						    (InhibitionConcentration);  // CB changed due to
					                                    // stimulance of growth
					                                    // for neg conc.
				}

				/* transfer partial derivatives to dfdc-Array of equation solver
				 */
				if (m_kr->grow)
					for (i = 0; i < n; i++)
						/* transfer der berechneten Ableitungen f�r die
						 * Bakteriengruppe */
						dfdc[BacteriaNumber][i + 1] += d2X_dtdS[i + 1];

				/* Berechnung der Ableitungen f�r die vom Bakteriellen Wachstum
				 * abh�ngigen Substanzen, unabh�ngig von maxkap */
				/* d2S_j / dt*dS_k = yield(j) * d2X/dt*dS_k */
				porosity1 = m_kr->GetReferenceVolume(BacteriaNumber - 1, node);
				for (i = 0; i < n; i++)
				{
					Yield = m_kr->ProductionStoch[i];
					if (fabs(Yield) > 1.E-30)
					{
						porosity2 = m_kr->GetReferenceVolume(i, node);
						for (j = 0; j < n; j++)
							if (fabs(d2X_dtdS[j + 1]) > 1.E-30)
								dfdc[i + 1][j + 1] += d2X_dtdS[j + 1] * Yield *
								                      porosity1 / porosity2;
					}
				}
			}  // Ende if BacteriaMass > 1e-40*/
		}      // Ende if type monod
	}          // Ende Schleife �ber nreactions*/

	/**********************************************************************************************/
	/* Berechnung der Ableitungen der Austauschprozesse */
	// calculate already occupied surfaces first
	if (m_krd->NumberLangmuir > 0)
	{
		// Initialise Surfaces for langmuir isotherms
		for (i = 0; i < m_krd->maxSurfaces; i++)
			occupiedSurface.push_back(0.0);

		for (r = 0; r < nreactions; r++)
		{
			m_kr = KinReact_vector[r];
			// CB new reaction switch for individual reactions
			if (m_kr->switched_off_node.size() > 0)
				if (m_kr->switched_off_node[node] == true) continue;
			if ((m_kr->getType().compare("exchange") == 0) &&
			    (m_kr->typeflag_exchange_langmuir))
			{
				Sp1 = m_kr->ex_species[0] + 1;
				surfaceID = m_kr->exSurfaceID;
				occupiedSurface[surfaceID] += c[Sp1];
			}
		}
	}  // if NumberLangmuir > 0

	/* Berechnung der Ableitungen der Austauschprozesse */
	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;

		if (m_kr->typeflag_exchange)
		{
			/* linearer Austausch mit kd */
			if (m_kr->typeflag_exchange_linear)
			{
				// Matrix
				Sp1 = m_kr->ex_species[0] + 1;
				porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
				density1 = m_kr->GetDensity(Sp1 - 1, node);
				// geloest
				Sp2 = m_kr->ex_species[1] + 1;
				porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);

				exch = m_kr->ex_param[0];
				kd = m_kr->ex_param[1];

				if (fabs(kd) < MKleinsteZahl)
				{
					// no kd, exchange between two species in solution
					dfdc[Sp1][Sp1] += -exch / porosity1;
					dfdc[Sp1][Sp2] += exch / porosity1;
					dfdc[Sp2][Sp1] += exch / porosity2;
					dfdc[Sp2][Sp2] += -exch / porosity2;
				}
				else
				{
					// with kd, exchange between matrix (mol/kg) and solution
					// (mol/l)
					foc = m_krd->node_foc[node];
					if (foc > MKleinsteZahl) kd = kd * foc;
					// else
					//  kd = 0;
					dfdc[Sp1][Sp1] += -exch;
					dfdc[Sp1][Sp2] += exch * kd;
					dfdc[Sp2][Sp1] += exch * porosity1 / porosity2 * density1;
					dfdc[Sp2][Sp2] +=
					    -exch * kd * porosity1 / porosity2 * density1;
				}
			}  // linear
			/* Langmuir Kinetik */
			if (m_kr->typeflag_exchange_langmuir)
			{
				// calculated already occupied surfaces above
				Sp1 = m_kr->ex_species[0] + 1;
				porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
				density1 = m_kr->GetDensity(Sp1 - 1, node);
				Sp2 = m_kr->ex_species[1] + 1;
				porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
				kadsorb = m_kr->ex_param[0];
				kdesorb = m_kr->ex_param[1];
				// SB_langmuir		surfaceID	= m_kr->exSurfaceID;
				// //Exchange.Langmuir[i].SurfaceID;
				totalSurface = m_krd->exSurface[m_kr->exSurfaceID];
				//      occupiedSurface is calculated just above
				adsorb = kadsorb * (totalSurface - occupiedSurface[surfaceID]);
				dfdc[Sp1][Sp1] += -kdesorb;
				dfdc[Sp2][Sp1] += kdesorb * porosity1 / porosity2 * density1;
				dfdc[Sp1][Sp2] += adsorb;
				dfdc[Sp2][Sp2] += -adsorb * porosity1 / porosity2 * density1;
				// additional derivatives due to occupied surface
				for (j = 0; j < nreactions; j++)
				{
					m_kr1 = KinReact_vector[j];
					if (m_kr1->getType().compare("exchange") == 0)
						if (m_kr1->typeflag_exchange_langmuir)
						{
							SpX = m_kr1->ex_species[0] + 1;
							surfaceID2 = m_kr1->exSurfaceID;
							if (surfaceID == surfaceID2)
							{
								dfdc[Sp1][SpX] += -kadsorb * c[Sp2];
								dfdc[Sp2][SpX] += kadsorb * c[Sp2] * porosity1 /
								                  porosity2 * density1;
							}
						}
				}
			}  // end if langmuir

			/* Freundlich Kinetik */
			if (m_kr->typeflag_exchange_freundlich)
			{
				Sp1 = m_kr->ex_species[0] + 1;
				porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
				density1 = m_kr->GetDensity(Sp1 - 1, node);
				Sp2 = m_kr->ex_species[1] + 1;
				porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
				exponent = m_kr->ex_param[2];
				parameter = m_kr->ex_param[1];
				exch = m_kr->ex_param[0];

				if (c[Sp2] > residual)
					// no linearisation required
					adsorb = exch * parameter * exponent *
					         pow(c[Sp2], (exponent - 1.0));
				else
					// linearisation required due to instability of c^x if
					// c<residual
					adsorb = exch * parameter * exponent *
					         pow(residual, (exponent - 1.0));

				dfdc[Sp1][Sp1] += -exch;
				dfdc[Sp2][Sp1] += exch * porosity1 / porosity2 * density1;
				dfdc[Sp1][Sp2] += adsorb;
				dfdc[Sp2][Sp2] += -adsorb * porosity1 / porosity2 * density1;
			}  // end freundlich
		}      // end if exchange
	}          // end loop over reactions

	//#ds
	/**********************************************************************************************/
	/* NAPL-dissolution */
	/**********************************************************************************************/
	for (r = 0; r < nreactions; r++)
	{
		m_kr = KinReact_vector[r];
		// CB new reaction switch for individual reactions
		if (m_kr->switched_off_node.size() > 0)
			if (m_kr->switched_off_node[node] == true) continue;

		if (m_kr->typeflag_napldissolution)
		{
			/* NAPL-L�sung */
			Sp1 = m_kr->ex_species[0] +
			      1;  // Exchange.Linear[i].Species1; should be NAPL
			//		porosity1	= m_kr->GetReferenceVolume(Sp1-1,node);
			blob = m_kr->blob_ID;
			m_kb = KinBlob_vector[blob];  // pointer to blob-properties set in
			                              // the reaction r
			Sp2 = m_kr->ex_species[1] +
			      1;  // Exchange.Linear[i].Species2; should be dissolved
			// CB this includes the saturation
			porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
			//#ds TODO	    saturation2 = ??
			saturation2 = 1.;

			/**************************/
			// CB 040808 Saturation for aqueous phase
			// Is this necessary? Saturation is included in GetReferenceVolume!!
			// I think it's wrong ,so I deactivate 09/2009
			// CRFProcess *m_pcs = NULL;
			// m_pcs = PCSGet("TWO_PHASE_FLOW");
			// if(m_pcs->pcs_type_number==0)
			//  m_pcs = pcs_vector[m_pcs->pcs_number+1]; // this is the
			//  saturation equation
			// int idxs1 = m_pcs->GetNodeValueIndex("SATURATION1"); // Sat of
			// water phase
			// saturation2 = m_pcs->GetNodeValue(node, idxs1);
			/**************************/
			exch = m_kb->Masstransfer_k * m_kb->current_Interfacial_area;
			/* Remark: In function derivs the dissolution velocity can be
			   reduced, if NAPL concentration falls below zero.
			   The partial derivatives should change, too, but this is not
			   considered here.
			   However, it works fine */

			//		dfdc[Sp1][Sp1] = 0     derivatives for NAPL-concentration
			// always zero
			//		dfdc[Sp2][Sp1] = 0

			// no dissolution or NAPL mass present
			if ((m_kr->current_Csat < c[Sp2]) || (c[Sp1] > MKleinsteZahl))
			{
				// d2CNAPL / dt dCmob = k*A
				// d2Cmob  / dt dCmob = -k*A/n/Sw
				dfdc[Sp1][Sp2] += exch;
				dfdc[Sp2][Sp2] += -exch / porosity2 / saturation2;
			}
		}  // NAPL-dissolution
	}      // loop over reactions r
	//#ds

	free_dvector(d2X_dtdS, 1, n);
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction class test output function
   Programing:
   05/2004 SB Implementation - adapted from OK rf_bc_new
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReact::TestWrite(void)
{
	int i, length, flag = 0;

	// Write Keyword
	cout << "888888888888888888888888888888888888888888888888888888888888888888"
	        "8888888888888888888888888888888888888" << endl;
	cout << " Test Output " << endl;
	cout << "888888888888888888888888888888888888888888888888888888888888888888"
	        "8888888888888888888888888888888888888" << endl;
	cout << "#REACTION" << endl;
	// Name of reaction
	cout << "$NAME" << endl
	     << name << endl;
	// Type of reaction
	cout << "$TYPE" << endl
	     << type << endl;
	// bacteria name
	cout << "$BACTERIANAME" << endl
	     << bacteria_name << endl;
	// ReactionEquation
	cout << "$EQUATION" << endl;
	for (i = 0; i < number_reactionpartner; i++)
	{
		if (stochmet[i] < 0.0)  // left side of equation
		{
			if (i == 0)
				cout << " " << fabs(stochmet[i]) << " " << reactionpartner[i];
			else
				cout << " + " << fabs(stochmet[i]) << " " << reactionpartner[i];
		}
		if (stochmet[i] > 0 && (flag > 0))  // remaining right hand side
			cout << " + " << fabs(stochmet[i]) << " " << reactionpartner[i];
		if (stochmet[i] > 0 &&
		    (flag == 0))  // " = " Sign and first term on right hand side
		{
			cout << " = " << fabs(stochmet[i]) << " " << reactionpartner[i];
			flag = 1;
		}
	}
	cout << endl;
	// Rateconstant and order
	cout << "$RATEKONSTANT" << endl
	     << rateconstant << "   " << rateorder << endl;
	cout << "$GROWTH" << endl
	     << grow << endl;
	// Monod terms
	cout << "$MONODTERMS" << endl
	     << number_monod << endl;
	for (i = 0; i < number_monod; i++)
		cout << monod[i]->species << "  " << monod[i]->concentration << "  "
		     << monod[i]->order << endl;
	// Inhibition terms
	cout << "$INHIBITIONTERMS" << endl
	     << number_inhibit << endl;
	for (i = 0; i < number_inhibit; i++)
		cout << inhibit[i]->species << "  " << inhibit[i]->concentration << "  "
		     << inhibit[i]->order << endl;
	// Production Terms
	cout << "$PRODUCTIONTERMS" << endl
	     << number_production << endl;
	for (i = 0; i < number_production; i++)
		cout << production[i]->species << "  " << production[i]->concentration
		     << "  " << production[i]->order << endl;
	// ProductionStochhelp Terms
	cout << "$PRODUCTIONSTOCH" << endl
	     << (int)ProdStochhelp.size() << endl;
	for (i = 0; i < (int)ProdStochhelp.size(); i++)
		cout << ProdStochhelp[i]->species << "  "
		     << ProdStochhelp[i]->concentration << endl;
	// exchange
	cout << "$EXCHANGE_PARAMETERS" << endl
	     << (int)ex_param.size() << endl;
	for (i = 0; i < (int)ex_param.size(); i++)
		cout << ex_param[i] << "  ";
	cout << endl;

	cout << endl;

	cout << "number_reactionpartner " << (int)number_reactionpartner << endl;
	cout << "bacteria_number " << (int)bacteria_number << endl;
	cout << "grow " << grow << endl;
	length = (int)ProductionStoch.size();
	cout << "length ProductionStoch: " << length << endl;
	for (i = 0; i < length; i++)
		cout << (int)ProductionStoch[i] << " ";
	cout << endl;

	length = (int)ex_species.size();
	cout << "length exSpecies: " << length << endl;
	for (i = 0; i < length; i++)
		cout << ex_species_names[i] << " " << ex_species[i] << endl;
	cout << endl;
	cout << " sorption type : " << exType << endl;

	// Test output
}

/**************************************************************************
   Reaction-Method:
   Task: returns true if NAPL dissolution is modeled
   required for calculation of NAPL Densities and Saturations in
   LOPPreTimeLoop_PCS and LOPTimeLoop_PCS
   Programing:
   08/2008 CB Implementation
**************************************************************************/
bool KNaplDissCheck(void)
{
	int j;
	CKinReact* m_kr = NULL;
	// OK411 CKinReactData *m_krd = NULL;
	int nreact;
	bool NAPLdiss = false;

	// check for NAPL dissolution
	nreact = (int)KinReact_vector.size();
	for (j = 0; j < nreact; j++)
	{
		m_kr = KinReact_vector[j];
		if (m_kr->getType().compare("NAPLdissolution") == 0)
		{
			NAPLdiss = true;
			break;
		}
	}
	return NAPLdiss;
}

/**************************************************************************
   Reaction-Method:
   Task: returns the Porevelocity of the mobile (water) phase in case of a NAPL
   dissolution model and TwoPhaseFlow; v at the node as a inverse distance
   weighted mean of the connecting elements velocities
   Programing:
   08/2008 CB Implementation
   10/2010 TF changed access to process type
   09/2011 TF changed access to coordinates of mesh node,
    - substituted access to mesh_element from pointer to direct access into the
vector
    - made the mesh node a const pointer
    - made the pointer to the mesh const, made the mesh itself const
    - substituted pow(x,2) by x*x
    - reduced scope of some loop variables
**************************************************************************/
double CKinReact::GetNodePoreVelocity(long node_number)
{
	CRFProcess* pcs(PCSGetFlow());
	CFEMesh const* const msh(fem_msh_vector[0]);  // SB: ToDo hart gesetzt

	long group;
	long elem;
	long idxVx, idxVy, idxVz, idxs1;
	double vel_nod[3], vel_ele[3];
	double distance, weight, sum_w(0);
	double PoreVel(0), poro(0), satu = 1.0;  // default
	double theta = pcs->m_num->ls_theta;

	// Get node saturation of mobile (water) phase
	if (pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
	{
		if (pcs->pcs_type_number == 0)
			// this is the saturation equation
			pcs = pcs_vector[pcs->pcs_number + 1];
		// Sat of water phase
		idxs1 = pcs->GetNodeValueIndex("SATURATION1");
		satu = pcs->GetNodeValue(node_number, idxs1);
	}
	else if (pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
	{
		// Sat of water phase
		idxs1 = pcs->GetNodeValueIndex("SATURATION1");
		satu = pcs->GetNodeValue(node_number, idxs1);
	}

	// initialize data structures
	for (size_t i = 0; i < 3; i++)
		vel_nod[i] = vel_ele[i] = 0;

	// Get node coordinates
	MeshLib::CNode const* const node(msh->nod_vector[node_number]);
	double const* const coord(node->getData());  // Coordinates(coord);
	// get the indices of velocity of flow process
	pcs = PCSGet("TWO_PHASE_FLOW");
	idxVx = pcs->GetElementValueIndex("VELOCITY1_X");
	idxVy = pcs->GetElementValueIndex("VELOCITY1_Y");
	idxVz = pcs->GetElementValueIndex("VELOCITY1_Z");

	for (size_t el = 0; el < node->getConnectedElementIDs().size(); el++)
	{
		distance = weight = 0;  // initialize for each connected element
		elem = node->getConnectedElementIDs()[el];

		// get the porosity of current element elem
		group = 0;  // SB todo group = m_ele->GetPatchIndex(); Todo CB
		poro =
		    (mmp_vector[group])->Porosity(elem, theta);  // CB Now provides also
		// heterogeneous porosity, model 11
		// get the velocity components of element elem and divide by porosity
		vel_ele[0] = pcs->GetElementValue(elem, idxVx) / poro;
		vel_ele[1] = pcs->GetElementValue(elem, idxVy) / poro;
		vel_ele[2] = pcs->GetElementValue(elem, idxVz) / poro;
		// calculate distance node <-> element center of gravity
		double const* grav_c(msh->ele_vector[elem]->GetGravityCenter());
		for (size_t i = 0; i < 3; i++)
			distance +=
			    (coord[i] - grav_c[i]) *
			    (coord[i] - grav_c[i]);  // pow((coord[i] - grav_c[i]), 2);
		// linear inverse distance weight = 1/(distance)
		distance =
		    sqrt(distance);  // for quadratic interpolation uncomment this line
		weight = (1 / distance);
		sum_w += weight;
		for (size_t i = 0; i < 3; i++)
			// the 3 velocity components
			vel_nod[i] += vel_ele[i] * weight;
	}
	// normalize weighted sum by sum_of_weights sum_w
	for (size_t i = 0; i < 3; i++)
		vel_nod[i] *= 1 / sum_w;
	// absolute value of velocity vector and divide by saturation
	for (size_t i = 0; i < 3; i++)
		PoreVel += vel_nod[i] * vel_nod[i];  // pow(vel_nod[i], 2);
	PoreVel = sqrt(PoreVel) / satu;

	return PoreVel;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactionDeactivation()                            */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Deactivates kinetic reaction calculation at individual nodes based on  */
/* evaluation of reaction rates of the previous time step in local        */
/* neighborhoods around nodes, or based on comparison to previous         */
/*   concentrations at the nod                                            */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 08/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void CKinReactData::ReactionDeactivation(long nonodes)
{
	long node;
	int react_t =
	    10;  // reactions are calculated at all nodes every react_t timesteps
	size_t sp;
	const size_t Number_of_Components(cp_vec.size());

	// reactivate all nodes every n time steps
	if (((aktueller_zeitschritt) % react_t) == 0)  // when evaluated before loop
		                                           // over nodes, i.e. prepare
		                                           // for this time step

		// if(((aktueller_zeitschritt+1) % react_t ) == 0) {  // when evaluated
		// after  loop over nodes, i.e. prepare for next time step
		for (node = 0; node < nonodes; node++)
		{
			ReactDeact[node] = false;
			if (ReactDeactMode != 3)  // for all timesteps, anyway prepare the
				                      // concentrationmatrix for the next time
				                      // step, i.e. save the current
				                      // concentrations after Transport
				for (sp = 0; sp < Number_of_Components; sp++)
					concentrationmatrix[node][sp] =
					    cp_vec[sp]->getProcess()->GetNodeValue(node,
					                                           sp_varind[sp]);
		}
	// CB Now check if node can be deactivated for the next time step
	else
	{
		// only for mode 1: first calculate for each node the sum of rates and
		// store in vector,
		// which then is accessed in next loop over nodes, to evaluate the sum
		// of sum of rates
		// over the neighbours
		if (ReactDeactMode == 1)
			for (node = 0; node < nonodes; node++)
			{
				React_dCdT[node] = 0;
				for (sp = 0; sp < Number_of_Components; sp++)
				{
					double maxi = 1;
					// maxi = max(conc_old[sp+1], conc_new[sp+1]);
					if (maxi < MKleinsteZahl)
						React_dCdT[node] += 0.0;
					else
					{
						// This is C of last time step after reactions, i.e. the
						// old time level for this time step
						double Concentration =
						    cp_vec[sp]->getProcess()->GetNodeValue(
						        node, (sp_varind[sp] - 1));
						// This is C of previous time step after transport only,
						// which was stored in matrix
						double Concentration_old =
						    concentrationmatrix[node][sp];
						React_dCdT[node] +=
						    fabs((Concentration_old - Concentration) / maxi) /
						    dt;  // normalized by current local concentration
						// and now prepare concentrationmatrix for next time
						// step, i.e. save current concentrations after
						// Transport
						concentrationmatrix[node][sp] =
						    cp_vec[sp]->getProcess()->GetNodeValue(
						        node, sp_varind[sp]);
					}
				}
				// cout << nod << " " << React_dCdT[nod] << endl;
			}

		// this is the check, if a node may be deactivated for this time step
		for (node = 0; node < nonodes; node++)
		{
			double sumReact_dCdT = 0;

			switch (ReactDeactMode)
			{
				case 1:  // loop over no of connected nodes and their respective
					     // neighbours
					for (size_t nn = 0; nn < ReactNeighborhood[node].size();
					     nn++)
					{
						int node_idx = ReactNeighborhood[node][nn];
						sumReact_dCdT += React_dCdT[node_idx];
					}
					break;
				case 2:  // compare with C after transport of last timestep,
					     // loop over all components
					for (sp = 0; sp < Number_of_Components; sp++)
					{
						// This is C of current time step after transport
						double Concentration =
						    cp_vec[sp]->getProcess()->GetNodeValue(
						        node, sp_varind[sp]);
						// This is C of previous time step after transport only
						double Concentration_old =
						    concentrationmatrix[node][sp];
						sumReact_dCdT +=
						    fabs(Concentration - Concentration_old);
						// and now prepare the concentrationmatrix for the next
						// time step, i.e. save the current concentrations after
						// Transport
						concentrationmatrix[node][sp] = Concentration;
					}
					break;
				case 3:  // compare with C after reaction of last timestep, loop
					     // over all components
					for (sp = 0; sp < Number_of_Components; sp++)
					{
						// This is C of current time step after transport
						double Concentration =
						    cp_vec[sp]->getProcess()->GetNodeValue(
						        node, sp_varind[sp]);
						// This is C of previous time step after transport &
						// reaction
						double Concentration_old =
						    cp_vec[sp]->getProcess()->GetNodeValue(
						        node, (sp_varind[sp] - 1));
						sumReact_dCdT +=
						    fabs(Concentration - Concentration_old);
					}
					break;
				default:
					break;
			}

			// check if deactivation criterion is met
			if (sumReact_dCdT < ReactDeactEpsilon)
				ReactDeact[node] = true;  // negligible change, deactivate the
			                              // node for the next time step
			else
				ReactDeact[node] = false;  // sufficient change, reactivate the
			                               // node for the next time step
		}
	}  // else

	// Resets the reaction rates vector
	if (ReactDeactMode == 1)
		for (node = 0; node < nonodes; node++)
			React_dCdT[node] = 0;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactDeactReset_dCdT()                            */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Sets the C after reaction of last time step as new C after reaction    */
/* for deactivated nodes */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 11/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void CKinReactData::ReactDeactSetOldReactionTerms(long nonodes)
{
	long node;
	int sp;
	int Number_of_Components = (int)cp_vec.size();
	double Concentration;

	for (node = 0; node < nonodes; node++)
	{
		if (ReactDeact[node] == true)
			for (sp = 0; sp < Number_of_Components; sp++)
			{
				// Get the C after reactions of last time step (old time level,
				// index = 0)
				Concentration = cp_vec[sp]->getProcess()->GetNodeValue(
				    node, sp_varind[sp] - 1);
				// Set this C as the new concentration after reaction
				cp_vec[sp]->getProcess()->SetNodeValue(node, sp_varind[sp],
				                                       Concentration);
			}
	}
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactDeactPlotFlagsToTec()                        */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Prints flags for reaction deactivation in tecplot format               */
/* for all nodes                                                          */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 11/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void CKinReactData::ReactDeactPlotFlagsToTec()
{
	CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
	std::string eleType;

	if (m_msh->getNumberOfLines() > 0) eleType = "QUADRILATERAL";
	if (m_msh->getNumberOfQuads() > 0) eleType = "QUADRILATERAL";
	if (m_msh->getNumberOfHexs() > 0) eleType = "BRICK";
	if (m_msh->getNumberOfTris() > 0) eleType = "QUADRILATERAL";
	if (m_msh->getNumberOfTets() > 0) eleType = "TETRAHEDRON";
	if (m_msh->getNumberOfPrisms() > 0) eleType = "BRICK";

	if (NumberReactions > 0)
	{
		std::string filename(FileName + "_Deactivated_nodes.tec");
		std::ofstream aus;
		if (aktueller_zeitschritt == 1)
			aus.open(filename.c_str());
		else
			aus.open(filename.c_str(), ios::app);

		const size_t nnodes(m_msh->nod_vector.size());
		const size_t nele(m_msh->ele_vector.size());

		aus << "VARIABLES = "
		    << "\"x\""
		    << " "
		    << "\"y\""
		    << " "
		    << "\"z\""
		    << "\"active\"" << endl;
		aus << "ZONE T="
		    << "\"aktueller_zeitschritt=" << aktueller_zeitschritt << "\"";
		aus << ", N=" << nnodes << ", E=" << nele
		    << " F=FEPOINT, ET=" << eleType << endl;

		for (size_t i = 0; i < nnodes; i++)
		{
			double const* coord =
			    (m_msh->nod_vector[i])->getData();  // Coordinates(coord);
			aus << coord[0] << " " << coord[1] << " " << coord[2] << " ";
			if (is_a_CCBC[i] == true)
				aus << 0 << endl;
			else if (ReactDeact[i] == true)
				aus << 0 << endl;
			else
				aus << 1 << endl;
		}
		for (size_t i = 0; i < nele; i++)
			m_msh->ele_vector[i]->WriteIndex_TEC(aus);

		aus.close();
	}
}

void CKinReactData::Aromaticum(long nonodes)
{
	long node;
	double conc, lambda = 0.0;  // OK411
	// int pcsindex = 0;
	CRFProcess* m_pcs = NULL;
	int varindex = 0;
	int nospec = (int)sp_varind.size();

	for (int sp = 0; sp < nospec; sp++)
		if (cp_vec[sp]->getProcess()->nod_val_name_vector[0].compare(
		        "Aromaticum") == 0)
		{
			// pcsindex = sp_pcsind[sp];
			m_pcs = cp_vec[sp]->getProcess();
			varindex = sp_varind[sp];
			break;
		}
	for (int sp = 0; sp < nospec; sp++)
		if (cp_vec[sp]->compname.compare("Aromaticum") == 0)
		{
			lambda = cp_vec[sp]->decay_model_values[0];
			break;
		}

	for (node = 0; node < nonodes; node++)
	{
		conc = m_pcs->GetNodeValue(node, varindex);
		conc = conc * exp(-lambda * dt);
		m_pcs->SetNodeValue(node, varindex, conc);
	}
}
