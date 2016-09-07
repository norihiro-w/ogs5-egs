/*
 * LinearFunctionData.h
 *
 *  Created on: Sep 1, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef LINEARFUNCTIONDATA_H_
#define LINEARFUNCTIONDATA_H_

#include <fstream>

#include "exprtk/exprtk.hpp"

#include "display.h"

#if 1
class ExprtTkFunction
{
	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double> expression_t;
	typedef exprtk::parser<double> parser_t;
	typedef exprtk::parser_error::type error_t;

public:
	explicit ExprtTkFunction(const std::string& str_expression);

	double getValue(double x, double y, double z) const;
	std::string getExpression() const { return _exp_str; }

private:
	std::string _exp_str;
	mutable double _x, _y, _z;
	symbol_table_t symbol_table;
	expression_t expression;
};

class LinearFunctionData
{
public:
	LinearFunctionData(std::ifstream& ins, int num_var = -1);
	~LinearFunctionData();

	double getValue(size_t dom_i, double x, double y, double z) const;
	double getValue(double x, double y, double z) const;
	size_t* getSubDomIndex() const { return (size_t*)&_subdom_index[0]; }

	std::string getExpression(size_t dom_i) const
	{
		for (size_t i = 0; i < _subdom_index.size(); i++)
			if (dom_i == _subdom_index[i]) return _subdom_f[i]->getExpression();
		return "";
	}

private:
	std::vector<size_t> _subdom_index;
	std::vector<ExprtTkFunction*> _subdom_f;
};

#else
/*!
   \class LinearFunctionData
   \brief Define a linear function for IC, BC and ST

   WW 24.08.2011
 */
class LinearFunctionData
{
public:
	/*!
	   \brief Constrcutor of class LinearFunctionData

	   \param ifstream &ins: file
	   \param num_var: number of data
	   \param sub_domain: if sub_domain

	   WW 24.08.2011
	 */
	LinearFunctionData(std::ifstream& ins, int num_var = -1);
	~LinearFunctionData();

	double getValue(size_t dom_i, double x, double y, double z) const;
	double getValue(double x, double y, double z) const;
	size_t* getSubDomIndex() const { return _subdom_index; }

private:
	size_t _ndata;
	size_t* _subdom_index;
	// Coefficents for linear distribution function
	// f = a0+b0*x+c0*y+d0*z
	double* _a0, *_b0, *_c0, *_d0;
};
#endif

#endif /* LINEARFUNCTIONDATA_H_ */
