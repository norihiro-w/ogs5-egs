/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
 * \file Test_Base.cpp
 * 29/4/2010 LB Initial implementation
 *
 * Tests for the Base directory
 */

// ** INCLUDES **
#include "gtest.h"

#include "swap.h"
#include "StringTools.h"

TEST(Base, SwapInt)
{
	int arg0 = 5;
	int arg1 = 10;
	BASELIB::swap(arg0, arg1);
	ASSERT_EQ(arg0, 10);
	ASSERT_EQ(arg1, 5);
}

TEST(Base, SwapDouble)
{
	double arg0 = 5.0;
	double arg1 = 10.0;
	BASELIB::swap(arg0, arg1);
	ASSERT_EQ(arg0, 10.0);
	ASSERT_EQ(arg1, 5.0);
}

TEST(Base, SwapString)
{
	std::string arg0 = "5";
	std::string arg1 = "10";
	BASELIB::swap(arg0, arg1);
	ASSERT_EQ(arg0, std::string("10"));
	ASSERT_EQ(arg1, std::string("5"));
}


std::vector<std::pair<std::string, std::string> > parse(std::string const& misc_setting)
{
	std::vector<std::pair<std::string, std::string> > vec_para;
	std::list<std::string> lst = splitString(misc_setting, ' ');
	for (std::list<std::string>::iterator itr = lst.begin();
	     itr != lst.end();
	     ++itr)
	{
		// key-value or only key
		std::string& str1 = *itr;
		std::string val = "";
		if (str1.find('-') == std::string::npos) continue;
		++itr;
		if (itr != lst.end())
		{
			std::string& str2 = *itr;
			if (str2[0] != '-')
				val = str2;
			else
				--itr;
		}
		vec_para.push_back(std::make_pair(str1, val));
		std::cout << "\t " << str1 << " = " << val.c_str() << "\n";
		if (itr == lst.end())
			break;
	}
	return vec_para;
}

TEST(PETSC, ParseArgs)
{
	std::vector<std::pair<std::string, std::string> >  res1(parse("-test1 -boomer_amg_maxitr 10 -test2 -boomer_amg_tol 1e-12 -test3"));
	ASSERT_EQ(5u, res1.size());
	ASSERT_EQ("-test1", res1[0].first);
	ASSERT_EQ("", res1[0].second);
	ASSERT_EQ("-boomer_amg_maxitr", res1[1].first);
	ASSERT_EQ("10", res1[1].second);
	ASSERT_EQ("-test2", res1[2].first);
	ASSERT_EQ("", res1[2].second);
	ASSERT_EQ("-boomer_amg_tol", res1[3].first);
	ASSERT_EQ("1e-12", res1[3].second);
	ASSERT_EQ("-test3", res1[4].first);
	ASSERT_EQ("", res1[4].second);
}

