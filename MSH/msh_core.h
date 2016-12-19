/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef msh_core_INC
#define msh_core_INC

#include <iostream>

namespace MeshLib
{

class CCore
{
public:
	CCore(size_t id);
	virtual ~CCore() {}

	size_t GetIndex() const { return index; }

	bool GetMark() const { return mark; }

	bool GetOrder() const { return quadratic; }

	bool Interior() const
	{
		if (boundary_type == 'I')
			return true;
		else
			return false;
	}

	void SetBoundaryType(char type) { boundary_type = type; }

	char GetBoundaryType() const { return boundary_type; }

	void SetOrder(bool order) { quadratic = order; }

	void SetMark(bool state) { mark = state; }

	void SetIndex(size_t lvalue) { index = lvalue; }

	virtual void Write(std::ostream& os = std::cout) const { os << "\n"; }

protected:
	size_t index;
	char boundary_type;
	bool mark;       // e.g. marked to be refined or active
	bool quadratic;  // higher order
};

}  // namespace MeshLib

#endif
