/* ===========================================================================
 * Copyright (c) 2016-2017 Giacomo Resta
 *
 * This file is part of TightBinding++.
 *
 * TightBinding++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TightBinding++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ===========================================================================
 */

#ifndef __TBPP_PYTHON_CONTEXT_H__
#define __TBPP_PYTHON_CONTEXT_H__

#include <Python.h>
#include <tbpp/context.h>

extern PyTypeObject ContextPyType;
#define ContextPy_Check(_v) PyObject_TypeCheck((_v), &ContextPyType)

typedef struct {
    PyObject_HEAD
    tbpp::ContextPtr ptr;
} ContextPy;

PyObject* context_wrap(tbpp::ContextPtr ptr);
tbpp::ContextPtr context_unwrap(PyObject* obj);

#endif /* __TBPP_PYTHON_CONTEXT_H__ */
