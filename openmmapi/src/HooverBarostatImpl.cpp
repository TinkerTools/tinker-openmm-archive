/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/HooverBarostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/Context.h"
#include "openmm/kernels.h"
#include "openmm/OpenMMException.h"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace OpenMM;
using std::vector;

const float BOLTZMANN = 1.380658e-23f; // (J/K)
const float AVOGADRO = 6.0221367e23f;
const float RGAS = BOLTZMANN*AVOGADRO; // (J/(mol K))
const float BOLTZ = RGAS/1000;         // (kJ/(mol K))

HooverBarostatImpl::HooverBarostatImpl(const HooverBarostat& owner) : owner(owner), step(0) {
}

void HooverBarostatImpl::initialize(ContextImpl& context) {
    if (!context.getSystem().usesPeriodicBoundaryConditions())
        throw OpenMMException("A barostat cannot be used with a non-periodic system");
    kernel = context.getPlatform().createKernel(HooverBarostatKernel::Name(), context);
    kernel.getAs<HooverBarostatKernel>().initialize(context.getSystem(), owner);
    counter=0;
    frequency= owner.getFrequency();
}

void HooverBarostatImpl::updateContextState(ContextImpl& context) {
    if(context.getStepUsesVirial()){
	        kernel.getAs<HooverBarostatKernel>().execute(context,owner);
    }
    counter+=1;
    if(counter%frequency==0 && counter>0){
	context.setStepUsesVirial(true);
    }else{
	context.setStepUsesVirial(false);
    }
}
std::map<std::string, double> HooverBarostatImpl::getDefaultParameters() {
    std::map<std::string, double> parameters;
    return parameters;
}

std::vector<std::string> HooverBarostatImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(HooverBarostatKernel::Name());
    return names;
}

