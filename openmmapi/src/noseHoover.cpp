/******************************************************************************
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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

#include "openmm/NoseHoover.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

NoseHooverIntegrator::NoseHooverIntegrator(double stepSize,double inputpres,double inputtemp, double inputtaupres, int inputnfree) {
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    temperature=inputtemp;
    pressure=inputpres;    
    nfree=inputnfree;
    taupres=inputtaupres;
}

void NoseHooverIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateNoseHooverKernel::Name(), contextRef);
    kernel.getAs<IntegrateNoseHooverKernel>().initialize(contextRef.getSystem(),*this);
}

void NoseHooverIntegrator::cleanup() {
    kernel = Kernel();
}

vector<string> NoseHooverIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateNoseHooverKernel::Name());
    return names;
}

double NoseHooverIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateNoseHooverKernel>().computeKineticEnergy(*context, *this);
}

void NoseHooverIntegrator::step(int steps) {
        context->calcForcesAndEnergy(true, false);
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    for (int i = 0; i < steps; ++i) {
        kernel.getAs<IntegrateNoseHooverKernel>().execute(*context, *this);
    }
}
double NoseHooverIntegrator::getTemperature() const{
	return temperature;
}
double NoseHooverIntegrator::getPressure() const{
           return pressure;
}int NoseHooverIntegrator::getNfree() const{
	return nfree;
}double NoseHooverIntegrator::gettaupres() const{
	 return taupres;
}
