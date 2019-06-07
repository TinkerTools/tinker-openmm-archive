#ifndef OPENMM_HooverBarostat_H_
#define OPENMM_HooverBarostat_H_

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

#include "Force.h"
#include <string>
#include "internal/windowsExport.h"

namespace OpenMM {


class OPENMM_EXPORT HooverBarostat : public Force {
public:
    /**
     * This is the name of the parameter which stores the current pressure acting on
     * the system (in bar).
     */
    /**
     * Create a HooverBarostat.
     *
     * @param defaultPressure     the default pressure acting on the system (in bar)
     * @param defaultTemperature  the default temperature at which the system is being maintained (in Kelvin)
     * @param frequency           the frequency at which Monte Carlo pressure changes should be attempted (in time steps)
     */
    HooverBarostat(double stepsize, double defaultTemperature, double pressure, int nfree, int frequency, double tau, int innersteps);
    /**
     * Get the frequency (in time steps) at which Barostat changes should be attempted.  If this is set to
     * 0, the barostat is disabled.
     */
    int getFrequency() const {
        return frequency;
    }
    /**
     * Set the frequency (in time steps) at which Barostat changes should be attempted.  If this is set to
     * 0, the barostat is disabled.
     */
    void setFrequency(int freq) {
        frequency = freq;
    }
    /**
     * Get the default temperature at which the system is being maintained, measured in Kelvin.
     */
    double getDefaultTemperature() const {
        return defaultTemperature;
    }
    double getPressure() const{
	return pressure;
    }
    /**
     * Set the default temperature at which the system is being maintained.  This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param temp     the system temperature, measured in Kelvin.
     */
    void setDefaultTemperature(double temp) {
        defaultTemperature = temp;
    }
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }double getStepSize() const{
	return stepsize; 
    }int getnfree() const{
	return nfree;
    }double gettau() const{
	return tau;
    }int getinnersteps() const{
	return innersteps;
    }
protected:
    ForceImpl* createImpl() const;
private:
    double pressure, defaultTemperature, stepsize, tau;
    int frequency,nfree,innersteps;
};

} // namespace OpenMM

#endif /*OPENMM_HooverBarostat_H_*/
