#ifndef OPENMM_AMOEBA_CT_FORCE_H_
#define OPENMM_AMOEBA_CT_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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

#include "openmm/Force.h"
#include "internal/windowsExportAmoeba.h"
#include <vector>

namespace OpenMM {

/**
 * This class implements a pairwise exponential potential used to model CT forces.
 *
 * To use it, create an AmoebaCTForce object then call addParticle() once for each particle.  Then call
 * addCTprByOldTypes() to use special CTpr parameters.  After particles have been added,
 * you can modify its force field parameters by calling setParticleParameters() and setCTprParameters().  This
 * will have no effect on Contexts that already exist unless you call updateParametersInContext().
 *
 * A unique feature of this class is that the interaction site for a particle does not need to be
 * exactly at the particle's location.  Instead, it can be placed a fraction of the distance from that
 * particle to another one.  This is typically done for hydrogens to place the interaction site slightly
 * closer to the parent atom.  The fraction is known as the "reduction factor", since it reduces the distance
 * from the parent atom to the interaction site.
 */

class OPENMM_EXPORT_AMOEBA AmoebaCTForce : public Force {
public:
    /**
     * This is an enumeration of the different methods that may be used for handling long range nonbonded forces.
     */
    enum NonbondedMethod {
        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.
         */
        CutoffPeriodic = 1,
    };

    /**
     * Create an Amoeba CTForce.
     */
    AmoebaCTForce();

    /**
     * Get the number of particles.
     */
    int getNumParticles() const {
        return parameters.size();
    }

    /**
     * Set the force field parameters for a CT particle.
     *
     * @param particleIndex   the particle index
     * @param CTprType        the CT class/type number
     * @param apre            CT apre
     * @param bexp            CT bexp
     * @param lambda          CT lambda
     */
    void setParticleParameters(int particleIndex, int CTprType, double apre, double bexp, double lambda);

    /**
     * Get the force field parameters for a CT particle.
     *
     * @param particleIndex        the particle index
     * @param[out] CTprType        the CT class/type number
     * @param[out] apre            CT apre
     * @param[out] bexp            CT bexp
     * @param[out] lambda          CT lambda
     */
    void getParticleParameters(int particleIndex, int& CTprType, double& apre, double& bexp, double& lambda) const;


    /**
     * Add the force field parameters for a CT particle.
     *
     * @param CTprType       the CT class/type number
     * @param apre           CT apre
     * @param bexp           CT bexp
     * @param lambda         CT lambda
     * @return index of added particle
     */
    int addParticle(int CTprType, double apre, double bexp, double lambda);

    /**
     * Compute the combined apres and bexps using combining rules.
     */
    void computeCombinedApreBexp();

    /**
     * Get the number of CT classes/types.
     */
    int getNumCTprTypes() const {
        return numCTprTypes;
    }

    /**
     * Set the number of CT classes/types.
     *
     * @param newNum the new number of CT classes/types
     */
    void setNumCTprTypes(int newNum);

    /**
     * Get the new CT class/type.
     *
     * @param oldType the old CT class/type
     * @return the new CT class/type number
     */
    int getNewCTprType(int oldType) const;

    /**
     * Get the old CT class/type number from a new CT class/type number.
     *
     * @param newType the new CT class/type number
     * @return the old CT class/type number
     */
    int getOldCTprType(int newType) const;

    /**
     * Set an old CT class/type number with a new CT class/type.
     *
     * @param newType the new CT class/type number
     * @param oldType the old CT class/type number
     */
    void setOldCTprType(int newType, int oldType);

    /**
     * Resize some internal storing variables.
     *
     * @param newSize the new number of CT classes/types.
     */
    void resize(int newSize);

    /**
     * Set CTpr parameters by old CT class/type numbers.
     *
     * @param oldtype1        the old CT class/type number 1
     * @param oldtype2        the old CT class/type number 2
     * @param combinedApre    combined apre
     * @param combinedBexp    combined bexp
     */
    void setCTprParametersByOldTypes(int oldtype1, int oldtype2, double combinedApre, double combinedBexp);

    /**
     * Add CTpr parameters by old CT class/type numbers.
     *
     * @param oldtype1        the old CT class/type number 1
     * @param oldtype2        the old CT class/type number 2
     * @param combinedApre    combined apre
     * @param combinedBexp    combined bexp
     * @return the indexed of the added pair
     */
    int addCTprByOldTypes(int oldtype1, int oldtype2, double combinedApre, double combinedBexp);

    /**
     * Set CTpr parameters by new CT class/type numbers.
     *
     * @param ntype1          the new CT class/type number 1
     * @param ntype2          the new CT class/type number 2
     * @param combinedApre    combined apre
     * @param combinedBexp    combined bexp
     */
    void setCTprParameters(int ntype1, int ntype2, double combinedApre, double combinedBexp);

    /**
     * Get CTpr parameters by new CT class/type numbers.
     *
     * @param ntype1          the new CT class/CT number 1
     * @param ntype2          the new CT class/CT number 2
     * @param combinedApre    combined apre
     * @param combinedBexp    combined bexp
     */
    void getCTprParameters(int ntype1, int ntype2, double& combinedApre, double& combinedBexp) const;

    /**
     * Add CTpr parameters by new CT class/type numbers.
     *
     * @param ntype1          the new CT class/type number 1
     * @param ntype2          the new CT class/type number 2
     * @param combinedApre    combined apre
     * @param combinedBexp combined bexp
     */
    int addCTpr(int ntype1, int ntype2, double combinedApre, double combinedBexp);

    /**
     * Set apre combining rule
     *
     * @param apreCombiningRule   apre combining rule:  'ARITHMETIC', 'GEOMETRIC'.
     */
    void setApreCombiningRule(const std::string& apreCombiningRule);

    /**
     * Get apre combining rule
     *
     * @return apreCombiningRule   apre combining rule:  'ARITHMETIC', 'GEOMETRIC'. 
     */
    const std::string& getApreCombiningRule(void) const;

    /**
     * Set bexp combining rule
     *
     * @param bexpCombiningRule   bexp combining rule:   'ARITHMETIC', 'GEOMETRIC'. 
     */
    void setBexpCombiningRule(const std::string& bexpCombiningRule);

    /**
     * Get bexp combining rule
     *
     * @return bexpCombiningRule   bexp combining rule:  'ARITHMETIC', 'GEOMETRIC'.
     */
    const std::string& getBexpCombiningRule(void) const;


    /**
     * Set exclusions for specified particle
     *
     * @param particleIndex particle index
     * @param exclusions vector of exclusions
     */
    void setParticleExclusions(int particleIndex, const std::vector<int>& exclusions);

    /**
     * Get exclusions for specified particle
     *
     * @param particleIndex   particle index
     * @param[out] exclusions vector of exclusions
     */
    void getParticleExclusions(int particleIndex, std::vector<int>& exclusions) const;

    /**
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @return the cutoff distance, measured in nm
     */

    double getCutoffDistance() const;
    
    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);

    /**
     * Set the cutoff distance.
     * 
     * @deprecated This method exists only for backward compatibility.  Use setCutoffDistance() instead.
     */
    void setCutoff(double cutoff);

    /**
     * Get the cutoff distance.
     * 
     * @deprecated This method exists only for backward compatibility.  Use getCutoffDistance() instead.
     */
    double getCutoff() const;

    /**
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;

    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * The only information this method updates is the values of per-particle parameters.  All other aspects of the Force
     * (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be changed by reinitializing the Context.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if nonbondedMethod uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return nonbondedMethod == AmoebaCTForce::CutoffPeriodic;
    }

protected:
    ForceImpl* createImpl() const;

private:
    class CTInfo;
    class CTprInfo;
    NonbondedMethod nonbondedMethod;
    double cutoff;
    int numCTprTypes;

    std::string apreCombiningRule;
    std::string bexpCombiningRule;

    std::vector< std::vector<int> > exclusions; // size = number of atoms
    std::vector<CTInfo> parameters; // size = number of atoms
    std::vector<CTprInfo> arguments; // size = (number of CT classes/types)^2
    std::vector<int> typeMap; // size = number of CT classes/types
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class AmoebaCTForce::CTInfo {
public:
    int CTprType;
    double apre, bexp, lambda;

    CTInfo()
        : CTprType(-1)
        , apre(0.0)
        , bexp(0.0)  
        , lambda(1.0) {}

    CTInfo(int CTprType, double apre, double bexp, double lambda)
        : CTprType(CTprType)
        , apre(apre)
        , bexp(bexp)
        , lambda(lambda) {}
};

/**
 * This is an internal class used to record information about CT pair.
 * @private
 */
class AmoebaCTForce::CTprInfo {
public:
    double combinedApre;
    double combinedBexp;

    CTprInfo()
        : combinedApre(0.0)
        , combinedBexp(0.0) {}

    CTprInfo(double combinedApre, double combinedBexp)
        : combinedApre(combinedApre)
        , combinedBexp(combinedBexp) {}
};
} // namespace OpenMM

#endif /*OPENMM_AMOEBA_CT_FORCE_H_*/

