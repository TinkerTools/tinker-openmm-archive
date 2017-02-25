
%pythonprepend OpenMM::Platform::registerPlatform(Platform *platform) %{
    if not platform.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomNonbondedForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomNonbondedForce::addInteractionGroup(const std::set< int > &set1, const std::set< int > &set2) %{
    set1 = list(set1)
    set2 = list(set2)
%}
%pythonprepend OpenMM::CustomNonbondedForce::setInteractionGroupParameters(int index, const std::set< int > &set1, const std::set< int > &set2) %{
    set1 = list(set1)
    set2 = list(set2)
%}
%pythonprepend OpenMM::CustomHbondForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CompoundIntegrator::addIntegrator(Integrator *integrator) %{
    if not integrator.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::System::setVirtualSite(int index, VirtualSite *virtualSite) %{
    if not virtualSite.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::System::addForce(Force *force) %{
    if not force.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomCompoundBondForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomManyParticleForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomGBForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}