// ParFormula.h
//
// David Adams
// July 2021
//
// Interface for a formula with N variables and M named parameters.
// Methods are provided to set the parameters and evaluate
// given variables using those settings.

#ifndef ParFormula_H
#define ParFormula_H

#include <string>
#include <vector>

class ParFormula {

public:

  using Value = double;   // Type for variables, parameters and evaluation
  using Values = std::vector<Value>;
  using Index = unsigned int;
  using Name = std::string;
  using Names = std::vector<Name>;

public:  // const methods

  // Formula name.
  virtual Name name() const { return ""; }

  // String representation of the formula.
  // Expected but not required that prameters are in square brackets
  // and variables x, y, z or x[0], x[1], ....
  // E.g. for parameters Par1 and Par2:   [Par1]*sin(x + [Par2])
  virtual Name formulaString() const =0;

  // Number of variables.
  virtual Index nvar() const =0;

  // Parameter count and names.
  virtual Index npar() const { return pars().size(); }
  virtual Names pars() const =0;

  // Return the set, unset and reset parameter names.
  virtual Names setPars() const =0;
  virtual Names unsetPars() const =0;
  virtual Names resetPars() const =0;

  // Return if all parameters are set.
  virtual bool ready() const { return setPars().size() == npar(); }

  // Default evaluation value.
  virtual Value defaultEval() const { return 0.0; }

  // Evaluate the formula.
  virtual double eval(const Values& vars) const =0;
  virtual double eval(Value var) const =0;

public:  // non-const methods

  // Set a parameter value.
  // Return 0 for success.
  virtual int setParValue(Name parnam, Value parval) =0;

  // Set the default evaluation value.
  // Return 0 for success.
  virtual int setDefaultEval(Value val) { return 1; }

  // Unset the parameters.
  virtual int unsetParValues() =0;

};

#endif
