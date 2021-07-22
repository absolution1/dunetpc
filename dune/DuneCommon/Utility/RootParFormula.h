// RootParFormula.h
//
// David Adams
// July 2021
//
// Implementation of ParFormula base on TFormula.

#ifndef RootParFormula_H
#define RootParFormula_H

#include "dune/DuneInterface/Utility/ParFormula.h"
#include "TFormula.h"

class RootParFormula : public ParFormula {

public:

  // Ctor from a string formula.
  RootParFormula(Name snam, Name sform);

  // Name.
  Name name() const override { return m_ptf->GetName(); }

  // Formula.
  Name formulaString() const override { return m_sform; }

  // Variable dimension.
  Index nvar() const override { return m_ptf->GetNdim(); }

  // Parameter counts.
  Index npar() const override { return m_parNames.size(); }
  Names pars() const override { return m_parNames; }
  Names setPars() const override;
  Names unsetPars() const override;
  Names resetPars() const override;

  // Return if a parameter appears in the equation.
  bool isPar(Name parnam) const override;

  // Default return.
  Value defaultEval() const override { return m_defval; }

  // Evaluate.
  double eval(const Values& vars) const override {
    return ready() && nvar() <= vars.size()  ? m_ptf->EvalPar(&vars[0], &m_parValues[0]) : m_defval;
  }
  double eval(Value var) const override {
    return ready() && nvar() <= 1 ? m_ptf->Eval(var) : m_defval;
  }

  // Set a parameter value.
  int setParValue(Name parnam, Value parval) override;

  // Set default return.
  int setDefaultEval(Value val) override { m_defval = val; return 0; }

  // Unset paramter values.
  int unsetParValues() override;

private:

  TFormula* m_ptf;
  Name m_sform;
  Names m_parNames;
  Values m_parValues;
  std::vector<Index> m_setCounts;
  Value m_defval;

};

#endif

