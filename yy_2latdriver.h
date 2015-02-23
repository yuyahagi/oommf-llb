/** FILE: yy_2latdriver.h                 -*-Mode: c++-*-
 *
 * Abstract YY_2LatDriver class for ferrimagnet simulation. It fills states 
 * and initiates steps and registers the command "Oxs_Run" with the Tcl
 * interpreter.
 *
 * Copyright (C) 2015 Yu Yahagi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _YY_2LATDRIVER
#define _YY_2LATDRIVER

#include <string>

#include "nb.h"
#include "ext.h"
#include "labelvalue.h"
#include "mesh.h"
#include "simstate.h"
#include "key.h"
#include "outputderiv.h"
#include "scalarfield.h"
#include "vectorfield.h"

#include "driver.h"

OC_USE_STRING;

/* End includes */

class YY_2LatEvolver; // Forward references
struct OxsRunEvent;

class YY_2LatDriver: public Oxs_Driver {
private:
  Oxs_ConstKey<Oxs_SimState> current_state1;
  Oxs_ConstKey<Oxs_SimState> current_state2;

  mutable Oxs_MeshValue<OC_REAL8m> Ms1, Ms2;  // Saturation magnetization
  mutable Oxs_MeshValue<OC_REAL8m> Ms1_inverse, Ms2_inverse;  // 1/Ms
  mutable Oxs_MeshValue<OC_REAL8m> Ms01, Ms02;  // Saturation magnetization
  mutable Oxs_MeshValue<OC_REAL8m> Ms01_inverse, Ms02_inverse;  // 1/Ms0
  Oxs_OwnedPointer<Oxs_VectorField> m01, m02; // Initial spin configuration

  // Additional outputs for sublattices
  Oxs_VectorFieldOutput<YY_2LatDriver> spin1_output;
  Oxs_VectorFieldOutput<YY_2LatDriver> spin2_output;
  void Fill__spin1_output(const Oxs_SimState&);
  void Fill__spin2_output(const Oxs_SimState&);
  Oxs_VectorFieldOutput<YY_2LatDriver> magnetization1_output;
  Oxs_VectorFieldOutput<YY_2LatDriver> magnetization2_output;
  void Fill__magnetization1_output(const Oxs_SimState&);
  void Fill__magnetization2_output(const Oxs_SimState&);
  Oxs_ScalarOutput<YY_2LatDriver> aveMx1_output, aveMx2_output;
  Oxs_ScalarOutput<YY_2LatDriver> aveMy1_output, aveMy2_output;
  Oxs_ScalarOutput<YY_2LatDriver> aveMz1_output, aveMz2_output;
  void Fill__aveMsub_output(const Oxs_SimState&);
  // Override output function for average M
  void Fill__aveM_output(const Oxs_SimState&);

  void UpdateSpinAngleData(
      const Oxs_SimState& state,
      const Oxs_SimState& state1,
      const Oxs_SimState& state2) const;

  // Internal "Run" interface.  The difference with the external
  // interface is that the internal version includes a stage_increment
  // parameter for use by the SetStage method.
  virtual void Run(vector<OxsRunEvent>& results,OC_INT4m stage_increment);
  // Overriding Oxs_Driver::Run().

  // Routines defined by child classes to detect stage and
  // problem end events.
  virtual OC_BOOL ChildIsStageDone(
      const Oxs_SimState& state,
      const Oxs_SimState& state1,
      const Oxs_SimState& state2) const =0;
  virtual OC_BOOL ChildIsRunDone(
      const Oxs_SimState& state,
      const Oxs_SimState& state1,
      const Oxs_SimState& state2) const =0;
  // Disable default functions for standard simulation
  virtual OC_BOOL ChildIsStageDone(
      const Oxs_SimState& state) const
  { return 0; }
  virtual OC_BOOL ChildIsRunDone(
      const Oxs_SimState& state) const
  { return 0; }

  // Disable copy constructor and assignment operator by declaring
  // them without defining them.
  YY_2LatDriver(const YY_2LatDriver&);
  YY_2LatDriver& operator=(const YY_2LatDriver&);

protected:
  YY_2LatDriver(const char* name,        // Child instance id
             Oxs_Director* newdtr,    // App director
	     const char* argstr);     // Args

  // Handler of two sublattice simulation states
  void CombineSubLattices(Oxs_SimState& istate1,
      Oxs_SimState& istate2) const;

  void SetStartValues(
      Oxs_Key<Oxs_SimState>& initial_state,
      Oxs_Key<Oxs_SimState>& initial_state1,
      Oxs_Key<Oxs_SimState>& initial_state2) const;

  // Disable SetStartValues for 1 lattice simulation for safety
  void SetStartValues(Oxs_SimState& istate) const
  { return; }
  void SetStartValues(Oxs_Key<Oxs_SimState>& initial_state) const
  { return; }

  virtual OC_BOOL Init();  // All children of YY_2LatDriver *must* call
  /// this function in their Init() routines.  The main purpose
  /// of this function is to initialize the current state.

public:

  virtual ~YY_2LatDriver();

  // Disable default GetInitialState()
  virtual Oxs_ConstKey<Oxs_SimState> GetInitialState() const
  { 
    Oxs_ConstKey<Oxs_SimState> dummy;
    return dummy;
  }

  virtual void GetInitialState(
      Oxs_ConstKey<Oxs_SimState>& state,
      Oxs_ConstKey<Oxs_SimState>& state1,
      Oxs_ConstKey<Oxs_SimState>& state2) =0;

  OC_BOOL IsStageDone(
      const Oxs_SimState& state,
      const Oxs_SimState& state1,
      const Oxs_SimState& state2) const;
  OC_BOOL IsRunDone(
      const Oxs_SimState& state,
      const Oxs_SimState& state1,
      const Oxs_SimState& state2) const;

  virtual  OC_BOOL
  Step(Oxs_ConstKey<Oxs_SimState> current_state,
      Oxs_ConstKey<Oxs_SimState> current_state1,
      Oxs_ConstKey<Oxs_SimState> current_state2,
       const Oxs_DriverStepInfo& step_info,
       Oxs_Key<Oxs_SimState>& next_state,
       Oxs_Key<Oxs_SimState>& next_state1,
       Oxs_Key<Oxs_SimState>& next_state2)=0;

  // External problem "Run" interface; called from director.
  virtual void Run(vector<OxsRunEvent>& results) {
    this->Run(results,1); // Call internal Run interface with default
    // stage increment (==1).
  }

  virtual OC_BOOL InitNewStage(
      Oxs_ConstKey<Oxs_SimState> state,
      Oxs_ConstKey<Oxs_SimState> state1,
      Oxs_ConstKey<Oxs_SimState> state2,
      Oxs_ConstKey<Oxs_SimState> prevstate,
      Oxs_ConstKey<Oxs_SimState> prevstate1,
      Oxs_ConstKey<Oxs_SimState> prevstate2) =0;
};

#endif // _YY_2LATDRIVER
