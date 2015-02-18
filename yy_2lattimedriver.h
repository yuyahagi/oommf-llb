/** FILE: yy_2lattimedriver.h            -*-Mode: c++-*-
 *
 * Concrete YY_2LatDriver class for ferrimagnet simulation
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

#ifndef _YY_2LATTIMEDRIVER
#define _YY_2LATTIMEDRIVER

#include <vector>

#include "simstate.h"
#include "key.h"
#include "yy_2latdriver.h"
#include "yy_2lattimeevolver.h"

OC_USE_STD_NAMESPACE;

/* End includes */

class YY_2LatTimeDriver:public YY_2LatDriver {
private:
  Oxs_OwnedPointer<YY_2LatTimeEvolver> evolver_obj; // Evolver basket
  Oxs_Key<YY_2LatTimeEvolver> evolver_key;

  vector<OC_REAL8m> stopping_time;  // Seconds
  vector<OC_REAL8m> stopping_dm_dt; // deg/ns

  Oxs_Output* max_dm_dt_obj_ptr; // Pointer to object providing
  /// max dm/dt data.  This is needed to determine StageDone events.

  // State-based outputs, maintained by the driver.  These are
  // conceptually public, but are specified private to force
  // clients to use the output_map interface in Oxs_Director.
#define OSO_DECL(name) \
void Fill__##name##_output(const Oxs_SimState&); \
Oxs_ScalarOutput<YY_2LatTimeDriver> name##_output
  OSO_DECL(last_timestep);
  OSO_DECL(simulation_time);
#undef OSO_DECL

  // Done checks, called by parent YY_2LatDriver::IsStageDone and
  // YY_2LatDriver::IsRunDone functions.
  virtual OC_BOOL ChildIsStageDone(const Oxs_SimState& state) const;
  virtual OC_BOOL ChildIsRunDone(const Oxs_SimState& state) const;

  // Disable copy constructor and assignment operator by declaring
  // them without defining them.
  YY_2LatTimeDriver(const YY_2LatTimeDriver&);
  YY_2LatTimeDriver& operator=(const YY_2LatTimeDriver&);

public:
  YY_2LatTimeDriver(const char* name,    // Child instance id
                 Oxs_Director* newdtr, // App director
                 const char* argstr);  // MIF input block parameters
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  virtual OC_BOOL Init();
  virtual ~YY_2LatTimeDriver();

  virtual void StageRequestCount(unsigned int& min,
				 unsigned int& max) const;
  // Number of stages wanted by driver

  // Generic interface
  virtual Oxs_ConstKey<Oxs_SimState> GetInitialState() const;

  // Use FillState* and FillNewStageState* routines inherited from
  // parent.

  virtual OC_BOOL InitNewStage(Oxs_ConstKey<Oxs_SimState> state,
                            Oxs_ConstKey<Oxs_SimState> prevstate);

  virtual  OC_BOOL
  Step(Oxs_ConstKey<Oxs_SimState> current_state,
       const Oxs_DriverStepInfo& step_info,
       Oxs_Key<Oxs_SimState>& next_state);
  //Step(Oxs_ConstKey<Oxs_SimState> current_state,
  //     const YY_2LatDriverStepInfo& step_info,
  //     Oxs_Key<Oxs_SimState>& next_state);
  // Returns true if step was successful, false if
  // unable to step as requested.

  // Time driver interface
  virtual void FillStateSupplemental(Oxs_SimState& work_state) const;
  /// FillStateSupplemental is called from time evolvers to adjust timestep.
};

#endif // _YY_2LATTIMEDRIVER
