/** FILE: yy_2lattimedriver.cc            -*-Mode: c++-*-
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

#include <string>

#include "nb.h"
#include "director.h"
#include "simstate.h"
#include "key.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

#include "yy_2lattimedriver.h"
#include "yy_2lattimeevolver.h"

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(YY_2LatTimeDriver);

/* End includes */

// Constructor
YY_2LatTimeDriver::YY_2LatTimeDriver(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : YY_2LatDriver(name,newdtr,argstr), max_dm_dt_obj_ptr(NULL)
{
  // Process arguments
  OXS_GET_INIT_EXT_OBJECT("evolver",YY_2LatTimeEvolver,evolver_obj);
  evolver_key.Set(evolver_obj.GetPtr());
  // Dependency lock on YY_2LatTimeEvolver object is
  // held until *this is destroyed.

  if(!HasInitValue("stopping_dm_dt")) {
    stopping_dm_dt.push_back(0.0); // Default is no control
  } else {
    GetGroupedRealListInitValue("stopping_dm_dt",stopping_dm_dt);
  }

  if(!HasInitValue("stopping_time")) {
    stopping_time.push_back(0.0); // Default is no control
  } else {
    GetGroupedRealListInitValue("stopping_time",stopping_time);
  }

  // Delete basename
  // If there is "basename" option, the following should take it.
  // If there is not, the following returns default value "oxs".
  String basename = GetStringInitValue("basename", "oxs");

  VerifyAllInitArgsUsed();

  last_timestep_output.Setup(
           this,InstanceName(),"Last time step","s",0,
	   &YY_2LatTimeDriver::Fill__last_timestep_output);
  simulation_time_output.Setup(
	   this,InstanceName(),"Simulation time","s",0,
	   &YY_2LatTimeDriver::Fill__simulation_time_output);

  last_timestep_output.Register(director,0);
  simulation_time_output.Register(director,0);

  // Reserve space for initial state (see GetInitialState() below)
  //director->ReserveSimulationStateRequest(3);
}

void YY_2LatTimeDriver::GetInitialState(
    Oxs_ConstKey<Oxs_SimState>& state,
    Oxs_ConstKey<Oxs_SimState>& state1,
    Oxs_ConstKey<Oxs_SimState>& state2)
{
  Oxs_Key<Oxs_SimState> tempstate;
  Oxs_Key<Oxs_SimState> tempstate1;
  Oxs_Key<Oxs_SimState> tempstate2;
  director->GetNewSimulationState(tempstate);
  director->GetNewSimulationState(tempstate1);
  director->GetNewSimulationState(tempstate2);

  SetStartValues(tempstate,tempstate1,tempstate2);

  // Set pointers to the sublattices
  Oxs_SimState& tempstate_ = tempstate.GetWriteReference();
  Oxs_SimState& tempstate1_ = tempstate1.GetWriteReference();
  Oxs_SimState& tempstate2_ = tempstate2.GetWriteReference();
  tempstate_.lattice1 = &tempstate1_;
  tempstate_.lattice2 = &tempstate2_;
  tempstate1_.total_lattice = &tempstate_;
  tempstate1_.lattice2 = &tempstate2_;
  tempstate2_.total_lattice = &tempstate_;
  tempstate2_.lattice1 = &tempstate1_;
  tempstate1_.lattice_type = Oxs_SimState::LATTICE1;
  tempstate2_.lattice_type = Oxs_SimState::LATTICE2;

  tempstate.GetReadReference();  // Release write lock.
  tempstate1.GetReadReference();  // Release write lock.
  tempstate2.GetReadReference();  // Release write lock.
  /// The read lock will be automatically released when the
  /// key "initial_state" is destroyed.

  state = tempstate;
  state1 = tempstate1;
  state2 = tempstate2;

  return;
}

OC_BOOL YY_2LatTimeDriver::Init()
{ 
  YY_2LatDriver::Init();  // Run init routine in parent.
  /// This will call YY_2LatTimeDriver::GetInitialState().

  // Get pointer to output object providing max dm/dt data
  const YY_2LatTimeEvolver* evolver = evolver_key.GetPtr();
  if(evolver==NULL) {
    throw Oxs_ExtError(this,"PROGRAMMING ERROR: No evolver found?");
  }
  String output_name = String(evolver->InstanceName());
  output_name += String(":Max dm/dt");
  max_dm_dt_obj_ptr
    =  director->FindOutputObjectExact(output_name.c_str());
  if(max_dm_dt_obj_ptr==NULL) {
    throw Oxs_ExtError(this,"Unable to identify unique"
                         " Max dm/dt output object");
  }

  return 1;
}

YY_2LatTimeDriver::~YY_2LatTimeDriver()
{}

void YY_2LatTimeDriver::StageRequestCount
(unsigned int& min,
 unsigned int& max) const
{ // Number of stages wanted by driver

  YY_2LatDriver::StageRequestCount(min,max);

  unsigned int count = static_cast<OC_UINT4m>(stopping_dm_dt.size());
  if(count>min) min=count;
  if(count>1 && count<max) max=count;
  // Treat length 1 lists as imposing no upper constraint.

  count =  static_cast<OC_UINT4m>(stopping_time.size());
  if(count>min) min=count;
  if(count>1 && count<max) max=count;
  // Treat length 1 lists as imposing no upper constraint.
}

OC_BOOL
YY_2LatTimeDriver::ChildIsStageDone(
    const Oxs_SimState& state, 
    const Oxs_SimState& state1, 
    const Oxs_SimState& state2) const
{
  OC_UINT4m stage_index = state.stage_number;

  // Stage time check
  OC_REAL8m stop_time=0.;
  if(stage_index >= stopping_time.size()) {
    stop_time = stopping_time[stopping_time.size()-1];
  } else {
    stop_time = stopping_time[stage_index];
  }
  if(stop_time>0.0
     && stop_time-state.stage_elapsed_time<=stop_time*OC_REAL8_EPSILON*2) {
    return 1; // Stage done
  }

  // TODO: check both sublattices?
  // dm_dt check
  Tcl_Interp* mif_interp = director->GetMifInterp();
  if(max_dm_dt_obj_ptr==NULL ||
     max_dm_dt_obj_ptr->Output(&state,mif_interp,0,NULL) != TCL_OK) {
    String msg=String("Unable to obtain Max dm/dt output: ");
    if(max_dm_dt_obj_ptr==NULL) {
      msg += String("PROGRAMMING ERROR: max_dm_dt_obj_ptr not set."
		    " Driver Init() probably not called.");
    } else {
      msg += String(Tcl_GetStringResult(mif_interp));
    }
    throw Oxs_ExtError(this,msg.c_str());
  }
  OC_BOOL err;
  OC_REAL8m max_dm_dt = Nb_Atof(Tcl_GetStringResult(mif_interp),err);
  if(err) {
    String msg=String("Error detected in StageDone method"
		      " --- Invalid Max dm/dt output: ");
    msg += String(Tcl_GetStringResult(mif_interp));
    throw Oxs_ExtError(this,msg.c_str());
  }
  OC_REAL8m stop_dm_dt=0.;
  if(stage_index >= stopping_dm_dt.size()) {
    stop_dm_dt = stopping_dm_dt[stopping_dm_dt.size()-1];
  } else {
    stop_dm_dt = stopping_dm_dt[stage_index];
  }
  if(stop_dm_dt>0.0 && max_dm_dt <= stop_dm_dt) {
    return 1; // Stage done
  }

  // If control gets here, then stage not done
  return 0;
}

OC_BOOL
YY_2LatTimeDriver::ChildIsRunDone(
    const Oxs_SimState& /* state */,
    const Oxs_SimState& /* state1 */,
    const Oxs_SimState& /* state2 */) const
{
  // No child-specific checks at this time...
  return 0; // Run not done
}

void YY_2LatTimeDriver::FillStateSupplemental(Oxs_SimState& work_state) const
{
  OC_REAL8m work_step = work_state.last_timestep;
  OC_REAL8m base_time = work_state.stage_elapsed_time - work_step;

  // Insure that step does not go past stage stopping time
  OC_UINT4m stop_index = work_state.stage_number;
  OC_REAL8m stop_value=0.0;
  if(stop_index >= stopping_time.size()) {
    stop_value = stopping_time[stopping_time.size()-1];
  } else {
    stop_value = stopping_time[stop_index];
  }
  if(stop_value>0.0) {
    OC_REAL8m timediff = stop_value-work_state.stage_elapsed_time;
    if(timediff<=0) { // Over step
      // In the degenerate case where dm_dt=0, work_step will be
      // large (==1) and work_state.stage_elapsed_time will also
      // be large.  In that case, timediff will be numerically
      // poor because stop_value << work_state.stage_elapsed_time.
      // Check for this, and adjust sums accordingly.
      if(work_step>stop_value) { // Degenerate case
        work_step -= work_state.stage_elapsed_time;
        work_step += stop_value;
      } else {                   // Normal case
        work_step += timediff;
      }
      if(work_step<=0.0) work_step = stop_value*OC_REAL8_EPSILON; // Safety
      work_state.last_timestep = work_step;
      work_state.stage_elapsed_time = stop_value;
    } else if(timediff < 2*stop_value*OC_REAL8_EPSILON) {
      // Under step, but close enough for government work
      work_state.last_timestep += timediff;
      work_state.stage_elapsed_time = stop_value;
    } else if(0.25*work_step>timediff) {
      // Getting close to stage boundary.  Foreshorten.
      OC_REAL8m tempstep = (3*work_step+timediff)*0.25;
      work_state.last_timestep = tempstep;
      work_state.stage_elapsed_time = base_time+tempstep;
    }
  }
}

OC_BOOL
YY_2LatTimeDriver::Step(
    Oxs_ConstKey<Oxs_SimState> base_state,
    Oxs_ConstKey<Oxs_SimState> base_state1,
    Oxs_ConstKey<Oxs_SimState> base_state2,
    const Oxs_DriverStepInfo& stepinfo,
    Oxs_Key<Oxs_SimState>& next_state,
    Oxs_Key<Oxs_SimState>& next_state1,
    Oxs_Key<Oxs_SimState>& next_state2)
{ // Returns true if step was successful, false if
  // unable to step as requested.

  // Put write lock on evolver in order to get a non-const
  // pointer.  Use a temporary variable, temp_key, so
  // write lock is automatically removed when temp_key
  // is destroyed.
  Oxs_Key<YY_2LatTimeEvolver> temp_key = evolver_key;
  YY_2LatTimeEvolver& evolver = temp_key.GetWriteReference();
  return evolver.Step(
      this,
      base_state,
      base_state1,
      base_state2,
      stepinfo,
      next_state,
      next_state1,
      next_state2);
}

OC_BOOL
YY_2LatTimeDriver::InitNewStage(
    Oxs_ConstKey<Oxs_SimState> state,
    Oxs_ConstKey<Oxs_SimState> state1,
    Oxs_ConstKey<Oxs_SimState> state2,
    Oxs_ConstKey<Oxs_SimState> prevstate,
    Oxs_ConstKey<Oxs_SimState> prevstate1,
    Oxs_ConstKey<Oxs_SimState> prevstate2)
{
  // Put write lock on evolver in order to get a non-const
  // pointer.  Use a temporary variable, temp_key, so
  // write lock is automatically removed when temp_key
  // is destroyed.
  Oxs_Key<YY_2LatTimeEvolver> temp_key = evolver_key;
  YY_2LatTimeEvolver& evolver = temp_key.GetWriteReference();
  // TODO: Default is NOP but check with the evolver.
  OC_BOOL result = evolver.InitNewStage(this,state,prevstate);
  OC_BOOL result1 = evolver.InitNewStage(this,state1,prevstate1);
  OC_BOOL result2 = evolver.InitNewStage(this,state2,prevstate2);
  return result & result1 & result2;
  //return evolver.InitNewStage(this,state,prevstate);  // Default is NOP
}


////////////////////////////////////////////////////////////////////////
// State-based outputs, maintained by the driver.  These are
// conceptually public, but are specified private to force
// clients to use the output_map interface in Oxs_Director.

#define OSO_FUNC(NAME) \
void YY_2LatTimeDriver::Fill__##NAME##_output(const Oxs_SimState& state) \
{ NAME##_output.cache.state_id=state.Id(); \
  NAME##_output.cache.value=state.NAME; }

OSO_FUNC(last_timestep)

void
YY_2LatTimeDriver::Fill__simulation_time_output(const Oxs_SimState& state)
{
  simulation_time_output.cache.state_id = state.Id();
  simulation_time_output.cache.value =
    state.stage_start_time + state.stage_elapsed_time;
}
