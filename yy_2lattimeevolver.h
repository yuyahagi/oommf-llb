/** FILE: yy_2lattimeevolver.h                 -*-Mode: c++-*-
 *
 * Abstract time evolver class for two lattices
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

#ifndef _YY_2LATTIMEEVOLVER
#define _YY_2LATTIMEEVOLVER

#include "evolver.h"
#include "output.h"

/* End includes */

class YY_2LatTimeDriver; // Forward references
//struct YY_2LatDriverStepInfo;
struct Oxs_DriverStepInfo;

class YY_2LatTimeEvolver:public Oxs_Evolver {
private:

  OC_UINT4m energy_calc_count; // Number of times GetEnergyDensity
  /// has been called in current problem run.

  Oxs_MeshValue<OC_REAL8m> temp_energy;     // Scratch space used by
  Oxs_MeshValue<ThreeVector> temp_field; // GetEnergyDensity().

  // Outputs maintained by this interface layer.  These are conceptually
  // public, but are specified private to force clients to use the
  // output_map interface.
  void UpdateEnergyOutputs(const Oxs_SimState&);
  void FillEnergyCalcCountOutput(const Oxs_SimState&);
  Oxs_ScalarOutput<YY_2LatTimeEvolver> total_energy_output;
  Oxs_ScalarFieldOutput<YY_2LatTimeEvolver> total_energy_density_output;
  Oxs_VectorFieldOutput<YY_2LatTimeEvolver> total_field_output;
  Oxs_ScalarOutput<YY_2LatTimeEvolver> energy_calc_count_output;

  // Disable copy constructor and assignment operator by declaring
  // them without defining them.
  YY_2LatTimeEvolver(const YY_2LatTimeEvolver&);
  YY_2LatTimeEvolver& operator=(const YY_2LatTimeEvolver&);

protected:

#if REPORT_TIME
  mutable Nb_StopWatch steponlytime;
#endif

  YY_2LatTimeEvolver(const char* name,      // Child instance id
                 Oxs_Director* newdtr);  // App director
  YY_2LatTimeEvolver(const char* name,
                 Oxs_Director* newdtr,
                 const char* argstr);      // MIF block argument string

  virtual OC_BOOL Init();  // All children of YY_2LatTimeEvolver *must*
  /// call this function in their Init() routines.  The main purpose
  /// of this function is to initialize output variables.

  void GetEnergyDensity(const Oxs_SimState& state,
                        Oxs_MeshValue<OC_REAL8m>& energy,
                        Oxs_MeshValue<ThreeVector>* mxH_req,
                        Oxs_MeshValue<ThreeVector>* H_req,
                        OC_REAL8m& pE_pt,
                        OC_REAL8m& total_E);

  void GetEnergyDensity(const Oxs_SimState& state,
			Oxs_MeshValue<OC_REAL8m>& energy,
			Oxs_MeshValue<ThreeVector>* mxH_req,
			Oxs_MeshValue<ThreeVector>* H_req,
			OC_REAL8m& pE_pt) {
    // This interface for backwards compatibility
    OC_REAL8m dummy_E;
    GetEnergyDensity(state,energy,mxH_req,H_req,pE_pt,dummy_E);
  }

public:
  virtual ~YY_2LatTimeEvolver();

  virtual OC_BOOL
  InitNewStage(const YY_2LatTimeDriver* /* driver */,
               Oxs_ConstKey<Oxs_SimState> /* state */,
               Oxs_ConstKey<Oxs_SimState> /* prevstate */) { return 1; }
  /// Default implementation is a NOP.  Children may override.
  /// NOTE: prevstate may be "INVALID".  Children should check
  ///       before use.

  virtual OC_BOOL
  Step(const YY_2LatTimeDriver* driver,
       Oxs_ConstKey<Oxs_SimState> current_state,
       //const YY_2LatDriverStepInfo& step_info,
       const Oxs_DriverStepInfo& step_info,
       Oxs_Key<Oxs_SimState>& next_state) = 0;
  // Returns true if step was successful, false if
  // unable to step as requested.  The evolver object
  // is responsible for calling driver->FillState()
  // and driver->FillStateSupplemental() to fill
  // next_state as needed.
};

#endif // _YY_2LATTIMEEVOLVER
