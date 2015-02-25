/** FILE: yy_2lattimeevolver.cc                 -*-Mode: c++-*-
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

#include "director.h"
#include "energy.h"

#include "yy_2lattimeevolver.h"

/* End includes */

// Constructors
YY_2LatTimeEvolver::YY_2LatTimeEvolver
(const char* name,     // Child instance id
 Oxs_Director* newdtr) // App director
  : Oxs_Evolver(name,newdtr), energy_calc_count(0)
{}

YY_2LatTimeEvolver::YY_2LatTimeEvolver
(const char* name,
 Oxs_Director* newdtr,
 const char* argstr)      // MIF block argument string
  : Oxs_Evolver(name,newdtr,argstr), energy_calc_count(0)
{
  total_energy_output.Setup(this,InstanceName(),
                            "Total energy","J",1,
                            &YY_2LatTimeEvolver::UpdateEnergyOutputs);
  total_energy_density_output.Setup(this,InstanceName(),
                           "Total energy density","J/m^3",1,
                           &YY_2LatTimeEvolver::UpdateEnergyOutputs);
  total_field_output.Setup(this,InstanceName(),
                           "Total field","A/m",1,
                           &YY_2LatTimeEvolver::UpdateEnergyOutputs);
  energy_calc_count_output.Setup(this,InstanceName(),
                           "Energy calc count","",0,
                           &YY_2LatTimeEvolver::FillEnergyCalcCountOutput);
  /// Note: MSVC++ 6.0 requires fully qualified member names

  total_energy_output.Register(director,-5);
  total_energy_density_output.Register(director,-5);
  total_field_output.Register(director,-5);
  energy_calc_count_output.Register(director,-5);

  // Eventually, caching should be handled by controlling Tcl script.
  // Until then, request caching of scalar energy output by default.
  total_energy_output.CacheRequestIncrement(1);

}


OC_BOOL YY_2LatTimeEvolver::Init()
{
#if REPORT_TIME
  Oc_TimeVal cpu,wall;
  steponlytime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"   Step-only    .....   %7.2f cpu /%7.2f wall,"
            " (%.1000s)\n",
            double(cpu),double(wall),InstanceName());
  }
  steponlytime.Reset();
#endif // REPORT_TIME

  energy_calc_count=0;

  // Release scratch space.
  temp_energy.Release();
  temp_field.Release();

  return Oxs_Evolver::Init();
}


YY_2LatTimeEvolver::~YY_2LatTimeEvolver()
{
#if REPORT_TIME
  Oc_TimeVal cpu,wall;
  steponlytime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"   Step-only    .....   %7.2f cpu /%7.2f wall,"
            " (%.1000s)\n",
            double(cpu),double(wall),InstanceName());
  }
#endif // REPORT_TIME
}


// GetEnergyDensity: Note that mxH is returned, as opposed to MxH.
// This relieves this routine from needing to know what Ms is, and saves
// an unneeded multiplication (since the evolver is just going to divide
// it back out again to calculate dm/dt (as opposed again to dM/dt)).
// The returned energy array is average energy density for the
// corresponding cell in J/m^3; mxH is in A/m, pE_pt (partial derivative
// of E with respect to t) is in J/s.  Any of mxH or H may be
// NULL, which disables assignment for that variable.
void YY_2LatTimeEvolver::GetEnergyDensity
(const Oxs_SimState& state,
 Oxs_MeshValue<OC_REAL8m>& energy,
 Oxs_MeshValue<ThreeVector>* mxH1_req,
 Oxs_MeshValue<ThreeVector>* mxH2_req,
 Oxs_MeshValue<ThreeVector>* H1_req,
 Oxs_MeshValue<ThreeVector>* H2_req,
 OC_REAL8m& pE_pt,
 OC_REAL8m& total_E)
{
  // Update call count
  ++energy_calc_count;

  // Extract simulation states for sublattices
  const Oxs_SimState& state1 = *(state.lattice1);
  const Oxs_SimState& state2 = *(state.lattice2);
  Oxs_MeshValue<OC_REAL8m> energy1, energy2;

  /// If field is requested by both H_req and total_field_output,
  /// then fill H_req first, and copy to field_cache at end.
  // lattice 1
  Oxs_MeshValue<ThreeVector>* H1_fill = H1_req;
  Oxs_MeshValue<ThreeVector>* field1_cache = NULL;
  if(total_field_output.GetCacheRequestCount()>0) {
    total_field_output.cache.state_id=0;
    field1_cache = &total_field_output.cache.value;
    if(H1_fill==NULL) H1_fill = field1_cache;
  }
  // lattice 2
  Oxs_MeshValue<ThreeVector>* H2_fill = H2_req;
  Oxs_MeshValue<ThreeVector>* field2_cache = NULL;
  if(total_field_output.GetCacheRequestCount()>0) {
    total_field_output.cache.state_id=0;
    field2_cache = &total_field_output.cache.value;
    if(H2_fill==NULL) H2_fill = field2_cache;
  }

  // Set up energy computation output data structure
  Oxs_ComputeEnergyData oced1(state);
  oced1.scratch_energy = &temp_energy;
  oced1.scratch_H      = &temp_field;
  oced1.energy_accum   = &energy1;
  oced1.H_accum        = H1_fill;
  oced1.mxH_accum      = mxH1_req;
  oced1.energy         = NULL;  // Required null
  oced1.H              = NULL;  // Required null
  oced1.mxH            = NULL;  // Required null

  Oxs_ComputeEnergyData oced2(state);
  oced2.scratch_energy = &temp_energy;
  oced2.scratch_H      = &temp_field;
  oced2.energy_accum   = &energy2;
  oced2.H_accum        = H2_fill;
  oced2.mxH_accum      = mxH2_req;
  oced2.energy         = NULL;  // Required null
  oced2.H              = NULL;  // Required null
  oced2.mxH            = NULL;  // Required null

  UpdateFixedSpinList(state1.mesh);
  UpdateFixedSpinList(state2.mesh);

  // TODO: Set fixed spin for each sublattice?
  Oxs_ComputeEnergyExtraData oceed(GetFixedSpinList(),0);

  // Compute total energy and torque
#if REPORT_TIME
  OC_BOOL sot_running = steponlytime.IsRunning();
  if(sot_running) {
    steponlytime.Stop();
  }
#endif // REPORT_TIME
  YY_2LatComputeEnergies(state,oced1,oced2,director->GetEnergyObjects(),oceed);
#if REPORT_TIME
  if(sot_running) {
    steponlytime.Start();
  }
#endif // REPORT_TIME

  // Add energy density for sublattices
  // TODO: Does this involve double counting?
  energy.AdjustSize(state.mesh);
  for(OC_INDEX i=0; i<state.mesh->Size(); i++) {
    energy[i] = energy1[i] + energy2[i];
  }

  if(total_energy_density_output.GetCacheRequestCount()>0) {
    // Energy density field output requested.  Copy results
    // to output cache.
    total_energy_density_output.cache.state_id=0;
    total_energy_density_output.cache.value = energy;
    total_energy_density_output.cache.state_id=state.Id();
  }

  if(field1_cache!=NULL) {
    if(field1_cache!=H1_fill) {
      // Field requested by both H_req and total_field_output,
      // so copy from H_req to field1_cache.
      *field1_cache = *H1_fill;
    }
    field1_cache=NULL; // Safety
    total_field_output.cache.state_id=state.Id();
  }

  if(field2_cache!=NULL) {
    if(field2_cache!=H2_fill) {
      // Field requested by both H_req and total_field_output,
      // so copy from H_req to field1_cache.
      *field2_cache = *H2_fill;
    }
    field2_cache=NULL; // Safety
    //total_field_output.cache.state_id=state.Id();
  }

  // Store total energy sum if output object total_energy_output
  // has cache enabled.
  if (total_energy_output.GetCacheRequestCount()>0) {
    total_energy_output.cache.state_id=0;
    total_energy_output.cache.value=oced1.energy_sum;
    total_energy_output.cache.state_id=state.Id();
  }

  pE_pt = oced1.pE_pt;  // Export pE_pt value
  total_E = oced1.energy_sum; // Export total energy
}

void YY_2LatTimeEvolver::UpdateEnergyOutputs(const Oxs_SimState& state)
{
  if(state.Id()==0) { // Safety
    return;
  }
  Oxs_MeshValue<OC_REAL8m> energy(state.mesh);
  OC_REAL8m pE_pt;
  GetEnergyDensity(state,energy,NULL,NULL,NULL,NULL,pE_pt);
}

void YY_2LatTimeEvolver::FillEnergyCalcCountOutput(const Oxs_SimState& state)
{
  energy_calc_count_output.cache.state_id=state.Id();
  energy_calc_count_output.cache.value
    = static_cast<OC_REAL8m>(energy_calc_count);
}
