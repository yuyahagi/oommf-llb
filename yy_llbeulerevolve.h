/** FILE: yy_llbeulerevolve.h                 -*-Mode: c++-*-
 *
 * Euler evolver class for Landau-Lifshitz-Bloch equation including thermal 
 * fluctuations. It is based on thetaevolve.cc and .h written by Oliver 
 * Lemcke released under GPLv2 license.
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

#ifndef _YY_LLBEULEREVOLVE
#define _YY_LLBEULEREVOLVE

#include "nb.h"

#include "timeevolver.h"
#include "key.h"
#include "tclcommand.h"
#include "output.h"
#include "scalarfield.h"

/* End includes */

#define DEFAULT_M_E_TOL 1e-4

class YY_LLBEulerEvolve:public Oxs_TimeEvolver {
private:
  mutable OC_UINT4m mesh_id;

  // =======================================================================
  // Stepsize control and error criteria.
  // =======================================================================
  // At this point, adaptive stepsize is not implemented for T > 0K.
  OC_REAL8m min_timestep;   // Seconds
  OC_REAL8m max_timestep;   // Seconds
  OC_REAL8m fixed_timestep; // Seconds -> min_timestep = max_timestep

  OC_REAL8m allowed_error_rate;
  OC_REAL8m allowed_absolute_step_error;
  OC_REAL8m allowed_relative_step_error;
  OC_REAL8m step_headroom;

  OC_REAL8m start_dm;

  const OC_UINT4m energy_accum_count_limit ;
  OC_UINT4m energy_accum_count;

  // Evolver control parameters
  OC_BOOL do_precess;
  OC_BOOL allow_signed_gamma;
  enum GammaStyle { GS_INVALID, GS_LL, GS_G };
  GammaStyle gamma_style; // Landau-Lifshitz or Gilbert
  OC_BOOL use_stochastic; // Include stochastic field

  // =======================================================================
  // Spatially variable coefficients
  // =======================================================================
  // Some of them are temperature dependent (e.g., alpha_t, alpha_l) and
  // have to be updated after changing temperature.
  Oxs_OwnedPointer<Oxs_ScalarField> gamma_init;
  mutable Oxs_MeshValue<OC_REAL8m> gamma;             // LL gyromagnetic ratio
  Oxs_OwnedPointer<Oxs_ScalarField> alpha_t_init;
  mutable Oxs_MeshValue<OC_REAL8m> alpha_t, alpha_t0; // transverse
  mutable Oxs_MeshValue<OC_REAL8m> alpha_l;           // longitudinal
  // alpha_t0 is the value at T = 0 K.

  void UpdateMeshArrays(const Oxs_SimState& state);

  // Parameters used for longitudinal susceptibility
  // Exchange parameter J = nJ_0 and atomistic magnetic moment mu, where
  // n is the number of neighboring atoms
  Oxs_OwnedPointer<Oxs_ScalarField> J_init, mu_init;
  Oxs_MeshValue<OC_REAL8m> J, mu;
  // Currie temperature in Kelvin, calculated from J and mu
  mutable Oxs_MeshValue<OC_REAL8m> Tc;

  // Members for calculating m_e, equilibrium spin polarization at
  // temperature T and chi_l, longitudinal susceptibility.
  mutable Oxs_MeshValue<OC_REAL8m> m_e, chi_l;
  void CalculateLongField(const Oxs_SimState& state,
      Oxs_MeshValue<ThreeVector>& longfield) const;
  // Langevin function and its derivative
  OC_REAL8m Langevin(OC_REAL8m x) const;
  OC_REAL8m LangevinDeriv(OC_REAL8m x) const;
  void Update_m_e_chi_l(OC_REAL8m tol) const;
  void Update_m_e_chi_l() const {
    return Update_m_e_chi_l(DEFAULT_M_E_TOL);
  }

  // =======================================================================
  // Caches and scratch spaces
  // =======================================================================
  // Data cached from last state
  OC_UINT4m energy_state_id;
  Oxs_MeshValue<OC_REAL8m> energy;
  OC_REAL8m next_timestep;

  // Ms MeshValue array for the new state. This is required in order to
  // store trial Ms values while keeping the original Ms.
  Oxs_MeshValue<OC_REAL8m> Ms_B, Ms_inverse_B;
  Oxs_MeshValue<OC_REAL8m> *Ms_ptr_A, *Ms_ptr_B;
  Oxs_MeshValue<OC_REAL8m> *Ms_inverse_ptr_A, *Ms_inverse_ptr_B;

  // Scratch space
  Oxs_MeshValue<OC_REAL8m> new_energy;
  Oxs_MeshValue<ThreeVector> new_dm_dt_t;
  Oxs_MeshValue<ThreeVector> new_dm_dt_l;

  Oxs_MeshValue<ThreeVector> total_field;

  // =======================================================================
  // Support for stage-varying temperature
  // =======================================================================
  const OC_REAL8m KBoltzmann;           // Boltzmann constant
  Oxs_OwnedPointer<Oxs_ScalarField> temperature_init;
  Oxs_MeshValue<OC_REAL8m> temperature; // in Kelvin
  Oxs_MeshValue<OC_REAL8m> kB_T;        // KBoltzmann*temperature
  // Make sure kB_T gets updated when temperature is changed.
  OC_UINT4m iteration_Tcalculated;
  // Stores iteration for whic the stochastic field is calculated.

  OC_BOOL has_tempscript;
  vector<Nb_TclCommandLineOption> tempscript_opts;
  Nb_TclCommand tempscript_cmd;
  // Stores last stage number to detect stage change
  OC_INDEX last_stage_number;
  // The following also updates kB_T.
  void UpdateStageTemperature(const Oxs_SimState& stage);

  // =======================================================================
  // Random functions and supports (for stochastic field)
  // =======================================================================
  OC_REAL8m Gaussian_Random(const OC_REAL8m muGaus,
      const OC_REAL8m  sigmaGaus);
  OC_BOOL gaus2_isset;    
  OC_REAL8m gaus2;        
  // returns a gaussian distributed random number 
  // with average muGaus und standard deviation sigmaGaus.
  // The algorithm here computes two values but returns only one at a time.
  // The second value (and a flag) is stored in gaus2 and is returned 
  // second time the method is called.

  // seed to initialize the generator with, can be any integer
  // (beware: -|n| is used in this method)
  OC_INT4m uniform_seed;  
  OC_BOOL has_uniform_seed;

  // constant part of the variance of the thermal field
  Oxs_MeshValue<OC_REAL8m> hFluctVarConst_t;
  Oxs_MeshValue<OC_REAL8m> hFluctVarConst_l;
  void FillHFluctConst(const Oxs_Mesh* mesh);

  // Current values of the thermal field
  // These values are stored in arrays because dm_dt is sometimes 
  // calculated more than once per iteration.
  Oxs_MeshValue<ThreeVector> hFluct_t;  // transverse
  Oxs_MeshValue<ThreeVector> hFluct_l;  // longitudinal

  void Calculate_dm_dt
  (const Oxs_SimState& state_,
   const Oxs_MeshValue<ThreeVector>& mxH_,
   const Oxs_MeshValue<ThreeVector>& total_field_,
   OC_REAL8m pE_pt_,
   Oxs_MeshValue<ThreeVector>& dm_dt_t_,
   Oxs_MeshValue<ThreeVector>& dm_dt_l_,
   OC_REAL8m& max_dm_dt_,
   OC_REAL8m& dE_dt_,
   OC_REAL8m& min_timestep_);
  /// Imports: state_, mxH_, pE_pt
  /// Exports: dm_dt_t_, dm_dt_l_, max_dm_dt_, dE_dt_

  // =======================================================================
  // Outputs
  // =======================================================================
  void UpdateDerivedOutputs(const Oxs_SimState&);
  Oxs_ScalarOutput<YY_LLBEulerEvolve> max_dm_dt_output;
  Oxs_ScalarOutput<YY_LLBEulerEvolve> dE_dt_output;
  Oxs_ScalarOutput<YY_LLBEulerEvolve> delta_E_output;
  Oxs_VectorFieldOutput<YY_LLBEulerEvolve> dm_dt_t_output;
  Oxs_VectorFieldOutput<YY_LLBEulerEvolve> dm_dt_l_output;
  Oxs_VectorFieldOutput<YY_LLBEulerEvolve> mxH_output;

  // =======================================================================
  // Support for verifying correct energy terms
  // =======================================================================
  static const String SupportedEnergyNameList[];
  OC_BOOL IsSupportedEnergy(const Oxs_Energy*) const;

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  virtual OC_BOOL Init();
  YY_LLBEulerEvolve(const char* name,     // Child instance id
     Oxs_Director* newdtr, // App director
     const char* argstr);  // MIF input block parameters
  virtual ~YY_LLBEulerEvolve();

  virtual  OC_BOOL
  Step(const Oxs_TimeDriver* driver,
       Oxs_ConstKey<Oxs_SimState> current_state,
       const Oxs_DriverStepInfo& step_info,
       Oxs_Key<Oxs_SimState>& next_state);
  // Returns true if step was successful, false if
  // unable to step as requested.
};

const String YY_LLBEulerEvolve::SupportedEnergyNameList[] = {
  "Oxs_UZeeman",
  "Oxs_FixedZeeman",
  "Oxs_ScriptUZeeman",
  "Oxs_StageZeeman",
  "Oxs_TransformZeeman",
  "Oxs_Demag",
  "YY_LLBExchange6Ngbr",
  "YY_LLBUniaxialAnisotropy"
};

/**
 * Notes on timestep control and error criteria
 * Same as in thetaevolve.h, mostly same as in eulerevolve.h except for
 * energy_accum_count and energy_accum_count_limit.
 *
 * At this point, adaptive stepsize is not implemented for T != 0K.
 *
 * Error-based step size control parameters. Each may be disabled
 * by setting to -1.  There is an additional step size control that
 * insures that energy is monotonically non-increasing (up to
 * estimated rounding error).
 *
 * OC_REAL8m allowed_error_rate;  
 * Step size is adjusted so that the estimated maximum error (across all 
 * spins) divided by the step size is smaller than this value.  The units
 * internally are radians per second, converted from the value specified in
 * the input MIF file, which is in deg/sec.
 * 
 * OC_REAL8m allowed_absolute_step_error; 
 * Similar to allowed_error_rate, but without the step size adjustment.
 * Internal units are radians; MIF input units are degrees.
 * 
 * OC_REAL8m allowed_relative_step_error; 
 * Step size is adjusted so that the estimated maximum error (across all 
 * spins) divided by [maximum dm/dt (across all spins) * step size] is 
 * smaller than this value.  This value is non-dimensional, representing the
 * allowed relative (proportional) error, presumably in (0,1).
 * 
 * OC_REAL8m step_headroom; 
 * The 3 control parameters above can be used to estimate the step size that
 * would just fit the control requirements.  Because this is only an 
 * estimate, if the step size is actually set to that value there is a good 
 * chance that the requirement will not be met.  So instead, we leave some 
 * headroom by setting the step size to the computed value multiplied by
 * step_headroom. This is a non-dimensional quantity, which should be in the
 * range (0,1).
 * 
 * OC_REAL8m start_dm;
 * The next timestep is based on the error from the last step.  If there 
 * is no last step (either because this is the first step, or because the 
 * last state handled by this routine is different from the incoming 
 * current_state), then timestep is calculated so that 
 * max_dm_dt * timestep = start_dm.
 *
 * const OC_UINT4m energy_accum_count_limit ;
 * OC_UINT4m energy_accum_count;
 * The total energy field in Oxs_SimState is computed by accumulating
 * the dE into the total energy from the previous state.  It order to
 * keep this value from becoming too inaccurate, we recalculate total
 * energy directly from the energy densities after each
 * energy_accum_count_limit passes.  NB: The code here assumes that a
 * single sequence of states is being fed to this routine.  If this
 * is not the case, then the accum count needs to be tied to the state
 * id.
 */


#undef DEFAULT_M_E_TOL
#endif // _YY_LLBEULEREVOLVE
