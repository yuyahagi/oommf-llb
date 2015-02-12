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

class YY_LLBEulerEvolve:public Oxs_TimeEvolver {
private:
  mutable OC_UINT4m mesh_id;     // Used by gamma and alpha meshvalues to
  void UpdateMeshArrays(const Oxs_Mesh* mesh);

  // Base step size control parameters
  OC_REAL8m min_timestep;   // Seconds
  OC_REAL8m max_timestep;   // Seconds
  OC_REAL8m fixed_timestep; // Seconds -> min_timestep = max_timestep

  // Error-based step size control parameters.  Each may be disabled
  // by setting to -1.  There is an additional step size control that
  // insures that energy is monotonically non-increasing (up to
  // estimated rounding error).
  OC_REAL8m allowed_error_rate;  /// Step size is adjusted so
  /// that the estimated maximum error (across all spins) divided
  /// by the step size is smaller than this value.  The units
  /// internally are radians per second, converted from the value
  /// specified in the input MIF file, which is in deg/sec.

  OC_REAL8m allowed_absolute_step_error; /// Similar to allowed_error_rate,
  /// but without the step size adjustment.  Internal units are
  /// radians; MIF input units are degrees.

  OC_REAL8m allowed_relative_step_error; /// Step size is adjusted so that
  /// the estimated maximum error (across all spins) divided by
  /// [maximum dm/dt (across all spins) * step size] is smaller than
  /// this value.  This value is non-dimensional, representing the
  /// allowed relative (proportional) error, presumably in (0,1).

  OC_REAL8m step_headroom; /// The 3 control parameters above can be
  /// used to estimate the step size that would just fit the control
  /// requirements.  Because this is only an estimate, if the step size
  /// is actually set to that value there is a good chance that the
  /// requirement will not be met.  So instead, we leave some headroom
  /// by setting the step size to the computed value multiplied by
  /// step_headroom.  This is a non-dimensional quantity, which should
  /// be in the range (0,1).

  // The total energy field in Oxs_SimState is computed by accumulating
  // the dE into the total energy from the previous state.  It order to
  // keep this value from becoming too inaccurate, we recalculate total
  // energy directly from the energy densities after each
  // energy_accum_count_limit passes.  NB: The code here assumes that a
  // single sequence of states is being fed to this routine.  If this
  // is not the case, then the accum count needs to be tied to the state
  // id.
  const OC_UINT4m energy_accum_count_limit ;
  OC_UINT4m energy_accum_count;

  // Spatially variable Landau-Lifshitz-Gilbert damping coef
  Oxs_OwnedPointer<Oxs_ScalarField> gamma_init;
  mutable Oxs_MeshValue<OC_REAL8m> gamma;       // LL gyromagnetic ratio
  Oxs_OwnedPointer<Oxs_ScalarField> alpha_t_init;
  mutable Oxs_MeshValue<OC_REAL8m> alpha_t, alpha_t0;       // LL damping coef, transverse
  mutable Oxs_MeshValue<OC_REAL8m> alpha_l;       // LL damping coef, longitudinal
  // alpha_t0 is the value at T = 0 K.
  OC_BOOL do_precess;  // If false, then do pure damping
  OC_BOOL allow_signed_gamma; // If false, then force gamma negative
  enum GammaStyle { GS_INVALID, GS_LL, GS_G }; // Landau-Lifshitz or Gilbert
  GammaStyle gamma_style;

  // The next timestep is based on the error from the last step.  If
  // there is no last step (either because this is the first step,
  // or because the last state handled by this routine is different
  // from the incoming current_state), then timestep is calculated
  // so that max_dm_dt * timestep = start_dm.
  OC_REAL8m start_dm;

  // Data cached from last state
  OC_UINT4m energy_state_id;
  Oxs_MeshValue<OC_REAL8m> energy;
  OC_REAL8m next_timestep;
  
  
  /**
  Variables and constants for the temperature dependant part of this evolver
  **/
  const OC_REAL8m KBoltzmann;            // Boltzmann constant
  OC_REAL8m kB_T;                        // thermal energy
  Oxs_MeshValue<OC_REAL8m> hFluctVarConst_t;  // constant part of the variance of the thermal field
  Oxs_MeshValue<OC_REAL8m> hFluctVarConst_l;
  Oxs_MeshValue<ThreeVector> hFluct_t; // current values of the thermal field
                                      // due to the fact, that dm_dt is sometimes calculated more
                                      // than one time per iteration these values are to be stored in an array
  Oxs_MeshValue<ThreeVector> hFluct_l; // longitudinal
  void FillHFluctConst(const Oxs_Mesh* mesh);
  void InitHFluct(const Oxs_Mesh* mesh);

  Oxs_OwnedPointer<Oxs_ScalarField> Tc_init;
  Oxs_MeshValue<OC_REAL8m> Tc;  // Currie temperature in Kelvin

  // constant part of the additional drift term that arises in stochastic caculus
  Oxs_MeshValue<OC_REAL8m> inducedDriftConst_t;
  Oxs_MeshValue<OC_REAL8m> inducedDriftConst_l;
  OC_REAL8m temperature;            // in Kelvin
  OC_UINT4m iteration_Tcalculated;  // keep in mind for which iteration the thermal field is already calculated
  OC_BOOL   ito_calculus;
  // use alternative interpretation of the stochastic Langevin equation
  // (in this case the drift term is omitted)                                    

  /**
  Support for stage-varying temperature
  **/
  void SetTemperature(const Oxs_Mesh* mesh_, OC_REAL8m newtemp);
  OC_BOOL has_tempscript;
  vector<Nb_TclCommandLineOption> tempscript_opts;
  Nb_TclCommand tempscript_cmd;
  OC_REAL8m GetStageTemp(OC_UINT4m stage) const;

  /**
  Variables for the Random distributions
  **/
  // seed to initialize the generator with, can be any integer
  // (beware: -|n| is used in this method)
  OC_INT4m uniform_seed;  
  OC_BOOL has_uniform_seed;

  // the here used gaussian distribution algorithm 
  // computes two values at a time, but only one is returned.
  // Therefor the second value (plus a flag) must be stored outside the method itself,
  // to be returned when the method is called for the second time
  OC_BOOL gaus2_isset;    
  OC_REAL8m gaus2;        

  /**
  Random functions (needed for temperature effects)
  **/
  OC_REAL8m Gaussian_Random (const OC_REAL8m muGaus, const OC_REAL8m  sigmaGaus);
  //returns a gaussian distributed random number 
  //with average muGaus und standard deviation sigmaGaus

  // Outputs
  void UpdateDerivedOutputs(const Oxs_SimState&);
  Oxs_ScalarOutput<YY_LLBEulerEvolve> max_dm_dt_output;
  Oxs_ScalarOutput<YY_LLBEulerEvolve> dE_dt_output;
  Oxs_ScalarOutput<YY_LLBEulerEvolve> delta_E_output;
  Oxs_VectorFieldOutput<YY_LLBEulerEvolve> dm_dt_t_output;
  Oxs_VectorFieldOutput<YY_LLBEulerEvolve> dm_dt_l_output;
  Oxs_VectorFieldOutput<YY_LLBEulerEvolve> mxH_output;

  // Scratch space
  Oxs_MeshValue<OC_REAL8m> new_energy;
  Oxs_MeshValue<ThreeVector> new_dm_dt_t;
  Oxs_MeshValue<ThreeVector> new_dm_dt_l;

  Oxs_MeshValue<ThreeVector> total_field;
  Oxs_MeshValue<OC_REAL8m> Ms0; // Ms at T = 0K.
  OC_BOOL isMs0Set;

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

#endif // _YY_LLBEULEREVOLVE
