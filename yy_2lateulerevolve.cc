/** FILE: yy_2lateulerevolve.cc                 -*-Mode: c++-*-
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

#include <math.h>

#include "nb.h"
#include "director.h"
#include "simstate.h"
#include "key.h"
#include "energy.h"    // Needed to make MSVC++ 5 happy
#include "meshvalue.h"
#include "rectangularmesh.h"
#include "scalarfield.h"

#include "yy_2lattimedriver.h"
#include "yy_2lateulerevolve.h"

// Oxs_Ext registration support
OXS_EXT_REGISTER(YY_2LatEulerEvolve);

/* End includes */

void YY_2LatEulerEvolve::UpdateStageTemperature(const Oxs_SimState& state)
{
  if(!has_tempscript) return;

  const Oxs_Mesh* mesh = state.mesh;
  const OC_REAL8m stage = state.stage_number;
  const OC_REAL8m size = mesh->Size();
  OC_INDEX index;
  if((index = tempscript_opts[0].position)>=0) { // stage
    tempscript_cmd.SetCommandArg(index,stage);
  }
  if((index = tempscript_opts[1].position)>=0) { // stage_time
    tempscript_cmd.SetCommandArg(index,state.stage_elapsed_time);
  }
  if((index = tempscript_opts[2].position)>=0) { // total_time
    tempscript_cmd.SetCommandArg(index,state.stage_start_time+state.stage_elapsed_time);
  }

  vector<String> params;
  tempscript_cmd.SaveInterpResult();
  tempscript_cmd.Eval();
  tempscript_cmd.GetResultList(params);
  tempscript_cmd.RestoreInterpResult();

  OXS_GET_EXT_OBJECT(params,Oxs_ScalarField,temperature_init);
  temperature_init->FillMeshValue(mesh,temperature);
  kB_T.AdjustSize(mesh);
  for(OC_INDEX i=0; i<size; i++) {
    kB_T[i] = KBoltzmann*temperature[i];
  }
}

// Constructor
YY_2LatEulerEvolve::YY_2LatEulerEvolve(
    const char* name,     // Child instance id
    Oxs_Director* newdtr, // App director
    const char* argstr)   // MIF input block parameters
    : YY_2LatTimeEvolver(name,newdtr,argstr),
    mesh_id(0), min_timestep(0.), max_timestep(1e-10),
    energy_accum_count_limit(25),
    energy_state_id(0),next_timestep(0.),
    KBoltzmann(1.38062e-23),
    iteration_hFluct1_calculated(0),
    iteration_hFluct2_calculated(0),
    has_tempscript(0),
    last_stage_number(0)
{
  // Process arguments
  // For now, it works with a fixed time step but there still are min_ and
  // max_timestep for future implementation of adaptive stepsize.
  fixed_timestep = GetRealInitValue("fixed_timestep",1e-16);
  min_timestep = max_timestep = fixed_timestep;
  if(max_timestep<=0.0) {
    char buf[4096];
    Oc_Snprintf(buf,sizeof(buf),
    "Invalid parameter value:"
    " Specified max time step is %g (should be >0.)",
    max_timestep);
    throw Oxs_Ext::Error(this,buf);
  }

  allowed_error_rate = GetRealInitValue("error_rate",-1);
  if(allowed_error_rate>0.0) {
    allowed_error_rate *= PI*1e9/180.; // Convert from deg/ns to rad/s
  }
  allowed_absolute_step_error
    = GetRealInitValue("absolute_step_error",0.2);
  if(allowed_absolute_step_error>0.0) {
    allowed_absolute_step_error *= PI/180.; // Convert from deg to rad
  }
  allowed_relative_step_error
    = GetRealInitValue("relative_step_error",0.2);

  step_headroom = GetRealInitValue("step_headroom",0.85);
  if(step_headroom<=0.) {
    throw Oxs_Ext::Error(this,"Invalid initialization detected:"
       " step_headroom value must be bigger than 0.");
  }

  if(HasInitValue("alpha_t1")) {
    OXS_GET_INIT_EXT_OBJECT("alpha_t1",Oxs_ScalarField,alpha_t1_init);
  } else {
    alpha_t1_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
                          (MakeNew("Oxs_UniformScalarField",director,
                                   "value 0.5")));
  }

  if(HasInitValue("alpha_t2")) {
    OXS_GET_INIT_EXT_OBJECT("alpha_t2",Oxs_ScalarField,alpha_t2_init);
  } else {
    alpha_t2_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
                          (MakeNew("Oxs_UniformScalarField",director,
                                   "value 0.5")));
  }

  // Flag to include stochastic field
  use_stochastic = GetRealInitValue("use_stochastic",0);

  // User may specify either gamma_G (Gilbert) or
  // gamma_LL (Landau-Lifshitz).  Code uses "gamma"
  // which is LL form.
  gamma1_style = GS_INVALID;
  if(HasInitValue("gamma_G1") && HasInitValue("gamma_LL1")) {
    throw Oxs_Ext::Error(this,"Invalid Specify block; "
       "both gamma_G1 and gamma_LL1 specified.");
  } else if(HasInitValue("gamma_G1")) {
    OXS_GET_INIT_EXT_OBJECT("gamma_G1",Oxs_ScalarField,gamma1_init);
    gamma1_style = GS_G;
  } else if(HasInitValue("gamma_LL1")) {
    OXS_GET_INIT_EXT_OBJECT("gamma_LL1",Oxs_ScalarField,gamma1_init);
    gamma1_style = GS_LL;
  } else {
    gamma1_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
                          (MakeNew("Oxs_UniformScalarField",director,
                                   "value 2.211e5")));
  }

  gamma2_style = GS_INVALID;
  if(HasInitValue("gamma_G2") && HasInitValue("gamma_LL2")) {
    throw Oxs_Ext::Error(this,"Invalid Specify block; "
       "both gamma_G2 and gamma_LL2 specified.");
  } else if(HasInitValue("gamma_G2")) {
    OXS_GET_INIT_EXT_OBJECT("gamma_G2",Oxs_ScalarField,gamma2_init);
    gamma2_style = GS_G;
  } else if(HasInitValue("gamma_LL2")) {
    OXS_GET_INIT_EXT_OBJECT("gamma_LL2",Oxs_ScalarField,gamma2_init);
    gamma2_style = GS_LL;
  } else {
    gamma2_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
                          (MakeNew("Oxs_UniformScalarField",director,
                                   "value 2.211e5")));
  }

  do_precess = GetIntInitValue("do_precess",1);

  start_dm = GetRealInitValue("start_dm",0.01);
  start_dm *= PI/180.; // Convert from deg to rad

  // Get time dependent multiplier to scale temperature
  if(HasInitValue("tempscript")) {
    has_tempscript=1;
    String cmdoptreq = GetStringInitValue("tempscript_args",
                                          "stage stage_time total_time");
    tempscript_opts.push_back(Nb_TclCommandLineOption("stage",1));
    tempscript_opts.push_back(Nb_TclCommandLineOption("stage_time",1));
    tempscript_opts.push_back(Nb_TclCommandLineOption("total_time",1));
    tempscript_cmd.SetBaseCommand(InstanceName(),
				  director->GetMifInterp(),
				  GetStringInitValue("tempscript"),
				  Nb_ParseTclCommandLineRequest(InstanceName(),
								 tempscript_opts,
								 cmdoptreq));
  } else {
    temperature_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
                          (MakeNew("Oxs_UniformScalarField",director,
                                   "value 0.0")));
  }

  // set temperature to zero to get an estimate for a reasonable stepsize
  // or use it for comparison (acts like eulerevolve with temperature=0K)
  if(!has_tempscript){ // That is, T = 0.
    min_timestep = 0.;    
    max_timestep = 1e-10; 
  }

  if(HasInitValue("uniform_seed")) {
    uniform_seed = GetIntInitValue("uniform_seed");
    has_uniform_seed = 1;
  } else {
    has_uniform_seed = 0;
  }

  gaus2_isset = 0;    //no gaussian random numbers calculated yet

  // Setup outputs
  max_dm_dt_output.Setup(this,InstanceName(),"Max dm/dt","deg/ns",0,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);
  dE_dt_output.Setup(this,InstanceName(),"dE/dt","J/s",0,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);
  delta_E_output.Setup(this,InstanceName(),"Delta E","J",0,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);
  dm_dt_t1_output.Setup(this,InstanceName(),"dm/dt (trans.)1","rad/s",1,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);
  dm_dt_l1_output.Setup(this,InstanceName(),"dm/dt (long.)1","rad/s",1,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);
  mxH1_output.Setup(this,InstanceName(),"mxH1","A/m",1,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);
  dm_dt_t2_output.Setup(this,InstanceName(),"dm/dt (trans.)2","rad/s",1,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);
  dm_dt_l2_output.Setup(this,InstanceName(),"dm/dt (long.)2","rad/s",1,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);
  mxH2_output.Setup(this,InstanceName(),"mxH2","A/m",1,
     &YY_2LatEulerEvolve::UpdateDerivedOutputs);

  VerifyAllInitArgsUsed();
}   // end Constructor

OC_BOOL YY_2LatEulerEvolve::Init()
{
  // Register outputs
  max_dm_dt_output.Register(director,-5);
  dE_dt_output.Register(director,-5);
  delta_E_output.Register(director,-5);
  dm_dt_t1_output.Register(director,-5);
  dm_dt_l1_output.Register(director,-5);
  mxH1_output.Register(director,-5);
  dm_dt_t2_output.Register(director,-5);
  dm_dt_l2_output.Register(director,-5);
  mxH2_output.Register(director,-5);

  // dm_dt and mxH output caches are used for intermediate storage,
  // so enable caching.
  dm_dt_t1_output.CacheRequestIncrement(1);
  dm_dt_l1_output.CacheRequestIncrement(1);
  mxH1_output.CacheRequestIncrement(1);
  dm_dt_t2_output.CacheRequestIncrement(1);
  dm_dt_l2_output.CacheRequestIncrement(1);
  mxH2_output.CacheRequestIncrement(1);

  alpha_t10.Release(); alpha_t1.Release(); alpha_l1.Release();
  alpha_t20.Release(); alpha_t2.Release(); alpha_l2.Release();
  gamma1.Release(); gamma2.Release();
  energy.Release();
  total_field1.Release();
  total_field2.Release();
  new_energy.Release();
  new_dm_dt_t1.Release();
  new_dm_dt_l1.Release();
  new_dm_dt_t2.Release();
  new_dm_dt_l2.Release();

  hFluct_t1.Release(); hFluct_l1.Release();
  hFluct_t2.Release(); hFluct_l2.Release();
  hFluctVarConst_t1.Release(); hFluctVarConst_l1.Release();
  hFluctVarConst_t2.Release(); hFluctVarConst_l2.Release();

  energy_state_id=0;   // Mark as invalid state
  next_timestep=0.;    // Dummy value
  energy_accum_count=energy_accum_count_limit; // Force cold count
  // on first pass

  // (Re)initialize random number generator
  if(has_uniform_seed) {
    Oc_Srand(uniform_seed); //initialize Random number generator
  } else {
    // Default seed value is time dependent
    Oc_Srand();
  }

  return YY_2LatTimeEvolver::Init();  // Initialize parent class.
  // Do this after child output registration so that
  // UpdateDerivedOutputs gets called before the parent
  // total_energy_output update function.
}

YY_2LatEulerEvolve::~YY_2LatEulerEvolve()
{}

void YY_2LatEulerEvolve::Calculate_dm_dt(
    const Oxs_SimState& state_,
    const Oxs_MeshValue<ThreeVector>& mxH_,
    const Oxs_MeshValue<ThreeVector>& total_field_,
    OC_REAL8m pE_pt_,
    Oxs_MeshValue<ThreeVector>& dm_dt_t_,
    Oxs_MeshValue<ThreeVector>& dm_dt_l_,
    OC_REAL8m& max_dm_dt_,
    OC_REAL8m& dE_dt_,
    OC_REAL8m& min_timestep_)
{
  // Imports: state_, mxH_, pE_pt
  // Exports: dm_dt_t_, dm_dt_l_, max_dm_dt_, dE_dt_
  const Oxs_Mesh* mesh_ = state_.mesh;
  const OC_INDEX size = mesh_->Size(); // Assume all imports are compatible
  const Oxs_MeshValue<OC_REAL8m>& Ms_ = *(state_.Ms);
  const Oxs_MeshValue<OC_REAL8m>& Ms_inverse_ = *(state_.Ms_inverse);
  const Oxs_MeshValue<OC_REAL8m>& Ms0_ = *(state_.Ms0);
  const Oxs_MeshValue<OC_REAL8m>& Ms0_inverse_ = *(state_.Ms0_inverse);
  const Oxs_MeshValue<ThreeVector>& spin_ = state_.spin;
  OC_UINT4m iteration_now = state_.iteration_count;
  ThreeVector scratch_t;
  ThreeVector scratch_l;
  ThreeVector dm_fluct_t;
  ThreeVector dm_fluct_l;
  OC_REAL8m hFluctSigma_t;
  OC_REAL8m hFluctSigma_l;
  dm_dt_t_.AdjustSize(mesh_);
  dm_dt_l_.AdjustSize(mesh_);
  OC_INDEX i;

  if(state_.lattice_type==Oxs_SimState::TOTAL) {
    throw Oxs_ExtError(this, "PROGRAMMING ERROR: YY_2LatEulerEvolve::"
        "Calculate_dm_dt() is called with a wrong type of simulation"
        " state.");
  }

  iteration_now++;
  // if not done, hFluct for first step may be calculated too often
  
  if(mesh_id != mesh_->Id() || !gamma1.CheckMesh(mesh_)) {
    // First go or mesh change detected
    alpha_t1_init->FillMeshValue(mesh_,alpha_t10);
    alpha_t2_init->FillMeshValue(mesh_,alpha_t20);
    gamma1_init->FillMeshValue(mesh_,gamma1);
    gamma2_init->FillMeshValue(mesh_,gamma2);
    if(!allow_signed_gamma) {
      for(i=0;i<size;++i) {
        gamma1[i] = fabs(gamma1[i]);
        gamma2[i] = fabs(gamma2[i]);
      }
    }
    if(gamma1_style == GS_G) { // Convert to LL form
      for(i=0;i<size;++i) {
        gamma1[i] /= (1+alpha_t10[i]*alpha_t10[i]);
      }
    }

    if(gamma1_style == GS_G) { // Convert to LL form
      for(i=0;i<size;++i) {
        gamma2[i] /= (1+alpha_t20[i]*alpha_t20[i]);
      }
    }

    // Prepare temperature-dependent mesh value arrays
    alpha_t1.AdjustSize(mesh_);
    alpha_t2.AdjustSize(mesh_);
    alpha_l1.AdjustSize(mesh_);
    alpha_l2.AdjustSize(mesh_);
    hFluctVarConst_t1.AdjustSize(mesh_);
    hFluctVarConst_t2.AdjustSize(mesh_);
    hFluctVarConst_l1.AdjustSize(mesh_);
    hFluctVarConst_l2.AdjustSize(mesh_);

    // Other mesh value arrays
    total_field1.AdjustSize(mesh_);
    total_field2.AdjustSize(mesh_);
    hFluct_t1.AdjustSize(mesh_);     
    hFluct_t2.AdjustSize(mesh_);     
    hFluct_l1.AdjustSize(mesh_);     
    hFluct_l2.AdjustSize(mesh_);     
    iteration_hFluct1_calculated = 0;     
    iteration_hFluct2_calculated = 0;     

    // Update stage-dependent temperature and temperature-dependent parameters
    UpdateStageTemperature(*(state_.total_lattice));
    UpdateMeshArrays(*(state_.total_lattice));

    // Set pointers for the temperature-dependent parameters in states.
    switch(state_.lattice_type) {
    case Oxs_SimState::LATTICE1:
      state_.T = &temperature;
      state_.lattice2->T = &temperature;
      break;
    case Oxs_SimState::LATTICE2:
      state_.lattice1->T = &temperature;
      state_.T = &temperature;
      break;
    default:
      // Program should not reach here.
      break;
    }
  }

  // Judge the type of lattice (sublattice1 or 2) and use corresponding
  // parameters for calculation
  const Oxs_MeshValue<OC_REAL8m>* Tc;
  const Oxs_MeshValue<OC_REAL8m>* alpha_t;
  const Oxs_MeshValue<OC_REAL8m>* alpha_l;
  const Oxs_MeshValue<OC_REAL8m>* gamma;
  const Oxs_MeshValue<OC_REAL8m>* m_e;
  const Oxs_MeshValue<OC_REAL8m>* chi_l;
  Oxs_MeshValue<OC_REAL8m>* hFluctVarConst_t;
  Oxs_MeshValue<OC_REAL8m>* hFluctVarConst_l;
  Oxs_MeshValue<ThreeVector>* hFluct_t;
  Oxs_MeshValue<ThreeVector>* hFluct_l;
  OC_UINT4m* iteration_hFluct_calculated;

  switch(state_.lattice_type) {
  case Oxs_SimState::LATTICE1:
    Tc = state_.Tc;
    alpha_t = &alpha_t1;
    alpha_l = &alpha_l1;
    gamma = &gamma1;
    m_e = state_.m_e;
    chi_l = state_.chi_l;
    hFluctVarConst_t = &hFluctVarConst_t1;
    hFluctVarConst_l = &hFluctVarConst_l1;
    hFluct_t = &hFluct_t1;
    hFluct_l = &hFluct_l1;
    iteration_hFluct_calculated = &iteration_hFluct1_calculated;
    break;
  case Oxs_SimState::LATTICE2:
    Tc = state_.Tc;
    alpha_t = &alpha_t2;
    alpha_l = &alpha_l2;
    gamma = &gamma2;
    m_e = state_.m_e;
    chi_l = state_.chi_l;
    hFluctVarConst_t = &hFluctVarConst_t2;
    hFluctVarConst_l = &hFluctVarConst_l2;
    hFluct_t = &hFluct_t2;
    hFluct_l = &hFluct_l2;
    iteration_hFluct_calculated = &iteration_hFluct2_calculated;
    break;
  default:
    // Program should not reach here.
    break;
  }

  if (use_stochastic && iteration_now > *iteration_hFluct_calculated) {
    // i.e. if thermal field is not calculated for this step
    for(i=0;i<size;i++){
      if(Ms_[i] != 0){
        // Only sqrt(delta_t) is multiplied for stochastic functions
        // opposed to dm_dt * delta_t for deterministic functions.
        // This is the standard deviation of the gaussian distribution
        // used to represent the thermal perturbations
        hFluctSigma_t = sqrt((*hFluctVarConst_t)[i] / fixed_timestep);
        hFluctSigma_l = sqrt((*hFluctVarConst_l)[i] / fixed_timestep);

        (*hFluct_t)[i].x = hFluctSigma_t*Gaussian_Random(0.0, 1.0);
        (*hFluct_t)[i].y = hFluctSigma_t*Gaussian_Random(0.0, 1.0);
        (*hFluct_t)[i].z = hFluctSigma_t*Gaussian_Random(0.0, 1.0);
        (*hFluct_l)[i].x = hFluctSigma_l*Gaussian_Random(0.0, 1.0);
        (*hFluct_l)[i].y = hFluctSigma_l*Gaussian_Random(0.0, 1.0);
        (*hFluct_l)[i].z = hFluctSigma_l*Gaussian_Random(0.0, 1.0);
      }
    }
  }

  for(i=0;i<size;i++) {
    if(Ms_[i]==0) {
      dm_dt_t_[i].Set(0.0,0.0,0.0);
      dm_dt_l_[i].Set(0.0,0.0,0.0);
    } else {
      OC_REAL8m cell_alpha_t = (*alpha_t)[i];
      OC_REAL8m cell_alpha_l = (*alpha_l)[i];
      OC_REAL8m cell_gamma = (*gamma)[i];
      OC_REAL8m cell_m_inverse = Ms0_[i]*Ms_inverse_[i];

      // deterministic part
      scratch_t = mxH_[i];
      scratch_t *= -cell_gamma; // -|gamma|*(mxH)

      if(do_precess) {
        dm_dt_t_[i]  = scratch_t;
        dm_dt_l_[i].Set(0.0,0.0,0.0);
      } else {
        dm_dt_t_[i].Set(0.0,0.0,0.0);
        dm_dt_l_[i].Set(0.0,0.0,0.0);
      }

      // Transverse damping term
      if(use_stochastic) {
        // Note: The stochastic field is NOT included in the first term of 
        // the LLB equation. See PRB 85, 014433 (2012). The second form of 
        // LLB is the above article is implemented here.
        dm_fluct_t = spin_[i] ^ (*hFluct_t)[i];  // cross product mxhFluct_t
        dm_fluct_t *= -cell_gamma;
        scratch_t += dm_fluct_t;  // -|gamma|*mx(H+hFluct_t)
      }
      scratch_t ^= spin_[i];
      // -|gamma|((mx(H+hFluct_t))xm) = |gamma|(mx(mx(H+hFluct_t)))
      scratch_t *= -cell_alpha_t*cell_m_inverse; // -|alpha*gamma|(mx(mx(H+hFluct_t)))
      dm_dt_t_[i] += scratch_t;

      // Longitudinal terms
      OC_REAL8m temp = spin_[i]*total_field_[i];
      temp *= -cell_gamma*cell_alpha_l;
      temp *= cell_m_inverse;
      scratch_l = temp*spin_[i];
      dm_dt_l_[i] += scratch_l;

      // Check for overshooting
      scratch_l = dm_dt_l_[i]*fixed_timestep;
      scratch_l += spin_[i];
      if( scratch_l*spin_[i]<0.0 ) {
        dm_dt_l_[i] = -1*spin_[i];
        dm_dt_l_[i].x /= fixed_timestep;
        dm_dt_l_[i].y /= fixed_timestep;
        dm_dt_l_[i].z /= fixed_timestep;
      }

      if(temperature[i] != 0 && use_stochastic) {
        // Longitudinal stochastic field parallel to spin
        dm_dt_l_[i] += (*hFluct_l)[i].x*spin_[i]*cell_m_inverse;
      }
    }
  }

  // now hFluct_t is definetely calculated for this iteration
  *iteration_hFluct_calculated = iteration_now;

  // Zero dm_dt at fixed spin sites
  UpdateFixedSpinList(mesh_);
  const OC_INDEX fixed_count = GetFixedSpinCount();
  for(OC_INDEX j=0;j<fixed_count;j++) {
    dm_dt_t_[GetFixedSpin(j)].Set(0.,0.,0.);
  }

  // Collect statistics
  OC_REAL8m max_dm_dt_sq=0.0;
  OC_REAL8m dE_dt_sum=0.0;
  OC_INDEX max_index=0;
  for(i=0;i<size;i++) {
    ThreeVector tempvec = dm_dt_t_[i];
    tempvec += dm_dt_l_[i];
    OC_REAL8m dm_dt_sq = tempvec.MagSq();
    if(dm_dt_sq>0.0) {
      dE_dt_sum += -1*MU0*fabs((*gamma)[i]*(*alpha_t)[i])
        *mxH_[i].MagSq() * Ms_[i] * mesh_->Volume(i);
      if(dm_dt_sq>max_dm_dt_sq) {
        max_dm_dt_sq=dm_dt_sq;
        max_index = i;
      }
    }
  }

  max_dm_dt_ = sqrt(max_dm_dt_sq);
  dE_dt_ = dE_dt_sum; // Transverse terms
  dE_dt_ += pE_pt_;
  // TODO: What about the longitudinal terms?
  /// The first term is (partial E/partial M)*dM/dt, the
  /// second term is (partial E/partial t)*dt/dt.  Note that,
  /// provided Ms_[i]>=0, that by constructions dE_dt_sum above
  /// is always non-negative, so dE_dt_ can only be made positive
  /// by positive pE_pt_.

  if(!has_tempscript) { // temperature == 0 at all cells
    // Get bound on smallest stepsize that would actually
    // change spin new_max_dm_dt_index:
    OC_REAL8m min_ratio = DBL_MAX/2.;
    if(fabs(dm_dt_t_[max_index].x)>=1.0 ||
       min_ratio*fabs(dm_dt_t_[max_index].x) > fabs(spin_[max_index].x)) {
     min_ratio = fabs(spin_[max_index].x/dm_dt_t_[max_index].x);
    }
   if(fabs(dm_dt_t_[max_index].y)>=1.0 ||
      min_ratio*fabs(dm_dt_t_[max_index].y) > fabs(spin_[max_index].y)) {
     min_ratio = fabs(spin_[max_index].y/dm_dt_t_[max_index].y);
   }
   if(fabs(dm_dt_t_[max_index].z)>=1.0 ||
      min_ratio*fabs(dm_dt_t_[max_index].z) > fabs(spin_[max_index].z)) {
      min_ratio = fabs(spin_[max_index].z/dm_dt_t_[max_index].z);
    }
  min_timestep_ = min_ratio * OC_REAL8_EPSILON;
  }
  else {min_timestep_ = fixed_timestep;}

  return;
} // end Calculate_dm_dt

OC_BOOL
YY_2LatEulerEvolve::Step(const YY_2LatTimeDriver* driver,
          Oxs_ConstKey<Oxs_SimState> current_state,
          Oxs_ConstKey<Oxs_SimState> current_state1,
          Oxs_ConstKey<Oxs_SimState> current_state2,
          const Oxs_DriverStepInfo& /* step_info */,
          Oxs_Key<Oxs_SimState>& next_state,
          Oxs_Key<Oxs_SimState>& next_state1,
          Oxs_Key<Oxs_SimState>& next_state2)
{
  const OC_REAL8m max_step_increase = 1.25;
  const OC_REAL8m max_step_decrease = 0.5;

  OC_INDEX size,i; // Mesh size and indexing variable

  const Oxs_SimState& cstate = current_state.GetReadReference();
  const Oxs_SimState& cstate1 = current_state1.GetReadReference();
  const Oxs_SimState& cstate2 = current_state2.GetReadReference();
  Oxs_SimState& workstate = next_state.GetWriteReference();
  Oxs_SimState& workstate1 = next_state1.GetWriteReference();
  Oxs_SimState& workstate2 = next_state2.GetWriteReference();
  driver->FillState(cstate,workstate);
  driver->FillState(cstate1,workstate1);
  driver->FillState(cstate2,workstate2);

  // Set pointers to the sublattice
  workstate.lattice1 = &workstate1;
  workstate.lattice2 = &workstate2;
  workstate1.total_lattice = &workstate;
  workstate1.lattice2 = &workstate2;
  workstate2.total_lattice = &workstate;
  workstate2.lattice1 = &workstate1;
  workstate1.lattice_type = Oxs_SimState::LATTICE1;
  workstate2.lattice_type = Oxs_SimState::LATTICE2;

  if(cstate.mesh->Id() != workstate.mesh->Id()) {
    throw Oxs_Ext::Error(this,
        "YY_2LatEulerEvolve::Step: Oxs_Mesh not fixed across steps.");
  }

  if(cstate.Id() != workstate.previous_state_id) {
    throw Oxs_Ext::Error(this,
        "YY_2LatEulerEvolve::Step: State continuity break detected.");
  }

  // Pull cached values out from cstate.
  // If cstate.Id() == energy_state_id, then cstate has been run
  // through either this method or UpdateDerivedOutputs.  Either
  // way, all derived state data should be stored in cstate,
  // except currently the "energy" mesh value array, which is
  // stored independently inside *this.  Eventually that should
  // probably be moved in some fashion into cstate too.
  if(energy_state_id != cstate.Id()) {
    // cached data out-of-date
    UpdateDerivedOutputs(cstate);
  }
  OC_BOOL cache_good = 1;
  OC_REAL8m max_dm_dt;
  OC_REAL8m dE_dt, delta_E, pE_pt;
  OC_REAL8m timestep_lower_bound;  // Smallest timestep that can actually
  /// change spin with max_dm_dt (due to OC_REAL8_EPSILON restrictions).
  /// The next timestep is based on the error from the last step.  If
  /// there is no last step (either because this is the first step,
  /// or because the last state handled by this routine is different
  /// from the incoming current_state), then timestep is calculated
  /// so that max_dm_dt * timestep = start_dm.

  cache_good &= cstate.GetDerivedData("Max dm/dt",max_dm_dt);
  cache_good &= cstate.GetDerivedData("dE/dt",dE_dt);
  cache_good &= cstate.GetDerivedData("Delta E",delta_E);
  cache_good &= cstate.GetDerivedData("pE/pt",pE_pt);
  cache_good &= cstate.GetDerivedData("Timestep lower bound",
              timestep_lower_bound);
  cache_good &= (energy_state_id == cstate.Id());
  cache_good &= (dm_dt_t1_output.cache.state_id == cstate.Id());
  cache_good &= (dm_dt_l1_output.cache.state_id == cstate.Id());
  cache_good &= (dm_dt_t2_output.cache.state_id == cstate.Id());
  cache_good &= (dm_dt_l2_output.cache.state_id == cstate.Id());

  if(!cache_good) {
    throw Oxs_Ext::Error(this,
       "YY_2LatEulerEvolve::Step: Invalid data cache.");
  }

  const Oxs_MeshValue<ThreeVector>& dm_dt_t1 = dm_dt_t1_output.cache.value;
  const Oxs_MeshValue<ThreeVector>& dm_dt_l1 = dm_dt_l1_output.cache.value;
  const Oxs_MeshValue<ThreeVector>& dm_dt_t2 = dm_dt_t2_output.cache.value;
  const Oxs_MeshValue<ThreeVector>& dm_dt_l2 = dm_dt_l2_output.cache.value;

  // Negotiate with driver over size of next step
  OC_REAL8m stepsize = next_timestep;

  if(stepsize<=0.0) {
    if(start_dm < sqrt(DBL_MAX/4) * max_dm_dt) {
      stepsize = start_dm / max_dm_dt;
    } else {
      stepsize = sqrt(DBL_MAX/4);
    }
  }
   OC_BOOL forcestep=0;
  // Insure step is not outside requested step bounds
  if(stepsize<min_timestep) {
    // the step has to be forced here,to make sure we don't produce
    // an infinite loop
    stepsize = min_timestep;
    forcestep = 1;
    }
  if(stepsize>max_timestep) stepsize = max_timestep;

  workstate.last_timestep=stepsize;
  workstate1.last_timestep=stepsize;
  workstate2.last_timestep=stepsize;
  if(stepsize<timestep_lower_bound) {
    workstate.last_timestep=timestep_lower_bound;
    workstate1.last_timestep=timestep_lower_bound;
    workstate2.last_timestep=timestep_lower_bound;
  }

  if(cstate.stage_number != last_stage_number) {
    // New stage
    last_stage_number = cstate.stage_number;
    workstate.stage_start_time = cstate.stage_start_time
                                + cstate.stage_elapsed_time;
    workstate.stage_elapsed_time = workstate.last_timestep;
    workstate1.stage_start_time = cstate.stage_start_time
                                + cstate.stage_elapsed_time;
    workstate1.stage_elapsed_time = workstate1.last_timestep;
    workstate2.stage_start_time = cstate.stage_start_time
                                + cstate.stage_elapsed_time;
    workstate2.stage_elapsed_time = workstate2.last_timestep;

    // Update stage-dependent temperature and temperature-dependent parameters
    UpdateStageTemperature(workstate);
    UpdateMeshArrays(workstate);
  } else {
    workstate.stage_start_time = cstate.stage_start_time;
    workstate.stage_elapsed_time = cstate.stage_elapsed_time
                                  + workstate.last_timestep;
    workstate1.stage_start_time = cstate.stage_start_time;
    workstate1.stage_elapsed_time = cstate.stage_elapsed_time
                                  + workstate1.last_timestep;
    workstate2.stage_start_time = cstate.stage_start_time;
    workstate2.stage_elapsed_time = cstate.stage_elapsed_time
                                  + workstate2.last_timestep;
  }
  workstate.iteration_count = cstate.iteration_count + 1;
  workstate.stage_iteration_count = cstate.stage_iteration_count + 1;
  driver->FillStateSupplemental(workstate);
  workstate1.iteration_count = cstate1.iteration_count + 1;
  workstate1.stage_iteration_count = cstate1.stage_iteration_count + 1;
  driver->FillStateSupplemental(workstate1);
  workstate2.iteration_count = cstate2.iteration_count + 1;
  workstate2.stage_iteration_count = cstate2.stage_iteration_count + 1;
  driver->FillStateSupplemental(workstate2);

  if(workstate.last_timestep>stepsize) {
    // Either driver wants to force this stepsize (in order to end stage 
    // exactly at boundary), or else suggested stepsize is smaller than
    // timestep_lower_bound.
    forcestep=1;
  }
  stepsize = workstate.last_timestep;

  // Put new spin configuration in next_state
  workstate.spin.AdjustSize(workstate.mesh); // Safety
  workstate1.spin.AdjustSize(workstate.mesh);
  workstate2.spin.AdjustSize(workstate.mesh);
  size = workstate.spin.Size();
  Oxs_MeshValue<OC_REAL8m>& wMs = *(workstate.Ms);
  Oxs_MeshValue<OC_REAL8m>& wMs_inverse = *(workstate.Ms_inverse);
  Oxs_MeshValue<OC_REAL8m>& wMs1 = *(workstate1.Ms);
  Oxs_MeshValue<OC_REAL8m>& wMs_inverse1 = *(workstate1.Ms_inverse);
  Oxs_MeshValue<OC_REAL8m>& wMs2 = *(workstate2.Ms);
  Oxs_MeshValue<OC_REAL8m>& wMs_inverse2 = *(workstate2.Ms_inverse);
  ThreeVector tempspin;

  // Sublattice 1
  for(i=0;i<size;++i) {
    // Transverse movement
    tempspin = dm_dt_t1[i];
    tempspin *= stepsize;

    // For improved accuracy, adjust step vector so that
    // to first order m0 + adjusted_step = v/|v| where
    // v = m0 + step.
    OC_REAL8m adj = 0.5 * tempspin.MagSq();
    tempspin -= adj*cstate1.spin[i];
    tempspin *= 1.0/(1.0+adj);
    tempspin += cstate1.spin[i];
    tempspin.MakeUnit();
    workstate1.spin[i] = tempspin;

    // Longitudinal movement
    tempspin = dm_dt_l1[i]*stepsize;
    tempspin += cstate1.spin[i];

    // Update Ms in the next state.
    // Both of wMs and wMs_inverse should be updated at the same time.
    OC_REAL8m Ms_temp = wMs1[i];
    wMs1[i] = sqrt(tempspin.MagSq())*Ms_temp;
    if(tempspin*cstate1.spin[i]<0.0) {  // Dot product
      // If spin overshoots to the opposite direction, keep Ms positive
      // and flip spin direction.
      workstate1.spin[i] *= -1;
    }
    if(wMs1[i] > (*cstate1.Ms0)[i]) {
      // Ms cannot be >Ms0.
      wMs1[i] = (*cstate1.Ms0)[i];
    }
    if(wMs1[i] != 0.0) {
      wMs_inverse1[i] = 1.0/wMs1[i];
    } else {
      wMs_inverse1[i] = 0.0;
    }
  }
  const Oxs_SimState& nstate1
    = next_state1.GetReadReference();  // Release write lock

  // Sublattice 2
  for(i=0;i<size;++i) {
    // Transverse movement
    tempspin = dm_dt_t2[i];
    tempspin *= stepsize;

    // For improved accuracy, adjust step vector so that
    // to first order m0 + adjusted_step = v/|v| where
    // v = m0 + step.
    OC_REAL8m adj = 0.5 * tempspin.MagSq();
    tempspin -= adj*cstate2.spin[i];
    tempspin *= 1.0/(1.0+adj);
    tempspin += cstate2.spin[i];
    tempspin.MakeUnit();
    workstate2.spin[i] = tempspin;

    // Longitudinal movement
    tempspin = dm_dt_l2[i]*stepsize;
    tempspin += cstate2.spin[i];

    // Update Ms in the next state.
    // Both of wMs2 and wMs_inverse2 should be updated at the same time.
    OC_REAL8m Ms_temp = wMs2[i];
    wMs2[i] = sqrt(tempspin.MagSq())*Ms_temp;
    if(tempspin*cstate2.spin[i]<0.0) {  // Dot product
      // If spin overshoots to the opposite direction, keep Ms positive
      // and flip spin direction.
      workstate2.spin[i] *= -1;
    }
    if(wMs2[i] > (*cstate2.Ms0)[i]) {
      // Ms cannot be >Ms0.
      wMs2[i] = (*cstate2.Ms0)[i];
    }
    if(wMs2[i] != 0.0) {
      wMs_inverse2[i] = 1.0/wMs2[i];
    } else {
      wMs_inverse2[i] = 0.0;
    }
  }
  const Oxs_SimState& nstate2
    = next_state2.GetReadReference();  // Release write lock

  // Total spin
  for(i=0; i<size; i++) {
    tempspin = wMs1[i]*workstate1.spin[i];
    tempspin += wMs2[i]*workstate2.spin[i];
    wMs[i] = sqrt(tempspin.MagSq());
    tempspin.MakeUnit();
    workstate.spin[i] = tempspin;
    if(wMs[i] != 0.0) {
      wMs_inverse[i] = 1.0/wMs[i];
    } else {
      wMs_inverse[i] = 0.0;
    }
  }
  const Oxs_SimState& nstate
    = next_state.GetReadReference();  // Release write lock

  //  Calculate delta E
  OC_REAL8m new_pE_pt1, new_pE_pt2;
  GetEnergyDensity(
      nstate,
      new_energy,
      &mxH1_output.cache.value,
      &mxH2_output.cache.value,
      &total_field1,
      &total_field2,
      new_pE_pt1);
  mxH1_output.cache.state_id=nstate.Id();
  mxH2_output.cache.state_id=nstate.Id();
  const Oxs_MeshValue<ThreeVector>& mxH1 = mxH1_output.cache.value;
  const Oxs_MeshValue<ThreeVector>& mxH2 = mxH2_output.cache.value;

  OC_REAL8m dE=0.0;
  OC_REAL8m var_dE=0.0;
  OC_REAL8m total_E=0.0;
  for(i=0;i<size;++i) {
    OC_REAL8m vol = nstate1.mesh->Volume(i);
    OC_REAL8m e = energy[i];
    total_E += e * vol;
    OC_REAL8m new_e = new_energy[i];
    dE += (new_e - e) * vol;
    var_dE += (new_e*new_e + e*e)*vol*vol;
  }
  var_dE *= 256*OC_REAL8_EPSILON*OC_REAL8_EPSILON/3.; // Variance, assuming
  /// error in each energy[i] term is independent, uniformly
  /// distributed, 0-mean, with range +/- 16*OC_REAL8_EPSILON*energy[i].
  /// It would probably be better to get an error estimate directly
  /// from each energy term.

  // Get error estimate.  See step size adjustment discussion in
  // MJD Notes II, p72 (18-Jan-2001).

  OC_REAL8m new_max_dm_dt, new_max_dm_dt2;
  OC_REAL8m new_dE_dt1, new_timestep_lower_bound;
  OC_REAL8m new_dE_dt2, new_timestep_lower_bound2;

  // For sublattice 1
  Calculate_dm_dt(
      nstate1, 
      mxH1, 
      total_field1, 
      new_pE_pt1, 
      new_dm_dt_t1,
      new_dm_dt_l1,
      new_max_dm_dt, 
      new_dE_dt1, 
      new_timestep_lower_bound);

  // For sublattice 2
  Calculate_dm_dt(
      nstate2,
      mxH2,
      total_field2,
      new_pE_pt2,
      new_dm_dt_t2,
      new_dm_dt_l2,
      new_max_dm_dt2,
      new_dE_dt2,
      new_timestep_lower_bound2);

  // TODO: check sublattice 2 as well
  OC_REAL8m max_error=0;
  for(i=0;i<size;++i) { 
    ThreeVector temp = dm_dt_t1[i] + dm_dt_l1[i];
    temp -= new_dm_dt_t1[i];
    temp -= new_dm_dt_l1[i];
    OC_REAL8m temp_error = temp.MagSq();
    if(temp_error>max_error) max_error = temp_error;
  }
  max_error = sqrt(max_error)/2.0; // Actual (local) error
  /// estimate is max_error * stepsize

  // Energy check control
#ifdef FOO
  OC_REAL8m expected_dE = 0.5 * (dE_dt+new_dE_dt1) * stepsize;
  OC_REAL8m dE_error = dE - expected_dE;
  OC_REAL8m max_allowed_dE = expected_dE + 0.25*fabs(expected_dE);
  max_allowed_dE += OC_REAL8_EPSILON*fabs(total_E);
  max_allowed_dE += 2*sqrt(var_dE);
#else
  OC_REAL8m max_allowed_dE = 0.5 * (pE_pt+new_pE_pt1) * stepsize
    + OC_MAX(OC_REAL8_EPSILON*fabs(total_E),2*sqrt(var_dE));
  /// The above says essentially that the spin adjustment can
  /// increase the energy by only as much as pE/pt allows; in
  /// the absence of pE/pt, the energy should decrease.  I
  /// think this may be problematic, if at the start of a stage
  /// the spins are near equilibrium, and the applied field is
  /// ramping up slowly.  In this case there won't be much "give"
  /// in the spin configuration with respect to pE/pm.  But I
  /// haven't seen an example of this yet, so we'll wait and see.
  /// -mjd, 27-July-2001.
#endif

  // Check step and adjust next_timestep.  The relative error
  // check is a bit fudged, because rather than limiting the
  // relative error uniformly across the sample, we limit it
  // only at the position that has the maximum absolute error
  // (i.e., max_error is max *absolute* error).  I haven't
  // tested to see if uniformly limiting relative error is
  // workable (it might be too restrictive for most purposes),
  // but the present setup seems to solve the problem of convergence
  // near equilibrium.  -mjd, 2001-02-23.
  // NOTE: Since all three error controls (error_rate,
  //  absolute_step_error, and relative_step_error) assume error
  //  grows linearly with step size, we can check up front to see
  //  which control is most restrictive, store that constraint in
  //  working_allowed_error, and then adjust the step size without
  //  regard to which control is being exercised.
  OC_REAL8m working_allowed_error
    = max_step_increase*max_error/step_headroom;
  if(allowed_error_rate>=0.
     && working_allowed_error>allowed_error_rate) {
    working_allowed_error=allowed_error_rate;
  }
  if(allowed_absolute_step_error>=0.
     && stepsize*working_allowed_error>allowed_absolute_step_error) {
    working_allowed_error=allowed_absolute_step_error/stepsize;
  }
  if(allowed_relative_step_error>=0.
     && working_allowed_error>allowed_relative_step_error*max_dm_dt) {
    working_allowed_error = allowed_relative_step_error * max_dm_dt;
  }
  if(!forcestep) {
    next_timestep=1.0;  // Size relative to current step
    if(max_error>working_allowed_error) {
      next_timestep = step_headroom*working_allowed_error/max_error;
    } else if(dE>max_allowed_dE) {
      // Energy check
      next_timestep=0.5;
    }
    if(next_timestep<1.0) {
      // Reject step
      if(next_timestep<max_step_decrease)
  next_timestep=max_step_decrease;
      next_timestep *= stepsize;
      return 0;
    }
  }

  // TODO: Report value for lattice2 too.
  // Otherwise, accept step.  Calculate next step using
  // estimate of step size that would just meet the error
  // restriction (with "headroom" safety margin).
  next_timestep = max_step_increase;
  if(next_timestep*max_error>step_headroom*working_allowed_error) {
    next_timestep = step_headroom*working_allowed_error/max_error;
  }
  if(next_timestep<max_step_decrease)
    next_timestep=max_step_decrease;
  next_timestep *= stepsize;
  if(!nstate.AddDerivedData("Timestep lower bound",
          new_timestep_lower_bound) ||
     !nstate.AddDerivedData("Max dm/dt",new_max_dm_dt) ||
     !nstate.AddDerivedData("dE/dt",new_dE_dt1) ||
     !nstate.AddDerivedData("Delta E",dE) ||
     !nstate.AddDerivedData("pE/pt",new_pE_pt1)) {
    throw Oxs_Ext::Error(this,
       "YY_2LatEulerEvolve::Step:"
       " Programming error; data cache already set.");
  }

  dm_dt_t1_output.cache.value.Swap(new_dm_dt_t1);
  dm_dt_l1_output.cache.value.Swap(new_dm_dt_l1);
  dm_dt_t1_output.cache.state_id = nstate.Id();
  dm_dt_l1_output.cache.state_id = nstate.Id();
  dm_dt_t2_output.cache.value.Swap(new_dm_dt_t2);
  dm_dt_l2_output.cache.value.Swap(new_dm_dt_l2);
  dm_dt_t2_output.cache.state_id = nstate.Id();
  dm_dt_l2_output.cache.state_id = nstate.Id();

  energy.Swap(new_energy);
  energy_state_id = nstate.Id();

  return 1;  // Good step
}   // end Step

// Call with the total_lattice state and it updates values for both
// sublattices.
void YY_2LatEulerEvolve::UpdateMeshArrays(const Oxs_SimState& state)
{
  mesh_id = 0; // Mark update in progress
  const Oxs_Mesh* mesh = state.mesh;
  const Oxs_MeshValue<OC_REAL8m>& Ms10 = *(state.lattice1->Ms0);
  const Oxs_MeshValue<OC_REAL8m>& Ms20 = *(state.lattice2->Ms0);
  const Oxs_MeshValue<OC_REAL8m>& Ms10_inverse = *(state.lattice1->Ms0_inverse);
  const Oxs_MeshValue<OC_REAL8m>& Ms20_inverse = *(state.lattice2->Ms0_inverse);
  const Oxs_MeshValue<OC_REAL8m>& Tc1 = *(state.lattice1->Tc);
  const Oxs_MeshValue<OC_REAL8m>& Tc2 = *(state.lattice2->Tc);
  const OC_INDEX size = mesh->Size();
  OC_INDEX i;
  // Note: Tc1 and Tc2 are functions of T, m_e1, and m_e2.

  for(i=0;i<size;i++) {
    if(temperature[i] > Tc1[i]) {
      alpha_t1[i] = 2./3.*alpha_t10[i];
      alpha_l1[i] = alpha_t1[i];
    } else {
      alpha_l1[i] = alpha_t10[i]*2*temperature[i]/(3*Tc1[i]);
      alpha_t1[i] = alpha_t10[i]*(1-temperature[i]/(3*Tc1[i]));
    }
    if(temperature[i] > Tc2[i]) {
      alpha_t2[i] = 2./3.*alpha_t20[i];
      alpha_l2[i] = alpha_t2[i];
    } else {
      alpha_l2[i] = alpha_t20[i]*2*temperature[i]/(3*Tc2[i]);
      alpha_t2[i] = alpha_t20[i]*(1-temperature[i]/(3*Tc2[i]));
    }

    // Update variance of stochastic field
    // h_fluctVarConst_t = 2*(alpha_t-alpha_l)*kB*T/(MU0*gamma*alpha_t^2*Ms0*Vol)
    // h_fluctVarConst_l = 2*(alpha_l-gamma)*kB*T/(MU0*Ms0*Vol)
    OC_REAL8m cell_alpha_t1 = fabs(alpha_t1[i]);
    OC_REAL8m cell_alpha_t2 = fabs(alpha_t2[i]);
    OC_REAL8m cell_alpha_l1 = fabs(alpha_l1[i]);
    OC_REAL8m cell_alpha_l2 = fabs(alpha_l2[i]);
    OC_REAL8m cell_gamma1 = fabs(gamma1[i]);
    OC_REAL8m cell_gamma2 = fabs(gamma2[i]);
    OC_REAL8m cell_vol = mesh->Volume(i);
    hFluctVarConst_t1[i] = 2*fabs(cell_alpha_t1-cell_alpha_l1);
    hFluctVarConst_t1[i] *= kB_T[i]*Ms10_inverse[i];
    hFluctVarConst_t1[i] /= MU0*cell_gamma1
      *cell_alpha_t1*cell_alpha_t1*cell_vol;
    hFluctVarConst_t2[i] = 2*fabs(cell_alpha_t2-cell_alpha_l2);
    hFluctVarConst_t2[i] *= kB_T[i]*Ms20_inverse[i];
    hFluctVarConst_t2[i] /= MU0*cell_gamma2
      *cell_alpha_t2*cell_alpha_t2*cell_vol;
    hFluctVarConst_l1[i] = 2*cell_alpha_l1*cell_gamma1;
    hFluctVarConst_l1[i] *= kB_T[i]*Ms10_inverse[i];
    hFluctVarConst_l1[i] /= MU0*cell_vol;
    hFluctVarConst_l2[i] = 2*cell_alpha_l2*cell_gamma2;
    hFluctVarConst_l2[i] *= kB_T[i]*Ms20_inverse[i];
    hFluctVarConst_l2[i] /= MU0*cell_vol;
  }

  mesh_id = mesh->Id();
}

void YY_2LatEulerEvolve::UpdateDerivedOutputs(const Oxs_SimState& state)
{ // This routine fills all the YY_2LatEulerEvolve Oxs_ScalarOutput's to
  // the appropriate value based on the import "state", and any of
  // Oxs_VectorOutput's that have CacheRequest enabled are filled.
  // It also makes sure all the expected WOO objects in state are
  // filled.
  max_dm_dt_output.cache.state_id
    = dE_dt_output.cache.state_id
    = delta_E_output.cache.state_id
    = 0;  // Mark change in progress

  // If temperature has not been set up, do so. It is required in exchange
  // energy calculation.
  if(state.lattice1->T==NULL) {
    UpdateStageTemperature(state);
    state.lattice1->T = &temperature;
    state.lattice2->T = &temperature;
  }

  OC_REAL8m dummy_value;
  if(!state.GetDerivedData("Max dm/dt",max_dm_dt_output.cache.value) ||
     !state.GetDerivedData("dE/dt",dE_dt_output.cache.value) ||
     !state.GetDerivedData("Delta E",delta_E_output.cache.value) ||
     !state.GetDerivedData("pE/pt",dummy_value) ||
     !state.GetDerivedData("Timestep lower bound",dummy_value) ||
     (dm_dt_t1_output.GetCacheRequestCount()>0
      && dm_dt_t1_output.cache.state_id != state.Id()) ||
     (dm_dt_l1_output.GetCacheRequestCount()>0
      && dm_dt_l1_output.cache.state_id != state.Id()) ||
     (mxH1_output.GetCacheRequestCount()>0
      && mxH1_output.cache.state_id != state.Id()) ||
     (dm_dt_t2_output.GetCacheRequestCount()>0
      && dm_dt_t2_output.cache.state_id != state.Id()) ||
     (dm_dt_l2_output.GetCacheRequestCount()>0
      && dm_dt_l2_output.cache.state_id != state.Id()) ||
     (mxH2_output.GetCacheRequestCount()>0
      && mxH2_output.cache.state_id != state.Id()) ) {

    // Missing at least some data, so calculate from scratch

    // Calculate H and mxH outputs
    Oxs_MeshValue<ThreeVector>& mxH1 = mxH1_output.cache.value;
    Oxs_MeshValue<ThreeVector>& mxH2 = mxH2_output.cache.value;
    OC_REAL8m pE_pt;
    GetEnergyDensity(
        state,
        energy,
        &mxH1,
        &mxH2,
        &total_field1,
        &total_field2,
        pE_pt);
    energy_state_id=state.Id();
    mxH1_output.cache.state_id=state.Id();
    mxH2_output.cache.state_id=state.Id();

    // TODO: How can I include pE_pt2?
    if(!state.GetDerivedData("pE/pt",dummy_value)) {
      state.AddDerivedData("pE/pt",pE_pt);
    }

    // Calculate dm/dt, Max dm/dt and dE/dt
    OC_REAL8m timestep_lower_bound;
    Oxs_MeshValue<ThreeVector>& dm_dt_t1
      = dm_dt_t1_output.cache.value;
    Oxs_MeshValue<ThreeVector>& dm_dt_l1
      = dm_dt_l1_output.cache.value;
    dm_dt_t1_output.cache.state_id=0;
    dm_dt_l1_output.cache.state_id=0;
    Calculate_dm_dt(
        *(state.lattice1),
        mxH1,
        total_field1,
        pE_pt,
        dm_dt_t1,
        dm_dt_l1,
        max_dm_dt_output.cache.value,
        dE_dt_output.cache.value,
        timestep_lower_bound);
    dm_dt_t1_output.cache.state_id=state.Id();
    dm_dt_l1_output.cache.state_id=state.Id();

    Oxs_MeshValue<ThreeVector>& dm_dt_t2
      = dm_dt_t2_output.cache.value;
    Oxs_MeshValue<ThreeVector>& dm_dt_l2
      = dm_dt_l2_output.cache.value;
    dm_dt_t2_output.cache.state_id=0;
    dm_dt_l2_output.cache.state_id=0;
    Calculate_dm_dt(
        *(state.lattice2),
        mxH2,
        total_field2,
        pE_pt,
        dm_dt_t2,
        dm_dt_l2,
        max_dm_dt_output.cache.value,
        dE_dt_output.cache.value,
        timestep_lower_bound);
    dm_dt_t2_output.cache.state_id=state.Id();
    dm_dt_l2_output.cache.state_id=state.Id();

    if(!state.GetDerivedData("Max dm/dt",dummy_value)) {
      state.AddDerivedData("Max dm/dt",max_dm_dt_output.cache.value);
    }
    if(!state.GetDerivedData("dE/dt",dummy_value)) {
      state.AddDerivedData("dE/dt",dE_dt_output.cache.value);
    }
    if(!state.GetDerivedData("Timestep lower bound",dummy_value)) {
      state.AddDerivedData("Timestep lower bound",
         timestep_lower_bound);
    }

    if(!state.GetDerivedData("Delta E",dummy_value)) {
      if(state.previous_state_id!=0 && state.stage_iteration_count>0) {
  // Strictly speaking, we should be able to create dE for
  // stage_iteration_count==0 for stages>0, but as a practical
  // matter we can't at present.  Should give this more thought.
  // -mjd, 27-July-2001
  throw Oxs_Ext::Error(this,
     "YY_2LatEulerEvolve::UpdateDerivedOutputs:"
     " Can't derive Delta E from single state.");
      }
      state.AddDerivedData("Delta E",0.0);
      dummy_value = 0.;
    }
    delta_E_output.cache.value=dummy_value;
  }

  max_dm_dt_output.cache.value*=(180e-9/PI);
  /// Convert from radians/second to deg/ns

  max_dm_dt_output.cache.state_id
    = dE_dt_output.cache.state_id
    = delta_E_output.cache.state_id
    = state.Id();
}   // end UpdateDerivedOutputs

OC_REAL8m YY_2LatEulerEvolve::Gaussian_Random(const OC_REAL8m muGaus,
    const OC_REAL8m sigmaGaus)
{
  // Box-Muller algorithm, see W.H. Press' "Numerical recipes" chapter7.2 
  // for details.
  OC_REAL8m R, gaus1, FAC;

  if (!gaus2_isset) {
    R = 1.;
    while (R >= 1.){
      gaus1 = 2. * Oc_UnifRand() - 1.;
      gaus2 = 2. * Oc_UnifRand() - 1.;
      R  = gaus1*gaus1 + gaus2*gaus2;
    }
    gaus2_isset = 1;
    FAC = sqrt(-2. * log(R) / R);
    gaus1 = gaus1 * FAC * sigmaGaus + muGaus;
    gaus2 = gaus2 * FAC * sigmaGaus + muGaus;
    return gaus1;
  }

  gaus2_isset = false;
  return gaus2;
}
