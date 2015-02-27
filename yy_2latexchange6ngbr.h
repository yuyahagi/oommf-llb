/** FILE: yy_2lat_util.h                 -*-Mode: c++-*-
 *
 * 6 neighbor exchange energy on rectangular mesh for 2-lattice
 * simulations. It is based on Oxs_Exchange6Ngbr class.
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

#ifndef _YY_2LATEXCHANGE6NGBR
#define _YY_2LATEXCHANGE6NGBR

#include "atlas.h"
#include "key.h"
#include "chunkenergy.h"
#include "energy.h"
#include "mesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "rectangularmesh.h"

/* End includes */

#define DEFAULT_M_E_TOL 1e-4

class YY_2LatExchange6Ngbr : public Oxs_ChunkEnergy {
private:
  enum ExchangeCoefType {
    A_UNKNOWN, A_TYPE, LEX_TYPE
  }  excoeftype;

  OC_INDEX coef_size;
  OC_REAL8m** coef1;
  OC_REAL8m** coef2;
  OC_REAL8m** coef12;
  mutable Oxs_Key<Oxs_Atlas> atlaskey;  
  Oxs_OwnedPointer<Oxs_Atlas> atlas;
  mutable Oxs_ThreadControl thread_control;
  mutable OC_UINT4m mesh_id;
  mutable Oxs_MeshValue<OC_INT4m> region_id;

  // Support for threaded maxang calculations
  mutable vector<OC_REAL8m> maxdot;

  void CalcEnergyA(const Oxs_SimState& state,
                   Oxs_ComputeEnergyDataThreaded& ocedt,
                   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
                   OC_INDEX node_start,OC_INDEX node_stop,
                   int threadnumber) const;
  void CalcEnergyLex(const Oxs_SimState& state,
                     Oxs_ComputeEnergyDataThreaded& ocedt,
                     Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
                     OC_INDEX node_start,OC_INDEX node_stop,
                     int threadnumber) const;

  // Parameters used for longitudinal susceptibility
  // Exchange parameter J = nJ_0 and atomistic magnetic moment mu, where
  // n is the number of neighboring atoms
  Oxs_OwnedPointer<Oxs_ScalarField> J01_init, J02_init;
  Oxs_OwnedPointer<Oxs_ScalarField> J012_init, J021_init;
  Oxs_OwnedPointer<Oxs_ScalarField> mu1_init, mu2_init;
  mutable Oxs_MeshValue<OC_REAL8m> J01, J02;
  mutable Oxs_MeshValue<OC_REAL8m> J012, J021;
  mutable Oxs_MeshValue<OC_REAL8m> mu1, mu2;
  // Currie temperature in Kelvin, calculated from J and mu
  mutable Oxs_MeshValue<OC_REAL8m> Tc1, Tc2;

  // Members for calculating m_e, equilibrium spin polarization at
  // temperature T and chi_l, longitudinal susceptibility.
  mutable OC_INDEX last_stage_number;
  mutable Oxs_MeshValue<OC_REAL8m> m_e1, m_e2;
  mutable Oxs_MeshValue<OC_REAL8m> chi_l1, chi_l2;
  // Langevin function and its derivative
  OC_REAL8m Langevin(OC_REAL8m x) const;
  OC_REAL8m LangevinDeriv(OC_REAL8m x) const;
  void Update_m_e_chi_l(const Oxs_SimState& state, OC_REAL8m tol) const;
  void Update_m_e_chi_l(const Oxs_SimState& state) const {
    return Update_m_e_chi_l(state, DEFAULT_M_E_TOL);
  }
  mutable Oxs_MeshValue<OC_REAL8m> G1, G2;

  // Supplied outputs, in addition to those provided by Oxs_Energy.
  Oxs_ScalarOutput<YY_2LatExchange6Ngbr> maxspinangle_output;
  Oxs_ScalarOutput<YY_2LatExchange6Ngbr> stage_maxspinangle_output;
  Oxs_ScalarOutput<YY_2LatExchange6Ngbr> run_maxspinangle_output;
  void UpdateDerivedOutputs(const Oxs_SimState& state);
  String MaxSpinAngleStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Max Spin Angle";
    return dummy_name;
  }
  String StageMaxSpinAngleStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Stage Max Spin Angle";
    return dummy_name;
  }
  String RunMaxSpinAngleStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Run Max Spin Angle";
    return dummy_name;
  }

protected:
  virtual void GetEnergy(const Oxs_SimState& state,
			 Oxs_EnergyData& oed) const {
    GetEnergyAlt(state,oed);
  }

  virtual void ComputeEnergy(const Oxs_SimState& state,
                             Oxs_ComputeEnergyData& oced) const {
    ComputeEnergyAlt(state,oced);
  }

  virtual void ComputeEnergyChunkInitialize
  (const Oxs_SimState& state,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   int number_of_threads) const;

  virtual void ComputeEnergyChunkFinalize
  (const Oxs_SimState& state,
   const Oxs_ComputeEnergyDataThreaded& ocedt,
   const Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   int number_of_threads) const;

  virtual void ComputeEnergyChunk(const Oxs_SimState& state,
                                  Oxs_ComputeEnergyDataThreaded& ocedt,
                                  Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
                                  OC_INDEX node_start,OC_INDEX node_stop,
                                  int threadnumber) const;

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  YY_2LatExchange6Ngbr(const char* name,     // Child instance id
		    Oxs_Director* newdtr, // App director
		    const char* argstr);  // MIF input block parameters
  virtual ~YY_2LatExchange6Ngbr();
  virtual OC_BOOL Init();
};


#undef DEFAULT_M_E_TOL
#endif // _YY_2LATEXCHANGE6NGBR
