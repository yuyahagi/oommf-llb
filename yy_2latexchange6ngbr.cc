/** FILE: yy_2lat_util.cc                 -*-Mode: c++-*-
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

#include <string>

#include "atlas.h"
#include "nb.h"
#include "key.h"
#include "director.h"
#include "mesh.h"
#include "meshvalue.h"
#include "oxswarn.h"
#include "simstate.h"
#include "threevector.h"
#include "rectangularmesh.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

#include "yy_2lat_util.h"
#include "yy_2latexchange6ngbr.h"

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(YY_2LatExchange6Ngbr);

/* End includes */

// Revision information, set via CVS keyword substitution
static const Oxs_WarningMessageRevisionInfo revision_info
  (__FILE__,
   "$Revision: 1.48 $",
   "$Date: 2012-09-07 02:44:30 $",
   "$Author: donahue $",
   "Michael J. Donahue (michael.donahue@nist.gov)");

// Constructor
YY_2LatExchange6Ngbr::YY_2LatExchange6Ngbr(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_ChunkEnergy(name,newdtr,argstr),
    coef_size(0), mesh_id(0),
    coef1(NULL), coef2(NULL), coef12(NULL),
    last_stage_number(-1),
    tol(1e-4), tolsq(1e-4)
{
  // Process arguments
  OXS_GET_INIT_EXT_OBJECT("atlas",Oxs_Atlas,atlas);
  atlaskey.Set(atlas.GetPtr());
  /// Dependency lock is held until *this is deleted.

  if(HasInitValue("J01")) {
    OXS_GET_INIT_EXT_OBJECT("J01",Oxs_ScalarField,J01_init);
  } else {
    throw Oxs_Ext::Error(this,"Exchange parameter J01 not specified.\n");
  }

  if(HasInitValue("J02")) {
    OXS_GET_INIT_EXT_OBJECT("J02",Oxs_ScalarField,J02_init);
  } else {
    throw Oxs_Ext::Error(this,"Exchange parameter J02 not specified.\n");
  }

  if(HasInitValue("J012")) {
    OXS_GET_INIT_EXT_OBJECT("J012",Oxs_ScalarField,J012_init);
  } else {
    throw Oxs_Ext::Error(this,"Exchange parameter J012 not specified.\n");
  }

  if(HasInitValue("J021")) {
    OXS_GET_INIT_EXT_OBJECT("J021",Oxs_ScalarField,J021_init);
  } else {
    throw Oxs_Ext::Error(this,"Exchange parameter J021 not specified.\n");
  }

  if(HasInitValue("atom_moment1")) {
    OXS_GET_INIT_EXT_OBJECT("atom_moment1",Oxs_ScalarField,mu1_init);
  } else {
    throw Oxs_Ext::Error(this, "Atomic magnetic moment atom_moment1"
        " is not specified.");
  }

  if(HasInitValue("atom_moment2")) {
    OXS_GET_INIT_EXT_OBJECT("atom_moment2",Oxs_ScalarField,mu2_init);
  } else {
    throw Oxs_Ext::Error(this, "Atomic magnetic moment atom_moment2"
        " is not specified.");
  }

  // Determine number of regions, and check that the
  // count lies within the allowed range.
  coef_size = atlas->GetRegionCount();
  if(coef_size<1) {
    String msg = String("Oxs_Atlas object ")
      + atlas->InstanceName()
      + String(" must contain at least one region.");

    throw Oxs_ExtError(msg.c_str());
  }
  const OC_INDEX sqrootbits = (sizeof(OC_INDEX)*8-1)/2;
  // Note assumption that OC_INDEX is a signed type.
  // Also, sizeof(OC_INDEX)*8-1 is odd, so sqrootbits is shy
  // by half a bit.  We adjust for this below by multiplying
  // by a rational approximation to sqrt(2).
  OC_INDEX coef_check = 1;
  for(OC_INDEX ibit=0;ibit<sqrootbits;++ibit) coef_check *= 2;
  if(sizeof(OC_INDEX)==1) {
    coef_check = 11;  // Table look-up <g>
  } else if(sizeof(OC_INDEX)<5) {
    // Use sqrt(2) > 239/169
    coef_check = (coef_check*239+168)/169 - 1;
  } else {
    // Use sqrt(2) > 275807/195025.  This is accurate for
    // sizeof(OC_INDEX) <= 9, but will round low above that.
    coef_check = (coef_check*275807+195024)/195025 - 1;
  }
  const OC_INDEX max_coef_size = coef_check;
  if(coef_size>max_coef_size) {
      char buf[4096];
      Oc_Snprintf(buf,sizeof(buf),
		  "Oxs_Atlas object %.1000s has too many regions: %lu"
                  " (max allowed: %lu)",
                  atlas->InstanceName(),
                  (unsigned long)coef_size,
                  (unsigned long)max_coef_size);
      throw Oxs_ExtError(this,buf);
  }

  // Determine coef matrix fill type
  String typestr;
  vector<String> params1;
  vector<String> params2;
  vector<String> params12;
  OC_REAL8m default_coef1 = 0.0;
  OC_REAL8m default_coef2 = 0.0;
  OC_REAL8m default_coef12 = 0.0;
  OC_BOOL has_A = HasInitValue("A1")
    || HasInitValue("A2") || HasInitValue("A12");
  OC_BOOL has_lex = HasInitValue("lex1")
    || HasInitValue("lex2") || HasInitValue("lex12");
  if(has_A && has_lex) {
    throw Oxs_ExtError(this,"Invalid exchange coefficient request:"
			 " both A and lex specified; only one should"
			 " be given.");
  } else if(has_lex) {
    excoeftype = LEX_TYPE;
    typestr = "lex";
    default_coef1 = GetRealInitValue("default_lex1",0.0);
    default_coef2 = GetRealInitValue("default_lex2",0.0);
    default_coef12 = GetRealInitValue("default_lex12",0.0);
    FindRequiredInitValue("lex1",params1);
    FindRequiredInitValue("lex2",params2);
    FindRequiredInitValue("lex12",params12);
  } else {
    excoeftype = A_TYPE;
    typestr = "A";
    default_coef1 = GetRealInitValue("default_A1",0.0);
    default_coef2 = GetRealInitValue("default_A2",0.0);
    default_coef12 = GetRealInitValue("default_A12",0.0);
    FindRequiredInitValue("A1",params1);
    FindRequiredInitValue("A2",params2);
    FindRequiredInitValue("A12",params12);
  }
  if( params1.empty() || params2.empty() || params12.empty() ) {
    String msg = String("Empty parameter list for key \"")
      + typestr + String("\"");
    throw Oxs_ExtError(this,msg);
  }
  if( params1.size()%3!=0 || params2.size()%3!=0 || params12.size()%3!=0 ) {
      char buf[4096];
      Oc_Snprintf(buf,sizeof(buf),
		  "Number of elements in %.80s sub-list must be"
		  " divisible by 3"
		  " (actual sub-list size: %u)",
		  typestr.c_str(),(unsigned int)params1.size(), ", ",
		  typestr.c_str(),(unsigned int)params2.size(), ", ",
		  typestr.c_str(),(unsigned int)params12.size());
      throw Oxs_ExtError(this,buf);
  }

  // Allocate A matrix.  Because raw pointers are used, a memory leak
  // would occur if an uncaught exception occurred beyond this point.
  // So, we put the whole bit in a try-catch block.
  // =======================================================================
  // Lattice 1
  // =======================================================================
  try {
    coef1 = new OC_REAL8m*[coef_size];
    coef1[0] = NULL; // Safety, in case next alloc throws an exception
    coef1[0] = new OC_REAL8m[coef_size*coef_size];
    OC_INDEX ic;
    for(ic=1;ic<coef_size;ic++) coef1[ic] = coef1[ic-1] + coef_size;

    // Fill A matrix
    for(ic=0;ic<coef_size*coef_size;ic++) coef1[0][ic] = default_coef1;
    for(OC_INT4m ip=0;static_cast<size_t>(ip)<params1.size();ip+=3) {
      OC_INT4m i1 = atlas->GetRegionId(params1[ip]);
      OC_INT4m i2 = atlas->GetRegionId(params1[ip+1]);
      if(i1<0 || i2<0) {
        // Unknown region(s) requested
        char buf[4096];
        char* cptr=buf;
        if(i1<0) {
          char item[96];  // Safety
          item[80]='\0';
          Oc_Snprintf(item,sizeof(item),"%.82s",params1[ip].c_str());
          if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
          Oc_Snprintf(buf,sizeof(buf),
                      "First entry in %.80s[%ld] sub-list, \"%.85s\","
                      " is not a known region in atlas \"%.1000s\".  ",
                      typestr.c_str(),long(ip/3),item,
                      atlas->InstanceName());
          cptr += strlen(buf);
        }
        if(i2<0) {
          char item[96];  // Safety
          item[80]='\0';
          Oc_Snprintf(item,sizeof(item),"%.82s",params1[ip+1].c_str());
          if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
          Oc_Snprintf(cptr,sizeof(buf)-(cptr-buf),
                      "Second entry in %.80s[%ld] sub-list, \"%.85s\","
                      " is not a known region in atlas \"%.1000s\".  ",
                      typestr.c_str(),long(ip/3),item,
                      atlas->InstanceName());
        }
        String msg = String(buf);
        msg += String("Known regions:");
        vector<String> regions;
        atlas->GetRegionList(regions);
        for(OC_INT4m j=0;static_cast<size_t>(j)<regions.size();++j) {
          msg += String(" \n");
          msg += regions[j];
        }
        throw Oxs_ExtError(this,msg);
      }
      OC_BOOL err;
      OC_REAL8m coef1pair = Nb_Atof(params1[ip+2].c_str(),err);
      if(err) {
        char item[96];  // Safety
        item[80]='\0';
        Oc_Snprintf(item,sizeof(item),"%.82s",params1[ip+2].c_str());
        if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
        char buf[4096];
        Oc_Snprintf(buf,sizeof(buf),
                    "Third entry in %.80s[%ld] sub-list, \"%.85s\","
                    " is not a valid floating point number.",
                    typestr.c_str(),long(ip/3),item);
        throw Oxs_ExtError(this,buf);
      }
      coef1[i1][i2]=coef1pair;
      coef1[i2][i1]=coef1pair; // coef should be symmetric
    }
    DeleteInitValue(typestr+"1");

  }
  catch(...) {
    if(coef_size>0 && coef1!=NULL) { // Release coef memory
      delete[] coef1[0]; delete[] coef1;
      coef1 = NULL;
      coef_size = 0; // Safety
    }
    throw;
  }

  // =======================================================================
  // Lattice 2
  // =======================================================================
  try {
    coef2 = new OC_REAL8m*[coef_size];
    coef2[0] = NULL; // Safety, in case next alloc throws an exception
    coef2[0] = new OC_REAL8m[coef_size*coef_size];
    OC_INDEX ic;
    for(ic=1;ic<coef_size;ic++) coef2[ic] = coef2[ic-1] + coef_size;

    // Fill A matrix
    for(ic=0;ic<coef_size*coef_size;ic++) coef2[0][ic] = default_coef2;
    for(OC_INT4m ip=0;static_cast<size_t>(ip)<params2.size();ip+=3) {
      OC_INT4m i1 = atlas->GetRegionId(params2[ip]);
      OC_INT4m i2 = atlas->GetRegionId(params2[ip+1]);
      if(i1<0 || i2<0) {
        // Unknown region(s) requested
        char buf[4096];
        char* cptr=buf;
        if(i1<0) {
          char item[96];  // Safety
          item[80]='\0';
          Oc_Snprintf(item,sizeof(item),"%.82s",params2[ip].c_str());
          if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
          Oc_Snprintf(buf,sizeof(buf),
                      "First entry in %.80s[%ld] sub-list, \"%.85s\","
                      " is not a known region in atlas \"%.1000s\".  ",
                      typestr.c_str(),long(ip/3),item,
                      atlas->InstanceName());
          cptr += strlen(buf);
        }
        if(i2<0) {
          char item[96];  // Safety
          item[80]='\0';
          Oc_Snprintf(item,sizeof(item),"%.82s",params2[ip+1].c_str());
          if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
          Oc_Snprintf(cptr,sizeof(buf)-(cptr-buf),
                      "Second entry in %.80s[%ld] sub-list, \"%.85s\","
                      " is not a known region in atlas \"%.1000s\".  ",
                      typestr.c_str(),long(ip/3),item,
                      atlas->InstanceName());
        }
        String msg = String(buf);
        msg += String("Known regions:");
        vector<String> regions;
        atlas->GetRegionList(regions);
        for(OC_INT4m j=0;static_cast<size_t>(j)<regions.size();++j) {
          msg += String(" \n");
          msg += regions[j];
        }
        throw Oxs_ExtError(this,msg);
      }
      OC_BOOL err;
      OC_REAL8m coef2pair = Nb_Atof(params2[ip+2].c_str(),err);
      if(err) {
        char item[96];  // Safety
        item[80]='\0';
        Oc_Snprintf(item,sizeof(item),"%.82s",params2[ip+2].c_str());
        if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
        char buf[4096];
        Oc_Snprintf(buf,sizeof(buf),
                    "Third entry in %.80s[%ld] sub-list, \"%.85s\","
                    " is not a valid floating point number.",
                    typestr.c_str(),long(ip/3),item);
        throw Oxs_ExtError(this,buf);
      }
      coef2[i1][i2]=coef2pair;
      coef2[i2][i1]=coef2pair; // coef should be symmetric
    }
    DeleteInitValue(typestr+"2");

  }
  catch(...) {
    if(coef_size>0 && coef2!=NULL) { // Release coef memory
      delete[] coef2[0]; delete[] coef2;
      coef2 = NULL;
      coef_size = 0; // Safety
    }
    throw;
  }

  // =======================================================================
  // Between lattice 1 and 2
  // =======================================================================
  try {
    coef12 = new OC_REAL8m*[coef_size];
    coef12[0] = NULL; // Safety, in case next alloc throws an exception
    coef12[0] = new OC_REAL8m[coef_size*coef_size];
    OC_INDEX ic;
    for(ic=1;ic<coef_size;ic++) coef12[ic] = coef12[ic-1] + coef_size;

    // Fill A matrix
    for(ic=0;ic<coef_size*coef_size;ic++) coef12[0][ic] = default_coef12;
    for(OC_INT4m ip=0;static_cast<size_t>(ip)<params12.size();ip+=3) {
      OC_INT4m i1 = atlas->GetRegionId(params12[ip]);
      OC_INT4m i2 = atlas->GetRegionId(params12[ip+1]);
      if(i1<0 || i2<0) {
        // Unknown region(s) requested
        char buf[4096];
        char* cptr=buf;
        if(i1<0) {
          char item[96];  // Safety
          item[80]='\0';
          Oc_Snprintf(item,sizeof(item),"%.82s",params12[ip].c_str());
          if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
          Oc_Snprintf(buf,sizeof(buf),
                      "First entry in %.80s[%ld] sub-list, \"%.85s\","
                      " is not a known region in atlas \"%.1000s\".  ",
                      typestr.c_str(),long(ip/3),item,
                      atlas->InstanceName());
          cptr += strlen(buf);
        }
        if(i2<0) {
          char item[96];  // Safety
          item[80]='\0';
          Oc_Snprintf(item,sizeof(item),"%.82s",params12[ip+1].c_str());
          if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
          Oc_Snprintf(cptr,sizeof(buf)-(cptr-buf),
                      "Second entry in %.80s[%ld] sub-list, \"%.85s\","
                      " is not a known region in atlas \"%.1000s\".  ",
                      typestr.c_str(),long(ip/3),item,
                      atlas->InstanceName());
        }
        String msg = String(buf);
        msg += String("Known regions:");
        vector<String> regions;
        atlas->GetRegionList(regions);
        for(OC_INT4m j=0;static_cast<size_t>(j)<regions.size();++j) {
          msg += String(" \n");
          msg += regions[j];
        }
        throw Oxs_ExtError(this,msg);
      }
      OC_BOOL err;
      OC_REAL8m coef12pair = Nb_Atof(params12[ip+2].c_str(),err);
      if(err) {
        char item[96];  // Safety
        item[80]='\0';
        Oc_Snprintf(item,sizeof(item),"%.82s",params12[ip+2].c_str());
        if(item[80]!='\0') strcpy(item+80,"..."); // Overflow
        char buf[4096];
        Oc_Snprintf(buf,sizeof(buf),
                    "Third entry in %.80s[%ld] sub-list, \"%.85s\","
                    " is not a valid floating point number.",
                    typestr.c_str(),long(ip/3),item);
        throw Oxs_ExtError(this,buf);
      }
      coef12[i1][i2]=coef12pair;
      coef12[i2][i1]=coef12pair; // coef should be symmetric
    }
    DeleteInitValue(typestr+"12");

    VerifyAllInitArgsUsed();

    // Setup outputs
    maxspinangle_output.Setup(this,InstanceName(),"Max Spin Ang","deg",1,
                              &YY_2LatExchange6Ngbr::UpdateDerivedOutputs);
    maxspinangle_output.Register(director,0);
    stage_maxspinangle_output.Setup(this,InstanceName(),"Stage Max Spin Ang","deg",1,
                                    &YY_2LatExchange6Ngbr::UpdateDerivedOutputs);
    stage_maxspinangle_output.Register(director,0);
    run_maxspinangle_output.Setup(this,InstanceName(),"Run Max Spin Ang","deg",1,
                                  &YY_2LatExchange6Ngbr::UpdateDerivedOutputs);
    run_maxspinangle_output.Register(director,0);
  }
  catch(...) {
    if(coef_size>0 && coef12!=NULL) { // Release coef memory
      delete[] coef12[0]; delete[] coef12;
      coef12 = NULL;
      coef_size = 0; // Safety
    }
    throw;
  }

} // End Constructor

YY_2LatExchange6Ngbr::~YY_2LatExchange6Ngbr()
{
  if(coef_size>0 && coef1!=NULL) {
    delete[] coef1[0];
    delete[] coef1;
    delete[] coef2[0];
    delete[] coef2;
    delete[] coef12[0];
    delete[] coef12;
  }
}

OC_BOOL YY_2LatExchange6Ngbr::Init()
{
  mesh_id = 0;
  region_id.Release();
  J01.Release(); J02.Release();
  J012.Release(); J021.Release();
  mu1.Release(); mu2.Release();
  Tc1.Release(); Tc2.Release();
  m_e1.Release(); m_e2.Release();
  chi_l1.Release(); chi_l2.Release();
  G1.Release(); G2.Release();
  Lambdai11.Release(); Lambdai12.Release();
  Lambdai21.Release(); Lambdai22.Release();
  return Oxs_Energy::Init();
}

void YY_2LatExchange6Ngbr::CalcEnergyA
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,
 OC_INDEX node_stop,
 int threadnumber
 ) const
{
  // Depending on the lattice type, specify the coefficients to be used.
  // Sucscript A corresponds to this lattice, B is the other.
  const Oxs_MeshValue<ThreeVector> *spinA, *spinB;
  OC_REAL8m** coefA;
  OC_REAL8m** coefB;
  const Oxs_MeshValue<OC_REAL8m> *MsA, *MsB;
  const Oxs_MeshValue<OC_REAL8m> *MsA_inverse, *MsB_inverse;
  const Oxs_MeshValue<OC_REAL8m> *Ms0A_inverse, *Ms0B_inverse;
  const Oxs_MeshValue<OC_REAL8m> *J0A, *J0AB;
  const Oxs_MeshValue<OC_REAL8m> *muA, *muB;
  const Oxs_MeshValue<OC_REAL8m> *m_eA, *m_eB;
  const Oxs_MeshValue<OC_REAL8m> *chi_lA, *chi_lB;
  const Oxs_MeshValue<OC_REAL8m> *GA, *GB;
  switch(state.lattice_type) {
  case Oxs_SimState::TOTAL:
    throw Oxs_ExtError(this, "Programming error: CalcEnergyA was called"
        " with a wrong type of simulation state with lattice_type ="
        " TOTAL.");
    break;
  case Oxs_SimState::LATTICE1:
    spinA = &(state.spin);
    spinB = &(state.lattice2->spin);
    MsA = state.Ms;
    MsB = state.lattice2->Ms;
    MsA_inverse = state.Ms_inverse;
    MsB_inverse = state.lattice2->Ms_inverse;
    Ms0A_inverse = state.Ms0_inverse;
    Ms0B_inverse = state.lattice2->Ms0_inverse;
    coefA = coef1; coefB = coef2;
    muA = &mu1; muB = &mu2;
    J0A = &J01;
    J0AB = &J012;
    m_eA = &m_e1; m_eB = &m_e2;
    chi_lA = &chi_l1; chi_lB = &chi_l2;
    GA = &G1; GB = &G2;
    break;
  case Oxs_SimState::LATTICE2:
    spinA = &(state.spin);
    spinB = &(state.lattice1->spin);
    MsA = state.Ms;
    MsB = state.lattice1->Ms;
    MsA_inverse = state.Ms_inverse;
    MsB_inverse = state.lattice1->Ms_inverse;
    Ms0A_inverse = state.Ms0_inverse;
    Ms0B_inverse = state.lattice1->Ms0_inverse;
    coefA = coef2;
    coefB = coef1;
    muA = &mu2; muB = &mu1;
    J0A = &J02;
    J0AB = &J021;
    m_eA = &m_e2; m_eB = &m_e1;
    chi_lA = &chi_l2; chi_lB = &chi_l1;
    GA = &G2; GB = &G1;
    break;
  }

  // Downcast mesh
  const Oxs_CommonRectangularMesh* mesh
    = dynamic_cast<const Oxs_CommonRectangularMesh*>(state.mesh);
  if(mesh==NULL) {
    String msg =
      String("Import mesh (\"")
      + String(state.mesh->InstanceName())
      + String("\") to YY_2LatExchange6Ngbr::GetEnergyA()"
             " routine of object \"") + String(InstanceName())
      + String("\" is not a rectangular mesh object.");
    throw Oxs_ExtError(msg.c_str());
  }

  // If periodic, collect data for distance determination
  // Periodic boundaries?
  int xperiodic=0, yperiodic=0, zperiodic=0;
  const Oxs_PeriodicRectangularMesh* pmesh
    = dynamic_cast<const Oxs_PeriodicRectangularMesh*>(mesh);
  if(pmesh!=NULL) {
    xperiodic = pmesh->IsPeriodicX();
    yperiodic = pmesh->IsPeriodicY();
    zperiodic = pmesh->IsPeriodicZ();
  }

  OC_INDEX xdim = mesh->DimX();
  OC_INDEX ydim = mesh->DimY();
  OC_INDEX zdim = mesh->DimZ();
  OC_INDEX xydim = xdim*ydim;
  OC_INDEX xyzdim = xdim*ydim*zdim;

  OC_REAL8m wgtx = -1.0/(mesh->EdgeLengthX()*mesh->EdgeLengthX());
  OC_REAL8m wgty = -1.0/(mesh->EdgeLengthY()*mesh->EdgeLengthY());
  OC_REAL8m wgtz = -1.0/(mesh->EdgeLengthZ()*mesh->EdgeLengthZ());

  OC_REAL8m hcoef_t = -2/MU0;

  Nb_Xpfloat energy_sum = 0;
  OC_REAL8m thread_maxdot = maxdot[threadnumber];
  // Note: For maxangle calculation, it suffices to check
  // spin[j]-spin[i] for j>i.

  OC_INDEX x,y,z;
  mesh->GetCoords(node_start,x,y,z);

  OC_INDEX i = node_start;
  while(i<node_stop) {
    OC_INDEX xstop = xdim;
    if(xdim-x>node_stop-i) xstop = x + (node_stop-i);
    while(x<xstop) {
      ThreeVector baseA = (*spinA)[i];
      ThreeVector baseB = (*spinB)[i];
      OC_REAL8m MsiA = (*MsA)[i];
      OC_REAL8m MsiiA = (*MsA_inverse)[i];
      if(0.0 == MsiiA) {
        if(ocedt.energy) (*ocedt.energy)[i] = 0.0;
        if(ocedt.H)      (*ocedt.H)[i].Set(0.,0.,0.);
        if(ocedt.mxH)    (*ocedt.mxH)[i].Set(0.,0.,0.);
        ++i;   ++x;
        continue;
      }
      OC_REAL8m* ArowA = coefA[region_id[i]];
      OC_REAL8m* ArowB = coefB[region_id[i]];
      OC_REAL8m* ArowAB = coef12[region_id[i]];
      OC_REAL8m miA = (*MsA)[i]*(*Ms0A_inverse)[i];
      ThreeVector sum(0.,0.,0.);
      ThreeVector sum_l(0.,0.,0.);

      {
        // Exchange with the other sublattice
        OC_REAL8m LambdaiAA = (1.0+(*GB)[i])/(*chi_lA)[i];
        OC_REAL8m LambdaiAB = fabs((*m_eB)[i])
          *fabs((*J0AB)[i])/((*m_eA)[i]*(*muA)[i]*MU0);
        ThreeVector PB = baseA^baseB;
        PB ^= baseA;
        PB *= (*MsB)[i]*(*Ms0B_inverse)[i]; // PB = -nA x (nA x mB)
        OC_REAL8m tauB = abs(baseB*baseA);  // Dot product
        tauB *= (*MsB)[i]*(*Ms0B_inverse)[i]; // |tauB| = mB dot nA

        if((*state.T)[i]!=0.0) {
          OC_REAL8m beta = 1.0/(KB*(*state.T)[i]);
          OC_REAL8m A_AA = beta*fabs((*J0A)[i]);
          OC_REAL8m A_AB = beta*fabs((*J0AB)[i]);
          OC_REAL8m B = Langevin(A_AA*miA+A_AB*tauB);
          OC_REAL8m dB = LangevinDeriv(A_AA*miA+A_AB*tauB);
          if( (*state.T)[i]<(*state.Tc)[i] || 1-B/miA>0 ) {
            // This condition prevents unstable bursts of spin polarizations
            // when m ~ 0.
            sum_l -= (1-B/miA)/(MU0*(*muA)[i]*beta*dB)*baseA;
          }
        }

        sum += (*J0AB)[i]/(*muA)[i]*PB;
        sum *= -0.5*MsiA;  // This will be divided by MsiA at the end.
      }

      if(z>0 || zperiodic) {
        OC_INDEX j = i-xydim;
        if(z==0) j += xyzdim;
        OC_REAL8m ApairA = ArowA[region_id[j]];
        if(ApairA!=0 && (*MsA_inverse)[j]!=0.0) {
          ThreeVector diffA = ((*MsA)[j]*(*Ms0A_inverse)[j]*(*spinA)[j] 
              - miA*baseA);
          OC_REAL8m dot = diffA.MagSq();
          sum += ApairA*wgtz*diffA;
          if(dot>thread_maxdot) thread_maxdot = dot;
        }
      }
      if(y>0 || yperiodic) {
        OC_INDEX j = i-xdim;
        if(y==0) j += xydim;
        OC_REAL8m ApairA = ArowA[region_id[j]];
        if(ApairA!=0.0 && (*MsA_inverse)[j]!=0.0) {
          ThreeVector diffA = ((*MsA)[j]*(*Ms0A_inverse)[j]*(*spinA)[j] 
              - miA*baseA);
          OC_REAL8m dot = diffA.MagSq();
          sum += ApairA*wgty*diffA;
          if(dot>thread_maxdot) thread_maxdot = dot;
        }
      }
      if(x>0 || xperiodic) {
        OC_INDEX j = i-1;
        if(x==0) j += xdim;
        OC_REAL8m ApairA = ArowA[region_id[j]];
        if(ApairA!=0.0 && (*MsA_inverse)[j]!=0.0) {
          ThreeVector diffA = ((*MsA)[j]*(*Ms0A_inverse)[j]*(*spinA)[j] 
              - miA*baseA);
          OC_REAL8m dot = diffA.MagSq();
          sum += ApairA*wgtx*diffA;
          if(dot>thread_maxdot) thread_maxdot = dot;
        }
      }
      if(x<xdim-1 || xperiodic) {
        OC_INDEX j = i+1;
        if(x==xdim-1) j -= xdim;
        OC_REAL8m ApairA = ArowA[region_id[j]];
        if((*MsA_inverse)[j]!=0.0) {
          sum += ApairA*wgtx*( (*MsA)[j]*(*Ms0A_inverse)[j]*(*spinA)[j] 
              - miA*baseA );
        }
      }
      if(y<ydim-1 || yperiodic) {
        OC_INDEX j = i+xdim;
        if(y==ydim-1) j -= xydim;
        OC_REAL8m ApairA = ArowA[region_id[j]];
        if((*MsA_inverse)[j]!=0.0) {
          sum += ApairA*wgty*( (*MsA)[j]*(*Ms0A_inverse)[j]*(*spinA)[j] 
              - miA*baseA );
        }
      }
      if(z<zdim-1 || zperiodic) {
        OC_INDEX j = i+xydim;
        if(z==zdim-1) j -= xyzdim;
        OC_REAL8m ApairA = ArowA[region_id[j]];
        if((*MsA_inverse)[j]!=0.0) {
          sum += ApairA*wgtz*( (*MsA)[j]*(*Ms0A_inverse)[j]*(*spinA)[j] 
              - miA*baseA );
        }
      }

      OC_REAL8m ei = baseA.x*sum.x + baseA.y*sum.y + baseA.z*sum.z;
      OC_REAL8m hmult = hcoef_t*MsiiA;
      sum.x *= hmult;  sum.y *= hmult;   sum.z *= hmult;
      sum += sum_l;
      OC_REAL8m tx = baseA.y*sum.z - baseA.z*sum.y;
      OC_REAL8m ty = baseA.z*sum.x - baseA.x*sum.z;
      OC_REAL8m tz = baseA.x*sum.y - baseA.y*sum.x;

      energy_sum += ei;
      if(ocedt.energy)       (*ocedt.energy)[i] = ei;
      if(ocedt.energy_accum) (*ocedt.energy_accum)[i] += ei;
      if(ocedt.H)       (*ocedt.H)[i] = sum;
      if(ocedt.H_accum) (*ocedt.H_accum)[i] += sum;
      if(ocedt.mxH)       (*ocedt.mxH)[i] = ThreeVector(tx,ty,tz);
      if(ocedt.mxH_accum) (*ocedt.mxH_accum)[i] += ThreeVector(tx,ty,tz);
      ++i;   ++x;
    }
    x=0;
    if((++y)>=ydim) {
      y=0;
      ++z;
    }
  }

  ocedtaux.energy_total_accum += energy_sum.GetValue() * mesh->Volume(0);
  /// All cells have same volume in an Oxs_RectangularMesh.

  maxdot[threadnumber] = thread_maxdot;
}


void YY_2LatExchange6Ngbr::ComputeEnergyChunkInitialize
(const Oxs_SimState& state, // One of the sublattice states
 Oxs_ComputeEnergyDataThreaded& /* ocedt */,
 Oxs_ComputeEnergyDataThreadedAux& /* ocedtaux */,
 int number_of_threads) const
{
  // Update temperature-related parameters at the new stage
  if(last_stage_number==-1) { // First go
    last_stage_number = state.stage_number;

    const Oxs_Mesh* mesh = state.mesh;
    const OC_INDEX size = mesh->Size();
    Tc1.AdjustSize(mesh);
    Tc2.AdjustSize(mesh);
    J01_init->FillMeshValue(mesh,J01);
    J02_init->FillMeshValue(mesh,J02);
    J012_init->FillMeshValue(mesh,J012);
    J021_init->FillMeshValue(mesh,J021);
    mu1_init->FillMeshValue(mesh,mu1);
    mu2_init->FillMeshValue(mesh,mu2);
    m_e1.AdjustSize(mesh);
    m_e2.AdjustSize(mesh);
    chi_l1.AdjustSize(mesh);
    chi_l2.AdjustSize(mesh);
    G1.AdjustSize(mesh);
    G2.AdjustSize(mesh);
    Lambdai11.AdjustSize(mesh);
    Lambdai12.AdjustSize(mesh);
    Lambdai21.AdjustSize(mesh);
    Lambdai22.AdjustSize(mesh);

    // Set pointers for the temperature-dependent parameters in states.
    switch(state.lattice_type) {
    case Oxs_SimState::LATTICE1:
      state.m_e = &m_e1;
      state.lattice2->m_e = &m_e2;
      state.chi_l = &chi_l1;
      state.lattice2->chi_l = &chi_l2;
      state.Tc = &Tc1;
      state.lattice2->Tc = &Tc2;
      break;
    case Oxs_SimState::LATTICE2:
      state.lattice1->m_e = &m_e1;
      state.m_e = &m_e2;
      state.lattice1->chi_l = &chi_l1;
      state.chi_l = &chi_l2;
      state.lattice1->Tc = &Tc1;
      state.Tc = &Tc2;
      break;
    default:
      // Program should not reach here.
      break;
    }

    // Tc is also updated after calculation of m_e, chi_l
    Update_m_e(*(state.total_lattice), 1e-4);
  }

  if(state.stage_number != last_stage_number) { // New stage
    last_stage_number = state.stage_number;
    Update_m_e(*(state.total_lattice), 1e-4);
  }

  // chi_l depends on instantaneous magnetization 
  // so calculate it at each step
  Update_chi_l(*(state.total_lattice));

  if(maxdot.size() != (vector<OC_REAL8m>::size_type)number_of_threads) {
    maxdot.resize(number_of_threads);
  }
  for(int i=0;i<number_of_threads;++i) {
    maxdot[i] = 0.0; // Minimum possible value for (m_i-m_j).MagSq()
  }
}


void YY_2LatExchange6Ngbr::ComputeEnergyChunkFinalize
(const Oxs_SimState& state,
 const Oxs_ComputeEnergyDataThreaded& /* ocedt */,
 const Oxs_ComputeEnergyDataThreadedAux& /* ocedtaux */,
 int number_of_threads) const
{
  // Set max angle data
  OC_REAL8m total_maxdot = 0.0;
  for(int i=0;i<number_of_threads;++i) {
    if(maxdot[i]>total_maxdot) total_maxdot = maxdot[i];
  }
  const OC_REAL8m arg = 0.5*Oc_Sqrt(total_maxdot);
  const OC_REAL8m maxang = ( arg >= 1.0 ? 180.0 : asin(arg)*(360.0/PI));

  OC_REAL8m dummy_value;
  String msa_name = MaxSpinAngleStateName();
  if(state.GetDerivedData(msa_name,dummy_value)) {
    // Ideally, energy values would never be computed more than once
    // for any one state, but in practice it seems inevitable that
    // such will occur on occasion.  For example, suppose a user
    // requests output on a state obtained by a stage crossing (as
    // opposed to a state obtained through a normal intrastage step);
    // a subsequent ::Step operation will re-compute the energies
    // because not all the information needed by the step transition
    // machinery is cached from an energy computation.  Even user
    // output requests on post-::Step states is problematic if some of
    // the requested output is not cached as part of the step
    // proceedings.  A warning is put into place below for debugging
    // purposes, but in general an error is raised only if results
    // from the recomputation are different than originally.
#ifndef NDEBUG
    //static Oxs_WarningMessage maxangleset(3);
    //maxangleset.Send(revision_info,OC_STRINGIFY(__LINE__),
    //                 "Programming error?"
    //                 " YY_2LatExchange6Ngbr max spin angle set twice.");
#endif
    // Max angle is computed by taking acos of the dot product
    // of neighboring spin vectors.  The relative error can be
    // quite large if the spins are nearly parallel.  The proper
    // error comparison is between the cos of the two values.
    // See NOTES VI, 6-Sep-2012, p71.
    OC_REAL8m diff = (dummy_value-maxang)*(PI/180.);
    OC_REAL8m sum  = (dummy_value+maxang)*(PI/180.);
    if(sum > PI ) sum = 2*PI - sum;
    if(fabs(diff*sum)>8*OC_REAL8_EPSILON) {
      char errbuf[1024];
      Oc_Snprintf(errbuf,sizeof(errbuf),
                  "Programming error:"
                  " YY_2LatExchange6Ngbr max spin angle set to"
                  " two different values;"
                  " orig val=%#.17g, new val=%#.17g",
                  dummy_value,maxang);
      throw Oxs_ExtError(this,errbuf);
    }
  } else {
    state.AddDerivedData(msa_name,maxang);
  }

  // Run and stage angle data depend on data from the previous state.
  // In the case that the energy (and hence max stage and run angle)
  // for the current state was computed previously, then the previous
  // state may have been dropped.  So, compute and save run and stage
  // angle data iff they are not already computed.

  // Check stage and run max angle data from previous state
  const Oxs_SimState* oldstate = NULL;
  OC_REAL8m stage_maxang = -1;
  OC_REAL8m run_maxang = -1;
  String smsa_name = StageMaxSpinAngleStateName();
  String rmsa_name = RunMaxSpinAngleStateName();
  if(state.previous_state_id &&
     (oldstate
      = director->FindExistingSimulationState(state.previous_state_id)) ) {
    if(oldstate->stage_number != state.stage_number) {
      stage_maxang = 0.0;
    } else {
      if(oldstate->GetDerivedData(smsa_name,dummy_value)) {
        stage_maxang = dummy_value;
      }
    }
    if(oldstate->GetDerivedData(rmsa_name,dummy_value)) {
      run_maxang = dummy_value;
    }
  }
  if(stage_maxang<maxang) stage_maxang = maxang;
  if(run_maxang<maxang)   run_maxang = maxang;

  // Stage max angle data
  if(!state.GetDerivedData(smsa_name,dummy_value)) {
    state.AddDerivedData(smsa_name,stage_maxang);
  }

  // Run max angle data
  if(!state.GetDerivedData(rmsa_name,dummy_value)) {
    state.AddDerivedData(rmsa_name,run_maxang);
  }
}


void YY_2LatExchange6Ngbr::ComputeEnergyChunk
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,
 OC_INDEX node_stop,
 int threadnumber
 ) const
{
#ifndef NDEBUG
  if(node_stop>state.mesh->Size() || node_start>node_stop) {
    throw Oxs_ExtError(this,"Programming error:"
                       " Invalid node_start/node_stop values");
  }
#endif

  const OC_INDEX size = state.mesh->Size();
  if(size<1) {
    return;
  }

  if(mesh_id !=  state.mesh->Id() || !atlaskey.SameState()) {
    // Setup region mapping.
    // NB: At a lower level, this may potentially involve calls back
    // into the Tcl interpreter.  Per Tcl spec, only the thread
    // originating the interpreter is allowed to make calls into it, so
    // only threadnumber == 0 can do this processing.  Any other thread
    // must block until that processing is complete.
    thread_control.Lock();
    if(Oxs_ThreadError::IsError()) {
      if(thread_control.count>0) {
        // Release a blocked thread
        thread_control.Notify();
      }
      thread_control.Unlock();
      return; // What else?
    }
    if(threadnumber != 0) {
      if(mesh_id != state.mesh->Id() || !atlaskey.SameState()) {
        // If above condition is false, then the main thread came
        // though and initialized everything between the time of
        // the previous check and this thread's acquiring of the
        // thread_control mutex; in which case, "never mind".
        // Otherwise:
        ++thread_control.count; // Multiple threads may progress to this
        /// point before the main thread (threadnumber == 0) grabs the
        /// thread_control mutex.  Keep track of how many, so that
        /// afterward they may be released, one by one.  (The main
        /// thread will Notify control_wait.cond once; after that
        /// as each waiting thread is released, the newly released
        /// thread sends a Notify to wake up the next one.
        thread_control.Wait(0);
        --thread_control.count;
        int condcheckerror=0;
        if(mesh_id !=  state.mesh->Id() || !atlaskey.SameState()) {
          // Error?
          condcheckerror=1;
          Oxs_ThreadPrintf(stderr,"Invalid condition in"
                           " YY_2LatExchange6Ngbr::ComputeEnergyChunk(),"
                           " thread number %d\n",threadnumber);
        }
        if(thread_control.count>0) {
          // Free a waiting thread.
          thread_control.Notify();
        }
        thread_control.Unlock();
        if(condcheckerror || Oxs_ThreadError::IsError()) {
          return; // What else?
        }
      } else {
        if(thread_control.count>0) {
          // Free a waiting thread.  (Actually, it can occur that the
          // thread_control will be grabbed by another thread that is
          // blocked at the first thread_control mutex Lock() call above
          // rather than on the ConditionWait, in which case this
          // ConditionNotify will be effectively lost.  But that is
          // okay, because then *that* thread will Notify when it
          // releases the mutex.)
          thread_control.Notify();
        }
        thread_control.Unlock();
      }
    } else {
      // Main thread (threadnumber == 0)
      try {
        region_id.AdjustSize(state.mesh);
        ThreeVector location;
        for(OC_INDEX i=0;i<size;i++) {
          state.mesh->Center(i,location);
          if((region_id[i] = atlas->GetRegionId(location))<0) {
            String msg = String("Import mesh to YY_2LatExchange6Ngbr::GetEnergy()"
                                " routine of object ")
              + String(InstanceName())
                + String(" has points outside atlas ")
              + String(atlas->InstanceName());
            throw Oxs_ExtError(msg.c_str());
          }
        }
        atlaskey.Set(atlas.GetPtr());
        mesh_id = state.mesh->Id();
      } catch(Oxs_ExtError& err) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String(err));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(String& serr) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(serr);
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(const char* cerr) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String(cerr));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(...) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String("Error in "
            "YY_2LatExchange6Ngbr::ComputeEnergyChunk"));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      }
      if(thread_control.count>0) {
        // Free a waiting thread.  (Actually, it can occur that the
        // thread_control will be grabbed by another thread that is
        // blocked at the first thread_control mutex Lock() call above
        // rather than on the ConditionWait, in which case this
        // ConditionNotify will be effectively lost.  But that is
        // okay, because then *that* thread will Notify when it
        // releases the mutex.)
        thread_control.Notify();
      }
      thread_control.Unlock();
    }
  }
  if(excoeftype == LEX_TYPE) {
    //CalcEnergyLex(state,ocedt,ocedtaux,node_start,node_stop,threadnumber);
  } else {
    CalcEnergyA(state,ocedt,ocedtaux,node_start,node_stop,threadnumber);
  }
}

void YY_2LatExchange6Ngbr::Update_m_e(
    const Oxs_SimState& state,  // Total lattice state
    OC_REAL8m tol_in = 1e-4) const
{
  // Solve for the equilibrium spin polarization m_e using 2 variable
  // Newton method.
  const OC_REAL8m size = state.mesh->Size();
  tol = fabs(tol_in);
  tolsq = tol_in*tol_in;

  Oxs_MeshValue<OC_REAL8m>& Ms1 = *(state.lattice1->Ms);
  Oxs_MeshValue<OC_REAL8m>& Ms2 = *(state.lattice2->Ms);
  Oxs_MeshValue<OC_REAL8m>& Ms01_inverse = *(state.lattice1->Ms0_inverse);
  Oxs_MeshValue<OC_REAL8m>& Ms02_inverse = *(state.lattice2->Ms0_inverse);

  OC_REAL8m A11, A12, A21, A22;
  OC_REAL8m x10, x20, y10, y20;
  OC_REAL8m x1, x2, y1, y2;
  OC_REAL8m dx1, dx2, dy1, dy2;
  OC_REAL8m dL1, dL2;
  OC_REAL8m J11, J12, J21, J22, det, deti;
  for(OC_INDEX i=0; i<size; i++) {
    const OC_REAL8m kB_T = KB*(*(state.lattice1->T))[i];
    if(kB_T == 0) {
      m_e1[i]=1.0;
      m_e2[i]=1.0;

      Tc1[i] = (J01[i]+fabs(J012[i]))/(3*KB);
      Tc2[i] = (J02[i]+fabs(J021[i]))/(3*KB);
      continue;
    }
    const OC_REAL8m beta = 1.0/kB_T;
    A11 = beta*J01[i];
    A12 = beta*fabs(J012[i]);
    A21 = beta*fabs(J021[i]);
    A22 = beta*J02[i];

    x1 = 0.8; x2 = 0.8;
    y1 = Langevin(A11*x1+A12*x2)-x1;
    y2 = Langevin(A21*x1+A22*x2)-x2;
    do {
      x10 = x1; x20 = x2; y10 = y1; y20 = y2;
      dL1 = LangevinDeriv(A11*x1+A12*x2);
      dL2 = LangevinDeriv(A21*x1+A22*x2);
      // Jacobian
      J11 = A11*dL1-1;
      J12 = A12*dL1;
      J21 = A21*dL2;
      J22 = A22*dL2-1;
      det = J11*J22-J12*J21;
      if(det == 0.0) {
        // No more change. Calculate parameters with current x1, x2.
        break;
      }
      deti = 1.0/det;
      dx1 = -deti*(J22*y1-J12*y2);
      dx2 = -deti*(-J21*y1+J11*y2);
      x1 = x10 + dx1;
      x2 = x20 + dx2;
      y1 = Langevin(A11*x1+A12*x2)-x1;
      y2 = Langevin(A21*x1+A22*x2)-x2;
    } while( dx1*dx1>tolsq || dx2*dx2>tolsq );
    m_e1[i] = x1>tol ? x1 : 0.0;
    m_e2[i] = x2>tol ? x2 : 0.0;

    OC_REAL8m m1 = Ms1[i]*Ms01_inverse[i];
    OC_REAL8m m2 = Ms2[i]*Ms02_inverse[i];

    if(m_e1[i]>tol && m_e2[i]>tol) {
      Tc1[i] = (J01[i]*m1+fabs(J012[i])*m2)/(3*KB*m1);
      Tc2[i] = (J02[i]*m2+fabs(J021[i])*m1)/(3*KB*m2);
    } else {
      // MFA for T>Tc. See Eq. (10) in Ostler, PRB 84, 024407 (2011).
      const OC_REAL8m a = J012[i]*J021[i]-J01[i]*J02[i];
      const OC_REAL8m b = J01[i]+J02[i];
      Tc1[i] = 2*a/( 3*KB*b*( sqrt(1+4*a/(b*b))-1 ) );
      Tc2[i] = 2*a/( 3*KB*b*( sqrt(1+4*a/(b*b))-1 ) );
    }
  }
}

OC_REAL8m YY_2LatExchange6Ngbr::Langevin(OC_REAL8m x) const
{
  OC_REAL8m temp = exp(2*x)+1;
  if(!Nb_IsFinite(temp)) return x>0 ? 1.0 : -1.0; // Large input
  temp /= temp-2; // temp = coth(x);
  return temp-1/x;
}

OC_REAL8m YY_2LatExchange6Ngbr::LangevinDeriv(OC_REAL8m x) const
{
  if(fabs(x)<OC_REAL4_EPSILON) return 1./3.;
  OC_REAL8m temp = sinh(x);
  if(!Nb_IsFinite(temp)) return 1.0/(x*x); // Large input
  return -1.0/(temp*temp)+1.0/(x*x);
}

void YY_2LatExchange6Ngbr::Update_chi_l(const Oxs_SimState& state) const
{
  const OC_REAL8m size = state.mesh->Size();
  const Oxs_MeshValue<OC_REAL8m>& Ms1 = *(state.lattice1->Ms);
  const Oxs_MeshValue<OC_REAL8m>& Ms2 = *(state.lattice2->Ms);
  const Oxs_MeshValue<OC_REAL8m>& Ms01_inverse = *(state.lattice1->Ms0_inverse);
  const Oxs_MeshValue<OC_REAL8m>& Ms02_inverse = *(state.lattice2->Ms0_inverse);

  OC_REAL8m kB_T, beta;
  OC_REAL8m A11, A12, A21, A22;
  OC_REAL8m m1, m2, dL1, dL2;

  for(OC_INDEX i=0; i<size; i++) {
    kB_T = KB*(*(state.lattice1->T))[i];
    if(kB_T == 0) {
      dL1 = dL2 = 1.0;
      G1[i] = ( fabs(J021[i])*fabs(J012[i]) - mu1[i]/mu2[i]*fabs(J021[i])*J02[i] )
        /( J01[i]*J02[i] - fabs(J021[i]*J012[i]) );
      G2[i] = ( fabs(J021[i])*fabs(J012[i]) - mu2[i]/mu1[i]*fabs(J012[i])*J01[i] )
        /( J01[i]*J02[i] - fabs(J021[i]*J012[i]) );
      chi_l1[i] = MU0*mu2[i]/fabs(J021[i])*G1[i];
      chi_l2[i] = MU0*mu1[i]/fabs(J012[i])*G2[i];

    } else {
      beta = 1.0/kB_T;
      A11 = beta*J01[i];
      A12 = beta*fabs(J012[i]);
      A21 = beta*fabs(J021[i]);
      A22 = beta*J02[i];
      m1 = Ms1[i]*Ms01_inverse[i];
      m2 = Ms2[i]*Ms02_inverse[i];
      dL1 = LangevinDeriv(A11*m1+A12*m2);
      dL2 = LangevinDeriv(A21*m1+A22*m2);

      G1[i] = ( A21*A12*dL1*dL2 + (mu1[i]/mu2[i])*A21*dL1*(1-A22*dL2) )
        /( (1-A11*dL1)*(1-A22*dL2) - A21*A12*dL1*dL2 );
      G2[i] = ( A21*A12*dL1*dL2 + (mu2[i]/mu1[i])*A12*dL2*(1-A11*dL1) )
        /( (1-A11*dL1)*(1-A22*dL2) - A21*A12*dL1*dL2 );

      chi_l1[i] = MU0*mu2[i]/fabs(J021[i])*G1[i];
      chi_l2[i] = MU0*mu1[i]/fabs(J012[i])*G2[i];
    }
  }
}

void YY_2LatExchange6Ngbr::UpdateDerivedOutputs(const Oxs_SimState& state)
{
  maxspinangle_output.cache.state_id
    = stage_maxspinangle_output.cache.state_id
    = run_maxspinangle_output.cache.state_id
    = 0;  // Mark change in progress

  String dummy_name = MaxSpinAngleStateName();
  if(!state.GetDerivedData(dummy_name.c_str(),
                           maxspinangle_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    maxspinangle_output.cache.value = -1.0;
  }

  String stage_dummy_name = StageMaxSpinAngleStateName();
  if(!state.GetDerivedData(stage_dummy_name.c_str(),
                           stage_maxspinangle_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    stage_maxspinangle_output.cache.value = -1.0;
  }

  String run_dummy_name = RunMaxSpinAngleStateName();
  if(!state.GetDerivedData(run_dummy_name.c_str(),
                           run_maxspinangle_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    run_maxspinangle_output.cache.value = -1.0;
  }

  maxspinangle_output.cache.state_id
    = stage_maxspinangle_output.cache.state_id
    = run_maxspinangle_output.cache.state_id
    = state.Id();
}
