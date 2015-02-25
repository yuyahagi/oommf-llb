/** FILE: yy_2lat_util.cc                 -*-Mode: c++-*-
 *
 * OOMMF 2 lattice extension utilities.
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

#include <assert.h>
#include <string>
#include "chunkenergy.h"
#include "energy.h"
#include "mesh.h"

#include "yy_2lat_util.h"

struct Oxs_ComputeEnergies_ChunkStruct {
public:
  Oxs_ChunkEnergy* energy;
  Oxs_ComputeEnergyDataThreaded ocedt;
  Oxs_ComputeEnergyDataThreadedAux ocedtaux;
  Oxs_ComputeEnergies_ChunkStruct()
    : energy(0) {}
};

class Oxs_ComputeEnergiesChunkThread : public Oxs_ThreadRunObj {
public:
  static Oxs_JobControl<ThreeVector> job_basket;
  /// job_basket is static, so only one "set" of this class is allowed.

  const Oxs_SimState* state;
  vector<Oxs_ComputeEnergies_ChunkStruct> energy_terms;

  Oxs_MeshValue<ThreeVector>* mxH;
  Oxs_MeshValue<ThreeVector>* mxH_accum;
  Oxs_MeshValue<ThreeVector>* mxHxm;
  const vector<OC_INDEX>* fixed_spins;
  OC_REAL8m max_mxH;

  OC_INDEX cache_blocksize;

  OC_BOOL accums_initialized;

  Oxs_ComputeEnergiesChunkThread()
    : state(0),
      mxH(0),mxH_accum(0),
      mxHxm(0), fixed_spins(0),
      max_mxH(0.0),
      cache_blocksize(0), accums_initialized(0) {}

  void Cmd(int threadnumber, void* data);

  static void Init(int thread_count,
                   const Oxs_StripedArray<ThreeVector>* arrblock) {
    job_basket.Init(thread_count,arrblock);
  }

  // Note: Default copy constructor and assignment operator,
  // and destructor.
};

void YY_2LatComputeEnergies(
    const Oxs_SimState& state,
    Oxs_ComputeEnergyData& oced1,
    Oxs_ComputeEnergyData& oced2,
    const vector<Oxs_Energy*>& energies,
    Oxs_ComputeEnergyExtraData& oceed)
{

  if(state.lattice_type != Oxs_SimState::TOTAL) {
    String msg = String("Programming error:"
        " YY_2LatComputeEnergies was called with non-TOTAL state.");
    throw Oxs_ExtError(msg);
  }

  if(state.Id()==0) {
    String msg = String("Programming error:"
                        " Invalid (unlocked) state detected"
                        " in YY_2LatComputeEnergies");
    throw Oxs_ExtError(msg);
  }

  if(oced1.scratch_energy==NULL || oced1.scratch_H==NULL) {
    // Bad input
    String msg = String("Oxs_ComputeEnergyData object in function"
                        " YY_2LatComputeEnergies"
                        " contains NULL scratch pointers.");
    throw Oxs_ExtError(msg);
  }

  if(oced1.energy != NULL || oced1.H != NULL || oced1.mxH != NULL) {
    String msg = String("Programming error in function"
                        " YY_2LatComputeEnergies:"
                        " non-NULL energy, H, and/or mxH imports.");
    throw Oxs_ExtError(msg);
  }

  const Oxs_SimState& state1 = *(state.lattice1);
  const Oxs_SimState& state2 = *(state.lattice2);

  const int thread_count = Oc_GetMaxThreadCount();

  if(oced1.energy_accum) {
    oced1.energy_accum->AdjustSize(state.mesh);
  }
  if(oced2.energy_accum) {
    oced2.energy_accum->AdjustSize(state.mesh);
  }
  if(oced1.H_accum) {
    oced1.H_accum->AdjustSize(state.mesh);
  }
  if(oced2.H_accum) {
    oced2.H_accum->AdjustSize(state.mesh);
  }
  if(oced1.mxH_accum) {
    oced1.mxH_accum->AdjustSize(state.mesh);
  }
  if(oced2.mxH_accum) {
    oced2.mxH_accum->AdjustSize(state.mesh);
  }
  if(oceed.mxHxm) {
    oceed.mxHxm->AdjustSize(state.mesh);
  }

  oced1.energy_sum = 0.0;
  oced2.energy_sum = 0.0;
  oced1.pE_pt = 0.0;
  oced2.pE_pt = 0.0;
  oceed.max_mxH = 0.0;

  if(energies.size() == 0) {
    // No energies.  Zero requested outputs and return.
    OC_INDEX size = state.mesh->Size();
    if(oced1.energy_accum) {
      for(OC_INDEX i=0; i<size; ++i) (*(oced1.energy_accum))[i] = 0.0;
    }
    if(oced2.energy_accum) {
      for(OC_INDEX i=0; i<size; ++i) (*(oced2.energy_accum))[i] = 0.0;
    }
    if(oced1.H_accum) {
      for(OC_INDEX i=0; i<size; ++i) (*(oced1.H_accum))[i] = ThreeVector(0.0,0.0,0.0);
    }
    if(oced2.H_accum) {
      for(OC_INDEX i=0; i<size; ++i) (*(oced2.H_accum))[i] = ThreeVector(0.0,0.0,0.0);
    }
    if(oced1.mxH_accum) {
      for(OC_INDEX i=0; i<size; ++i) (*(oced1.mxH_accum))[i] = ThreeVector(0.0,0.0,0.0);
    }
    if(oced2.mxH_accum) {
      for(OC_INDEX i=0; i<size; ++i) (*(oced2.mxH_accum))[i] = ThreeVector(0.0,0.0,0.0);
    }
    if(oceed.mxHxm) {
      for(OC_INDEX i=0; i<size; ++i) (*(oceed.mxHxm))[i] = ThreeVector(0.0,0.0,0.0);
    }
    return;
  }

  if(oced1.mxH_accum==0 && oceed.mxHxm!=0) {
    // Hack mxHxm into mxH_accum.  We can identify this situation
    // by checking mxH_accum == mxHxm, and undo at the end.  Also
    // The Oxs_ComputeEnergiesChunkThread objects know about this
    // and respond appropriately.
    oced1.mxH_accum = oceed.mxHxm;
  }
  if(oced2.mxH_accum==0 && oceed.mxHxm!=0) {
    oced2.mxH_accum = oceed.mxHxm;
  }


  vector<Oxs_ComputeEnergies_ChunkStruct> chunk1, chunk2;
  vector<Oxs_Energy*> nonchunk1, nonchunk2;

  // Initialize those parts of ChunkStruct that are independent
  // of any particular energy term.
  Oxs_ComputeEnergies_ChunkStruct foo1, foo2;
  foo1.ocedt.state_id = state.Id();
  foo1.ocedt.scratch_energy = oced1.scratch_energy;
  foo1.ocedt.scratch_H      = oced1.scratch_H;
  foo1.ocedt.energy_accum   = oced1.energy_accum;
  foo1.ocedt.H_accum        = oced1.H_accum;
  foo1.ocedt.mxH_accum      = oced1.mxH_accum;
  foo2.ocedt.state_id = state.Id();
  foo2.ocedt.scratch_energy = oced2.scratch_energy;
  foo2.ocedt.scratch_H      = oced2.scratch_H;
  foo2.ocedt.energy_accum   = oced2.energy_accum;
  foo2.ocedt.H_accum        = oced2.H_accum;
  foo2.ocedt.mxH_accum      = oced2.mxH_accum;
  for(vector<Oxs_Energy*>::const_iterator it = energies.begin();
      it != energies.end() ; ++it ) {
    Oxs_ChunkEnergy* ceptr =
      dynamic_cast<Oxs_ChunkEnergy*>(*it);
    if(ceptr != NULL) {
      // Set up and initialize chunk energy structures
      foo1.energy = ceptr;
      foo2.energy = ceptr;
      if(ceptr->energy_density_output.GetCacheRequestCount()>0) {
        ceptr->energy_density_output.cache.state_id=0;
        foo1.ocedt.energy = &(ceptr->energy_density_output.cache.value);
        foo2.ocedt.energy = &(ceptr->energy_density_output.cache.value);
        foo1.ocedt.energy->AdjustSize(state.mesh);
        foo2.ocedt.energy->AdjustSize(state.mesh);
      }
      if(ceptr->field_output.GetCacheRequestCount()>0) {
        ceptr->field_output.cache.state_id=0;
        foo1.ocedt.H = &(ceptr->field_output.cache.value);
        foo2.ocedt.H = &(ceptr->field_output.cache.value);
        foo1.ocedt.H->AdjustSize(state.mesh);
        foo2.ocedt.H->AdjustSize(state.mesh);
      }
      chunk1.push_back(foo1);
      chunk2.push_back(foo2);
    } else {
      nonchunk1.push_back(*it);
      nonchunk2.push_back(*it);
    }
  }

  // The "accum" elements are initialized on the first pass by
  // moving each accum pointer to the corresponding non-accum member.
  // After filling by the first energy term, the pointers are moved
  // back to the accum member.  This way we avoid a pass through
  // memory storing zeros, and a pass through memory loading zeros.
  // Zero load/stores are cheap in the chunk memory case, because
  // in that case the load/stores are just to and from cache, but
  // we prefer here to run non-chunk energies first so that we
  // can compute mxHxm and max |mxH| on the backsize of the chunk
  // energy runs (where m and mxH are in cache and so don't have
  // to be loaded).  Create a boolean to track initialization.
  OC_BOOL accums_initialized = 0;

  // Non-chunk energies //////////////////////////////////////
  for(vector<Oxs_Energy*>::const_iterator ncit = nonchunk1.begin();
      ncit != nonchunk1.end() ; ++ncit ) {
    Oxs_Energy& eterm = *(*ncit);  // Convenience

#if REPORT_TIME
    eterm.energytime.Start();
#endif // REPORT_TIME

    Oxs_ComputeEnergyData term_oced1(state1);
    term_oced1.scratch_energy = oced1.scratch_energy;
    term_oced1.scratch_H      = oced1.scratch_H;
    term_oced1.energy_accum = oced1.energy_accum;
    term_oced1.H_accum      = oced1.H_accum;
    term_oced1.mxH_accum    = oced1.mxH_accum;

    Oxs_ComputeEnergyData term_oced2(state2);
    term_oced2.scratch_energy = oced2.scratch_energy;
    term_oced2.scratch_H      = oced2.scratch_H;
    term_oced2.energy_accum = oced2.energy_accum;
    term_oced2.H_accum      = oced2.H_accum;
    term_oced2.mxH_accum    = oced2.mxH_accum;

    if(eterm.energy_density_output.GetCacheRequestCount()>0) {
      eterm.energy_density_output.cache.state_id=0;
      term_oced1.energy = &(eterm.energy_density_output.cache.value);
      term_oced2.energy = &(eterm.energy_density_output.cache.value);
      term_oced1.energy->AdjustSize(state.mesh);
      term_oced2.energy->AdjustSize(state.mesh);
    }

    if(eterm.field_output.GetCacheRequestCount()>0) {
      eterm.field_output.cache.state_id=0;
      term_oced1.H = &(eterm.field_output.cache.value);
      term_oced2.H = &(eterm.field_output.cache.value);
      term_oced1.H->AdjustSize(state.mesh);
      term_oced2.H->AdjustSize(state.mesh);
    }

    if(!accums_initialized) {
      // Initialize by filling
      term_oced1.energy_accum = 0;
      term_oced2.energy_accum = 0;
      term_oced1.H_accum = 0;
      term_oced2.H_accum = 0;
      term_oced1.mxH_accum = 0;
      term_oced2.mxH_accum = 0;
      if(term_oced1.energy == 0) term_oced1.energy = oced1.energy_accum;
      if(term_oced2.energy == 0) term_oced2.energy = oced2.energy_accum;
      if(term_oced1.H == 0)      term_oced1.H      = oced1.H_accum;
      if(term_oced2.H == 0)      term_oced2.H      = oced2.H_accum;
      if(term_oced1.mxH == 0)    term_oced1.mxH    = oced1.mxH_accum;
      if(term_oced2.mxH == 0)    term_oced2.mxH    = oced2.mxH_accum;
    }

    ++(eterm.calc_count);
    eterm.ComputeEnergy(state1,term_oced1);
    eterm.ComputeEnergy(state2,term_oced2);

    if(eterm.field_output.GetCacheRequestCount()>0) {
      eterm.field_output.cache.state_id=state.Id();
    }
    if(eterm.energy_density_output.GetCacheRequestCount()>0) {
      eterm.energy_density_output.cache.state_id=state.Id();
    }
    if(eterm.energy_sum_output.GetCacheRequestCount()>0) {
      eterm.energy_sum_output.cache.value=term_oced1.energy_sum;
      eterm.energy_sum_output.cache.state_id=state.Id();
    }

    if(!accums_initialized) {
      // If output buffer spaced was used instead of accum space, then
      // copy from output buffer to accum space.  This hurts from a
      // memory bandwidth perspective, but is rather hard to avoid.
      // (Options: Do accum initialization in chunk-energy branch,
      // but that hurts with respect to mxHxm and max |mxH| computations.
      // Or one could have the ComputeEnergy class fill more than one
      // array with the non-accum output (say, via a parameter that
      // says to set to accum rather than add to accum), but that is
      // rather awkward.  Instead, we assume that if the user wants
      // high speed then he won't enable term energy or H outputs.)
      if(oced1.energy_accum && term_oced1.energy != oced1.energy_accum) {
        *(oced1.energy_accum) = *(term_oced1.energy);
      }
      if(oced2.energy_accum && term_oced2.energy != oced2.energy_accum) {
        *(oced2.energy_accum) = *(term_oced2.energy);
      }
      if(oced1.H_accum      && term_oced1.H      != oced1.H_accum) {
        *(oced1.H_accum) = *(term_oced1.H);
      }
      if(oced2.H_accum      && term_oced2.H      != oced2.H_accum) {
        *(oced2.H_accum) = *(term_oced2.H);
      }
      if(oced1.mxH_accum    && term_oced1.mxH    != oced1.mxH_accum) {
        *(oced1.mxH_accum) = *(term_oced1.mxH);
      }
      if(oced2.mxH_accum    && term_oced2.mxH    != oced2.mxH_accum) {
        *(oced2.mxH_accum) = *(term_oced2.mxH);
      }
      accums_initialized = 1;
    }

    oced1.energy_sum += term_oced1.energy_sum;
    oced2.energy_sum += term_oced2.energy_sum;
    oced1.pE_pt += term_oced1.pE_pt;
    oced2.pE_pt += term_oced2.pE_pt;

#if REPORT_TIME
    eterm.energytime.Stop();
#endif // REPORT_TIME
  }


  // Chunk energies ///////////////////////////////////////////
#if REPORT_TIME
  Oxs_ChunkEnergy::chunktime.Start();
#endif

  // Compute cache_blocksize
  const OC_INDEX meshsize = state.mesh->Size();
  const OC_INDEX cache_size = 1024 * 1024;  // Should come from
  /// platform file or perhaps sysconf().

  const OC_INDEX recsize = sizeof(ThreeVector) + sizeof(OC_REAL8m);
  /// May want to query individual energies for this.

#define FUDGE 8
  OC_INDEX tcblocksize = (cache_size>FUDGE*recsize ?
                        cache_size/(FUDGE*recsize) : 1);
  if(thread_count*tcblocksize>meshsize) {
    tcblocksize = meshsize/thread_count;
  }
  if(0 == tcblocksize) {
    tcblocksize = 1;    // Safety
  } else if(0 != tcblocksize%16) {
    tcblocksize += 16 - (tcblocksize%16);  // Make multiple of 16
  }
  const OC_INDEX cache_blocksize = tcblocksize;

  // Thread control
  static Oxs_ThreadTree threadtree;

// =========================================================================
// Lattice 1
// =========================================================================
  Oxs_ComputeEnergiesChunkThread::Init(thread_count,
                                       state.spin.GetArrayBlock());

  vector<Oxs_ComputeEnergiesChunkThread> chunk_thread;
  chunk_thread.resize(thread_count);
  chunk_thread[0].state     = &state1;
  chunk_thread[0].energy_terms = chunk1; // Make copies.
  chunk_thread[0].mxH       = oced1.mxH;
  chunk_thread[0].mxH_accum = oced1.mxH_accum;
  chunk_thread[0].mxHxm     = oceed.mxHxm;
  chunk_thread[0].fixed_spins = oceed.fixed_spin_list;
  chunk_thread[0].cache_blocksize = cache_blocksize;
  chunk_thread[0].accums_initialized = accums_initialized;


  // Initialize chunk energy computations
  for(vector<Oxs_ComputeEnergies_ChunkStruct>::iterator it
        = chunk1.begin(); it != chunk1.end() ; ++it ) {
    Oxs_ChunkEnergy& eterm = *(it->energy);  // For code clarity
    Oxs_ComputeEnergyDataThreaded& ocedt = it->ocedt;
    Oxs_ComputeEnergyDataThreadedAux& ocedtaux = it->ocedtaux;
    eterm.ComputeEnergyChunkInitialize(state1,ocedt,ocedtaux,
                                       thread_count);
  }

  for(int ithread=1;ithread<thread_count;++ithread) {
    chunk_thread[ithread] = chunk_thread[0];
    threadtree.Launch(chunk_thread[ithread],0);
  }
  threadtree.LaunchRoot(chunk_thread[0],0);

  // Note: If chunk.size()>0, then we are guaranteed that accums are
  // initialized.  If accums_initialized is ever needed someplace
  // downstream, then uncomment the following line:
  // if(chunk.size()>0) accums_initialized = 1;

  // Finalize chunk energy computations
  for(OC_INDEX ei=0;static_cast<size_t>(ei)<chunk1.size();++ei) {

    Oxs_ChunkEnergy& eterm = *(chunk1[ei].energy);  // Convenience
    const Oxs_ComputeEnergyDataThreaded& ocedt = chunk1[ei].ocedt;
    const Oxs_ComputeEnergyDataThreadedAux& ocedtaux = chunk1[ei].ocedtaux;

    eterm.ComputeEnergyChunkFinalize(state1,ocedt,ocedtaux,
                                     thread_count);

    ++(eterm.calc_count);

    // For each energy term, loop though all threads and sum
    // energy and pE_pt contributions.
    OC_REAL8m pE_pt_term = chunk1[ei].ocedtaux.pE_pt_accum;
    for(int ithread=0;ithread<thread_count;++ithread) {
      pE_pt_term
        += chunk_thread[ithread].energy_terms[ei].ocedtaux.pE_pt_accum;
    }
    oced1.pE_pt += pE_pt_term;

    OC_REAL8m energy_term = chunk1[ei].ocedtaux.energy_total_accum;
    for(int ithread=0;ithread<thread_count;++ithread) {
      energy_term
        += chunk_thread[ithread].energy_terms[ei].ocedtaux.energy_total_accum;
    }
    oced1.energy_sum += energy_term;

    if(eterm.energy_sum_output.GetCacheRequestCount()>0) {
      eterm.energy_sum_output.cache.value=energy_term;
      eterm.energy_sum_output.cache.state_id=state.Id();
    }

    if(eterm.field_output.GetCacheRequestCount()>0) {
      eterm.field_output.cache.state_id=state.Id();
    }

    if(eterm.energy_density_output.GetCacheRequestCount()>0) {
      eterm.energy_density_output.cache.state_id=state.Id();
    }

#if REPORT_TIME
    Nb_StopWatch bar;
    bar.ThreadAccum(chunk1[ei].ocedtaux.energytime);
    for(int ithread=0;ithread<thread_count;++ithread) {
      bar.ThreadAccum
        (chunk_thread[ithread].energy_terms[ei].ocedtaux.energytime);

    }
    eterm.energytime.Accum(bar);
#endif // REPORT_TIME
  }

  if(oceed.mxHxm!=0 && oced1.mxH_accum == oceed.mxHxm) {
    // Undo mxHxm hack
    oced1.mxH_accum = 0;
  }

  oceed.max_mxH = 0.0;
  for(vector<Oxs_ComputeEnergiesChunkThread>::const_iterator cect
        = chunk_thread.begin(); cect != chunk_thread.end() ; ++cect ) {
    if(cect->max_mxH > oceed.max_mxH) oceed.max_mxH = cect->max_mxH;
  }

#if REPORT_TIME
  Oxs_ChunkEnergy::chunktime.Stop();
#endif

// =========================================================================
// Lattice 2
// =========================================================================
  Oxs_ComputeEnergiesChunkThread::Init(thread_count,
                                       state.spin.GetArrayBlock());

  chunk_thread.resize(thread_count);
  chunk_thread[0].state     = &state2;
  chunk_thread[0].energy_terms = chunk2; // Make copies.
  chunk_thread[0].mxH       = oced2.mxH;
  chunk_thread[0].mxH_accum = oced2.mxH_accum;
  chunk_thread[0].mxHxm     = oceed.mxHxm;
  chunk_thread[0].fixed_spins = oceed.fixed_spin_list;
  chunk_thread[0].cache_blocksize = cache_blocksize;
  chunk_thread[0].accums_initialized = accums_initialized;

  // Initialize chunk energy computations
  for(vector<Oxs_ComputeEnergies_ChunkStruct>::iterator it
        = chunk2.begin(); it != chunk2.end() ; ++it ) {
    Oxs_ChunkEnergy& eterm = *(it->energy);  // For code clarity
    Oxs_ComputeEnergyDataThreaded& ocedt = it->ocedt;
    Oxs_ComputeEnergyDataThreadedAux& ocedtaux = it->ocedtaux;
    eterm.ComputeEnergyChunkInitialize(state2,ocedt,ocedtaux,
                                       thread_count);
  }

  for(int ithread=1;ithread<thread_count;++ithread) {
    chunk_thread[ithread] = chunk_thread[0];
    threadtree.Launch(chunk_thread[ithread],0);
  }
  threadtree.LaunchRoot(chunk_thread[0],0);

  // Note: If chunk.size()>0, then we are guaranteed that accums are
  // initialized.  If accums_initialized is ever needed someplace
  // downstream, then uncomment the following line:
  // if(chunk.size()>0) accums_initialized = 1;

  // Finalize chunk energy computations
  for(OC_INDEX ei=0;static_cast<size_t>(ei)<chunk2.size();++ei) {

    Oxs_ChunkEnergy& eterm = *(chunk2[ei].energy);  // Convenience
    const Oxs_ComputeEnergyDataThreaded& ocedt = chunk2[ei].ocedt;
    const Oxs_ComputeEnergyDataThreadedAux& ocedtaux = chunk2[ei].ocedtaux;

    eterm.ComputeEnergyChunkFinalize(state2,ocedt,ocedtaux,
                                     thread_count);

    ++(eterm.calc_count);

    // For each energy term, loop though all threads and sum
    // energy and pE_pt contributions.
    OC_REAL8m pE_pt_term = chunk2[ei].ocedtaux.pE_pt_accum;
    for(int ithread=0;ithread<thread_count;++ithread) {
      pE_pt_term
        += chunk_thread[ithread].energy_terms[ei].ocedtaux.pE_pt_accum;
    }
    oced2.pE_pt += pE_pt_term;

    OC_REAL8m energy_term = chunk2[ei].ocedtaux.energy_total_accum;
    for(int ithread=0;ithread<thread_count;++ithread) {
      energy_term
        += chunk_thread[ithread].energy_terms[ei].ocedtaux.energy_total_accum;
    }
    oced2.energy_sum += energy_term;

    if(eterm.energy_sum_output.GetCacheRequestCount()>0) {
      eterm.energy_sum_output.cache.value=energy_term;
      eterm.energy_sum_output.cache.state_id=state.Id();
    }

    if(eterm.field_output.GetCacheRequestCount()>0) {
      eterm.field_output.cache.state_id=state.Id();
    }

    if(eterm.energy_density_output.GetCacheRequestCount()>0) {
      eterm.energy_density_output.cache.state_id=state.Id();
    }

#if REPORT_TIME
    Nb_StopWatch bar;
    bar.ThreadAccum(chunk2[ei].ocedtaux.energytime);
    for(int ithread=0;ithread<thread_count;++ithread) {
      bar.ThreadAccum
        (chunk_thread[ithread].energy_terms[ei].ocedtaux.energytime);

    }
    eterm.energytime.Accum(bar);
#endif // REPORT_TIME
  }

  if(oceed.mxHxm!=0 && oced2.mxH_accum == oceed.mxHxm) {
    // Undo mxHxm hack
    oced2.mxH_accum = 0;
  }

  oceed.max_mxH = 0.0;
  for(vector<Oxs_ComputeEnergiesChunkThread>::const_iterator cect
        = chunk_thread.begin(); cect != chunk_thread.end() ; ++cect ) {
    if(cect->max_mxH > oceed.max_mxH) oceed.max_mxH = cect->max_mxH;
  }

#if REPORT_TIME
  Oxs_ChunkEnergy::chunktime.Stop();
#endif

}
