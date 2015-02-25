/** FILE: yy_2lat_util.h                 -*-Mode: c++-*-
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

#ifndef _YY_2LAT_UTIL
#define _YY_2LAT_UTIL

void YY_2LatComputeEnergies(
    const Oxs_SimState& state,  // the "total" lattice
    Oxs_ComputeEnergyData& oced1,
    Oxs_ComputeEnergyData& oced2,
    const vector<Oxs_Energy*>& energies,
    Oxs_ComputeEnergyExtraData& oceed);
  // Compute sums of energies, fields, and/or torques for all energies
  // in "energies" import.  On entry, oced.energy_accum, oced.H_accum,
  // and oced.mxH_accum should be set or null as desired.
  // oced.scratch_energy and oced.scratch_H must be non-null.
  // oced.energy, oced.H and oced.mxH *must* be *null* on entry.  This
  // routine does not fill these fields, but rather the accumulated
  // values are collected as necessary in oced.*_accum entries.
  // (However, pointers to energy_accum, H_accum, and mxH_accum may
  // be temporarily swapped to energy, H, and mxH for initialization
  // purposes.  This is transparent to the Oxs_ComputeEnergies caller,
  // but will exercise the non-accum portions of callees.)
  //   This routine handles outputs, energy calculation counts, and
  // timers appropriately.
  //   Those "energies" members that are actually Oxs_ChunkEnergies will
  // use the ComputeEnergyChunk interface in a collated fashion to help
  // minimize memory bandwidth usage.  On threaded OOMMF builds, these
  // calls will be run in parallel and load balanced.  Also, the number
  // of threads launched will not exceed the number of chunks.  This is
  // to insure that the main thread (threadnumber == 0) has an
  // opportunity to run for initialization purposes in the
  // Oxs_ChunkEnergy::ComputeEnergyChunk() function.  (Oxs_ChunkEnergy
  // classes that make (or may make) call into the Tcl interpreter must
  // use threadnumber == 0 for those calls, as per Tcl specs.  So if
  // all threads with threadnumber != 0 block on ComputeEnergyChunk()
  // entry, then the threadnumber == 0 is guaranteed at least one call
  // into ComputeEnergyChunk().
  //
  // Update May-2009: The now preferred initialization method is to
  // use ComputeEnergyChunkInitialize.  The guarantee that threadnumber
  // 0 will always run is honored for backward compatibility, but new
  // code should use ComputeEnergyChunkInitialize instead.
  //
  //    Data in Oxs_ComputeEnergyExtraData are filled in on the backside
  // of the chunk compute code.  These results could be computed by the
  // client, but doing it here gives improved cache locality.

#endif  // _YY_2LAT_UTIL
