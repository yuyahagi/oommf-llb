OOMMF-LLB: Landau-Lifshitz-Bloch evolver for ferrimagnets
=========================================================

oommf-llb is an [OOMMF (Object-Oriented Micromagnetic Framework)](http://math.nist.gov/oommf/) extension to simulate the magnetization dynamics of ferrimagnets (two sublattices of spins coupled via exchange energy) using micromagnetic Landau-Lifshitz-Bloch (LLB) equation. The LLB equation was proposed as a valid equation of motion of magnetic spins at elevated temperature, even beyond the Currie temperature.


Features
--------

* Micromagnetic LLB evolver, valid at elevated temperature
* Stochastic field to include randomness
* Coupled 2-lattice spin systems (ferrimagnets)
* Temperature dependent interlattice exchange coupling
* Demagnetization field calculated with the total magnetization
* Uniaxial anisotropy can be specified for each sublattice


Installation
------------

You need the OOMMF source code (recommended v1.2a5 or above) and an building environment. See [OOMMF User's Guide](http://math.nist.gov/oommf/doc/userguide12a5/userguide/Basic_Installation.html).

### Caution ###

In order to circumvent the fundamental restriction of conventional micromagnetics (constant Ms), oommf-llb replaces several common Oxs (OOMMF eXtensible Solver) classes. This is done so with no intention to alter the behavior of the other Oxs classes but you are recommended to copy the entire OOMMF directory to another location and install the oommf-llb extension to the new copy. The extension does not change the files outside this directory. Alternatively, you can keep the original source files to be replaced in installation.

### Manual installation ###

1. Copy source codes (\*.cc and \*.h) in this directory to (OOMMF_DIR)/app/oxs/local/ where (OOMMF_DIR) is the directory OOMMF is installed in your system.

2. You also need to copy files in base_src/ and ext_src/ to (OOMMF_DIR)/app/oxs/base/ and (OOMMF_DIR)/app/oxs/ext/, respectively. These files will replace the original OOMMF source files so keep the original as backup or 

3. Go to (OOMMF_DIR) and run pimake as
    
    ```tclsh oommf.tcl pimake```

### Makefile configuration ###

If your system runs `make` and standard \*nix commands, you can automate the installation process. Makefile is configured for Mac OSX, on which this extension was developed. For the other environments, you need to modify the first part of Makefile. Note that Makefiles are also in subdirectories base_src/ and ext_src/ and you need to modify all of them.

1. Copy all the extension codes to (OOMMF_DIR)/app/oxs/contrib/llb/. Makefiles assume this relative directory location to the main OOMMF directory. You can also work from other directories if you modify Makefiles accordingly.

2. Edit Makefile in this directory, as well as in base_src/ and ext_src/. Replace 'darwin' to your platform name, e.g., 'wintel' for Windows and 'cygtel' for Cygwin environments. If you do not know your platform name, run `tclsh oommf.tcl pimake +platform` to check. You may also need to change the Tcl shell name from `tclsh` to `tclsh86` or `tclsh8.6` depending on your configuration.

3. Run `make` from (OOMMF_DIR)/app/oxs/contrib/llb/. This will copy all the extension codes to OOMMF directories and run `pimake` automatically. It also saves the original OOMMF sources to base_src/original/ and ext_src/original/ too.


Principles
----------

### Stochastic LLB micromagnetics ###

The stochastic evolver in this extension module is based on eq. (5) in [Evans et al, Phys. Rev. B 85, 014433 (2012)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.014433), where the new stochastic form of the LLB equation was shown to be consistent with Boltzmann distribution at higher temperature, unlike previously proposed form (eq. (1)). 

### Coupled LLB equation for ferrimagnets ###

An equation of motion for exchange-coupled spin sublattices (ferrimagnets) is derived in [Atxitia et al, Phys. Rev. B 86, 104414 (2012)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.104414). Each sublattice obeys the LLB equation with additional interlattice exchange terms in Heff, derived with mean field approximation. In oommf-llb, this is calculated in YY_2LatExchange6Ngbr. Note that, this formalism does not allow T > Tc. 

In the above-mentioned article, a single-macrospin model is assumed with no spatial variation. oommf-llb also incorporates the other spatially-varying energy contributions as in standard micromagnetic framework. These energy terms include the intralattice exchange and demagnetization field, as well as the ability to specify nonuniform parameters and temperature across the simulation mesh.


Oxs classes
-----------

### Single-lattice LLB ###

#### YY_LLBEulerEvolve ####

    Specify YY_LLBEulerEvolve:name {
        do_precess      < 0 | 1 >
        gamma_LL1       < value | scalarfield_spec >
        gamma_LL2       < value | scalarfield_spec >
        alpha_t1        < value | scalarfield_spec >
        alpha_t2        < value | scalarfield_spec >
        J1              < value | scalarfield_spec >
        J2              < value | scalarfield_spec >
        atom_moment1    < value | scalarfield_spec >
        atom_moment2    < value | scalarfield_spec >
        fixed_timestep  value
        tempscript      Tcl_script
        tempscript_args { args_request }
        # args_request is a subset of { stage stage_time total_time }
        uniform_seed    value
        use_stochastic  < 0 | 1 >
    }

### Two-lattice LLB ###

For two-lattice simulations, you need to use YY\_2LatEulerEvolve, YY\_2LatTimeDriver, and other YY\_2Lat\* energy terms, replacing Oxs\_EulerEvolve, Oxs\_TimeDriver, and other Oxs\_\* energy terms. The simulation can still run with inproper Oxs\_\* energy terms but with false results.

At this point, temperature above Tc is not supported. YY_2LatExchange6Ngbr throws and exception when T > Tc. A workaround is to simulate the two separate systems using YY_LLBEulerEvolve.

#### YY_2LatEulerEvolve ####

    Specify YY_2LatEulerEvolve {
        do_precess      < 0 | 1 >
        gamma_LL1       < value | scalarfield_spec >
        gamma_LL2       < value | scalarfield_spec >
        alpha_t1        < value | scalarfield_spec >
        alpha_t2        < value | scalarfield_spec >
        J1              < value | scalarfield_spec >
        J2              < value | scalarfield_spec >
        atom_moment1    < value | scalarfield_spec >
        atom_moment2    < value | scalarfield_spec >
        fixed_timestep  value
        tempscript      Tcl_script
        tempscript_args { args_request }
        # args_request is a subset of { stage stage_time total_time }
        uniform_seed    value
        use_stochastic  < 0 | 1 >
    }

#### YY_2LatTimeDriver ####

    Specify YY_2LatTimeDriver:name {
        basename              base_name
        evolver               evolver_name
        mesh                  mesh_name
        stopping_time         { stopping_time_list }
        stage_count           value
        Ms                    scalarfield_spec
        m0                    scalarfield_spec
        Ms1                   scalarfield_spec
        m01                   scalarfield_spec
        Ms2                   scalarfield_spec
        m02                   scalarfield_spec
        normalize_aveM_output < 0 | 1 >
    }

#### YY_2LatExchange6Ngbr ####

    Specify YY_2LatExchange6Ngbr {
        default_A    value
        atlas        atlas_spec
        A1 {
            region-1 region-1 A11
            region-1 region-2 A12
            ...
            region-m region-n Amn
        }
        A2 {
            region-1 region-1 A11
            region-1 region-2 A12
            ...
            region-m region-n Amn
        }
        J01          scalarfield_spec
        J02          scalarfield_spec
        J012         scalarfield_spec
        J021         scalarfield_spec
        atom_moment1 scalarfield_spec
        atom_moment2 scalarfield_spec
    }

#### YY_2LatUniaxialAnisotropy ####

    Specify YY_2LatUniaxialAnisotropy {
        K11   value
        axis1 { ex ey ez }
        K12   value
        axis2 { ex ey ez }
    }

#### YY_2LatDemag ####

    Specify YY_2LatDemag {}

Programmer's guide
------------------

### Modified common Oxs classes ###

#### Oxs_SimState ####

The following lines contain the `diff` of the Oxs_SimState class definition. Pointers to Ms and Ms_inverse `const` meshvalues were replaced with pointers with non-`const` mesh values in order to allow change in magnitude of magnetization. Accordingly, other Oxs classes (Oxs_Driver, Oxs_TimeDriver, etc) were modified to handle Ms and Ms_inverse with different type. In case the energy terms require Ms at T = 0K, Ms0 and Ms0_inverse is added.

Additional pointers to temperature T, Currie temperature Tc, equilibrium magnetization polarization m_e, and longitudinal susceptibility chi_l are included. Memory assignment to these variables are done in YY_LLBEulerEvolve for 1-lattice simulations, or in YY_2LatExchange6Ngbr for 2-lattice simulations.

For two-lattice simulations, three Oxs_SimState instances are used to express the simulation state at a given time, sublattice 1, 2, and the total magnetization as the vectorial sum of the two sublattices. The three states are connected with each other via the pointers, lattice1, lattice2, and total_lattice so that YY\_2Lat\* energy terms can refer to the other sublattice in Heff calculation. lattice_type indicates the type of a given Oxs_SimState instance and the pointer to the same type is kept NULL (e.g., the simulation state of the sublattice 1 has `state1.lattice1 == NULL`, `state1.lattice2 == &state2`, and `state1.total_lattice == &state`). These variables must be properly set when new simulation states are generated, first at the beginning of the simulation and then at each time step.

    -  const Oxs_MeshValue<OC_REAL8m>* Ms;
    -  const Oxs_MeshValue<OC_REAL8m>* Ms_inverse; // 1/Ms
    +  Oxs_MeshValue<OC_REAL8m>* Ms;
    +  Oxs_MeshValue<OC_REAL8m>* Ms_inverse; // 1/Ms
    +
    +  // Ms value at T = 0K, used in LLB simulations
    +  Oxs_MeshValue<OC_REAL8m>* Ms0;
    +  Oxs_MeshValue<OC_REAL8m>* Ms0_inverse;
    +
    +  // For calculation of the temperature-dependent parameters (e.g., 
    +  // anisotropy and exchange), include a pointer to temperature and related
    +  // parameters.
    +  mutable Oxs_MeshValue<OC_REAL8m> const* T;
    +  mutable Oxs_MeshValue<OC_REAL8m> const* Tc;
    +  mutable Oxs_MeshValue<OC_REAL8m> const* m_e;
    +  mutable Oxs_MeshValue<OC_REAL8m> const* chi_l;
    +
    +  // For 2 lattice simulation, pointer to the other sublattice
    +  enum LatticeType { TOTAL, LATTICE1, LATTICE2 } lattice_type;
    +  // Default: TOTAL for standard simulations
    +  Oxs_SimState* total_lattice;
    +  Oxs_SimState* lattice1;
    +  Oxs_SimState* lattice2;

#### Oxs_Driver and Oxs_TimeDriver ####

Several functions have been modified in order to handle non-const Ms mesh value and additional pointers to the temperature-dependent parameters (T, Tc, m_e, and chi_l). Also, an abstract YY\_2LatDriver is added as a friend class.

#### Oxs_Energy and Oxs_ChunkEnergy ####

A friend function YY_2LatComputeEnergies() is added to Oxs_Energy and Oxs_ChunkEnergy classes, in order to allow for energy calculation for two sublattices.


License
-------
This projected is released under the terms of the GNU General Public License version 3. See LICENSE.


Contact
-------

Submit an issue at <https://github.com/yuyahagi/oommf-llb/issues>

Yu Yahagi: <yuyahagi2@gmail.com>
