/****************************************************************
 *
 * force_elstat.c: Routine used for calculating pair/monopole/dipole
 *     forces/energies in various interpolation schemes.
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
 *
 ****************************************************************
 *
 * This file is part of potfit.
 *
 * potfit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * potfit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#if !defined(COULOMB)
#error force_elstat_csh.c compiled without COULOMB support
#endif

#include "potfit.h"

#include "force.h"
#include "functions.h"
#include "memory.h"
#include "potential_input.h"
#include "potential_output.h"
#include "splines.h"
#include "utils.h"

/****************************************************************
 *
 *  init_forces
 *      called after all parameters and potentials are read
 *      additional assignments and initializations can be made here
 *
 ****************************************************************/

void init_force(int is_worker)
{
#if defined(COULOMB)
  if (g_pot.apot_table.sw_kappa)  // FIXME is sw_kappa really correct here?
    init_tails(g_pot.apot_table.dp_kappa[0]);
#endif  // COULOMB
}

/****************************************************************
 *
 *  compute forces using pair potentials with spline interpolation
 *
 *  returns sum of squares of differences between calculated and reference
 *     values
 *
 *  arguments: *xi - pointer to short-range potential
 *             *forces - pointer to forces calculated from potential
 *             flag - used for special tasks
 *
 * When using the mpi-parallelized version of potfit, all processes but the
 * root process jump into this function immediately after initialization and
 * stay in here for an infinite loop, to exit only when a certain flag value
 * is passed from process 0. When a set of forces needs to be calculated,
 * the root process enters the function with a flag value of 0, broadcasts
 * the current potential table xi and the flag value to the other processes,
 * thus initiating a force calculation. Whereas the root process returns with
 * the result, the other processes stay in the loop. If the root process is
 * called with flag value 1, all processes exit the function without
 * calculating the forces.
 * If anything changes about the potential beyond the values of the parameters,
 * e.g. the location of the sampling points, these changes have to be broadcast
 * from rank 0 process to the higher ranked processes. This is done when the
 * root process is called with flag value 2. Then a potsync function call is
 * initiated by all processes to get the new potential from root.
 *
 * xi_opt is the array storing the potential parameters (usually it is the
 *     g_pot.opt_pot.table - part of the struct g_pot.opt_pot, but it can also
 *be
 *     modified from the current potential.
 *
 * forces is the array storing the deviations from the reference data, not
 *     only for forces, but also for energies, stresses or dummy constraints
 *     (if applicable).
 *
 * flag is an integer controlling the behaviour of calc_forces_pair.
 *    flag == 1 will cause all processes to exit calc_forces_pair after
 *             calculation of forces.
 *    flag == 2 will cause all processes to perform a potsync (i.e. broadcast
 *             any changed potential parameters from process 0 to the others)
 *             before calculation of forces
 *    all other values will cause a set of forces to be calculated. The root
 *             process will return with the sum of squares of the forces,
 *             while all other processes remain in the function, waiting for
 *             the next communication initiating another force calculation
 *             loop
 *
 ****************************************************************/

double calc_forces(double* xi_opt, double* forces, int flag)
{
  double tmpsum = 0.0;
  double sum = 0.0;
  int first = 0;
  int col = 0;
  int ne = 0;
  int size = 0;
  int i = flag;
  double* xi = NULL;
  apot_table_t* apt = &g_pot.apot_table;
  double charge[g_param.ntypes];
  double sum_charges;
  double dp_kappa;


  switch (g_pot.format_type) {
    case POTENTIAL_FORMAT_UNKNOWN:
      break;
    case POTENTIAL_FORMAT_ANALYTIC:
      xi = g_pot.calc_pot.table;
      break;
    case POTENTIAL_FORMAT_TABULATED_EQ_DIST:
    case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST:
      xi = xi_opt;
      break;
  }

  ne = g_pot.apot_table.total_ne_par;
  size = apt->number;

  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.0; /* sum of squares of local process */

#if defined(APOT) && !defined(MPI)
    if (g_pot.format_type == POTENTIAL_FORMAT_ANALYTIC) {
      apot_check_params(xi_opt);
      update_calc_table(xi_opt, xi, 0);
    }
#endif  // APOT && !MPI

#if defined(MPI)
/* exchange potential and flag value */
#if !defined(APOT)
    MPI_Bcast(xi, g_pot.calc_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif  // !APOT
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (flag == 1)
      break; /* Exception: flag 1 means clean up */

#if defined(APOT)
    if (g_mpi.myid == 0)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, g_calc.ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (g_pot.format_type == POTENTIAL_FORMAT_ANALYTIC)
      update_calc_table(xi_opt, xi, 0);
#else   // APOT
    /* if flag==2 then the potential parameters have changed -> sync */
    if (flag == 2)
      potsync();
#endif  // APOT
#endif  // MPI

    /* local arrays for electrostatic parameters */
    sum_charges = 0;
    for (i = 0; i < g_param.ntypes - 1; i++) {
      if (xi_opt[2 * size + ne + i]) {
        charge[i] = xi_opt[2 * size + ne + i];
        sum_charges += apt->ratio[i] * charge[i];
      } else {
        charge[i] = 0.0;
      }
    }
    apt->last_charge = -sum_charges / apt->ratio[g_param.ntypes - 1];
    charge[g_param.ntypes - 1] = apt->last_charge;
    if (xi_opt[2 * size + ne + g_param.ntypes - 1]) {
      dp_kappa = xi_opt[2 * size + ne + g_param.ntypes - 1];
    } else {
      dp_kappa = 0.0;
    }

    /* init second derivatives for splines */
    for (col = 0; col < g_calc.paircol; col++) {
      first = g_pot.calc_pot.first[col];
      switch (g_pot.format_type) {
        case POTENTIAL_FORMAT_UNKNOWN:
          error(1, "Unknown potential format detected! (%s:%d)\n", __FILE__,
                __LINE__);
        case POTENTIAL_FORMAT_ANALYTIC:
        case POTENTIAL_FORMAT_TABULATED_EQ_DIST: {
          spline_ed(g_pot.calc_pot.step[col], xi + first,
                    g_pot.calc_pot.last[col] - first + 1, *(xi + first - 2),
                    0.0, g_pot.calc_pot.d2tab + first);
          break;
        }
        case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST: {
          spline_ne(g_pot.calc_pot.xcoord + first, xi + first,
                    g_pot.calc_pot.last[col] - first + 1, *(xi + first - 2),
                    0.0, g_pot.calc_pot.d2tab + first);
        }
      }
    }

#if !defined(MPI)
    g_mpi.myconf = g_config.nconf;
#endif  // MPI

    /* region containing loop over configurations,
       also OMP-parallelized region */
    {
      int self;
      vector tmp_force;
      int h, j, type1, type2, uf;
#if defined(STRESS)
      int us, stresses;
#endif  // STRESS
      int n_i, n_j;
      double fnval, grad, fnval_tail, grad_tail, grad_i, grad_j;
      atom_t* atom;
      neigh_t* neigh;

#if defined(DEBUG)
      double coulener, csener, vdwener, angener;

      printf("Debug information for Core-Shell potential: \n");
      printf("Config   Ecoul    Evdw     Ecs     Eang \n") ;
#endif // DEBUG 
      /* loop over configurations: M A I N LOOP CONTAINING ALL ATOM-LOOPS */
      for (h = g_mpi.firstconf; h < g_mpi.firstconf + g_mpi.myconf; h++) {
        uf = g_config.conf_uf[h - g_mpi.firstconf];
#if defined(STRESS)
        us = g_config.conf_us[h - g_mpi.firstconf];
#endif  // STRESS
        /* reset energies and stresses */
        forces[g_calc.energy_p + h] = 0.0;
#if defined(STRESS)
        stresses = g_calc.stress_p + 6 * h;
        for (i = 0; i < 6; i++)
          forces[stresses + i] = 0.0;
#endif  // STRESS


        /* F I R S T LOOP OVER ATOMS: reset forces, dipoles */
        for (i = 0; i < g_config.inconf[h]; i++) { /* atoms */
          n_i = 3 * (g_config.cnfstart[h] + i);
          if (uf) {
            forces[n_i + 0] = -g_config.force_0[n_i + 0];
            forces[n_i + 1] = -g_config.force_0[n_i + 1];
            forces[n_i + 2] = -g_config.force_0[n_i + 2];
          } else {
            forces[n_i + 0] = 0.0;
            forces[n_i + 1] = 0.0;
            forces[n_i + 2] = 0.0;
          }
        } /* end F I R S T LOOP */

        /* S E C O N D loop: calculate short-range and monopole forces,
           calculate static field- and dipole-contributions */
        for (i = 0; i < g_config.inconf[h]; i++) { /* atoms */
          atom =
              g_config.conf_atoms + i + g_config.cnfstart[h] - g_mpi.firstatom;
          type1 = atom->type;
          n_i = 3 * (g_config.cnfstart[h] + i);
          for (j = 0; j < atom->num_neigh; j++) { /* neighbors */
            neigh = atom->neigh + j;
            type2 = neigh->type;
            col = neigh->col[0];

            /* updating tail-functions - only necessary with variing kappa */
            if (!apt->sw_kappa)
	    elstat_dsf(neigh->r, dp_kappa, &neigh->fnval_el,
                           &neigh->grad_el, &neigh->ggrad_el);
            /* In small cells, an atom might interact with itself */
            self = (neigh->nr == i + g_config.cnfstart[h]) ? 1 : 0;

            /* calculate short-range forces */
            if (neigh->r < g_pot.calc_pot.end[col]) {
              if (uf) {
                fnval = splint_comb_dir(&g_pot.calc_pot, xi, neigh->slot[0],
                                        neigh->shift[0], neigh->step[0], &grad);
              } else {
                fnval = splint_dir(&g_pot.calc_pot, xi, neigh->slot[0],
                                   neigh->shift[0], neigh->step[0]);
              }

              /* avoid double counting if atom is interacting with a
                 copy of itself */
              if (self) {
                fnval *= 0.5;
                grad *= 0.5;
              }
              forces[g_calc.energy_p + h] += fnval;

              if (uf) {
                tmp_force.x = neigh->dist_r.x * grad;
                tmp_force.y = neigh->dist_r.y * grad;
                tmp_force.z = neigh->dist_r.z * grad;
                forces[n_i + 0] += tmp_force.x;
                forces[n_i + 1] += tmp_force.y;
                forces[n_i + 2] += tmp_force.z;
                /* actio = reactio */
                n_j = 3 * neigh->nr;
                forces[n_j + 0] -= tmp_force.x;
                forces[n_j + 1] -= tmp_force.y;
                forces[n_j + 2] -= tmp_force.z;

#if defined(STRESS)
                /* calculate pair stresses */
                if (us) {
                  forces[stresses + 0] -= neigh->dist.x * tmp_force.x;
                  forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
                  forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
                  forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
                  forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
                  forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
                }
#endif  // STRESS
              }
            }

            /* calculate monopole forces */
            if (neigh->r < g_config.dp_cut &&
                (charge[type1] || charge[type2])) {
              fnval_tail = neigh->fnval_el;
              grad_tail = neigh->grad_el;

              grad_i = charge[type2] * grad_tail;
              if (type1 == type2) {
                grad_j = grad_i;
              } else {
                grad_j = charge[type1] * grad_tail;
              }
              fnval = charge[type1] * charge[type2] * fnval_tail;
              grad = charge[type1] * grad_i;

	      /* supress coulomb contribution if pair is a core-shell one */
              if ( (int) (g_pot.apot_table.cweight[col]) == 0 ) {  /* coreshell pair */
                 if (neigh->r < g_pot.calc_pot.end[col]) {  
                     fnval -= DP_EPS * charge[type1] * charge[type2] / neigh->r;
                     grad=0;
                 }
              }

              if (self) {
                grad_i *= 0.5;
                grad_j *= 0.5;
                fnval *= 0.5;
                grad *= 0.5;
              }

              forces[g_calc.energy_p + h] += fnval;

              if (uf) {
                tmp_force.x = neigh->dist.x * grad;
                tmp_force.y = neigh->dist.y * grad;
                tmp_force.z = neigh->dist.z * grad;
                forces[n_i + 0] += tmp_force.x;
                forces[n_i + 1] += tmp_force.y;
                forces[n_i + 2] += tmp_force.z;
                /* actio = reactio */
                n_j = 3 * neigh->nr;
                forces[n_j + 0] -= tmp_force.x;
                forces[n_j + 1] -= tmp_force.y;
                forces[n_j + 2] -= tmp_force.z;
#if defined(STRESS)
                /* calculate coulomb stresses */
                if (us) {
                  forces[stresses + 0] -= neigh->dist.x * tmp_force.x;
                  forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
                  forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
                  forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
                  forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
                  forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
                }
#endif  // STRESS
              }
            }

          } /* loop over neighbours */
        }   /* end S E C O N D loop over atoms */


        /* F I F T H  loop: self energy contributions and sum-up force
         * contributions */
        double qq;
	double fnval_cut, gtail_cut, ggrad_cut;
        elstat_value(g_config.dp_cut, dp_kappa, &fnval_cut, &gtail_cut, &ggrad_cut);
        for (i = 0; i < g_config.inconf[h]; i++) { /* atoms */
          atom =
              g_config.conf_atoms + i + g_config.cnfstart[h] - g_mpi.firstatom;
          type1 = atom->type;
          n_i = 3 * (g_config.cnfstart[h] + i);

          /* self energy contributions */
          if (charge[type1]) {
            qq = charge[type1] * charge[type1];
	    fnval = qq * ( DP_EPS * dp_kappa / sqrt(M_PI) +
	       (fnval_cut - gtail_cut * g_config.dp_cut * g_config.dp_cut )*0.5 );
            forces[g_calc.energy_p + h] -= fnval;
          }

          /* sum-up: whole force contributions flow into tmpsum */
          if (uf) {
#if defined(FWEIGHT)
            /* Weigh by absolute value of force */
            forces[n_i + 0] /= FORCE_EPS + atom->absforce;
            forces[n_i + 1] /= FORCE_EPS + atom->absforce;
            forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif  // FWEIGHT
#if defined(CONTRIB)
            if (atom->contrib)
#endif  // CONTRIB
              tmpsum += g_config.conf_weight[h] *
                        (dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) +
                         dsquare(forces[n_i + 2]));
          }
   //   printf(" at: %d          f.x %f  f.x %f  f.x %f \n " ,  i, forces[n_i + 0],forces[n_i + 1],forces[n_i + 2]);
        } /* end F I F T H loop over atoms */

        /* whole energy contributions flow into tmpsum */
        forces[g_calc.energy_p + h] /= (double)g_config.inconf[h];
        forces[g_calc.energy_p + h] -= g_config.force_0[g_calc.energy_p + h];
        tmpsum += g_config.conf_weight[h] * g_param.eweight *
                  dsquare(forces[g_calc.energy_p + h]);

#if defined(STRESS)
        /* whole stress contributions flow into tmpsum */
        if (uf && us) {
          for (i = 0; i < 6; i++) {
            forces[stresses + i] /= g_config.conf_vol[h - g_mpi.firstconf];
            forces[stresses + i] -= g_config.force_0[stresses + i];
            tmpsum += g_config.conf_weight[h] * g_param.sweight *
                      dsquare(forces[stresses + i]);
          }
        }
#endif  // STRESS
      } /* end M A I N loop over configurations */
    }   /* parallel region */

/* dummy constraints (global) */
#if defined(APOT)
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (g_mpi.myid == 0) {
      tmpsum += apot_punish(xi_opt, forces);
    }
#endif  // APOT

#if defined(MPI)
    /* reduce global sum */
    sum = 0.0;
    MPI_Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    if (g_mpi.myid == 0) { /* root node already has data in place */
      /* forces */
      MPI_Gatherv(MPI_IN_PLACE, g_mpi.myatoms, g_mpi.MPI_VECTOR, forces,
                  g_mpi.atom_len, g_mpi.atom_dist, g_mpi.MPI_VECTOR, 0,
                  MPI_COMM_WORLD);
      /* energies */
      MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, MPI_DOUBLE,
                  forces + g_calc.energy_p, g_mpi.conf_len, g_mpi.conf_dist,
                  MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined(STRESS)
      /* stresses */
      MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, g_mpi.MPI_STENS,
                  forces + g_calc.stress_p, g_mpi.conf_len, g_mpi.conf_dist,
                  g_mpi.MPI_STENS, 0, MPI_COMM_WORLD);
#endif  // STRESS
    } else {
      /* forces */
      MPI_Gatherv(forces + g_mpi.firstatom * 3, g_mpi.myatoms, g_mpi.MPI_VECTOR,
                  forces, g_mpi.atom_len, g_mpi.atom_dist, g_mpi.MPI_VECTOR, 0,
                  MPI_COMM_WORLD);
      /* energies */
      MPI_Gatherv(forces + g_calc.energy_p + g_mpi.firstconf, g_mpi.myconf,
                  MPI_DOUBLE, forces + g_calc.energy_p, g_mpi.conf_len,
                  g_mpi.conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined(STRESS)
      /* stresses */
      MPI_Gatherv(forces + g_calc.stress_p + 6 * g_mpi.firstconf, g_mpi.myconf,
                  g_mpi.MPI_STENS, forces + g_calc.stress_p, g_mpi.conf_len,
                  g_mpi.conf_dist, g_mpi.MPI_STENS, 0, MPI_COMM_WORLD);
#endif  // STRESS
    }
#else
    sum = tmpsum; /* global sum = local sum  */
#endif  // MPI

    /* root process exits this function now */
    if (g_mpi.myid == 0) {
      g_calc.fcalls++; /* Increase function call counter */
      if (isnan(sum)) {
#if defined(DEBUG)
        printf("\n--> Force is nan! <--\n\n");
#endif  // DEBUG
        return 10e10;
      } else
        return sum;
    }
  }

  /* once a non-root process arrives here, all is done. */
  return -1.0;
}
