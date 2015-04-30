
#include <stdlib.h>
#include <mpi.h>
#include "az_aztec.h" 
#include <stdio.h>
#include <math.h> 
/* Common Block Declarations */

union {
    struct {
	int nx;
    } _1;
    struct {
	int n;
    } _2;
} global_;

#define global_1 (global_._1)
#define global_2 (global_._2)

/* Global variables */

static double *buf_r_xm, *buf_s_xm, *buf_s_ym, *buf_s_zm,
                  *buf_s_xp, *buf_s_yp, *buf_s_zp, *buf_r_xp, 
		  *buf_r_ym, *buf_r_yp, *buf_r_zm, *buf_r_zp;

/* *************************************************************** */
/*     PROGRAM  CGMN/CARP-CG */
/* --------------------------------------------------------------- */
/*      Implementation of CGMN/CARP-CG with setup of examples 1 - 28. */

/*      Authors: Dan Gordon, Dept. of Computer Science, University of Haifa, */
/*               Email: gordon@cs.haifa.ac.il */
/*               Rachel Gordon, Dept. of Aerospace Engineering, */
/*               The Technion-Israel Insitute of Technology, */
/*               Email: rgordon@tx.technion.ac.il */
/*      Last update: June 28, 2011 */
/* --------------------------------------------------------------- */

/* Main program */ int main(int argc, char **argv)
/* Main program */ /*int MAIN(void) */
{

    double d1, d2;
    int ipr;
    static int example_no, read_update, kcg, order, model_no, icg;
    static double pi, freq, kn, kn1, kn2, kn3;   
    static double c_min, c_max;
    static double xmin, ymin, xmax, ymax, zmin, zmax;
    static double *x, *y, *cxy, *ccxy;
    float  x1, y1, c1;
    static double hx, hy,  t21, t1, t2;
    static int i, j, ij, i1, i2, i3, n, n0, jk, ik, is, ll, ny, nz, i_start;
    static int ierr, ikcg, iend, ndim, dim, nend, iinc, ii;
    static int *index_s_xp, *index_s_xm, *index_s_yp, *index_s_ym,
           *index_s_zp, *index_s_zm, 
	   *index_r_xp, *index_r_xm, *index_r_yp, *index_r_ym,
	   *index_r_zp, *index_r_zm;
    static int nproc, me, me_new, comm_cart;
    static int n_internal, n_external, n_update_internal, n_update_external, 
               avg_nonzero_per_row, total_nz;
    int *bindx, *update_inv, *iem, *iema;
    int *data_org , *external, *update_index, *extern_index, *update;      
    static double *val, *err, *rhsm, *xnew, *ssum;
    static double *pnewl_loc, *qnewl_loc, *snewl_loc, *xnewl_loc, 
           *xnewl_loc_init, *rhsm_loc;

    static int  rank_dest_x, rank_dest_y, direction_x, direction_y, 
	    direction_z, rank_dest_z;
    static int coords[3], dims[3], disp, ndims, kmax; 
    static int reorder, periods[3];
    static int ixp, ixm, iyp, iym, izp, izm;
    static double alpha, beta, sum, sum1, sum2, s2, s2p, s2s, pq, pqs, s2ps, 
           res0, rxmax, rmax, residmax, sumr; 
    static double sumrt, sumerr, sumerr1, sumerr2, err_l2;
    
    static double cons, omeg1;   

    static int nrow, nrow_half, nxst, nyst, nzst, nxend, nyend, nzend;
    static double *qnewl, *xnewl, *ecoef, *resid, *resid1;
    static int nlocal, nx_end, ny_end, nxloc, nyloc, nzloc, nx_st, ny_st;
    static int rank_source_x, rank_source_y, rank_source_z; 
    static int istart, nxproc, nyproc, nzproc;
    static int n_border, n_update, n_update_border, n_rec_length, n_neigh, n_id_neigh;
    static int n_total_send, n_send, n_send_length, n_send_list;
    int proc_config[AZ_PROC_SIZE];

    extern /* Subroutine */ int transf_2d(double *, double *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *), 
            transf_3d(double *, 
	    double *, int *, int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *);
    extern /* Subroutine */ int example_1d(int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *), 
	    example_2d(int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *, int *, int *), 
            example_3d(int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *),
	    bound_1d(int *, int *, int *, int *, int *, 
	    int *), 
            bound_2d(int *, int *, int *, int *,
	    int *, int *, int *, int *, int *, int *),
	    bound_3d(int *, int *, int *, int *, int *
	    , int *, int *, int *, int *, int *, int *
	    , int *, int *, int *);
    extern /* Subroutine */ int example_2d_img_all(int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *),
            example_2d_all(int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *), 
            example_2d_img(int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *),
	    example_3d_img(int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *);

    extern /* Subroutine */ int example1_err(double *, double *, 
	    int *, int *, int *, double *, int *),
            example2_err(double *, double *, 
	    int *, int *, int *, double *, int *),
            example3_err(double *, double *, 
	    int *, int *, int *, double *, int *),
            example4_err(double *, double *, 
	    int *, int *, int *, double *, int *), 
            example5_err(double *, double *, int *,
	    int *, int *, double *, int *), 
            example6_err(double *, double *, int *,
	    int *, int *, double *, int *), 
            example7_err(double *, 
	    double *, int *, int *, int *, double *, int *), 
	    example8_err(double *, double *, int *, int *, 
	    int *, double *, int *),
            example9_err(double *, double *, 
	    int *, int *, double *, int *), 
            example10_err(double *,
	    double *, int *, double *, int *), 
            example11_err(
	    double *, double *, int *, int *, double *, int *), 
	    example12_err(double *, double *, int *, int *, 
	    int *, double *, int *), 
            example13_err(double *, double *, 
	    int *, int *, int *, double *, int *), 
            example14_err(double *, double *, int *, int *, 
	    int *, double *, int *), 
            example15_err(double *, double *, int *, 
	    int *, double *, int *), 
            example16_err(double *,
	    double *, int *, int *, double *), 
	    example17_err(double *, double *, int *, int *, 
	    double *), 
            example18_err(double *, double *, int *, int *, 
	    int *, double *), 
            example19_err(double *, double *, int *, int *, 
	    int *, double *),
            example20_err(double *, double *, int *, int *, 
	    double *,int *), 
            example21_err(double *, double *, int *, int *,
	    double *),
            example22_err(double *, double *, int *, int *, 
	    int *, double *), 
            example23_err(double *, double *, int *, int *, 
	    double *, int *), 
            example24_err(double *, double *, int *, int *,
	    double *), 
            example25_err(double *, double *, int *, int *, 
	    int *, double *), 
	    example26_err(double *, double *, int *, int *, 
	    double *, double *, double *, double *, 
	    double *, int *, int *), 
            example27_err(double *, double *, int *, int *, 
	    double *, double *, double *, double *, 
            double *, int *, int *), 
	    example28_err(double *, double *, int *, int *, 
	    double *, double *, int *, int *);

    extern /* Subroutine */ int create_matrix_row_5pt_example21_2rows(
	    int *, int *, double *, int *, int *, 
	    double *, int *, int *, int *, int *, 
	    double *, double *, double *),
            create_matrix_row_9pt_example23_2rows(
	    int *, int *, double *, int *, int *, 
	    double *, int *, int *, int *, int *, 
	    double *, double *, double *, int *),
            create_matrix_row_5pt_example24_2rows(
	    int *, int *, double *, int *, int *, 
	    double *, int *, int *, int *, int *, 
	    double *, double *, double *, double *, 
	    double *, int *, int *),
            create_matrix_row_5_9_pt_example24_2rows(
            int *, int *, double *, int *, int *, double *,
            int *, int *, int *, int *, double *,
            double *, double *, double *, double *,
            int *, int *, int *, double *),
            create_matrix_row_7pt_example25_2rows(
	    int *, int *, double *, int *, int *, 
	    double *, int *, int *, int *, int *, int 
	    *, int *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, int *, int *, int *),
  /*  extern */ /* Subroutine */ 
            create_matrix_row_9pt_example26_2rows(int *, int *,
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *, double *, double *, 
	    double *, int *), 
            create_matrix_row_9pt_example27_2rows(
	    int *, int *, double *, int *, int *, 
	    double *, int *, int *, int *, int *, 
	    double *, double *, double *, int *); 

    extern /* Subroutine */ int 
	    create_matrix_row_7pt_example1(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *, int *, int *), 
	    create_matrix_row_7pt_example2(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *, int *, int *), 
	    create_matrix_row_7pt_example3(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *, int *, int *), 
	    create_matrix_row_7pt_example4(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *, int *, int *), 
	    create_matrix_row_7pt_example5(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *, int *, int *), 
	    create_matrix_row_7pt_example6(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *, int *, int *), 
	    create_matrix_row_7pt_example7(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *, int *, int *), 
	    create_matrix_row_7pt_example8(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *, int *, int *), 
	    create_matrix_row_5pt_example9(int *, int *, double 
	    *, int *, int *, double *, int *, int *, 
	    int *, int *),
  /*  extern */ /* Subroutine *//* int */
            create_matrix_row_5pt_example10(int *, 
	    int *, double *, int *, int *, double *, 
	    int *, int *), create_matrix_row_5pt_example11(int *
	    , int *, double *, int *, int *, double *, 
	    int *, int *, int *, int *),
            create_matrix_row_7pt_example12(int *, 
	    int *, double *, int *, int *, double *, 
	    int *, int *, int *, int *, int *, int *),
	     create_matrix_row_7pt_example13(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *, int *, int *),
            create_matrix_row_7pt_example14(int *, 
	    int *, double *, int *, int *, double *, 
	    int *, int *, int *, int *, int *, int *),
	    create_matrix_row_5pt_example15(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *), 
	    create_matrix_row_5pt_example16(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *), 
	    create_matrix_row_5pt_example17(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *), 
	    create_matrix_row_7pt_example18(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *, int *, int *), 
	    create_matrix_row_7pt_example19(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *, int *, int *), 
	    create_matrix_row_5pt_example20(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *, double *), 
	    create_matrix_row_7pt_example22(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *, int *, int *), 
	    create_matrix_row_9pt_example28(int *, int *, 
	    double *, int *, int *, double *, int *, 
	    int *, int *, int *, double *, int *); 

#define TRUE_ (1)
#define FALSE_ (0)
    FILE *f4;
    FILE *f7;
    FILE *f9;
    FILE *fp;

    MPI_Init ( &argc, &argv );
    MPI_Comm_size(MPI_COMM_WORLD, &nproc );
    MPI_Comm_rank(MPI_COMM_WORLD, &me );
    f7 = fopen("input.dat","r");
    if(f7==0) {printf("\n error opening input.dat"); exit;} 
    f9 = fopen("error.out1","w");
    if(f9==0) {printf("\n error opening error.out"); exit;}

/*   read input data */

    fscanf(f7, " %d %d %d ", &example_no, &global_1.nx, &read_update);
    fscanf(f7, " %d %d %d %d %d ", &nxproc, &nyproc, &nzproc, &kcg, &ierr);
    fscanf(f7, "%lf  %lf ", &residmax, &omeg1); 

    if (example_no == 20) {
        fscanf(f7, "%lf", &kn1);
    }
    if (example_no == 21) {
        fscanf(f7, "%lf  %lf %lf", &kn1, &kn2, &kn3);
    }
    if (example_no == 23) {
        fscanf(f7, "%lf  %lf %lf %d", &kn1, &kn2, &kn3, &order);
    }
    if (example_no == 26 || example_no == 27) {
        fscanf(f7, " %lf  %lf %lf %d", &kn1, &kn2, &kn3, &order);
    }
    if (example_no == 28) {
        fscanf(f7, "%lf  %d", &kn1, &order);
    }

/*  Description of the input data: */
/* ------------------------------------------------------------------------- */
/*     example_no - problem number */
/*     NX - the number of grid points in the x direction. */
/*          we assume NX=NY=NZ for all test cases except examples_no=24,25 */
/*     read_udate - Parameter that defines how the equations are distributed */
/*                  among the processors: */
/*                = 1 in case that  1*1*NZproc topology of processors */
/*                  is used in a 3D case or 1*NYproc*1 in a 2D case. */
/*                = 3 in case that NXproc*NYproc*NZproc processors */
/*                  topology is used in a 3D case, or NXproc*NYproc*1 */
/*                  in a 2D case,  where NXproc,NYproc,NZproc >= 1. In */
/*                  this case the indices of the equations of each */
/*                  processor are read from an external file. */
/*    NXproc,NYproc,NZproc - Number of sub-domains in X,Y,Z direction respectively. */
/*                           Each sub-domain is allocated to a different processor. */
/*                         (the total number of processors Np=NXproc*NYproc*NZproc) */
/*    KCG - the number of CARP-CG iterations. */
/*    ResidMax - the L2 norm of the residual required. */
/*    OMEG1 -  the relaxation parameter:  0 < OMEG1 < 2. */
/*    Kn1,Kn2,Kn3 - values of K for Helmholtz equation with three layers. */
/*    ierr = 0 - don't calculate the error during the iterative process
                 for examples having an analytic solution.     
           = 1 - calculate the error during the iterative process.  */

/* Print out the input data */

    if (me == 0) {
        printf("\n %s ","---------------------------------------------------------------");
        printf("\n    example_no= %d \n", example_no);
        printf("\n    NX= %d   ResidMax= %e   OMEG1= %e ", global_1.nx, residmax, omeg1);
        printf("\n    read_update= %d   No. of CG iterations= %d   ierr= %d \n",
                read_update,kcg,ierr);
        printf("\n    NXproc= %d    NYproc= %d     NZproc= %d",nxproc, nyproc, nzproc);

	if (example_no == 20) {
            printf("\n    Helmholtz equation:  K1= %5.1f", kn1);
	}
	if (example_no == 28) {
            printf("\n    Helmholtz equation:  K1= %5.1f    scheme order= %d",
               kn1, order);
	}
	if (example_no == 21) {
            printf("\n    Helmholtz equation:  K1= %5.1f    K2= %5.1f    K3=%5.1f",
               kn1, kn2, kn3);
	}
       if (example_no == 23 || example_no == 26 || example_no == 27) {
      printf("\n    Helmholtz equation:  K1= %5.1f    K2= %5.1f",kn1, kn2);
      printf( "    K3=%5.1f    scheme order= %d", kn3, order);
	}
        printf("\n %s ","---------------------------------------------------------------");
    }

    if (nxproc * nyproc * nzproc != nproc) {
	printf("\n  Nproc is incorrect, Nproc= %d", nproc);
        exit (0);
    }

/*  example 24 - Marmousi problem: */
/* -------------------------------- */
    if (example_no == 24) {

        fscanf(f7, "%lf %lf %lf %lf %lf %d %d %d %d", 
        &xmin, &xmax, &ymin, &ymax, &freq, &model_no, &global_1.nx, &ny, &order); 

/*	nrow = 2*(global_1.nx + 2) * (ny + 2); */
	nrow = 2*(global_1.nx + 4) * (ny + 4);
        nrow_half= 0.5*nrow;
        x = (double *) calloc(nrow_half + 1, sizeof(double));
        y = (double *) calloc(nrow_half + 1, sizeof(double));
        cxy = (double *) calloc(nrow_half + 1, sizeof(double));
        ccxy = (double *) calloc(nrow_half + 1, sizeof(double));


	if (me == 0) {
            printf("\n         Marmousi problem \n");
            printf("\n  Xmin= %f   Xmax= %f   Ymin= %f   Ymax= %f",
            xmin, xmax, ymin, ymax); 
            printf("\n  freq= %f  model_no= %d   NX= %d   NY= %d   scheme order=%d\n",
            freq, model_no, global_1.nx, ny, order);
	}
/* endif me=0 */
	if (model_no == 1) {
/*       Marmousi-1 model problem */
/*       open(4,FILE='/home/people/rgordon/Aztec.dir/marm2D.dir/marm2D/ */
/*    1data/marm_1.dat') */
            f4 = fopen("../Aztec.dir/marm2D.dir/marm2D/data/marm_1.dat","r");
            if(f4==0) {printf("\n error opening marm_1.dat"); exit;} 
	}
	if (model_no == 2) {
/*       Marmousi-2 model problem */
            f4 = fopen("../Aztec.dir/marm2D.dir/marm2D/data/marm_2.dat","r");
            if(f4==0) {printf("\n error opening marm_2.dat"); exit;} 
	}
	if (model_no == 3) {
/*       Marmousi-3 model problem */
            f4 = fopen("../Aztec.dir/marm2D.dir/marm2D/data/marm_3.dat","r");
            if(f4==0) {printf("\n error opening marm_3.dat"); exit;} 
	}
	if (model_no == 4) {
/*       Marmousi-4 model problem */
            f4 = fopen("../Aztec.dir/marm2D.dir/marm2D/data/marm_4.dat","r");
            if(f4==0) {printf("\n error opening marm_4.dat"); exit;} 
	}

/*  Set CXY() for each grid point(X();Y( )) */
/*  note: the value of c(x,y) is not available on the boundary and external grid */
/*        points. It will be taken as that of the adjacent grid point. */
	ij = 0;
	i1 = ny;
	for (j = 1; j <= i1; ++j) {
	    i2 = global_1.nx;
	    for (i = 1; i <= i2; ++i) {
		++ij;
		fscanf(f4," %f %f %f", &x1, &y1, &c1);
		x[ij - 1] = x1;
		y[ij - 1] = y1;
		ccxy[ij - 1]= c1;
/*              write(6,*) IJ,X(IJ),Y(IJ),CCXY(IJ) */
	    }
	}
	if (me == 0) {
            printf("\n  No. of grid points of input data for c(x,y) = %d", ij);
	}

/*  The first row (boundary + external grid points): */
        cxy[0] = ccxy[0];
        cxy[1] = ccxy[0];
/*	i1 = global_1.nx;  */
        is = 1;
	for (ii = 1; ii <= global_1.nx; ++ii) {
            ++is; 
	    cxy[is] = ccxy[ii-1];
	}
        ++is;
        cxy[is] = cxy[is-1];
        ++is;
        cxy[is] = cxy[is-1];

/*  define the second row:  */
        for (ii = 0; ii <= global_1.nx + 3; ++ii) {
            cxy[ii+global_1.nx+4]=cxy[ii];
        }

/*  The next NY rows: */
        is = 2*global_1.nx + 7;
	ij = 0;
        i1 = ny;
	for (j = 1; j <= i1; ++j) {
	    ++is;
	    cxy[is] = ccxy[ij];
	    ++is;
            cxy[is] = cxy[is-1]; 
/*          i2 = global_1.nx;
	    for (ii = 1; ii <= i2; ++ii) { */
	    for (ii = 1; ii <= global_1.nx; ++ii) {
		++is;
		++ij;
		cxy[is] = ccxy[ij - 1];
	    }
	    ++is;
	    cxy[is] = cxy[is - 1];
	    ++is;
	    cxy[is] = cxy[is - 1];
	}

/*  the (NY+3)th row (boundary points): */
        ++is;    
        i_start=is;
        ij -= global_1.nx;
        cxy[is] = ccxy[ij];
        ++is;
        cxy[is] = cxy[is - 1];
        i2 = global_1.nx;
        for (ii = 1; ii <= i2; ++ii) {
            ++is; 
            ++ij;
            cxy[is] = ccxy[ij-1];
        }
        ++is;
        cxy[is] = cxy[is - 1];
        ++is;
        cxy[is] = cxy[is - 1];

/*  The last (NY+4)th row (external points): */

        for (ii = 0; ii <= global_1.nx + 3; ++ii) {
            ++is;
            cxy[is] = cxy[i_start + ii];
        }

/*  Calculate C_max and C_min */

	if (me == 0) {
	    c_max = cxy[0];
	    c_min = cxy[0];
	    i1 = (global_1.nx + 4) * (ny + 4) - 1;
	    for (i = 1; i <= i1; ++i) {
		if (cxy[i] > c_max) {
		    c_max = cxy[i];
		}
		if (cxy[i] < c_min) {
		    c_min = cxy[i];
		}
	    }
            printf("\n    C_max= %f    C_min= %f ", c_max, c_min);
            printf("\n %s ","---------------------------------------------------------------");
	} /*  endif me=0 */

    }  /* end of example24 Marmousi problem  definition   */


/*  example 25 - 3D Helmholtz problem:  */
/* ------------------------------------ */
    if (example_no == 25) {
        fscanf(f7,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d", 
        &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, &kn1, &kn2, &kn3, &global_1.nx, &ny, &nz);

        if(me == 0){
        printf("\n  %s","    ");
        printf("\n              3D Helmholtz equation \n");
        printf("\n  Xmin= %f   Xmax= %f   Ymin= %f   Ymax= %f   Zmin= %f   Zmax= %f",
        xmin, xmax, ymin, ymax, zmin, zmax); 
        printf("\n    K1= %5.1f    K2= %5.1f    K3=%5.1f", kn1, kn2, kn3);
        printf("\n    NX= %d   NY= %d    NZ= %d \n", global_1.nx, ny, nz);
        printf("\n %s ","---------------------------------------------------------------");
       }  /* endif: if me  = 0 printout  */ 
    } /*  End reading data of example 25 */

/*  Set ndim - specifies the problem type. */

    ndim = 3;
    if (example_no == 10) {
	ndim = 1;
    }
    if (example_no == 9 || example_no == 11) {
	ndim = 2;
    }
    if (example_no == 15 || example_no == 16 || example_no == 20) {
	ndim = 2;
    }
    if (example_no == 17 || example_no == 28) {
	ndim = 4;
    }
/*  if (example_no == 21 || example_no == 24) { */
    if (example_no == 21 ) {
	ndim = 5;
    }
    if (example_no == 25) {
	ndim = 6;
    }
    if (example_no == 26) {
	ndim = 7;
    }
    if (example_no == 27) {
	ndim = 8;
    }
    if (example_no == 23) {
	ndim = 9;
    }
    if (example_no == 24) {
	ndim = 10;
    }

/* Set dim - the problem dimension (1D/2D/3D) */

    if(ndim == 3 || ndim == 6 ) dim = 3;
    if(ndim == 2 || ndim == 4 || ndim == 5 || ndim == 7 || ndim == 8 || ndim == 9 || ndim == 10) dim = 2;
    if(ndim == 1) dim = 1;

    if (example_no == 10 && nxproc > 1) {
	if (me == 0) {
            printf("\n  %s","    ");
            printf("\n  %s","The current version of the code can not solve example 10 on more than 1 processor");
	}
        exit (0) ; 
/*        stop */
    }

    if (example_no != 24 && example_no != 25) {
	ny = global_1.nx;
	nz = global_1.nx;
    }
    if (ndim == 10) {   /*  problem 24: Marmousi problem */
	n = (global_1.nx + 4 << 1) * (ny + 4); 
    }
    if (ndim == 9) {
	n = (global_1.nx + 4 << 1) * (ny + 3);
    }
/* problem 23 */
    if (ndim == 8) {
	n = (global_1.nx + 2 << 1) * (ny + 3);
    }
/* problem 27 */
    if (ndim == 7) {
	n = (global_1.nx + 3 << 1) * (ny + 2);
    }
/* problem 26 */
    if (ndim == 6) {
	n = (global_1.nx + 2 << 1) * (ny + 2) * (nz + 2);
    }
    if (ndim == 5) {
	n = (global_1.nx + 2 << 1) * (ny + 2);
    }
    if (ndim == 4) {
	n = (global_1.nx + 2) * (ny + 2);
    }
    if (ndim == 3) {
	n = global_1.nx * ny * nz;
    }
    if (ndim == 2) {
	n = global_1.nx * ny;
    }
    if (ndim == 1) {
	n = global_1.nx;
    }

/*     3D problems    */
/* ------------------ */
   if (dim == 3) {
      ndims = 3;
      dims[0] = nxproc;
      dims[1] = nyproc;
      dims[2] = nzproc;
      periods[0] = FALSE_;
      periods[1] = FALSE_;
      periods[2] = FALSE_;
      reorder = TRUE_;

/* Create a new communicator 'comm_cart' containing cartesian */
/* topology information */

      MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,
         &comm_cart);


/* Get my rank in new communicator */

      MPI_Comm_rank(comm_cart, &me_new );

/*    write(NOUT,*) 'my rank in new cartesian communicator:',me_new */

/* Get my cartesian coordinates */

      MPI_Cart_get(comm_cart, ndims, dims, periods, coords);

/* Get shifted source and destination ranks in each direction */
/*  for send/receive operation (i.e. get the rank of your neighbors) */

      direction_x = 0;
/* Coordinate dimension of send/rec 1-st dim */
      direction_y = 1;
      direction_z = 2;
      disp = 1;
/* displacement */
      MPI_Cart_shift( comm_cart, direction_x , disp , 
                      &rank_source_x, &rank_dest_x );
      disp = 1;
      MPI_Cart_shift( comm_cart,direction_y ,disp , 
                      &rank_source_y, &rank_dest_y );
      disp = 1;
      MPI_Cart_shift( comm_cart,direction_z ,disp , 
                      &rank_source_z, &rank_dest_z );

      printf("\n  %s","    ");
      printf("\n  me= %d  rank_source_x= %d   rank_dest_x= %d   rank_source_y= %d   rank_dest_y= %d",
               me,rank_source_x,rank_dest_x,rank_source_y,rank_dest_y);
      printf("   rank_source_z= %d   rank_dest_z= %d", rank_source_z,rank_dest_z);

    }

/*     2D problems    */
/* ------------------ */
    if (dim == 2)  {
       ndims = 2;
       dims[0] = nxproc;
       dims[1] = nyproc;
       periods[0] = FALSE_;
       periods[1] = FALSE_;
       reorder = TRUE_;

/* Create a new communicator 'comm_cart' containing cartesian */
/* topology information */

       MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,
         &comm_cart); 


/*   Get my rank in new communicator */

       MPI_Comm_rank(comm_cart, &me_new );
/*     write(NOUT,*) 'my rank in new cartesian communicator:',me_new */

/*   Get my cartesian coordinates */

       MPI_Cart_get(comm_cart, ndims, dims, periods,  coords);  

/*  Get shifted source and destination ranks in each direction       */
/*  for send/receive operation (i.e. get the rank of your neighbors) */

       direction_x = 0;
/* coordinate dimension of send/rec 1-st dim */
       direction_y = 1;
       disp = 1;
/* displacement */
       MPI_Cart_shift( comm_cart, direction_x ,disp ,   
                      &rank_source_x, &rank_dest_x );
       disp = 1;
       MPI_Cart_shift( comm_cart, direction_y ,disp ,    
                      &rank_source_y, &rank_dest_y );
 
       printf("\n  %s","    ");
       printf("\n  me= %d   rank_source_x= %d   rank_dest_x= %d   rank_source_y= %d   rank_dest_y= %d",
           me,rank_source_x,rank_dest_x,rank_source_y,rank_dest_y);
    }


/* Define the starting and ending grid point indices handled       */
/* by the current processor. */
/* --------------------------------------------------------------- */
/*   3D problems  */
/* -------------- */
    if (ndim == 3) {
       bound_3d(&global_1.nx, &ny, &nz, &nlocal, &nxst, &nxend, &nyst, &
		nyend, &nzst, &nzend, &nxproc, &nyproc, &nzproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d    NZst= %d   Nlocal= %d   NXend= %d   NYend= %d   NZend= %d",
              me, nxst, nyst, nzst, nlocal, nxend, nyend, nzend);
    }
    if (ndim == 6) {
       i1 = global_1.nx + 2;
       i2 = ny + 2;
       i3 = nz + 2;
       bound_3d(&i1, &i2, &i3, &nlocal, &nxst, &nxend, &nyst, &nyend,
		 &nzst, &nzend, &nxproc, &nyproc, &nzproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d    NZst= %d   Nlocal= %d   NXend= %d   NYend= %d   NZend= %d",
              me, nxst, nyst, nzst, nlocal, nxend, nyend, nzend);
    }

/*   2D problems  */
/* -------------- */
    if (ndim == 2) {
       bound_2d(&global_1.nx, &ny, &nlocal, &nxst, &nxend, &nyst, &nyend, &
		nxproc, &nyproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d   Nlocal= %d   NXend= %d   NYend= %d",
              me, nxst, nyst, nlocal, nxend, nyend);
    }
    if (ndim == 4) {
       i1 = global_1.nx + 2;
       i2 = ny + 2;
       bound_2d(&i1, &i2, &nlocal, &nxst, &nxend, &nyst, &nyend, &
		nxproc, &nyproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d   Nlocal= %d   NXend= %d   NYend= %d",
              me, nxst, nyst, nlocal, nxend, nyend);
    }
    if (ndim == 5) {
       i1 = global_1.nx + 2;
       i2 = ny + 2;
       bound_2d(&i1, &i2, &nlocal, &nxst, &nxend, &nyst, &nyend, &
		nxproc, &nyproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d   Nlocal= %d   NXend= %d   NYend= %d",
              me, nxst, nyst, nlocal, nxend, nyend);
    }
    if (ndim == 7) {
       i1 = global_1.nx + 3;
       i2 = ny + 2;
       bound_2d(&i1, &i2, &nlocal, &nxst, &nxend, &nyst, &nyend, &
		nxproc, &nyproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d   Nlocal= %d   NXend= %d   NYend= %d",
              me, nxst, nyst, nlocal, nxend, nyend);
    }
    if (ndim == 8) {
       i1 = global_1.nx + 2;
       i2 = ny + 3;
       bound_2d(&i1, &i2, &nlocal, &nxst, &nxend, &nyst, &nyend, &
		nxproc, &nyproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d   Nlocal= %d   NXend= %d   NYend= %d",
              me, nxst, nyst, nlocal, nxend, nyend);
    }
    if (ndim == 9) {
       i1 = global_1.nx + 4;
       i2 = ny + 3;
       bound_2d(&i1, &i2, &nlocal, &nxst, &nxend, &nyst, &nyend, &
		nxproc, &nyproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d   Nlocal= %d   NXend= %d   NYend= %d",
              me, nxst, nyst, nlocal, nxend, nyend);
    }
    if (ndim == 10) {
       i1 = global_1.nx + 4;
       i2 = ny + 4;
       bound_2d(&i1, &i2, &nlocal, &nxst, &nxend, &nyst, &nyend, &
		nxproc, &nyproc, &me_new);
       printf("\n  me= %d   NXst= %d   NYst= %d   Nlocal= %d   NXend= %d   NYend= %d",
              me, nxst, nyst, nlocal, nxend, nyend);
    }

/*   1D problems   */
/* --------------- */
    if (ndim == 1) {
       bound_1d(&global_1.nx, &nlocal, &nxst, &nxend, &nxproc, &me_new);
       printf("\n  me= %d   NXst= %d   Nlocal= %d   NXend= %d",
              me, nxst, nlocal, nxend);
    }


    AZ_set_proc_config(proc_config, comm_cart);     
/*  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);  */  
/*  me = proc_config[21];  */ 
    me = proc_config[AZ_node]; 

/*   1D / 2D /3D problems */
/* ------------------------ */
    if (ndim == 10) {
       nrow = 2*(global_1.nx + 4 ) * (ny + 4);
       if(order == 4 || order == 6) avg_nonzero_per_row = 2*9;
       if(order == 2 )  avg_nonzero_per_row = 2*5;
    }
    if (ndim == 9) {
/*     nrow = (global_1.nx + 4 << 1) * (ny + 3); */
       nrow = 2*(global_1.nx + 4 ) * (ny + 3);
       if(order == 4 || order == 6) avg_nonzero_per_row = 2*9;
       if(order == 2 )  avg_nonzero_per_row = 2*5;
    }
    if (ndim == 8) {
/*     nrow = (global_1.nx + 2 << 1) * (ny + 3); */
       nrow = 2*(global_1.nx + 2 ) * (ny + 3);
       if(order == 4 || order == 6) avg_nonzero_per_row = 2*9;
       if(order == 2 )  avg_nonzero_per_row = 2*5;
    }
    if (ndim == 7) {
/*     nrow = (global_1.nx + 3 << 1) * (ny + 2); */
       nrow = 2*(global_1.nx + 3 ) * (ny + 2);
       if(order == 4 || order == 6) avg_nonzero_per_row = 2*9;
       if(order == 2 )  avg_nonzero_per_row = 2*5;
    }
    if (ndim == 6) {
/*     nrow = (global_1.nx + 2 << 1) * (ny + 2) * (nz + 2); */
       nrow = 2*(global_1.nx + 2 ) * (ny + 2) * (nz + 2);
       avg_nonzero_per_row = 2*5;
    }
    if (ndim == 5) {
/*     nrow = (global_1.nx + 2 << 1) * (ny + 2); */
       nrow = 2*(global_1.nx + 2 ) * (ny + 2);
       avg_nonzero_per_row = 2*5;
    }
    if (ndim == 4) {
       nrow = (global_1.nx + 2) * (ny + 2);
       avg_nonzero_per_row = 5;
       if(example_no == 28) {
          if(order == 4 || order == 6) avg_nonzero_per_row = 9;
          if(order == 2 )  avg_nonzero_per_row = 5;
        }
    }
    if (ndim == 3) {
       nrow = global_1.nx * ny * nz;
       avg_nonzero_per_row = 7;
    }
    if (ndim == 2) {
       nrow = global_1.nx * ny;
       avg_nonzero_per_row = 5;
    }
    if (ndim == 1) {
       nrow = global_1.nx;
       avg_nonzero_per_row = 3;
    }

/* Define the equations handeled by the current processor */

    if (read_update == 1) {
	if (ndim == 1) {
	    AZ_read_update(&n_update, &update, proc_config,  nrow, 1, 0);
		}
	/* endif for NDIM=1 */

	if (ndim == 2) {
	    n_update = (nxend - nxst + 1) * (nyend - nyst + 1);
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nyst - 1) * global_1.nx + nxst - 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		update[i] = i + n0;
	    }
	}
        /* endif for NDIM=2 */

	if (ndim == 3) {
	    n_update = (nxend - nxst + 1) * (nyend - nyst + 1) * (nzend - 
		    nzst + 1);
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nzst - 1) * ny * global_1.nx + nxst - 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		update[i] = i + n0;
	    }
	}
	/* endif for NDIM=3 */

	if (ndim == 4) {
	    n_update = (nxend - nxst + 1) * (nyend - nyst + 1);
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nyst - 1) * (global_1.nx + 2) + nxst - 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
            update[i] = i + n0;
	    }
	}
	/* endif for NDIM=4 */

	if (ndim == 5) {
	    n_update = 2*(nxend - nxst + 1 ) * (nyend - nyst + 1);
/*	    n_update = (nxend - nxst + 1 << 1) * (nyend - nyst + 1); */
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nyst - 1) * (global_1.nx + 2) + nxst - 1 << 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		update[i] = i + n0;
	    }
	}
	/* endif for NDIM=5 */

	if (ndim == 6) {
/*   Note:  the domain is sub-devided into subdomains */
/*          in the z direction */
/*	    n_update = (nxend - nxst + 1 << 1) * (nyend - nyst + 1) * (
	    nzend - nzst + 1); */
	    n_update = 2*(nxend - nxst + 1 ) * (nyend - nyst + 1) * (
		    nzend - nzst + 1);
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nzst - 1) * (ny + 2) * (global_1.nx + 2) + nxst - 1 << 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
  	    update[i] = i + n0;
	    }
	/*           nrow = 2*(NX+2)*(NY+2)*(NZ+2) */
	}
	/* endif for NDIM=6 */

	if (ndim == 7) {
	    n_update = 2*(nxend - nxst + 1 ) * (nyend - nyst + 1);
/*	    n_update = (nxend - nxst + 1 << 1) * (nyend - nyst + 1); */
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nyst - 1) * (global_1.nx + 3) + nxst - 1 << 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		update[i] = i + n0;
	    }
	}
	/* endif for NDIM=7 */

	if (ndim == 8) {
/*	    n_update = (nxend - nxst + 1 << 1) * (nyend - nyst + 1); */
	    n_update = 2*(nxend - nxst + 1 ) * (nyend - nyst + 1);
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nyst - 1) * (global_1.nx + 2) + nxst - 1 << 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
	        update[i] = i + n0;
	    }
	}
	/* endif for NDIM=8 */

	if (ndim == 9) {
	    n_update = 2*(nxend - nxst + 1 ) * (nyend - nyst + 1);
/*          n_update = (nxend - nxst + 1 << 1) * (nyend - nyst + 1); */
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nyst - 1) * (global_1.nx + 4) + nxst - 1 << 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		update[i] = i + n0;
	    }
	}
	/* endif for NDIM=9 */

	if (ndim == 10) {
	    n_update = 2*(nxend - nxst + 1 ) * (nyend - nyst + 1);
            update     = (int *) AZ_allocate((n_update)*sizeof(int));
	    n0 = (nyst - 1) * (global_1.nx + 4) + nxst - 1 << 1;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		update[i] = i + n0;
	    }
	}
	/* endif for NDIM=10 */
    }
	/*  endif for  read_update=1 */

    if (read_update == 2) {
	AZ_read_update(&n_update, &update, proc_config,  nrow, 1, 2);
    }

    if (read_update == 3) {

/*  This corresponds to the case AZ_file of AZTEC (see page 53 of thr User's Guide Ver.2.1). 
    Each processor reads in the (global) indices of the equations assigned to it. */ 

    switch (me) {
	case 0:  goto L430;
	case 1:  goto L431;
	case 2:  goto L432;
	case 3:  goto L433;
	case 4:  goto L434;
	case 5:  goto L435;
	case 6:  goto L436;
	case 7:  goto L437;
	case 8:  goto L438;
	case 9:  goto L439;
	case 10:  goto L440;
	case 11:  goto L441;
	case 12:  goto L442;
	case 13:  goto L443;
	case 14:  goto L444;
	case 15:  goto L445;
    }
L430:    
        fp= fopen(".update_30","r");
        if(fp==0) {printf("\n error opening .update_30"); exit;} 
        goto L450;
L431:    
        fp = fopen(".update_31","r");
        if(fp==0) {printf("\n error opening .update_31"); exit;} 
        goto L450;
L432:    
        fp = fopen(".update_32","r");
        if(fp ==0) {printf("\n error opening .update_32"); exit;} 
        goto L450;
L433:    
        fp  = fopen(".update_33","r");
        if(fp ==0) {printf("\n error opening .update_33"); exit;} 
        goto L450;
L434:    
        fp  = fopen(".update_34","r");
        if(fp ==0) {printf("\n error opening .update_34"); exit;} 
        goto L450;
L435:    
        fp  = fopen(".update_35","r");
        if(fp  =0) {printf("\n error opening .update_35"); exit;} 
        goto L450;
L436:    
        fp  = fopen(".update_36","r");
        if(fp ==0) {printf("\n error opening .update_36"); exit;} 
        goto L450;
L437:    
        fp  = fopen(".update_37","r");
        if(fp ==0) {printf("\n error opening .update_37"); exit;} 
        goto L450;
L438:    
        fp  = fopen(".update_38","r");
        if(fp ==0) {printf("\n error opening .update_38"); exit;} 
        goto L450;
L439:    
        fp  = fopen(".update_39","r");
        if(fp ==0) {printf("\n error opening .update_39"); exit;} 
        goto L450;
L440:    
        fp  = fopen(".update_40","r");
        if(fp ==0) {printf("\n error opening .update_40"); exit;} 
        goto L450;
L441:    
        fp  = fopen(".update_41","r");
        if(fp ==0) {printf("\n error opening .update_41"); exit;} 
        goto L450;
L442:    
        fp  = fopen(".update_42","r");
        if(fp ==0) {printf("\n error opening .update_42"); exit;} 
        goto L450;
L443:    
        fp  = fopen(".update_43","r");
        if(fp ==0) {printf("\n error opening .update_43"); exit;} 
L444:    
        fp  = fopen(".update_44","r");
        if(fp ==0) {printf("\n error opening .update_44"); exit;} 
        goto L450;
L445:    
        fp  = fopen(".update_45","r");
        if(fp ==0) {printf("\n error opening .update_45"); exit;} 
        goto L450;
L450:   

/*	ll = me + 30; */
	fscanf (fp, " %d", &n_update);
        i1 = n_update - 1;
	for (i = 0; i <= i1; ++i) {
		fscanf(fp,"%d",&update[i]);
	}
    }
    printf("\n  me= %d   N_update= %d    nrow = %d", me, n_update, nrow);

/*  Allocate memory  */
	total_nz = n_update*avg_nonzero_per_row + 1;
        bindx = (int *) AZ_allocate(total_nz*sizeof(int));
        val   = (double *) AZ_allocate(total_nz*sizeof(double));   
        rhsm = (double *) calloc(nrow+1, sizeof(double));
        if((val == NULL) && (total_nz != 0)) {
           printf( "\n Error:  Not enough space to create matrix\n");
           printf( "\n   Try reducing the variable 'avg_nonzero_per_row'\n"); 
           exit(1);
        }
        if(dim == 3)  {
           n_send = (global_1.nx +2)*(global_1.nx +2);
           if(example_no ==25) n_send = 2*(global_1.nx +2)*(global_1.nx +2);
        }
	if(dim == 2) {
           n_send = global_1.nx +2;
	   if(example_no == 21 )  n_send = 2*(global_1.nx+2);
           if(example_no == 23 ) n_send = 2*(global_1.nx+4);
	   if(example_no == 24 )  n_send = 2*(global_1.nx+4);
	   if(example_no == 26 || example_no ==27) n_send = 2*(global_1.nx+3);
           }

        index_s_xp  = (int *) calloc(n_send  , sizeof(int));
        index_s_xm  = (int *) calloc(n_send  , sizeof(int));
        index_s_yp  = (int *) calloc(n_send  , sizeof(int));
        index_s_ym  = (int *) calloc(n_send  , sizeof(int));
        index_r_xp  = (int *) calloc(n_send  , sizeof(int));
        index_r_xm  = (int *) calloc(n_send  , sizeof(int));
        index_r_yp  = (int *) calloc(n_send  , sizeof(int));
        index_r_ym  = (int *) calloc(n_send  , sizeof(int));

/* more allocations: */

	buf_r_xm  = (double *) calloc(n_send  , sizeof(double));
	buf_s_xm  = (double *) calloc(n_send  , sizeof(double));
    	buf_r_ym  = (double *) calloc(n_send  , sizeof(double));
        buf_s_ym  = (double *) calloc(n_send  , sizeof(double));
        buf_r_xp  = (double *) calloc(n_send  , sizeof(double));
        buf_s_xp  = (double *) calloc(n_send  , sizeof(double));
        buf_r_yp  = (double *) calloc(n_send  , sizeof(double));
        buf_s_yp  = (double *) calloc(n_send  , sizeof(double));

    if(dim == 3) {
        index_s_zp  = (int *) calloc(n_send  , sizeof(int));
        index_s_zm  = (int *) calloc(n_send  , sizeof(int));
        index_r_zp  = (int *) calloc(n_send  , sizeof(int));
        index_r_zm  = (int *) calloc(n_send  , sizeof(int));
        buf_r_zm  = (double *) calloc(n_send  , sizeof(double));
        buf_s_zm  = (double *) calloc(n_send  , sizeof(double));
        buf_r_zp  = (double *) calloc(n_send  , sizeof(double));
        buf_s_zp  = (double *) calloc(n_send  , sizeof(double));
    }

/*   Define the indices of the external boundry grid points */
/*   which communicate with the neigboring processors */

    switch (ndim) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L8;
	case 9:  goto L9;
	case 10:  goto L10;
    }
L1:
    example_1d(&nxst, &nxend, &global_1.nx, index_s_xm, index_s_xp, 
	    index_r_xm, index_r_xp, &ixp, &ixm);
    goto L60;
L2:
    example_2d(&nxst, &nxend, &nyst, &nyend, &global_1.nx, &ny, 
	    index_s_xm, index_s_xp, index_s_ym, index_s_yp, 
	    index_r_xm, index_r_xp, index_r_ym, index_r_yp, &ixp, &
	    ixm, &iyp, &iym);
    goto L60;
L3:
    example_3d(&nxst, &nxend, &nyst, &nyend, &nzst, &nzend, &global_1.nx, &
	    ny, &nz, index_s_xm, index_s_xp, index_s_ym, index_s_yp, 
	    index_s_zm, index_s_zp, index_r_xm, index_r_xp, 
	    index_r_ym, index_r_yp, index_r_zm, index_r_zp, &ixp, &
	    ixm, &iyp, &iym, &izp, &izm);
    goto L60;
L4:
    if (example_no != 28) {
	i1 = global_1.nx + 2;
	i2 = ny + 2;
	example_2d(&nxst, &nxend, &nyst, &nyend, &i1, &i2, index_s_xm,
		 index_s_xp, index_s_ym, index_s_yp, index_r_xm, 
		index_r_xp, index_r_ym, index_r_yp, &ixp, &ixm, &iyp, 
		&iym);
    }
    if (example_no == 28) {
	if (order == 2) {
	    nx_st = 1;
	    nx_end = global_1.nx + 2;
	    ny_st = 1;
	    ny_end = ny + 2;
	}
	if (order == 4 || order == 6) {
	    nx_st = 0;
	    nx_end = global_1.nx + 3;
	    ny_st = 0;
	    ny_end = ny + 3;
	}
	i1 = global_1.nx + 2;
	i2 = ny + 2;
	example_2d_all(&nxst, &nxend, &nyst, &nyend, &i1, &i2, 
		index_s_xm, index_s_xp, index_s_ym, index_s_yp, 
		index_r_xm, index_r_xp, index_r_ym, index_r_yp, &ixp, 
		&ixm, &iyp, &iym, &nx_st, &nx_end, &ny_st, 
		&ny_end);
    }
    printf("\n  me= %d   ixm= %d   ixp= %d   iym= %d   iyp= %d", 
		    me, ixm, ixp, iym, iyp);
    goto L60;
L5:
    example_2d_img(&nxst, &nxend, &nyst, &nyend, &global_1.nx, &ny, 
	    index_s_xm, index_s_xp, index_s_ym, index_s_yp, 
	    index_r_xm, index_r_xp, index_r_ym, index_r_yp, &ixp, 
	    &ixm, &iyp, &iym);
    printf("\n  me= %d   ixm= %d   ixp= %d   iym= %d   iyp= %d", 
		    me, ixm, ixp, iym, iyp);
    goto L60;
L6:
    example_3d_img(&nxst, &nxend, &nyst, &nyend, &nzst, &nzend, &
	    global_1.nx, &ny, &nz, index_s_xm, index_s_xp, index_s_ym, 
	    index_s_yp, index_s_zm, index_s_zp, index_r_xm, 
	    index_r_xp, index_r_ym, index_r_yp, index_r_zm, 
	    index_r_zp, &ixp, &ixm, &iyp, &iym, &izp, &izm);
    printf("\n  me= %d   ixm= %d   ixp= %d   iym= %d   iyp= %d   izm= %d   izp= %d", 
		    me, ixm, ixp, iym, iyp, izm, izp);
    goto L60;
L7:
    if (order == 2) {
	nx_st = 1;
	nx_end = global_1.nx + 3;
	ny_st = 1;
	ny_end = ny + 2;
    }
    if (order == 4 || order == 6) {
	nx_st = 0;
	nx_end = global_1.nx + 4;
	ny_st = 0;
	ny_end = ny + 3;
    }
    i1 = global_1.nx + 1;
    example_2d_img_all(&nxst, &nxend, &nyst, &nyend, &i1, &ny, 
	    index_s_xm, index_s_xp, index_s_ym, index_s_yp, 
	    index_r_xm, index_r_xp, index_r_ym, index_r_yp, &ixp, &
	    ixm, &iyp, &iym, &nx_st, &nx_end, &ny_st, 
	    &ny_end);
    printf("\n  me= %d   ixm= %d   ixp= %d   iym= %d   iyp= %d", 
		    me, ixm, ixp, iym, iyp);
    goto L60;
L8:
    if (order == 2) {
	nx_st = 1;
	nx_end = global_1.nx + 2;
	ny_st = 1;
	ny_end = ny + 3;
    }
    if (order == 4 || order == 6) {
	nx_st = 0;
	nx_end = global_1.nx + 3;
	ny_st = 0;
	ny_end = ny + 4;
    }
    i1 = ny + 1;
    example_2d_img_all(&nxst, &nxend, &nyst, &nyend, &global_1.nx, &i1, 
	    index_s_xm, index_s_xp, index_s_ym, index_s_yp, 
	    index_r_xm, index_r_xp, index_r_ym, index_r_yp, &ixp, 
	    &ixm, &iyp, &iym, &nx_st, &nx_end, &ny_st, 
	    &ny_end);
    printf("\n  me= %d   ixm= %d   ixp= %d   iym= %d   iyp= %d", 
		    me, ixm, ixp, iym, iyp);
    goto L60;
L9:
    if (order == 2) {
	nx_st = 1;
	nx_end = global_1.nx + 4;
	ny_st = 1;
	ny_end = ny + 3;
    }
    if (order == 4 || order == 6) {
	nx_st = 0;
	nx_end = global_1.nx + 5;
	ny_st = 0;
	ny_end = ny + 4;
    }
    i1 = global_1.nx + 2;
    i2 = ny + 1;
    example_2d_img_all(&nxst, &nxend, &nyst, &nyend, &i1, &i2, 
	    index_s_xm, index_s_xp, index_s_ym, index_s_yp, 
	    index_r_xm, index_r_xp, index_r_ym, index_r_yp, &ixp, 
	    &ixm, &iyp, &iym, &nx_st, &nx_end, &ny_st, 
	    &ny_end);
    printf("\n  me= %d   ixm= %d   ixp= %d   iym= %d   iyp= %d", 
		    me, ixm, ixp, iym, iyp);
    goto L60;
L10:
    if (order == 2) {
	nx_st = 1;
	nx_end = global_1.nx + 4;
	ny_st = 1;
	ny_end = ny + 4;
    }
    if (order == 4 || order == 6) {
	nx_st = 0;
	nx_end = global_1.nx + 5;
	ny_st = 0;
	ny_end = ny + 5;
    }
    i1 = global_1.nx + 2;
    i2 = ny + 2;
    example_2d_img_all(&nxst, &nxend, &nyst, &nyend, &i1, &i2, 
	    index_s_xm, index_s_xp, index_s_ym, index_s_yp, 
	    index_r_xm, index_r_xp, index_r_ym, index_r_yp, &ixp, 
	    &ixm, &iyp, &iym, &nx_st, &nx_end, &ny_st, 
	    &ny_end);
    printf("\n  me= %d   ixm= %d   ixp= %d   iym= %d   iyp= %d", 
		    me, ixm, ixp, iym, iyp);
L60:

    nxloc = nxend - nxst + 1;
    if (ndim != 1) {
	nyloc = nyend - nyst + 1;
    }
    if (ndim == 3 || ndim == 6) {
	nzloc = nzend - nzst + 1;
    }
    bindx[0] = n_update + 1;
    nend = n_update;

/*  Define and save the matrix following 'AZTEC' sparse matrix format */
/* -------------------------------------------------------------------- */
    if (example_no == 21) {

/*  example 21 - 2D Helmholtz problem */

	i1 = n_update - 1;
	for (i = 0; i <= i1; i += 2) {
	    create_matrix_row_5pt_example21_2rows(&update[i], &i, val, 
		    bindx, &kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc, &kn1, 
		    &kn2, &kn3);
	}
	goto L290;
    }
/* endif example21 */

    if (example_no == 23) {
	i1 = n_update - 1;
	for (i = 0; i <= i1; i += 2) {
	    create_matrix_row_9pt_example23_2rows(&update[i], &i, val, 
		    bindx, &kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc, &kn1, 
		    &kn2, &kn3, &order);
	}
	goto L290;
    }
/* endif example23 */

    if (example_no == 24) {
/*  example 24 - Marmousi problem */
/*	pi = 3.14159265359793; */
/*      pi = atan(1.) * 4.; */
/*	hx = (xmax - xmin) / (global_1.nx - 1); */
/*	hy = (ymax - ymin) / (ny - 1); */
        i1 = n_update - 1;
	for (i = 0; i <= i1; i += 2) {
	    create_matrix_row_5_9_pt_example24_2rows(&update[i], &i, val, 
		    bindx, &kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc, cxy, 
		    &xmin, &xmax, &ymin, &ymax, &global_1.nx, &ny,
                    &order,&freq);
	}
	goto L290;
    }
/* endif example24 */

    if (example_no == 25) {
/*  example 25 - 3D Helmholtz problem */
	i1 = n_update - 1;
	for (i = 0; i <= i1; i += 2) {
	    create_matrix_row_7pt_example25_2rows(&update[i], &i, val, 
		    bindx, &kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, 
		    &nzloc, &kn1, &kn2, &kn3, &xmin, &xmax, &ymin, &ymax, 
		    &zmin, &zmax, &global_1.nx, &ny, &nz);
	}
	goto L290;
    }
/* endif create matrix for example25 */

    if (example_no == 26) {
	i1 = n_update - 1;
	for (i = 0; i <= i1; i += 2) {
	    create_matrix_row_9pt_example26_2rows(&update[i], &i, val, 
		    bindx, &kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc, &kn1, 
		    &kn2, &kn3, &order);
	}
	goto L290;
    }
/* endif example26 */

    if (example_no == 27) {
	i1 = n_update - 1;
	for (i = 0; i <= i1; i += 2) {
	    create_matrix_row_9pt_example27_2rows(&update[i], &i, val, 
		    bindx, &kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc, &kn1, 
		    &kn2, &kn3, &order);
	}
	goto L290;
    }
/* endif example27 */

    i1 = nend - 1;
    for (i = 0; i <= i1; ++i) {
	switch (example_no) {
	    case 1:  goto L31;
	    case 2:  goto L32;
	    case 3:  goto L33;
	    case 4:  goto L34;
	    case 5:  goto L35;
	    case 6:  goto L36;
	    case 7:  goto L37;
	    case 8:  goto L38;
	    case 9:  goto L39;
	    case 10:  goto L40;
	    case 11:  goto L41;
	    case 12:  goto L42;
	    case 13:  goto L43;
	    case 14:  goto L44;
	    case 15:  goto L45;
	    case 16:  goto L46;
	    case 17:  goto L47;
	    case 18:  goto L48;
	    case 19:  goto L49;
	    case 20:  goto L50;
	    case 21:  goto L51;
	    case 22:  goto L52;
	    case 23:  goto L290;
	    case 24:  goto L290;
	    case 25:  goto L290;
	    case 26:  goto L290;
	    case 27:  goto L290;
	    case 28:  goto L53;
	}
L31:
	create_matrix_row_7pt_example1(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L32:
	create_matrix_row_7pt_example2(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L33:
	create_matrix_row_7pt_example3(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L34:
	create_matrix_row_7pt_example4(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L35:
	create_matrix_row_7pt_example5(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L36:
	create_matrix_row_7pt_example6(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L37:
	create_matrix_row_7pt_example7(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L38:
	create_matrix_row_7pt_example8(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L39:
	create_matrix_row_5pt_example9(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc);
	goto L250;
L40:
	create_matrix_row_5pt_example10(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nxloc);
	goto L250;
L41:
	create_matrix_row_5pt_example11(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc);
	goto L250;
L42:
	create_matrix_row_7pt_example12(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L43:
	create_matrix_row_7pt_example13(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L44:
	create_matrix_row_7pt_example14(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L45:
	create_matrix_row_5pt_example15(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc);
	goto L250;
L46:
	create_matrix_row_5pt_example16(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc);
	goto L250;
L47:
	create_matrix_row_5pt_example17(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc);
	goto L250;
L48:
	create_matrix_row_7pt_example18(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L49:
	create_matrix_row_7pt_example19(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L50:
	create_matrix_row_5pt_example20(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc, &kn1);
	goto L250;
L51:
	goto L250;
L52:
	create_matrix_row_7pt_example22(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nzst, &nxloc, &nyloc, &nzloc);
	goto L250;
L53:
	create_matrix_row_9pt_example28(&update[i], &i, val, bindx, &
		kmax, rhsm, &nxst, &nyst, &nxloc, &nyloc, &kn1, &order);
L250:
	;
    }

L290:
    printf("\n  me= %d   Kmax= %d   update(i)= %d   NXloc= %d   NYloc= %d", 
		    me, kmax, update[nend-1], nxloc, nyloc);

/*  Define the averaging weights */
    iem  = (int *) calloc(nrow , sizeof(int));
    iema = (int *) calloc(nrow , sizeof(int));
    ecoef = (double *) calloc(nrow , sizeof(double));

/*  i1 = n - 1;
    for (i = 0; i <= i1; ++i) {
	iem[i] = 0;
    } */
    i1 = n_update - 1;
    for (i = 0; i <= i1; ++i) {
	iem[update[i]] = 1;
	i2 = bindx[i + 1] - 1;
	for (j = bindx[i]; j <= i2; ++j) {
	    iem[bindx[j]] = 1;
	}
    }

    MPI_Allreduce(iem, iema, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

    i1 = n - 1;
    for (i = 0; i <= i1; ++i) {
	ecoef[i] = 1. / iema[i];
/*        if(me.eq.0) */
/*    1   write(NOUT,*) 'me=',me,'  I=',I, ' IEMA(I)=',IEMA(I), */
/*    1   ' ECOEF(I)',ECOEF(I) */
    }

/*   Tranform the matrix to a local matrix and define various variables */
/*   to conform with AZTEC */

    AZ_transform(proc_config, &external, bindx, val, update, 
	    &update_index, &extern_index, &data_org, n_update, 0, 
	    0, 0, 0, AZ_MSR_MATRIX);

    n_update_internal = data_org[1];
    n_update_border = data_org[2];
    n_update_external = data_org[3];
    n_neigh = data_org[7];
    n_total_send = data_org[8];
    n_id_neigh = data_org[12];
    n_rec_length = data_org[262];
    n_send_length = data_org[512];
    n_send_list = data_org[762];

/*     write(NOUT,*)'me=',me,'  N_update_internal=',N_update_internal */
/*    1          ,'  N_update_border=',N_update_border */
/*    2          ,'  N_update_external=',N_update_external */
/*    3          ,'  N_neigh=',N_neigh */
/*    4          ,'  N_total_send=', N_total_send */
/*    5          ,'  N_id_neigh=',N_id_neigh */
/*    6          ,'  N_rec_length=',N_rec_length */
/*    7          ,'  N_send_length=',N_send_length */
/*    8          ,'  N_send_list=',N_send_list */


    printf("\n  me= %d   N_update_internal= %d   N_update_border= %d",
            me, n_update_internal,n_update_border);
    printf("   N_update_external= %d   N_neigh= %d", 
            n_update_external, n_neigh);
    printf("\n  N_total_send= %d   N_id_neigh= %d   N_rec_length= %d   N_send_length= %d   N_send_list= %d", 
            n_total_send, n_id_neigh, n_rec_length, n_send_length, n_send_list);

/*  Allocate memory  */

        update_inv = (int *) calloc(nrow , sizeof(int));
        xnew = (double *) calloc(nrow , sizeof(double));
        xnewl = (double *) calloc(nrow , sizeof(double));
        err  = (double *) calloc(nrow , sizeof(double));
        qnewl= (double *) calloc(nrow , sizeof(double));
        resid= (double *) calloc(nrow , sizeof(double));
        resid1= (double *) calloc(nrow , sizeof(double));
        xnewl_loc = (double *) calloc(n_update + data_org[AZ_N_external] , sizeof(double));
        xnewl_loc_init = (double *) calloc(n_update + data_org[AZ_N_external] , sizeof(double));
        rhsm_loc = (double *) calloc(n_update + data_org[AZ_N_external] , sizeof(double));
        pnewl_loc = (double *) calloc(n_update + data_org[AZ_N_external] , sizeof(double));
        qnewl_loc = (double *) calloc(n_update + data_org[AZ_N_external] , sizeof(double));
        snewl_loc = (double *) calloc(n_update + data_org[AZ_N_external] , sizeof(double));   
        ssum  = (double *) calloc(n_update  , sizeof(double));


/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

/*  Calculate the initial residual - RESIDUAL(0) */

    sumr = 0.;
    i1 = n_update - 1;
    for (i = 0; i <= i1; ++i) {
	d1 = rhsm[update[i]];
	sumr += d1 * d1;
    }
    MPI_Allreduce(&sumr, &sumrt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
    sumrt = sqrt(sumrt);
    printf("\n  RESIDUAL(0)= %.8f", sumrt);
    res0 = sumrt;

/*  Initiate vectors */
/*  i1 = n - 1;
    for (i = 0; i <= i1; ++i) {
	xnewl[i] = 0.;
    }  */

/*  i1 = n_update - 1;
    for (i = 0; i <= i1; ++i) {
	rhsm_loc[update_index[i]] = rhsm[update[i]];
    } */

    n_internal = n_update_internal;
    n_border = n_update_border;
    n_external = n_update_external;

    i1 = n_update - 1;
    for (i = 0; i <= i1; ++i) {
	update_inv[update_index[i]] = update[i];
	rhsm_loc[update_index[i]] = rhsm[update[i]];
/*  Define xnewl_loc[i] in case that xnewl_loc # 0. */
/*	xnewl_loc[i] = 0.;*/
    }

    t1 = MPI_Wtime();
    printf("\n  %s","   ");

/*  Define a normalization factor */

    i1 = n_update - 1;
    for (i = 0; i <= i1; ++i) {
	sum = val[i] * val[i];
	i2 = bindx[i + 1] - 1;
	for (j = bindx[i]; j <= i2; ++j) {
	    d1 = val[j];
	    sum += d1 * d1;
	}
	ssum[i] = omeg1 / sum;
/*      SSUM(I) = OMEG1 */
    }


/*  Start the initial forth and back sweeps of CGMN/CARP-CG */
/*  (required by CGMN/CARP-CG initial vectors S0( ),P0( ) ) */
/* -------------------------------------------------------------------- */
/*  Initiate the vector */
    i1 = n_update - 1;
    for (i = 0; i <= i1; ++i) {
	xnewl_loc_init[i] = xnewl_loc[i];
    }
    i1 = n_update_external - 1;
    for (i = 0; i <= i1; ++i) {
	xnewl_loc_init[extern_index[i]] = xnewl_loc[extern_index[
		i]];
    }
    for (icg = 1; icg <= 2; ++icg) {
	if (icg == 1) {
	    istart = 0;
	    iend = n_update - 1;
	    iinc = 1;
	} else {
	    istart = n_update - 1;
	    iend = 0;
	    iinc = -1;
	}
	i1 = iend;
	i2 = iinc;
	for (i = istart; i2 < 0 ? i >= i1 : i <= i1; i += i2) 
		{
	    sum = val[i] * xnewl_loc[i] - rhsm_loc[i];
	    i3 = bindx[i + 1] - 1;
	    for (j = bindx[i]; j <= i3; ++j) {
		sum += val[j] * xnewl_loc[bindx[j]];
	    }
	    cons = -sum * ssum[i];
	    xnewl_loc[i] += cons * val[i];
	    i3 = bindx[i + 1] - 1;
	    for (j = bindx[i]; j <= i3; ++j) {
		xnewl_loc[bindx[j]] += cons * val[j];
	    }
	}

/*  Transfer information: */

	if (nproc > 1) {
	    i2 = n_update - 1;
	    for (i = 0; i <= i2; ++i) {
		xnewl[update_inv[i]] = xnewl_loc[i];
	    }
	    i2 = n_update_external - 1;
	    for (i = 0; i <= i2; ++i) {
		xnewl[external[i]] = xnewl_loc[extern_index[i]];
	    }
	    if (dim == 3 ) {
		transf_3d(xnewl, ecoef, update_inv, &n_update, &
			comm_cart, &me, &rank_source_x, &rank_source_y, 
			&rank_source_z, &rank_dest_x, &rank_dest_y, &
			rank_dest_z, index_s_xm, index_s_xp, 
			index_s_ym, index_s_yp, index_s_zm, 
			index_s_zp, index_r_xm, index_r_xp, 
			index_r_ym, index_r_yp, index_r_zm, 
			index_r_zp, &ixp, &ixm, &iyp, &iym, &izp, &izm, &
			ndim);
	    } else {
		if (dim == 2 ) {
                transf_2d(xnewl, ecoef, update_inv, &n_update, &
	                comm_cart, &me, &rank_source_x, &
		        rank_source_y, &rank_dest_x, &rank_dest_y, 
		        index_s_xm, index_s_xp, index_s_ym, 
		        index_s_yp, index_r_xm, index_r_xp, 
		        index_r_ym, index_r_yp, &ixp, &ixm, &iyp, &
		        iym, &ndim);
		}
	    }
	    i2 = n_internal + n_border - 1;
	    for (i = n_internal; i <= i2; ++i) {
		xnewl_loc[i] = xnewl[update_inv[i]];
	    }
	    i2 = n_update_external - 1;
	    for (i = 0; i <= i2; ++i) {
		xnewl_loc[extern_index[i]] = xnewl[external[i]];
	    }
	}
/*  end if   (for: Nproc.gt.1) */
    }
/*  end of the first ICG (forth and back sweeps) loop. */

/*  Start the  main forth and back sweeps of CGMN/CARP-CG: */
/* --------------------------------------------------------- */
    i2 = kcg;
    for (ikcg = 1; ikcg <= i2; ++ikcg) {
	if (ikcg == 1) {

/*  Define the initial values of P( ), S( ) and X( ): */
/* ---------------------------------------------------- */
	    s2 = 0.;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		pnewl_loc[i] = xnewl_loc[i] - xnewl_loc_init[i];
		snewl_loc[i] = pnewl_loc[i];
		d1 = snewl_loc[i];
		s2 += d1 * d1;
		qnewl_loc[i] = pnewl_loc[i];
		xnewl_loc[i] = xnewl_loc_init[i];
	    }

	    i1 = n_update_external - 1;
	    for (i = 0; i <= i1; ++i) {
		pnewl_loc[extern_index[i]] = xnewl_loc[extern_index[
			i]] - xnewl_loc_init[extern_index[i]];
		snewl_loc[extern_index[i]] = pnewl_loc[extern_index[
			i]];
		qnewl_loc[extern_index[i]] = pnewl_loc[extern_index[
			i]];
		xnewl_loc[extern_index[i]] = 0.;
	    }
            MPI_Allreduce(&s2, &s2s, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
                      MPI_COMM_WORLD);

	} else {

/*   For IKCG >1:  update S( ) and P( ), initiate Q( ) */
/*   and calculate S**2 and BETA: */
/* ------------------------------------------------------- */

/* for IKCG >1 */
	    s2p = 0.;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		snewl_loc[i] -= alpha * qnewl_loc[i];
		d1 = snewl_loc[i];
		s2p += d1 * d1;
	    }
	    i1 = n_update_external - 1;
	    for (i = 0; i <= i1; ++i) {
		snewl_loc[extern_index[i]] -= alpha * qnewl_loc[
			extern_index[i]];
	    }

/*   Sum up S**2 from all processors: */

            MPI_Allreduce(&s2p, &s2ps, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
                      MPI_COMM_WORLD);

/*            write(6,*)'  IKCG=',IKCG,' me=',me,' *BETA= ',BETA */
/*    1                ,'  S2P= ', S2P,'  S2S=',S2PS */
	    beta = s2ps / s2s;
	    s2s = s2ps;

/*   Update P( ) and initiate Q( ) : */

	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		pnewl_loc[i] = snewl_loc[i] + beta * pnewl_loc[i];
		qnewl_loc[i] = pnewl_loc[i];
/*              write(IP,1011) IKCG,i,BETA,QNEWL_LOC(i) */
	    }
	    i1 = n_update_external - 1;
	    for (i = 0; i <= i1; ++i) {
		pnewl_loc[extern_index[i]] = snewl_loc[extern_index[
			i]] + beta * pnewl_loc[extern_index[i]];
		qnewl_loc[extern_index[i]] = pnewl_loc[extern_index[
			i]];
/*          write(IP,1012) IKCG,i,BETA,QNEWL_LOC(extern_index(i)) */
	    }
	}

/* Ccalculate Q( ):  do the forth sweep */

	i1 = n_update - 1;
	for (i = 0; i <= i1; ++i) {
	    sum = val[i] * qnewl_loc[i];
	    i3 = bindx[i + 1] - 1;
	    for (j = bindx[i]; j <= i3; ++j) {
		sum += val[j] * qnewl_loc[bindx[j]];
	    }
	    cons = -sum * ssum[i];
	    qnewl_loc[i] += cons * val[i];
	    i3 = bindx[i + 1] - 1;
	    for (j = bindx[i]; j <= i3; ++j) {
		qnewl_loc[bindx[j]] += cons * val[j];
	    }
	}

/*  Transfer information (for Nproc>1): */

	if (nproc > 1) {
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		qnewl[update_inv[i]] = qnewl_loc[i];
	    }
	    i1 = n_update_external - 1;
	    for (i = 0; i <= i1; ++i) {
		qnewl[external[i]] = qnewl_loc[extern_index[i]];
	    }
	    if (dim == 3 ) {
		transf_3d(qnewl, ecoef, update_inv, &n_update, &
			comm_cart, &me, &rank_source_x, &rank_source_y, 
			&rank_source_z, &rank_dest_x, &rank_dest_y, &
			rank_dest_z, index_s_xm, index_s_xp, 
			index_s_ym, index_s_yp, index_s_zm, 
			index_s_zp, index_r_xm, index_r_xp, 
			index_r_ym, index_r_yp, index_r_zm, 
			index_r_zp, &ixp, &ixm, &iyp, &iym, &izp, &izm, &
			ndim);
	    } else {
		if (dim == 2 ) {
		transf_2d(qnewl, ecoef, update_inv, &n_update, &
			comm_cart, &me, &rank_source_x, &
		        rank_source_y, &rank_dest_x, &rank_dest_y, 
			index_s_xm, index_s_xp, index_s_ym, 
			index_s_yp, index_r_xm, index_r_xp, 
			index_r_ym, index_r_yp, &ixp, &ixm, &iyp, &
			iym, &ndim);
		}
	    }
	    i1 = n_update - 1;
	    for (i = n_internal; i <= i1; ++i) {
		qnewl_loc[i] = qnewl[update_inv[i]];
	    }
	    i1 = n_update_external - 1;
	    for (i = 0; i <= i1; ++i) {
		qnewl_loc[extern_index[i]] = qnewl[external[i]];
	    }
	}

/*  Do the back sweep */

	for (i = n_update - 1; i >= 0; --i) {
	    sum = val[i] * qnewl_loc[i];
	    i1 = bindx[i + 1] - 1;
	    for (j = bindx[i]; j <= i1; ++j) {
		sum += val[j] * qnewl_loc[bindx[j]];
	    }
	    cons = -sum * ssum[i];
	    qnewl_loc[i] += cons * val[i];
	    i1 = bindx[i + 1] - 1;
	    for (j = bindx[i]; j <= i1; ++j) {
		qnewl_loc[bindx[j]] += cons * val[j];
	    }
	}

/*  Transfer information (for Nproc>1): */

	if (nproc > 1) {
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		qnewl[update_inv[i]] = qnewl_loc[i];
	    }
	    i1 = n_update_external - 1;
	    for (i = 0; i <= i1; ++i) {
		qnewl[external[i]] = qnewl_loc[extern_index[i]];
	    }
            if (dim == 3 ) {
		transf_3d(qnewl, ecoef, update_inv, &n_update, &
			comm_cart, &me, &rank_source_x, &rank_source_y, 
			&rank_source_z, &rank_dest_x, &rank_dest_y, &
			rank_dest_z, index_s_xm, index_s_xp, 
			index_s_ym, index_s_yp, index_s_zm, 
			index_s_zp, index_r_xm, index_r_xp, 
			index_r_ym, index_r_yp, index_r_zm, 
			index_r_zp, &ixp, &ixm, &iyp, &iym, &izp, &izm, &
			ndim);
	    } else {
	    if (dim == 2 ) {
		transf_2d(qnewl, ecoef, update_inv, &n_update, &
		        comm_cart, &me, &rank_source_x, &
			rank_source_y, &rank_dest_x, &rank_dest_y, 
			index_s_xm, index_s_xp, index_s_ym, 
			index_s_yp, index_r_xm, index_r_xp, 
			index_r_ym, index_r_yp, &ixp, &ixm, &iyp, &
			iym, &ndim);
		}
	    }
	    i1 = n_update - 1;
	    for (i = n_internal; i <= i1; ++i) {
		qnewl_loc[i] = qnewl[update_inv[i]];
	    }
	    i1 = n_update_external - 1;
	    for (i = 0; i <= i1; ++i) {
		qnewl_loc[extern_index[i]] = qnewl[external[i]];
	    }
	}

	i1 = n_update - 1;
	for (i = 0; i <= i1; ++i) {
	    qnewl_loc[i] = pnewl_loc[i] - qnewl_loc[i];
	}
	i1 = n_update_external - 1;
	for (i = 0; i <= i1; ++i) {
	    qnewl_loc[extern_index[i]] = pnewl_loc[extern_index[i]
		    ] - qnewl_loc[extern_index[i]];
	}

/*   Calculate ALPHA: */

	pq = 0.;
	i1 = n_update - 1;
	for (i = 0; i <= i1; ++i) {
	    pq += pnewl_loc[i] * qnewl_loc[i];
	}
        MPI_Allreduce(&pq, &pqs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
                      MPI_COMM_WORLD);
	alpha = s2s / pqs;

/*  Update  X( ): */

	i1 = n_update - 1;
	for (i = 0; i <= i1; ++i) {
	    xnewl_loc[i] += alpha * pnewl_loc[i];
	}
	i1 = n_update_external - 1;
	for (i = 0; i <= i1; ++i) {
	    xnewl_loc[extern_index[i]] += alpha * pnewl_loc[
		    extern_index[i]];
	}


/*  Calculate the residual: */
/* ---------------------------- */
/*      if (ikcg / 1 * 1 == ikcg) { */
/*	if (ikcg / 5 * 5 == ikcg) { */
	if (ikcg / 10 * 10 == ikcg) {
	    sumr = 0.;
	    sumr = 0.;
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		sum = val[i] * xnewl_loc[i] - rhsm_loc[i];
		i3 = bindx[i + 1] - 1;
		for (j = bindx[i]; j <= i3; ++j) {
		    sum += val[j] * xnewl_loc[bindx[j]];
		}
		d1 = sum;
		sumr += d1 * d1;
	    }
            MPI_Allreduce(&sumr, &sumrt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
                      MPI_COMM_WORLD);
	    sumrt = sqrt(sumrt);
	    if (me == 0) {
                printf("\n  IKCG= %d   RESIDUAL= %.9g     RESIDUAL/Res0= %.9g",
                ikcg, sumrt, d1 = sumrt / res0);
	    }
            if(ierr == 1 ) {   /* Calculate the l2 error and errmax  */

/*  Get XNEW( ) - the global vector of solution */
	    i1 = n - 1;
	    for (i = 0; i <= i1; ++i) {
		xnewl[i] = 0.;
	    }
	    i1 = n_update - 1;
	    for (i = 0; i <= i1; ++i) {
		xnewl[update[i]] = xnewl_loc[update_index[i]];
	    }

            MPI_Reduce(xnewl, xnew, n, MPI_DOUBLE_PRECISION, MPI_SUM, 
                      0, MPI_COMM_WORLD);

	    if (me == 0) {
		ipr=0;    
		switch (example_no) {
		    case 1:  goto L311;
		    case 2:  goto L312;
		    case 3:  goto L313;
		    case 4:  goto L314;
		    case 5:  goto L315;
		    case 6:  goto L316;
		    case 7:  goto L317;
		    case 8:  goto L318;
		    case 9:  goto L319;
		    case 10:  goto L320;
		    case 11:  goto L321;
		    case 12:  goto L322;
		    case 13:  goto L323;
		    case 14:  goto L324;
		    case 15:  goto L325;
		    case 16:  goto L380;
		    case 17:  goto L380;
		    case 18:  goto L380;
		    case 19:  goto L380;
		    case 20:  goto L330;
		    case 21:  goto L380;
		    case 22:  goto L380;
		    case 23:  goto L380;
		    case 24:  goto L380;
		    case 25:  goto L380;
		    case 26:  goto L336;
		    case 27:  goto L337;
		    case 28:  goto L338;
		}
L311:
		example1_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 1:\n");
		}
		goto L370;
L312:
		example2_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 2:\n");
		}
		goto L370;
L313:
		example3_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 3:\n");
		}
		goto L370;
L314:
		example4_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 4:\n");
		}
		goto L370;
L315:
		example5_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 5:\n");
		}
		goto L370;
L316:
		example6_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 6:\n");
		}
		goto L370;
L317:
		example7_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 7:\n");
		}
		goto L370;
L318:
		example8_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 8:\n");
		}
		goto L370;
L319:
		example9_err(xnew, err, &global_1.nx, &ny, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 9:\n");
		}
		goto L370;
L320:
		example10_err(xnew, err, &global_1.nx, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 10:\n");
		}
		goto L370;
L321:
		example11_err(xnew, err, &global_1.nx, &ny, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 11:\n");
		}
		goto L370;
L322:
		example12_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 12:\n");
		}
		goto L370;
L323:
		example13_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 13:\n");
		}
		goto L370;
L324:
		example14_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 14:\n");
		}
		goto L370;
L325:
		example15_err(xnew, err, &global_1.nx, &ny, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 15:\n");
		}
		goto L370;
L330:
		example20_err(xnew, err, &global_1.nx, &ny, &sumerr, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 20:");
		    fprintf(f9,"\n    from Erlangga, Vuik, Oosterlee, 2004 paper: \n");
		}
		goto L370;
L336:
		example26_err(xnew, err, &global_1.nx, &ny, &sumerr1, &
			sumerr2, &kn1, &kn2, &kn3, &order, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 26:");
		    fprintf(f9,"\n    from Erlangga & Turkel, 2011 paper \n");
		}
		goto L370;
L337:
		example27_err(xnew, err, &global_1.nx, &ny, &sumerr1, &
			sumerr2, &kn1, &kn2, &kn3, &order, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 27:");
		    fprintf(f9,"\n    from Erlangga & Turkel, 2011 paper \n");
		}
		goto L370;
L338:
		example28_err(xnew, err, &global_1.nx, &ny, &sumerr, &kn1, &order, &ipr);
		if (ikcg == 5 || ikcg == 10) {
                    fprintf(f9,"\n    Results for example 28:\n");
		}
L370:

/*  Calculate ERR_L2 */

		sum = 0.;
		i1 = n - 1;
		for (i = 0; i <= i1; ++i) {
		    d1 = err[i];
		    sum += d1 * d1;
		}
		if (example_no != 26 && example_no != 27) {
		    err_l2 = sqrt(sum) / sumerr;
		}
		if (example_no == 26 || example_no == 27) {
		    d1 = sumerr1;
		    d2 = sumerr2;
		    err_l2 = sqrt(sum) / sqrt(d1 * d1 + d2 * d2);
		}

/*  Calculate the maximum error */

		jk = 0;
		rmax = fabs(err[0]);
		i1 = n - 1;
		for (i = 1; i <= i1; ++i) {
		    if ((d1 = err[i], fabs(d1)) > rmax) {
			jk = i;
			rmax = (d1 = err[i], fabs(d1));
		    }
		}
  		fprintf(f9," I= %d    ERR_L2= %.9g    Err_max= %.9g \n", ikcg, err_l2, rmax);    

	    }  /* endif me = 0  */
            }  /* endif ierr=1   */
L380:
	    if (sumrt < residmax) {
		goto L150;
	    }
	}
    }
/*  End of CGMN/CARP-CG main loop */
L150:

/*  Get XNEW( ) - the global vector of solution */
    i2 = n - 1;
    for (i = 0; i <= i2; ++i) {
	xnewl[i] = 0.;
    }
    i2 = n_update - 1;
    for (i = 0; i <= i2; ++i) {
	xnewl[update[i]] = xnewl_loc[update_index[i]];
    }

    MPI_Reduce(xnewl, xnew,  n, MPI_DOUBLE_PRECISION, MPI_SUM, 
                      0, MPI_COMM_WORLD);

/*    get the runtime */

    t2 = MPI_Wtime();
    t21 = t2 - t1;
/*  printf("\n  %s","   "); */
    printf("\n  me= %d   TIME= %.4f   T1= %.4f   T2= %.4f", me, t21, t1, t2);

/*  Get the final L2 norm of the residual: */

    i2 = n - 1;
    for (i = 0; i <= i2; ++i) {
	resid1[i] = 0.;
    }
    sumr = 0.;
    i2 = n_update - 1;
    for (i = 0; i <= i2; ++i) {
	sum = val[i] * xnewl_loc[i] - rhsm_loc[i];
	i1 = bindx[i + 1] - 1;
	for (j = bindx[i]; j <= i1; ++j) {
	    sum += val[j] * xnewl_loc[bindx[j]];
	}
	d1 = sum;
	sumr += d1 * d1;
	resid1[update_inv[i]] = -sum;
    }
    MPI_Allreduce(&sumr, &sumrt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
                      MPI_COMM_WORLD);

/*  Get RESID( ) - the global residual vector: */

    MPI_Reduce(resid1, resid,  n, MPI_DOUBLE_PRECISION, MPI_SUM, 
                      0, MPI_COMM_WORLD);
    sumrt = sqrt(sumrt);


/*  Calculate the L2 norm of the error */

    if (me == 0) {

/*  Calculate the maximum residual: */

	ik = 0;
	rxmax = resid[0];
	i2 = n - 1;
	for (i = 1; i <= i2; ++i) {
	    if (resid[i] > rxmax) {
		ik = i;
		rxmax = resid[i];
	    }
	}

/*  Calculate the L2 norm of the error and the maximum error when */
/*  the analytic solution is known, and print out information */
/*  about the test case: */

        ipr= 1;
	switch (example_no) {
	    case 1:  goto L11;
	    case 2:  goto L12;
	    case 3:  goto L13;
	    case 4:  goto L14;
	    case 5:  goto L15;
	    case 6:  goto L16;
	    case 7:  goto L17;
	    case 8:  goto L18;
	    case 9:  goto L19;
	    case 10:  goto L20;
	    case 11:  goto L21;
	    case 12:  goto L22;
	    case 13:  goto L23;
	    case 14:  goto L24;
	    case 15:  goto L25;
	    case 16:  goto L26;
	    case 17:  goto L27;
	    case 18:  goto L28;
	    case 19:  goto L29;
	    case 20:  goto L30;
	    case 21:  goto L131;
	    case 22:  goto L132;
	    case 23:  goto L133;
	    case 24:  goto L134;
	    case 25:  goto L135;
	    case 26:  goto L136;
	    case 27:  goto L137;
	    case 28:  goto L138;
	}
L11:
	example1_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L12:
	example2_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L13:
	example3_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L14:
	example4_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L15:
	example5_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L16:
	example6_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L17:
	example7_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L18:
	example8_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L19:
	example9_err(xnew, err, &global_1.nx, &ny, &sumerr, &ipr);
	goto L70;
L20:
	example10_err(xnew, err, &global_1.nx, &sumerr, &ipr);
	goto L70;
L21:
	example11_err(xnew, err, &global_1.nx, &ny, &sumerr, &ipr);
	goto L70;
L22:
	example12_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L23:
	example13_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L24:
	example14_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr, &ipr);
	goto L70;
L25:
	example15_err(xnew, err, &global_1.nx, &ny, &sumerr, &ipr);
	goto L70;
L26:
	example16_err(xnew, err, &global_1.nx, &ny, &sumerr);
	goto L70;
L27:
	example17_err(xnew, err, &global_1.nx, &ny, &sumerr);
	goto L70;
L28:
	example18_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr);
	goto L70;
L29:
	example19_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr);
	goto L70;
L30:
	example20_err(xnew, err, &global_1.nx, &ny, &sumerr,&ipr);
	goto L70;
L131:
	example21_err(xnew, err, &global_1.nx, &ny, &sumerr);
	goto L70;
L132:
	example22_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr);
	goto L70;
L133:
	example23_err(xnew, err, &global_1.nx, &ny, &sumerr, &order);
	goto L70;
L134:
	example24_err(xnew, err, &global_1.nx, &ny, &sumerr);
	goto L70;
L135:
	example25_err(xnew, err, &global_1.nx, &ny, &nz, &sumerr);
	goto L70;
L136:
	example26_err(xnew, err, &global_1.nx, &ny, &sumerr1, &sumerr2, &
		kn1, &kn2, &kn3, &order, &ipr);
	goto L70;
L137:
	example27_err(xnew, err, &global_1.nx, &ny, &sumerr1, &sumerr2, &
		kn1, &kn2, &kn3, &order, &ipr);
	goto L70;
L138:
	example28_err(xnew, err, &global_1.nx, &ny, &sumerr, &kn1, &order, &ipr);

/* Note: the true (analytic) solution is unknown for   
         examples 16,17,18,19,21,22,23,24,25 */

L70:

/*  Calculate the L2 norm of the error and the maximum error */
/*  when the true solution is known */

	if (example_no < 16 || example_no == 20 || example_no == 26 || 
		example_no == 27 || example_no == 28) {
	    if (example_no == 26 || example_no == 27) {
		goto L605;
	    }
	    sum = 0.;
	    i2 = n - 1;
	    for (i = 0; i <= i2; ++i) {
		d1 = err[i];
		sum += d1 * d1;
	    }
	    err_l2 = sqrt(sum) / sumerr;
	    goto L610;
L605:
	    sum1 = 0.;
	    sum2 = 0.;
	    i2 = n - 1;
	    for (i = 0; i <= i2; i += 2) {
		d1 = err[i];
		sum1 += d1 * d1;
		d1 = err[i + 1];
		sum2 += d1 * d1;
	    }
	    d1 = sumerr1;
	    d2 = sumerr2;
	    err_l2 = sqrt(sum1 + sum2) / sqrt(d1 * d1 + d2 * d2);
	    printf("\n  ERR_L2_real = %.8f   ERR_L2_imag = %.8f", sqrt(sum1)/sumerr1, sqrt(sum2)/sumerr2);
L610:

/*  Calculate the maximum error */

	    jk = 0;
	    rmax = fabs(err[0]);
	    i2 = n - 1;
	    for (i = 1; i <= i2; ++i) {
		if ((d1 = err[i], fabs(d1)) > rmax) {
		    jk = i;
		    rmax = (d1 = err[i], fabs(d1));
		}
	    }
	}

	printf("\n  RESIDUAL= %.10g ", sumrt);
	if (example_no < 16 || example_no == 20 || example_no == 26 || 
		example_no == 27 || example_no == 28) {
	    printf("\n  ERR_L2= %.10g ", err_l2);
	}
	printf("\n  Res_max= %.10g   Res_max_Index= %d", rxmax, ik);
	if (example_no < 16 || example_no == 20 || example_no == 26 || 
		example_no == 27) {
            printf("\n  Err_max= %.10g   Err_max_Index= %d", rmax, jk);
	}
    }
/* endif  me=0 */
    
    MPI_Finalize(); 
    exit (0);
    return 0; 
} /* MAIN */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int bound_3d(int *nx, int *ny, int *nz, 
	int *nlocal, int *nxst, int *nxend, int *nyst, 
	int *nyend, int *nzst, int *nzend, int *nxproc, 
	int *nyproc, int *nzproc, int *me)
{
    static int n, i1, i2, i3, i4, nremx, nremy, nremz, nlocalx, nlocaly, 
	    nlocalz, nyzproc;

/*  Purpose: Define the starting and ending grid point indices handled   
    by the current processor, for a 3D problem. */

    n = *nx * *ny * *nz;
    nlocalx = *nx / *nxproc;
    nlocaly = *ny / *nyproc;
    nlocalz = *nz / *nzproc;
    nyzproc = *nyproc * *nzproc;
    i1 = *me / nyzproc;
    i2 = *me % nyzproc;
    i3 = i2 / *nzproc;
    i4 = i2 % *nzproc;
    *nxst = nlocalx * i1 + 1;
    *nxend = *nxst + nlocalx - 1;
    *nyst = nlocaly * i3 + 1;
    *nyend = *nyst + nlocaly - 1;
    *nzst = nlocalz * i4 + 1;
    *nzend = *nzst + nlocalz - 1;
    nremx = *nx % *nxproc;
    if (nremx != 0) {
	*nxst += min(nremx,i1);
	*nxend = *nxst + nlocalx - 1;
	if (i1 < nremx) {
	    ++(*nxend);
	}
    }
    nremy = *ny % *nyproc;
    if (nremy != 0) {
	*nyst += min(nremy,i3);
	*nyend = *nyst + nlocaly - 1;
	if (i3 < nremy) {
	    ++(*nyend);
	}
    }
    nremz = *nz % *nzproc;
    if (nremz != 0) {
	*nzst += min(nremz,i4);
	*nzend = *nzst + nlocalz - 1;
	if (i4 < nremz) {
	    ++(*nzend);
	}
    }
    *nlocal = (*nzend - *nzst + 1) * (*nyend - *nyst + 1) * (*nxend - *nxst + 
	    1);
    return 0;
} /* bound_3d */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int bound_2d(int *nx, int *ny, int *nlocal, 
	int *nxst, int *nxend, int *nyst, int *nyend, int 
	*nxproc, int *nyproc, int *me)
{
    static int n, i1, i2, nremx, nremy, nlocalx, nlocaly;

/*  Purpose: Define the starting and ending grid point indices handled   
    by the current processor, for a 2D problem. */

    n = *nx * *ny;
    nlocalx = *nx / *nxproc;
    nlocaly = *ny / *nyproc;
    i1 = *me / *nyproc;
    i2 = *me % *nyproc;
    *nxst = nlocalx * i1 + 1;
    *nxend = *nxst + nlocalx - 1;
    *nyst = nlocaly * i2 + 1;
    *nyend = *nyst + nlocaly - 1;
    nremx = *nx % *nxproc;
    if (nremx != 0) {
	*nxst += min(nremx,i1);
	*nxend = *nxst + nlocalx - 1;
	if (i1 < nremx) {
	    ++(*nxend);
	}
    }
    nremy = *ny % *nyproc;
    if (nremy != 0) {
	*nyst += min(nremy,i2);
	*nyend = *nyst + nlocaly - 1;
	if (i2 < nremy) {
	    ++(*nyend);
	}
    }
    *nlocal = (*nyend - *nyst + 1) * (*nxend - *nxst + 1);
    return 0;
} /* bound_2d */

/* ******************************************************************** */
/* ******************************************************************** */
/* Subroutine */ int bound_1d(int *nx, int *nlocal, int *nxst, 
	int *nxend, int *nxproc, int *me)
{
    static int n, i1, nremx, nlocalx;

/*  Purpose: Define the starting and ending grid point indices handled   
    by the current processor, for a 1D problem. */

    n = *nx;
    nlocalx = *nx / *nxproc;
    i1 = *me;
    *nxst = nlocalx * i1 + 1;
    *nxend = *nxst + nlocalx - 1;
    nremx = *nx % *nxproc;
    if (nremx != 0) {
	*nxst += min(nremx,i1);
	*nxend = *nxst + nlocalx - 1;
	if (i1 < nremx) {
	    ++(*nxend);
	}
    }
    *nlocal = *nxend - *nxst + 1;
    return 0;
} /* bound_1d */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example_3d(int *ist, int *iend, int *jst, 
	int *jend, int *kst, int *kend, int *nx, int *ny, 
	int *nz, int *index_s_xm, int *index_s_xp, int *
	index_s_ym, int *index_s_yp, int *index_s_zm, int *
	index_s_zp, int *index_r_xm, int *index_r_xp, int *
	index_r_ym, int *index_r_yp, int *index_r_zm, int *
	index_r_zp, int *ixp, int *ixm, int *iyp, int *iym, 
	int *izp, int *izm)
{
    /* System generated locals */
    int i1, i2, i3;

    /* Local variables */
    static int nxylocal, i, j, k, ii1, ii2, ii3, ii4, jcoef_ii_1, 
	    jcoef_ii_2, jcoef_ii_3, jcoef_ii_4, jcoef_ii_5, 
	    jcoef_ii_6, jcoef_ii_7, ii, igl, nxy, nxlocal, nylocal;

/*  Purpose:  Determine the indices of the external boundary grid points   
    comunicating with the neighboring  subdomains (3D, 7pt problems) */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --index_r_zp;
    --index_r_zm;
    --index_r_yp;
    --index_r_ym;
    --index_r_xp;
    --index_r_xm;
    --index_s_zp;
    --index_s_zm;
    --index_s_yp;
    --index_s_ym;
    --index_s_xp;
    --index_s_xm;

    /* Function Body */
    nxy = *nx * *ny;
    nxlocal = *iend - *ist + 1;
    nylocal = *jend - *jst + 1;
    nxylocal = nxlocal * nylocal;
    *ixm = 0;
    *ixp = 0;
    *iym = 0;
    *iyp = 0;
    *izm = 0;
    *izp = 0;
    ii = 0;
    i1 = *kend;
    for (k = *kst; k <= i1; ++k) {
	i2 = *jend;
	for (j = *jst; j <= i2; ++j) {
	    i3 = *iend;
	    for (i = *ist; i <= i3; ++i) {
		++ii;
		ii1 = (ii - 1) / nxylocal;
		ii2 = (ii - 1) % nxylocal;
		ii3 = ii2 / nxlocal;
		ii4 = ii2 % nxlocal;
		igl = (*kst - 1 + ii1) * nxy + *nx * (*jst - 1 + ii3) + *ist + 
			ii4;
		jcoef_ii_1 = igl;
		jcoef_ii_2 = igl + 1;
		jcoef_ii_3 = igl + *nx;
		jcoef_ii_4 = igl - 1;
		jcoef_ii_5 = igl - *nx;
		jcoef_ii_6 = igl + nxy;
		jcoef_ii_7 = igl - nxy;

		if (i == *ist && i != 1) {
		    ++(*ixm);
		    index_s_xm[*ixm] = jcoef_ii_4;
		    index_r_xm[*ixm] = jcoef_ii_1;
		}
		if (i == *iend && i != *nx) {
		    ++(*ixp);
		    index_s_xp[*ixp] = jcoef_ii_2;
		    index_r_xp[*ixp] = jcoef_ii_1;
		}
		if (j == *jst && j != 1) {
		    ++(*iym);
		    index_s_ym[*iym] = jcoef_ii_5;
		    index_r_ym[*iym] = jcoef_ii_1;
		}
		if (j == *jend && j != *ny) {
		    ++(*iyp);
		    index_s_yp[*iyp] = jcoef_ii_3;
		    index_r_yp[*iyp] = jcoef_ii_1;
		}
		if (k == *kst && k != 1) {
		    ++(*izm);
		    index_s_zm[*izm] = jcoef_ii_7;
		    index_r_zm[*izm] = jcoef_ii_1;
		}
		if (k == *kend && k != *nz) {
		    ++(*izp);
		    index_s_zp[*izp] = jcoef_ii_6;
		    index_r_zp[*izp] = jcoef_ii_1;
		}
	    }
	}
    }
    return 0;
} /* example_3d */

/* ********************************************************************* */
/* ******************************************************************** */

/* Subroutine */ int example_2d(int *ist, int *iend, int *jst, 
	int *jend, int *nx, int *ny, int *index_s_xm, 
	int *index_s_xp, int *index_s_ym, int *index_s_yp, 
	int *index_r_xm, int *index_r_xp, int *index_r_ym, 
	int *index_r_yp, int *ixp, int *ixm, int *iyp, 
	int *iym)
{
    /* System generated locals */
    int i1, i2;

    /* Local variables */
    static int nxylocal, i, j, ii1, ii2, ii3, ii4, jcoef_ii_1, 
	    jcoef_ii_2, jcoef_ii_3, jcoef_ii_4, jcoef_ii_5, ii, igl, 
	    nxy, nxlocal, nylocal;

/* Purpose: Determine the indices of the external boundary grid points   
   comunicating with the neighboring  subdomains (2D, 5pt problems) */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --index_r_yp;
    --index_r_ym;
    --index_r_xp;
    --index_r_xm;
    --index_s_yp;
    --index_s_ym;
    --index_s_xp;
    --index_s_xm;

    /* Function Body */
    nxy = *nx * *ny;
    nxlocal = *iend - *ist + 1;
    nylocal = *jend - *jst + 1;
    nxylocal = nxlocal * nylocal;
    *ixm = 0;
    *ixp = 0;
    *iym = 0;
    *iyp = 0;
    ii = 0;
    i1 = *jend;
    for (j = *jst; j <= i1; ++j) {
	i2 = *iend;
	for (i = *ist; i <= i2; ++i) {
	    ++ii;
	    ii1 = (ii - 1) / nxylocal;
	    ii2 = (ii - 1) % nxylocal;
	    ii3 = ii2 / nxlocal;
	    ii4 = ii2 % nxlocal;
	    igl = *nx * (*jst - 1 + ii3) + *ist + ii4;
	    jcoef_ii_1 = igl;
	    jcoef_ii_2 = igl + 1;
	    jcoef_ii_3 = igl + *nx;
	    jcoef_ii_4 = igl - 1;
	    jcoef_ii_5 = igl - *nx;

	    if (i == *ist && i != 1) {
		++(*ixm);
		index_s_xm[*ixm] = jcoef_ii_4;
		index_r_xm[*ixm] = jcoef_ii_1;
	    }
	    if (i == *iend && i != *nx) {
		++(*ixp);
		index_s_xp[*ixp] = jcoef_ii_2;
		index_r_xp[*ixp] = jcoef_ii_1;
	    }
	    if (j == *jst && j != 1) {
		++(*iym);
		index_s_ym[*iym] = jcoef_ii_5;
		index_r_ym[*iym] = jcoef_ii_1;
	    }
	    if (j == *jend && j != *ny) {
		++(*iyp);
		index_s_yp[*iyp] = jcoef_ii_3;
		index_r_yp[*iyp] = jcoef_ii_1;
	    }
	}
    }
    return 0;
} /* example_2d */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example_2d_all(int *ist, int *iend, int *
	jst, int *jend, int *nx, int *ny, int *index_s_xm, 
	int *index_s_xp, int *index_s_ym, int *index_s_yp, 
	int *index_r_xm, int *index_r_xp, int *index_r_ym, 
	int *index_r_yp, int *ixp, int *ixm, int *iyp, 
	int *iym, int *nx_st, int *nx_end, 
	int *ny_st, int *ny_end)
{
    /* System generated locals */
    int i1, i2;

    /* Local variables */
    static int nxylocal, i, j, ii1, ii2, ii3, ii4, jcoef_ii_1, 
	    jcoef_ii_2, jcoef_ii_3, jcoef_ii_4, jcoef_ii_5, ii, igl, 
	    nxy, nxlocal, nylocal;

/* Purpose: Determine the indices of the external boundary grid points   
   comunicating with the neighboring  subdomains (2D, 5pt problems) */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --index_r_yp;
    --index_r_ym;
    --index_r_xp;
    --index_r_xm;
    --index_s_yp;
    --index_s_ym;
    --index_s_xp;
    --index_s_xm;

    /* Function Body */
    nxy = *nx * *ny;
    nxlocal = *iend - *ist + 1;
    nylocal = *jend - *jst + 1;
    nxylocal = nxlocal * nylocal;
    *ixm = 0;
    *ixp = 0;
    *iym = 0;
    *iyp = 0;
    ii = 0;
/*     write(6,*)' Jst=',Jst,'  Jend=',Jend,'  Ist=',Ist,'  Iend=',Iend */
    i1 = *jend;
    for (j = *jst; j <= i1; ++j) {
	i2 = *iend;
	for (i = *ist; i <= i2; ++i) {
	    ++ii;
	    ii1 = (ii - 1) / nxylocal;
	    ii2 = (ii - 1) % nxylocal;
	    ii3 = ii2 / nxlocal;
	    ii4 = ii2 % nxlocal;
	    igl = *nx * (*jst - 1 + ii3) + *ist + ii4;
	    jcoef_ii_1 = igl;
	    jcoef_ii_2 = igl + 1;
	    jcoef_ii_3 = igl + *nx;
	    jcoef_ii_4 = igl - 1;
	    jcoef_ii_5 = igl - *nx;
/*           write(6,*)'J=',j,'  I=',i,'  IGL=',IGL */

	    if (i == *ist && i != 1 && (j > *ny_st && j < *ny_end)) {
		++(*ixm);
		index_s_xm[*ixm] = jcoef_ii_4;
		index_r_xm[*ixm] = jcoef_ii_1;
/*       write(6,*)' ixm=',ixm,'  Is_m =',JCOEF_II_4,'  Ir_m=',JCOEF_II_1 */
	    }
	    if (i == *iend && i != *nx && (j > *ny_st && j < *ny_end))
		     {
		++(*ixp);
		index_s_xp[*ixp] = jcoef_ii_2;
		index_r_xp[*ixp] = jcoef_ii_1;
/*       write(6,*)' ixp=',ixp,'  Is_p =',JCOEF_II_2,'  Ir_p=',JCOEF_II_1 */
	    }
	    if (j == *jst && j != 1 && (i > *nx_st && i < *nx_end)) {
		++(*iym);
		index_s_ym[*iym] = jcoef_ii_5;
		index_r_ym[*iym] = jcoef_ii_1;
/*       write(6,*)' iym=',iym,'  Is_m =',JCOEF_II_5,'  Ir_m=',JCOEF_II_1 */
	    }
	    if (j == *jend && j != *ny && (i > *nx_st && i < *nx_end))
		     {
		++(*iyp);
		index_s_yp[*iyp] = jcoef_ii_3;
		index_r_yp[*iyp] = jcoef_ii_1;
/*       write(6,*)' iyp=',iyp,'  Is_p =',JCOEF_II_3,'  Ir_p=',JCOEF_II_1 */
	    }
	}
    }
    return 0;
} /* example_2d_all */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example_3d_img(int *ist, int *iend, int *
	jst, int *jend, int *kst, int *kend, int *nx, int 
	*ny, int *nz, int *index_s_xm, int *index_s_xp, 
	int *index_s_ym, int *index_s_yp, int *index_s_zm, 
	int *index_s_zp, int *index_r_xm, int *index_r_xp, 
	int *index_r_ym, int *index_r_yp, int *index_r_zm, 
	int *index_r_zp, int *ixp, int *ixm, int *iyp, 
	int *iym, int *izp, int *izm)
{
    /* System generated locals */
    int i1, i2, i3;

    /* Local variables */
    static int nxylocal, i, j, k, ii1, ii2, ii3, ii4, jcoef_ii_1, 
	    jcoef_ii_2, jcoef_ii_3, jcoef_ii_4, jcoef_ii_5, 
	    jcoef_ii_6, jcoef_ii_7, ii, igl, nxy, nxy2, nxlocal, nylocal;

/*  Purpose:  Determine the indices of the external boundary grid points   
    comunicating with the neighboring subdomains (3D, 7pt problems) */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --index_r_zp;
    --index_r_zm;
    --index_r_yp;
    --index_r_ym;
    --index_r_xp;
    --index_r_xm;
    --index_s_zp;
    --index_s_zm;
    --index_s_yp;
    --index_s_ym;
    --index_s_xp;
    --index_s_xm;

    /* Function Body */
    nxy = (*nx + 2) * (*ny + 2);
    nxy2 = nxy << 1;
    nxlocal = *iend - *ist + 1;
    nylocal = *jend - *jst + 1;
    nxylocal = nxlocal * nylocal;
    *ixm = 0;
    *ixp = 0;
    *iym = 0;
    *iyp = 0;
    *izm = 0;
    *izp = 0;
    ii = 0;
    i1 = *kend;
    for (k = *kst; k <= i1; ++k) {
	i2 = *jend;
	for (j = *jst; j <= i2; ++j) {
	    i3 = *iend;
	    for (i = *ist; i <= i3; ++i) {
		++ii;
		ii1 = (ii - 1) / nxylocal;
		ii2 = (ii - 1) % nxylocal;
		ii3 = ii2 / nxlocal;
		ii4 = ii2 % nxlocal;
		igl = ((*kst - 1 + ii1) * nxy + (*nx + 2) * (*jst - 1 + ii3) + *
			ist + ii4 << 1) - 1;
		jcoef_ii_1 = igl;
		jcoef_ii_2 = igl + 2;
		jcoef_ii_3 = igl + (*nx + 2 << 1);
		jcoef_ii_4 = igl - 2;
		jcoef_ii_5 = igl - (*nx + 2 << 1);
		jcoef_ii_6 = igl + nxy2;
		jcoef_ii_7 = igl - nxy2;

		if (i == *ist && i > 1 && (j > 1 && j < *ny + 2) && (k > 
			1 && k < *nz + 2)) {
		    ++(*ixm);
		    index_s_xm[*ixm] = jcoef_ii_4;
		    index_r_xm[*ixm] = jcoef_ii_1;
		    ++(*ixm);
		    index_s_xm[*ixm] = jcoef_ii_4 + 1;
		    index_r_xm[*ixm] = jcoef_ii_1 + 1;
		}
		if (i == *iend && i < *nx + 2 && (j > 1 && j < *ny + 2) &&
			 (k > 1 && k < *nz + 2)) {
		    ++(*ixp);
		    index_s_xp[*ixp] = jcoef_ii_2;
		    index_r_xp[*ixp] = jcoef_ii_1;
		    ++(*ixp);
		    index_s_xp[*ixp] = jcoef_ii_2 + 1;
		    index_r_xp[*ixp] = jcoef_ii_1 + 1;
		}
		if (j == *jst && j > 1 && (i > 1 && i < *nx + 2) && (k > 
			1 && k < *nz + 2)) {
		    ++(*iym);
		    index_s_ym[*iym] = jcoef_ii_5;
		    index_r_ym[*iym] = jcoef_ii_1;
		    ++(*iym);
		    index_s_ym[*iym] = jcoef_ii_5 + 1;
		    index_r_ym[*iym] = jcoef_ii_1 + 1;
		}
		if (j == *jend && j < *ny + 2 && (i > 1 && i < *nx + 2) &&
			 (k > 1 && k < *nz + 2)) {
		    ++(*iyp);
		    index_s_yp[*iyp] = jcoef_ii_3;
		    index_r_yp[*iyp] = jcoef_ii_1;
		    ++(*iyp);
		    index_s_yp[*iyp] = jcoef_ii_3 + 1;
		    index_r_yp[*iyp] = jcoef_ii_1 + 1;
		}
		if (k == *kst && k > 1 && (i > 1 && i < *nx + 2) && (j > 
			1 && j < *ny + 2)) {
		    ++(*izm);
		    index_s_zm[*izm] = jcoef_ii_7;
		    index_r_zm[*izm] = jcoef_ii_1;
		    ++(*izm);
		    index_s_zm[*izm] = jcoef_ii_7 + 1;
		    index_r_zm[*izm] = jcoef_ii_1 + 1;
		}
		if (k == *kend && k < *nz && (i > 1 && i < *nx + 2) && (j 
			> 1 && j < *ny + 2)) {
		    ++(*izp);
		    index_s_zp[*izp] = jcoef_ii_6;
		    index_r_zp[*izp] = jcoef_ii_1;
		    ++(*izp);
		    index_s_zp[*izp] = jcoef_ii_6 + 1;
		    index_r_zp[*izp] = jcoef_ii_1 + 1;
		}
	    }
	}
    }
    return 0;
} /* example_3d_img */

/* ********************************************************************* */
/* ******************************************************************** */

/* Subroutine */ int example_2d_img(int *ist, int *iend, int *
	jst, int *jend, int *nx, int *ny, int *index_s_xm, 
	int *index_s_xp, int *index_s_ym, int *index_s_yp, 
	int *index_r_xm, int *index_r_xp, int *index_r_ym, 
	int *index_r_yp, int *ixp, int *ixm, int *iyp, 
	int *iym)
{
    /* System generated locals */
    int i1, i2;

    /* Local variables */
    static int nxylocal, i, j, ii2, ii3, ii4, jcoef_ii_1, jcoef_ii_2, 
	    jcoef_ii_3, jcoef_ii_4, jcoef_ii_5, ii, igl, nxlocal, 
	    nylocal;

/* Purpose: determine the indices of the external boundary grid points   
   comunicating with the neighboring subdomains (2D, 5pt problems)   
   for problems with an imaginary solution.   
   This routine follows Heroux paper using: (x1;y1, x2:y2,...,   ,xn;yn) */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --index_r_yp;
    --index_r_ym;
    --index_r_xp;
    --index_r_xm;
    --index_s_yp;
    --index_s_ym;
    --index_s_xp;
    --index_s_xm;

    /* Function Body */
    nxlocal = *iend - *ist + 1;
    nylocal = *jend - *jst + 1;
    nxylocal = nxlocal * nylocal;
    *ixm = 0;
    *ixp = 0;
    *iym = 0;
    *iyp = 0;
    ii = 0;
    i1 = *jend;
    for (j = *jst; j <= i1; ++j) {
	i2 = *iend;
	for (i = *ist; i <= i2; ++i) {
	    ++ii;
	    ii2 = ii - 1;
	    ii3 = ii2 / nxlocal;
	    ii4 = ii2 % nxlocal;
	    igl = ((*nx + 2) * (*jst - 1 + ii3) + *ist + ii4 << 1) - 1;
	    jcoef_ii_1 = igl;
	    jcoef_ii_2 = igl + 2;
	    jcoef_ii_3 = igl + (*nx + 2 << 1);
	    jcoef_ii_4 = igl - 2;
	    jcoef_ii_5 = igl - (*nx + 2 << 1);

	    if (i == *ist && i > 1 && (j > 1 && j < *ny + 2)) {
		++(*ixm);
		index_s_xm[*ixm] = jcoef_ii_4;
		index_r_xm[*ixm] = jcoef_ii_1;
		++(*ixm);
		index_s_xm[*ixm] = jcoef_ii_4 + 1;
		index_r_xm[*ixm] = jcoef_ii_1 + 1;
	    }
	    if (i == *iend && i < *nx + 2 && (j > 1 && j < *ny + 2)) {
		++(*ixp);
		index_s_xp[*ixp] = jcoef_ii_2;
		index_r_xp[*ixp] = jcoef_ii_1;
		++(*ixp);
		index_s_xp[*ixp] = jcoef_ii_2 + 1;
		index_r_xp[*ixp] = jcoef_ii_1 + 1;
	    }
	    if (j == *jst && j > 1 && (i > 1 && i < *nx + 2)) {
		++(*iym);
		index_s_ym[*iym] = jcoef_ii_5;
		index_r_ym[*iym] = jcoef_ii_1;
		++(*iym);
		index_s_ym[*iym] = jcoef_ii_5 + 1;
		index_r_ym[*iym] = jcoef_ii_1 + 1;
	    }
	    if (j == *jend && j < *ny + 2 && (i > 1 && i < *nx + 2)) {
		++(*iyp);
		index_s_yp[*iyp] = jcoef_ii_3;
		index_r_yp[*iyp] = jcoef_ii_1;
		++(*iyp);
		index_s_yp[*iyp] = jcoef_ii_3 + 1;
		index_r_yp[*iyp] = jcoef_ii_1 + 1;
	    }
	}
    }
    return 0;
} /* example_2d_img */

/* ********************************************************************* */
/* ******************************************************************** */

/* Subroutine */ int example_2d_img_all(int *ist, int *iend, 
	int *jst, int *jend, int *nx, int *ny, int *
	index_s_xm, int *index_s_xp, int *index_s_ym, int *
	index_s_yp, int *index_r_xm, int *index_r_xp, int *
	index_r_ym, int *index_r_yp, int *ixp, int *ixm, 
	int *iyp, int *iym, int *nx_st, int *
	nx_end, int *ny_st, int *ny_end)
{
/* System generated locals */
    int i1, i2;

/* Local variables */
    static int nxylocal, i, j, ii2, ii3, ii4, jcoef_ii_1, jcoef_ii_2, 
	    jcoef_ii_3, jcoef_ii_4, jcoef_ii_5, ii, igl, nxlocal, 
	    nylocal;

/* Purpose:  Determine the indices of the external boundary grid points   
   comunicating with the neighboring subdomains (2D, 5 or 9 point schemes)   
   for problems with an imaginary solution.   
   This routine follows Heroux paper using: (x1;y1, x2:y2,...,   ,xn;yn) */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --index_r_yp;
    --index_r_ym;
    --index_r_xp;
    --index_r_xm;
    --index_s_yp;
    --index_s_ym;
    --index_s_xp;
    --index_s_xm;

    /* Function Body */
    nxlocal = *iend - *ist + 1;
    nylocal = *jend - *jst + 1;
    nxylocal = nxlocal * nylocal;
    *ixm = 0;
    *ixp = 0;
    *iym = 0;
    *iyp = 0;
    ii = 0;
    i1 = *jend;
    for (j = *jst; j <= i1; ++j) {
	i2 = *iend;
	for (i = *ist; i <= i2; ++i) {
	    ++ii;
	    ii2 = ii - 1;
	    ii3 = ii2 / nxlocal;
	    ii4 = ii2 % nxlocal;
	    igl = ((*nx + 2) * (*jst - 1 + ii3) + *ist + ii4 << 1) - 1;
	    jcoef_ii_1 = igl;
	    jcoef_ii_2 = igl + 2;
	    jcoef_ii_3 = igl + (*nx + 2 << 1);
	    jcoef_ii_4 = igl - 2;
	    jcoef_ii_5 = igl - (*nx + 2 << 1);

	    if (i == *ist && i > 1 && (j > *ny_st && j < *ny_end)) {
		++(*ixm);
		index_s_xm[*ixm] = jcoef_ii_4;
		index_r_xm[*ixm] = jcoef_ii_1;
		++(*ixm);
		index_s_xm[*ixm] = jcoef_ii_4 + 1;
		index_r_xm[*ixm] = jcoef_ii_1 + 1;
	    }
	    if (i == *iend && i < *nx + 2 && (j > *ny_st && j < *
		    ny_end)) {
		++(*ixp);
		index_s_xp[*ixp] = jcoef_ii_2;
		index_r_xp[*ixp] = jcoef_ii_1;
		++(*ixp);
		index_s_xp[*ixp] = jcoef_ii_2 + 1;
		index_r_xp[*ixp] = jcoef_ii_1 + 1;
	    }
	    if (j == *jst && j > 1 && (i > *nx_st && i < *nx_end)) {
		++(*iym);
		index_s_ym[*iym] = jcoef_ii_5;
		index_r_ym[*iym] = jcoef_ii_1;
		++(*iym);
		index_s_ym[*iym] = jcoef_ii_5 + 1;
		index_r_ym[*iym] = jcoef_ii_1 + 1;
	    }
	    if (j == *jend && j < *ny + 2 && (i > *nx_st && i < *
		    nx_end)) {
		++(*iyp);
		index_s_yp[*iyp] = jcoef_ii_3;
		index_r_yp[*iyp] = jcoef_ii_1;
		++(*iyp);
		index_s_yp[*iyp] = jcoef_ii_3 + 1;
		index_r_yp[*iyp] = jcoef_ii_1 + 1;
	    }
	}
    }
    return 0;
} /* example_2d_img_all */

/* ********************************************************************* */
/* ******************************************************************** */

/* Subroutine */ int example_1d(int *ist, int *iend, int *nx, 
	int *index_s_xm, int *index_s_xp, int *index_r_xm, 
	int *index_r_xp, int *ixp, int *ixm)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    static int i, ii4, jcoef_ii_1, jcoef_ii_2, jcoef_ii_4, ii, igl,
	     nxlocal;

/*   Purpose: Determine the indices of the external boundary grid points   
            comunicating with the neighboring subdomains (1D, 3pt problems) */
/* -------------------------------------------------------------------------- */

    /* Parameter adjustments */
    --index_r_xp;
    --index_r_xm;
    --index_s_xp;
    --index_s_xm;

    /* Function Body */
    nxlocal = *iend - *ist + 1;
    *ixm = 0;
    *ixp = 0;
    ii = 0;
    i1 = *iend;
    for (i = *ist; i <= i1; ++i) {
	++ii;
	ii4 = (ii - 1) % nxlocal;
	igl = *ist + ii4;
	jcoef_ii_1 = igl;
	jcoef_ii_2 = igl + 1;
	jcoef_ii_4 = igl - 1;
	if (i == *ist && i != 1) {
	    ++(*ixm);
	    index_s_xm[*ixm] = jcoef_ii_4;
	    index_r_xm[*ixm] = jcoef_ii_1;
	}
	if (i == *iend && i != *nx) {
	    ++(*ixp);
	    index_s_xp[*ixp] = jcoef_ii_2;
	    index_r_xp[*ixp] = jcoef_ii_1;
	}
    }
    return 0;
} /* example_1d */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example1_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u, x, y, z;
    static int ii, nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 1:   

    Uxx + Uyy + Uzz +1000*Ux = F */
/* ------------------------------------------------ */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 1:");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  Uxx + Uyy + Uzz +1000*Ux = F ");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U = X*Y*Z*(1.-X)*(1.-Y)*(1.-Z)");
    printf("\n  %s"," -----------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:    

    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		x = (double) i * h;
		y = (double) j * h;
		z = (double) k * h;
		u = x * y * z * (1. - x) * (1. - y) * (1. - z);
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example1_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example2_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u, x, y, z;
    static int ii, nxy, nout;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 2:   
    Uxx + Uyy + Uzz +1000*dexp(X*Y*Z)*(Ux + Uy - Uz) = F */
/* ----------------------------------------------------------- */
    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 2:");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  Uxx + Uyy + Uzz +1000*dexp(X*Y*Z)*(Ux + Uy - Uz) = F ");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U = X + Y + Z");
    printf("\n  %s"," -------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:

    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		x = (double) i * h;
		y = (double) j * h;
		z = (double) k * h;
		u = x + y + z;
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example2_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example3_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u, x, y, z;
    static int ii;
    static double pi;
    static int nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 3:   
    Uxx + Uyy + Uzz +100*X*Ux - Y*Uy + Z*Uz + 100*(X+Y+Z)*U/(X*Y*Z) = F */
/* ----------------------------------------------------------------------- */
    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 3:");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  Uxx+Uyy+Uzz+100*X*Ux-Y*Uy+Z*Uz+100*(X+Y+Z)*U/(X*Y*Z) = F");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U = dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z)");
    printf("\n  %s"," --------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:
    
/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		x = (double) i * h;
		y = (double) j * h;
		z = (double) k * h;
		u = exp(x * y * z) * sin(pi * x) * sin(pi * y) * sin(pi * z);
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example3_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example4_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u, x, y, z;
    static int ii;
    static double pi;
    static int nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 4:   
    Uxx + Uyy + Uzz - 10**5*X**2*(Ux + Uy + Uz) = F */
/* ---------------------------------------------------- */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 4:");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  Uxx + Uyy + Uzz -10**5*X**2*(Ux + Uy + Uz) = F ");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U = dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z)");
    printf("\n  %s"," -----------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");   
L80:    

/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		x = (double) i * h;
		y = (double) j * h;
		z = (double) k * h;
		u = exp(x * y * z) * sin(pi * x) * sin(pi * y) * sin(pi * z);
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example4_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example5_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u, x, y, z;
    static int ii;
    static double pi;
    static int nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for EXAMPLE 5:   
    Uxx + Uyy + Uzz -1000*(1+X**2)*Ux + 100*(Uy + Uz) = F */
/* ------------------------------------------------------------ */
    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 5:");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  Uxx + Uyy + Uzz -1000*(1+X**2)*Ux + 100*(Uy + Uz) = F ");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U=dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z) ");
    printf("\n  %s"," ----------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:    

/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		x = (double) i * h;
		y = (double) j * h;
		z = (double) k * h;
		u = exp(x * y * z) * sin(pi * x) * sin(pi * y) * sin(pi * z);
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example5_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example6_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u, x, y, z;
    static int ii;
    static double pi;
    static int nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for EXAMPLE 6:   
    Uxx + Uyy + Uzz -1000*((1.-2*X)*Ux+(1-2*Y)*Uy+(1-2*Z)*Uz) = F */
/* ------------------------------------------------------------------ */
    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 6:");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  Uxx+Uyy+Uzz-1000*((1.-2*X)*Ux+(1-2*Y)*Uy+(1-2*Z)*Uz) = F ");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U = dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z)");
    printf("\n  %s"," --------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");   
L80:    
             
/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		x = (double) i * h;
		y = (double) j * h;
		z = (double) k * h;
		u = exp(x * y * z) * sin(pi * x) * sin(pi * y) * sin(pi * z);
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example6_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example7_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u, x, y, z;
    static int ii;
    static double pi;
    static int nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 7:   
    Uxx + Uyy + Uzz - 1000*X**2*Ux +1000*U = F   
    Where F is calculated for:   
    U=exp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z) */
/* ------------------------------------------------------ */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 7:");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  Uxx + Uyy + Uzz -1000*X**2*Ux + 1000*U = F ");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U = dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z)");
    printf("\n  %s"," -----------------------------------------------------");
    printf("\n  %s","   ");   
    printf("\n  %s","   ");
L80:    

/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		x = (double) i * h;
		y = (double) j * h;
		z = (double) k * h;
		u = exp(x * y * z) * sin(pi * x) * sin(pi * y) * sin(pi * z);
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example7_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example8_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u;
    static int ii, nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 8   
    (example F3D taken from Saad's book page 90):   

    -Uxx - Uyy - Uzz + 10*dexp(X*Y)*Ux +10*dexp(-X*Y)*Uy   
    +10*U*[Y*dexp(X*Y) - X*dexp(-X*Y)]= F   

    Where F is calculated from:   
         F=A*U ;    U=1. */
/* ------------------------------------------------------ */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 8:");
    printf("\n  %s","  (Saad's example F3D page 90).");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -Uxx - Uyy - Uzz + 10*dexp(X*Y)*Ux +");
    printf("\n  %s","  10*dexp(-X*Y)*Uy+10*U*[Y*dexp(X*Y)-X*dexp(-X*Y)] = F");
    printf("\n  %s","    ");
    printf("\n  %s","  Where F is calculated from: F=A*U ; for:");
    printf("\n  %s","  U=1. on the inner grid points;  U=0 on the boundaries.");
    printf("\n  %s","  --------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:    
    
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		u = 1.;
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example8_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example9_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2;
    double d1;

    /* Local variables */
    static double h;
    static int i, j;
    static double u, x, y;
    static int ii, nxy;
    static double sumerr1;
    double con; 

/*  Purpose:   
    Calculate the error for example 9   
    (taken from Kelley's book (1995) pages 26-27):   
               -D(a(x,y)Du) = F   

    Where:        a(x,y)= cos(x);   
    and F is calculated from:   
            U=10.*x*y*(1-x)*(1-y)*dexp(x**4.5) */
/* ------------------------------------------------ */
    
    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 9:");
    printf("\n  %s","  from Kelley's book.  "); 
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -D(a(x,y)DU) = F;  a(x,y)= cos(x) ");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U = 10*x*y*(1-x)*(1-y)*dexp(x**4.5)");
    printf("\n  %s"," -----------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:
    con=4.5;
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *ny;
    for (j = 1; j <= i1; ++j) {
	i2 = *nx;
	for (i = 1; i <= i2; ++i) {
	    x = (double) i * h;
	    y = (double) j * h;
	    u = x * 10. * y * (1. - x) * (1. - y) * exp(pow(x, con));
	    d1 = u;
	    sumerr1 += d1 * d1;
	    err[ii] = xnew[ii] - u;
/*           write(NOUT,100) II,U, XNEW(II),ERR(II) */
	    ++ii;
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example9_err */

/* ******************************************************************** */
/* ******************************************************************** */
/* Subroutine */ int example10_err(double *xnew, double *err, 
	int *nx, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double h;
    static int i;
    static double u, x;
    static int ii;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 10 (1D)   
    (taken from  Ph.D Dissertation, July 2000 page 89   
    example GK4.16):    

    The solution is:   
            U=-0.05*X**2 +0.05*X */
/* ---------------------------------------------------- */
    
    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 10:");
    printf("\n  %s","  1D problem from the Ph.D dissertation,");
    printf("\n  %s","  July 2000 example GK4.16, page 89.");
    printf("\n  %s","   "); 
    printf("\n  %s","  Domain: [0 1];");
    printf("\n  %s","  The analytic solution is:");
    printf("\n  %s","   U = -0.05*X**2 +0.05*X");
    printf("\n  %s","   ------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:    

    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nx;
    for (i = 1; i <= i1; ++i) {
	x = (double) i * h;
	d1 = x;
	u = d1 * d1 * -.05 + x * .05;
	d1 = u;
	sumerr1 += d1 * d1;
	err[ii] = xnew[ii] - u;
/*        write(NOUT,100) II,U, XNEW(II),ERR(II) */
	++ii;
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example10_err */

/* ******************************************************************** */
/* ******************************************************************** */
/* Subroutine */ int example11_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2;
    double d1;

    /* Local variables */
    static double h;
    static int i, j;
    static double u, x, y;
    static int ii, nxy;
    static double sumerr1;

/*  Purpose:  Calculate the error for   
    example 11 (problem F2DA taken from Saad's book pages 89-90):   
    where F is calculated from:   
         F=A*U ;   U=1.  on the inner grid points;  U=0 on the boundaries */
/* ---------------------------------------------------------------------- */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 11:");
    printf("\n  %s","  from Saad's book, example F2DA (pages 89-90).");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -Uxx -Uyy +D[(10*(x+y)*U]/Dx + D[(10*(x-y)*U]/Dy = F");
    printf("\n  %s","    ");
    printf("\n  %s","  Where F is calculated from: F=A*U ; for: ");
    printf("\n  %s","  U=1. on the inner grid points;  U=0 on the boundaries.");
    printf("\n  %s"," --------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *ny;
    for (j = 1; j <= i1; ++j) {
	i2 = *nx;
	for (i = 1; i <= i2; ++i) {
	    x = (double) i * h;
	    y = (double) j * h;
	    u = 1.;
	    d1 = u;
	    sumerr1 += d1 * d1;
	    err[ii] = xnew[ii] - u;
/*           write(NOUT,100) II,U, XNEW(II),ERR(II) */
	    ++ii;
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example11_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example12_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u, x, y, z;
    static int ii;
    static double pi;
    static int nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 12:   
      Uxx + Uyy + Uzz  = F */
/* -------------------------------------- */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 12:");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  Uxx + Uyy + Uzz = F ");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from the analytic solution:");
    printf("\n  %s","  U = dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z)");
    printf("\n  %s"," -----------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:    

/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		x = (double) i * h;
		y = (double) j * h;
		z = (double) k * h;
		u = exp(x * y * z) * sin(pi * x) * sin(pi * y) * sin(pi * 
			z);
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example12_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example13_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u;
    static int ii;
    static double sumerr1;

/*  Purpose:  Calculate the error for   
    example 13 (F2DA taken from Saad's book pages 89-90    
    extended to 3D), where F is calculated from:   
         F=A*U      U=1. */
/* ------------------------------------------------------ */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 13:");
    printf("\n  %s","  (Saad's example F2DA pages 89-90 extended to 3D).");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -Uxx -Uyy -Uzz +D[(10*(x+y)*U]/Dx + D[(10*(x-y)*U]/Dy = F ");
    printf("\n  %s","   ");
    printf("\n  %s","  Where F is calculated from: F=A*U for:; ");
    printf("\n  %s","  U = 1. on the inner grid points; U = 0. on the boundaries.");
    printf("\n  %s"," --------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:    

    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		u = 1.;
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example13_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example14_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2, i3;
    double d1;

    /* Local variables */
    static double h;
    static int i, j, k;
    static double u;
    static int ii, nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 14 = example 8  modified
    (1000 instead of 10)   
    (example F3D taken from Saad's book page 90):   

    -Uxx - Uyy - Uzz + 1000*dexp(X*Y)*Ux +1000*dexp(-X*Y)*Uy   
    +1000*U*[Y*dexp(X*Y) - X*dexp(-X*Y)]= F   

    where F is calculated from:   
         F=A*U ;    U=1. */
/* ----------------------------------------------------------- */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 14:");
    printf("\n  %s","  Similar to example8 but with 1000 instead of 10.");
    printf("\n  %s","  (Saad's book; example F3D page 90):");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -Uxx - Uyy - Uzz + 1000*dexp(X*Y)*Ux +1000*dexp(-X*Y)*Uy+");
    printf("\n  %s","  1000*U*[Y*dexp(X*Y) - X*dexp(-X*Y)] = F ");
    printf("\n  %s","    ");
    printf("\n  %s","  Where F is calculated from: F=A*U ; for:");
    printf("\n  %s","  U = 1. on the inner grid points; U = 0. on the boundaries");
    printf("\n  %s"," ---------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:    

    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *nz;
    for (k = 1; k <= i1; ++k) {
	i2 = *ny;
	for (j = 1; j <= i2; ++j) {
	    i3 = *nx;
	    for (i = 1; i <= i3; ++i) {
		u = 1.;
		d1 = u;
		sumerr1 += d1 * d1;
		err[ii] = xnew[ii] - u;
/*              write(NOUT,100) II,U, XNEW(II),ERR(II) */
		++ii;
	    }
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example14_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example15_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr, int *ipr)
{
    /* System generated locals */
    int i1, i2;
    double d1;

    /* Local variables */
    static double h;
    static int i, j;
    static double u;
    static int ii, nxy;
    static double sumerr1;

/*  Purpose:   
    Calculate the error for example 15.   
    This is problem F2DB from Saad's book "Iterative Methods   
    for Sparse Linear Systems" (page 90):   

    Domain:  [0 1]*[0 1];   

     -(a*Ux)x - (b*Uy)y + (10*(x+y)*U)x + (10*(x-y)*U)y =F   

    Where:   
      for:   0.25< x,y < 0.75   a=b=1000.   
      else:                     a=b=1.   

      F is calculated from:   
      F=A*U ;  U=1. on all inner grid points.   

      With Dirichlet boundary condition: U=0. on all boundaries.   

      The analytic solution is: U=1 */
/* ------------------------------------------------------------------ */

    if(*ipr == 0) goto L80;
    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 15:");
    printf("\n  %s","  Similar to example8 but with 1000 instead of 10.");
    printf("\n  %s","  (Saad's book; problem F2DB page 90):");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -(a*Ux)x - (b*Uy)y +(10*(x+y)*U)x + (10*(x-y)*U)y = F");
    printf("\n  %s","    ");
    printf("\n  %s","  Where: ");
    printf("\n  %s","      for: 0.25< x,y < 0.75  a=b=1000. ");
    printf("\n  %s","      else:                  a=b=1.");
    printf("\n  %s","    ");
    printf("\n  %s","  F is calculated from: F=A*U; for:  ");
    printf("\n  %s","  U = 1. on the inner grid points,  U = 0. on the boundaries. ");
    printf("\n  %s"," ---------------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:    

    nxy = *nx * *ny;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *ny;
    for (j = 1; j <= i1; ++j) {
	i2 = *nx;
	for (i = 1; i <= i2; ++i) {
	    u = 1.;
	    d1 = u;
	    sumerr1 += d1 * d1;
	    err[ii] = xnew[ii] - u;
/*           write(NOUT,100) II,U, XNEW(II),ERR(II) */
	    ++ii;
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example15_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int example16_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr)
{

/*  Purpose:   
    Calculate the error for example 16.   
    This is example #4 of Van Der Vorst (1992) BI-CGSTAB paper;   

    Domain:  [0 1]*[0 1];   

              -(a*Ux)x - (a*Uy)y  + b*Ux  =F   
    where:   
      a is discontinuous; see specification in the paper   
      b = 2*exp(2*(x**2+y**2))   

     and F = 0 evrywhere except   
         F = 100 in the square [0.45 0.55]*[0.45 0.55]   

    With Dirchlet boundary conditions:   
    U=0. on y=1. and U=1. on all other boundaries.   

    The analytic solution of this problem is unknown. */
/* -------------------------------------------------------------------------- */


    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 16:");
    printf("\n  %s","  - example #4 from Van Der Vorst BI-CGSTAB (1992) paper.");
    printf("\n  %s","    "); 
    printf("\n  %s","  Domain: [0 1]*[0 1];");
    printf("\n  %s","    "); 
    printf("\n  %s","  -(a*Ux)x - (a*Uy)y + 2*exp(2*(x**2+y**2))*Ux = F  ");
    printf("\n  %s","    ");
    printf("\n  %s","  a is discontinuous, see specification in Van Der Vorst (1992) paper.");
    printf("\n  %s","    ");
    printf("\n  %s","  F = 0. evrywhere except,  ");
    printf("\n  %s","  F = 100 in the square [0.45 0.55]*[0.45 0.55]");  
    printf("\n  %s","   ");
    printf("\n  %s","  With Dirichlet boundary condition: ");
    printf("\n  %s","  U=0 on y=1 and  U=1. on all other boundaries.");
    printf("\n  %s"," ----------------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");

    return 0;
} /* example16_err */

/* ******************************************************************** */
/* ******************************************************************* */
/* Subroutine */ int example17_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr)
{

/*  Purpose:   
    Calculate the error for example 17.   
    This is example #2  of Van Der Vorst (1992) paper.   

    Domain:  [0 1]*[0 1];   

        -(d*Ux)x - (d*Uy)y  = 1   

    Where:   
      for:  0.10<  x,y < 0.90   d=1000.   
      else:                     d=1.   

     With Dirichlet boundary condition: u=0 on y=0   
     and Neumann b. c. Du/Dn=0 on the other boundaries.   

    The analytic solution of this problem is unknown. */
/* ---------------------------------------------------------------------- */

    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 17:");
    printf("\n  %s","  - example #2 from Van Der Vorst BI-CGSTAB (1992) paper.");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -(d*Ux)x - (d*Uy)y =F  ");
    printf("\n  %s","    ");
    printf("\n  %s","  Where:");
    printf("\n  %s","     d=1000.  for:  0.10< x,y < 0.90  ");
    printf("\n  %s","     d=1.     else;                   ");
    printf("\n  %s","  F = 1. ");
    printf("\n  %s","   ");
    printf("\n  %s","  With Dirichlet boundary condition: U=0 on y=0.");
    printf("\n  %s","  and Neumann b.c. DU/Dn=0 on all other boundaries.");
    printf("\n  %s"," -------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");

    return 0;
} /* example17_err */

/* ********************************************************************* */
/* ******************************************************************** */
/* Subroutine */ int example18_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr)
{

/*  Purpose:   
    Calculate the error for example 18.   
    The problem is taken from  I. G. Graham and M. J. Hagger paper,   
    SIAM J. Sci. Comp., 1999, pp. 2041-2066, "Unstructured additive   
    Schwarz-CG method for elliptic problems with highly discontinuous   
    coefficients". Our version enables also the addition of convection   
    terms.   

    Domain:  [0 1]*[0 1]*[0 1];   

     -(a*Ux)x - (b*Uy)y - (c*Uz)z + d*Ux + e*Uy + f*Uz + g*U =F   

    Where:   
      d=e=f=prescribed ;  g=0;   

    and:   
    for:   0.3333< x,y,z < 0.6666  a=b=c (=1.E+04)   
    else:                          a=b=c (=1.E+00)   

    F=0.;   
    With Dirichlet boundary conditions:   
        U=1. on Y=0  and U=0 on all other boundaries.   

    The analytic solution of this problem is unknown. */
/* ------------------------------------------------------------------------ */

    printf("\n  %s","    ");
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 18:");
    printf("\n  %s","    ");
    printf("\n  %s","  from Graham&Hagger  SIAM J. Sci. Comp. 1999 paper.");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -(a*Ux)x - (b*Uy)y - (c*Uz)z + d*Ux + e*Uy + f*Uz + g*U =F");
    printf("\n  %s","    ");
    printf("\n  %s","  Where:");
    printf("\n  %s","     F = 0");
    printf("\n  %s","     g=1;  d=e=f=prescribed;");
    printf("\n  %s","   and");
    printf("\n  %s","     a=b=c =1.E+04  for:  0.3333< x,y,z < 0.6666");
    printf("\n  %s","     a=b=c =1.E+00  else;");
    printf("\n  %s","   ");
    printf("\n  %s","  With Dirichlet boundary condition:");
    printf("\n  %s","     U=1 on y=0 and U=0 on all other boundaries.");
    printf("\n  %s","  -------------------------------------------------------------");
    printf("\n  %s","   ");
    printf("\n  %s","   ");

    return 0;
} /* example18_err */

/* ********************************************************************* */
/* ******************************************************************** */
/* Subroutine */ int example19_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr)
{

/*  Purpose:   
    Calculate the error for example 19.   
    This is the two domains model problem of L. Gerardo-Giorda,   
    P. A. Tallec and F. Nataf 2004 paper on advection-diffusion with   
    discontinuous coefficients.   

    Domain:  [0 1]*[0 1]*[0 1];   

     -(v(x)*Ux)x-(v(x)*Uy)y-(v(x)*Uz)z+d*Ux+e*Uy+f*Uz+g*U =F   

    Where:   
      g=1.;   
      d=e=f=prescribed;   
    and:   
      v(x) is discontinuous:   
      for:  0. <  x <= 0.5  v(x)=v1  (1.E-01)   
      for:  0.5<  x < 1.0   v(x)=v2  (1.E-05)   

    F=0.   
    With Dirichlet boundary conditions:   
        U=1. on Y=0  and U=0 on all other boundaries.   

    The analytic solution of this problem is unknown. */
/* ----------------------------------------------------------------------- */

    printf("\n  %s","    "); 
    printf("\n  %s","    ");
    printf("\n  %s","  Results for example 19:");
    printf("\n  %s","  The model problem of L. Gerardo-Giorda, P. A. Tallec ");
    printf("\n  %s","  and F. Nataf, 2004 paper on convection-diffusion");
    printf("\n  %s","  with discontinuous coefficients.");
    printf("\n  %s","    ");
    printf("\n  %s","  Domain: [0 1]*[0 1]*[0 1];");
    printf("\n  %s","    ");
    printf("\n  %s","  -(a*Ux)x - (b*Uy)y - (c*Uz)z + d*Ux + e*Uy + fUz + g*U =F");
    printf("\n  %s","    ");
    printf("\n  %s","  Where:");
    printf("\n  %s","     F=0.;");
    printf("\n  %s","     g=1;  d=e=f=prescribed;");
    printf("\n  %s","  and:");
    printf("\n  %s","     a=b=c=v1 =1.E-01   for: 0.0 < x <=0.5");
    printf("\n  %s","     a=b=c=v2 =1.E-05   for: 0.5 < x < 1.0");
    printf("\n  %s","    ");
    printf("\n  %s","  With Dirichlet boundary condition:");
    printf("\n  %s","  U=1 on y=0 and U=0 on all other boundaries.");
    printf("\n  %s","  -------------------------------------------------------------");
    printf("\n  %s","    ");
    printf("\n  %s","    ");


    return 0;
} /* example19_err */

/* ******************************************************************** */
/* ******************************************************************** */
/* Subroutine */ int example20_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr,int *ipr)
{
    /* System generated locals */
    int i1, i2;
    double d1;

    /* Local variables */
    static double h;
    static int i, j;
    static double u, x, y;
    static int ii;
    static double pi, sumerr1;

/*  Purpose: calculate the error for the 2D Helmholtz problem.   
    Example20  taken from Erlangga, Vuik, Oosterlee, 2004 paper, in   
    applied Numerical Mathematics. Example #1.   

    Domain: [0 1]*[0 1]   

        Uxx + Uyy + K**2*U = (K**2 -5*PI**2)*sin(PI*X)*sin(2*PI*Y)   

    Dirichlet boundary conditions:  u=0 on the boundaries.   

    The analytic solution is:   
      u=sin(PI*X)*sin(2*PI*Y)   on the inner grid points, */
/* ----------------------------------------------------------------------- */
    if(*ipr == 0) goto L100;
    printf("\n  %s", "    ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Results for example 20: ");
    printf("\n  %s","   Similar to the example in Erlangga, Vuik, Oosterlee 2004 paper.");
    printf("\n  %s", "    ");
    printf("\n  %s","   Domain: [0 1]*[0 1]  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Uxx+Uyy+K**2*U=(K**2-5*PI**2)*sin(PI*X)*sin(2*PI*Y) ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Dirichlet boundary condition: U=0 on on all boundaries ");
    printf("\n  %s", "    ");
    printf("\n  %s","   The analytic solution is: U=sin(PI*X)*sin(2*PI*Y)");
    printf("\n  %s","----------------------------------------------------------" );
    printf("\n  %s","   ");
L100:
/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    h = 1. / (double) (*nx + 1);
    sumerr1 = 0.;
    ii = 0;
    i1 = *ny;
    for (j = 1; j <= i1; ++j) {
	i2 = *nx;
	for (i = 1; i <= i2; ++i) {
	    x = (double) i * h;
	    y = (double) j * h;
	    u = sin(pi * x) * sin(pi * 2. * y);
	    d1 = u;
	    sumerr1 += d1 * d1;
	    err[ii] = xnew[ii] - u;
/*           write(NOUT,100) II,U, XNEW(II),ERR(II) */
	    ++ii;
	}
    }
    *sumerr = sqrt(sumerr1);

    return 0;
} /* example20_err */

/* ********************************************************************* */
/* ******************************************************************** */
/* Subroutine */ int example21_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr)
{

/*  Purpose: calculate the error for a 2D Helmholtz problem.   
    Example21 is similar to example #2 from Erlangga,   
    Vuik, Oosterlee, 2004 paper, in applied Numerical Mathematics.   

    Domain: [0 1]*[0 1];   

         Uxx + Uyy + K**2*U = f   

    Where: f=0   

    With Dirichlet boundary condition on the bottom: u=0 on y=0, except at   
    the center point:   
            U=1   for x=0.5; y=0   
            U=0   for x#0.5; y=0   

    and Neumann b.c. Du/Dn -i*K*U=0 on all other boundaries. */
/* ----------------------------------------------------------------------- */

    printf("\n  %s", "    ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Results for example 21: ");
    printf("\n  %s","   Similar to example #2 in Erlangga, Vuik, Oosterlee 2004 paper.");
    printf("\n  %s", "    ");
    printf("\n  %s","   Domain: [0 1]*[0 1]  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Uxx+Uyy+K**2*U= F  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Where:  F=0.           ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Dirichlet boundary condition: U=0 on y=0, ");
    printf("\n  %s","   except:  U=1 at (x=0.5 , y=0).               ");
    printf("\n  %s","   Neumann b.c. DU/Dn -i*K*U=0 on all other boundaries");
    printf("\n  %s", "    ");
    printf("\n  %s","--------------------------------------------------------------------" );
    printf("\n  %s","   ");

    return 0;
} /* example21_err */

/* ********************************************************************* */
/* ******************************************************************** */
/* Subroutine */ int example22_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr)
{

/*  Purpose:   
    Calculate the error for example 22.   
    This is the two domains model problem of L. Gerardo-Giorda,   
    P. A. Tallec and F. Nataf 2004 paper on advection-diffusion with   
    discontinuous coefficients.   

    Domain:  [-0.5 0.5]*[-0.5 0.5]*[0 1];   

     -(v(x)*Ux)x-(v(x)*Uy)y-(v(x)*Uz)z+d*Ux+e*Uy+f*Uz+g*U =F   

    Where:   
      g=1.;   
      d=e=f=prescribed;   
    and:   
      v(x) is discontinuous:   
    The domain is subdevided into 8 subdomains 2*2*2   
      v1=1.E-1; v2=v4=v5=v7=1.E-6;   v3=1.E-2;  v6=1.E-3; v8=1.E-4   

    F=0.   
    With Dirichlet boundary conditions:   
        U=1. on Y=0  and U=0 on all other boundaries.   

    The analytic solution of this problem is unknown. */
/* ----------------------------------------------------------------------- */
    printf("\n  %s", "    ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Results for example 22:");
    printf("\n  %s","   This is the two domains model problem of L. Gerardo-Giorda,");
    printf("\n  %s","   P. A. Tallec and F. Nataf 2004 paper on advection-diffusion");
    printf("\n  %s","   with discontinuous coefficients. ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Domain: [-0.5 0.5]*[-0.5 0.5]*[0 1]");
    printf("\n  %s", "    ");
    printf("\n  %s", "  -(v(x)*Ux)x-(v(x)*Uy)y-(v(x)*Uz)z+d*Ux+e*Uy+f*Uz+g*U =F");
    printf("\n  %s", "  Where:");
    printf("\n  %s", "    F = 0.");
    printf("\n  %s", "    g=1.;");
    printf("\n  %s", "    d=e=f=prescribed;");
    printf("\n  %s", "  and:");
    printf("\n  %s", "    v(x) is discontinuous:");
    printf("\n  %s", "  The domain is subdevided into 8 subdomains 2*2*2");
    printf("\n  %s", "    v1=1.E-1; v2=v4=v5=v7=1.E-6;   v3=1.E-2;  v6=1.E-3; v8=1.E-4");
    printf("\n  %s", "    ");
    printf("\n  %s", "  With Dirichlet boundary condition:");
    printf("\n  %s", "      U=1. on Y=0  and U=0 on all other boundaries.");
    printf("\n  %s", "    ");
    printf("\n  %s", "  The analytic solution of this problem is unknown. ");
    printf("\n  %s", "  -------------------------------------------------------------------");
    printf("\n  %s", "    ");
    printf("\n  %s", "    ");

    return 0;
} /* example22_err */

/* ********************************************************************* */
/* ******************************************************************** */
/* Subroutine */ int example23_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr, int *order)
{

/*  Purpose: calculate the error for a 2D Helmholtz problem.   
    Example23 is similar to example #2 from Erlangga,   
    Vuik, Oosterlee, 2004 paper, in applied Numerical Mathematics.   
    The equations are solved using a fouth order difference scheme.   

    Domain: [0 1]*[0 1];   

            Uxx + Uyy + K**2*U = f   

     where: f=0   

     With Dirichlet boundary condition on the bottom: u=0 on y=0, except   
     at the center point:   
              U=1   for x=0.5; y=0   
              U=0   for x#0.5; y=0   

     and Neumann b.c. Du/Dn -i*K*U=0 on all other boundaries. */
/* ----------------------------------------------------------------------- */
    printf("\n  %s", "    ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Results for example 23:");
    printf("\n  %s","   Similar to example #2 in Erlangga, Vuik & Oosterlee 2004 paper.");
    printf("\n  %s", "    ");
    if (*order == 2) {
    printf("\n  %s","   Using a 2nd order finite difference scheme.     ");
    }
    if (*order == 4) {
    printf("\n  %s","   Using a 4th order finite difference scheme.     ");
    }
    if (*order == 6) {
    printf("\n  %s","   Using a 6th order finite difference scheme.     ");
    }
    printf("\n  %s", "    ");
    printf("\n  %s","   Domain: [-1/2 1/2]*[0 1]  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Uxx+Uyy+K**2*U= f  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Where:  f=0.           ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Dirichlet boundary condition: U=0 on y=0, except: ");
    printf("\n  %s","   U=1. at (x=0.5, y=0.)");
    printf("\n  %s","   Neumann b.c.:  DU/Dn -i*K*U=0 on all other boundaries.");
    printf("\n  %s","   " );
    printf("\n  %s","------------------------------------------------------------------" );
    printf("\n  %s","   ");

    return 0;
} /* example23_err */

/* ******************************************************************** */
/* ******************************************************************** */
/* Subroutine */ int example24_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr)
{

/*  Purpose: calculate the error for the 2D problem   
    Example24  taken from O. Schenk  2009 paper, in   
    SIAM J. Scientific Computing   
    Marmpusi problem:   

    Domain: [Xmin-hx  Xmax+hx]*[Ymin-hy Ymax+hy];   

            Uxx + Uyy + K**2*U = F   
            K=2*PI*freq/C   
    Neumann b.c. Du/Dn -i*K*U=gi on all  boundaries.   
    and:   
              f=1   for x=0.5, y=0   
              f=0    else */
/* ----------------------------------------------------------------------- */

    printf("\n  %s", "    "); 
    printf("\n  %s", "    ");
    printf("\n  %s","   Results for example 24:");
    printf("\n  %s","   Marmousi problem from O. Schenk 2009 paper. ");
    printf("\n  %s","    ");
    printf("\n  %s","   Uxx+Uyy+K**2*U= f");
    printf("\n  %s","    ");
    printf("\n  %s","   Where:  f=1/(hx*Hy)   for: x=0.5(Xmin+Xmax), ");
    printf("\n  %s","           f=0           else");
    printf("\n  %s","    ");
    printf("\n  %s","   Neumann b.c.: DU/Dn -i*K*U=0 on all boundaries");
    printf("\n  %s","    ");
    printf("\n  %s","  ----------------------------------------------------");
    printf("\n  %s", "    ");

    return 0;
} /* example24_err */

/* ************************************************************************ */
/* ******************************************************************** */
/* Subroutine */ int example25_err(double *xnew, double *err, 
	int *nx, int *ny, int *nz, double *sumerr)
{

/*  Purpose: calculate the error for the 3D problem   
    Example25  is similar to the 3D problem in  O. Schenk  2009 paper, in   
    SIAM J. Scientific Computing   

    Domain: [Xmin-hx  Xmax+hx]*[Ymin-hy Ymax+hy]*[Zmin-hz Zmax+hz];   

            Uxx + Uyy +Uzz + K**2*U = F   
            K=2*PI*freq/C   
    Neumann b.c. Du/Dn -i*K*U=gi on all  boundaries.   
    and:   
              f=1   for x=0.5*(Xmin+Xmax), y=0, z=0.5*(Zmin+Zmax)   
              f=0    else */
/* ----------------------------------------------------------------------- */
    printf("\n  %s", "    ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Results for example 25:");
    printf("\n  %s","   3D Helmholtz problem from O. Schenk 2009 paper.");
    printf("\n  %s","    ");
    printf("\n  %s","   Uxx+Uyy+Uzz+K**2*U= f");
    printf("\n  %s","    ");
    printf("\n  %s","   Where:   ");
    printf("\n  %s","   f = 1/(hx*Hy*Hz) for: x=0.5*(Xmin+Xmax), y=0, z=0.5*(Zmin+Zmax); ");
    printf("\n  %s","   f = 0. else      ");
    printf("\n  %s","  "); 
    printf("\n  %s","   Neumann b.c.:  DU/Dn -i*K*U=0 on all boundaries.");
    printf("\n  %s","   ------------------------------------------------------------------ ");
    printf("\n  %s","  "); 

    return 0;
} /* example25_err */

/* ************************************************************************ */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_7pt_example1(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1, d2; 

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, z;
    static int k1;
    static double pi;
    static int icn;
    static double sum, sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D EXAMPLE1 operator on an n x n xn cube.   

    Example 1:   
      Uxx + Uyy + Uzz +1000*Ux = F   

    Where F is calculated for:   
      U = X*Y*Z*(1.-X)*(1.-Y)*(1.-Z) */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    h = 1. / (double) (global_2.n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    sum1 = 0.;
    icn = 0;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % global_2.n != global_2.n - 1) {
	val[k] = -1. - h * 500.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % global_2.n != 0) {
	val[k] = h * 500. - 1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n;
    if (*row / global_2.n % global_2.n != global_2.n - 1) {
	val[k] = -1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n;
    if (*row / global_2.n % global_2.n != 0) {
	val[k] = -1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n * global_2.n;
    if (bindx[k] < global_2.n * global_2.n * global_2.n) {
	val[k] = -1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n * global_2.n;
    if (bindx[k] >= 0) {
	val[k] = -1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6.;
    f = (y * y * z + y * (z * z) + x * x * z + x * (
	    z * z) + x * x * y + x * (y * y) - y * z - x *
	     z - x * y) * 2. + (y * z + x * 2. * (y * y) * z + 
	    x * 2. * y * (z * z) + y * y * (z * z) - x * 
	    2. * y * z - y * y * z - y * (z * z) - x * 
	    2. * (y * y) * (z * z)) * 1.e3 - y * y * 
	    2. * (z * z) - x * x * 2. * (z * z) - 
	    x * x * 2. * (y * y);

    rhs[*row] = -f * h * h;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum1 += d1 * d1;
    sum = sqrt(sum1);
    rhs[*row] /= sum;
    val[*location] /= sum;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[*kmax - i] /= sum;
    }

    return 0;
} /* create_matrix_row_7pt_example1 */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int create_matrix_row_7pt_example2(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, z;
    static int k1, icn;
    static double sum, sum1;

/*   Purpose:   
     Add one row to an MSR matrix corresponding to a 7pt discrete   
     approximation to the 3D EXAMPLE2 operator on an n x n xn cube.   

     Example 2:   

     Uxx + Uyy + Uzz +1000*dexp(X*Y*Z)*(Ux + Uy - Uz) = F   

     Where F is calculated for:   
     U = X + Y + Z */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    h = 1. / (double) (global_2.n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;
    f = exp(x * y * z) * 1000.;
    rhs[*row] = -f * h * h;

    sum1 = 0.;
    icn = 0;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % global_2.n == global_2.n - 1) {
	rhs[*row] -= (-1. - exp(x * y * z) * 500. * h) * (y + 1. + z);
    }
    if (*row % global_2.n != global_2.n - 1) {
	val[k] = -1. - exp(x * y * z) * 500. * h;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % global_2.n == 0) {
	rhs[*row] -= (exp(x * y * z) * 500. * h - 1.) * (y + 0. + z);
    }
    if (*row % global_2.n != 0) {
	val[k] = exp(x * y * z) * 500. * h - 1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n;
    if (*row / global_2.n % global_2.n == global_2.n - 1) {
	rhs[*row] -= (-1. - exp(x * y * z) * 500. * h) * (x + 1. + z);
    }
    if (*row / global_2.n % global_2.n != global_2.n - 1) {
	val[k] = -1. - exp(x * y * z) * 500. * h;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n;
    if (*row / global_2.n % global_2.n == 0) {
	rhs[*row] -= (exp(x * y * z) * 500. * h - 1.) * (x + 0. + z);
    }
    if (*row / global_2.n % global_2.n != 0) {
	val[k] = exp(x * y * z) * 500. * h - 1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n * global_2.n;
    if (bindx[k] >= global_2.n * global_2.n * global_2.n) {
	rhs[*row] -= (exp(x * y * z) * 500. * h - 1.) * (x + y + 1.);
    }
    if (bindx[k] < global_2.n * global_2.n * global_2.n) {
	val[k] = exp(x * y * z) * 500. * h - 1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n * global_2.n;
    if (bindx[k] < 0) {
	rhs[*row] -= (-1. - exp(x * y * z) * 500. * h) * (x + y + 0.);
    }
    if (bindx[k] >= 0) {
	val[k] = -1. - exp(x * y * z) * 500. * h;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6.;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum1 += d1 * d1;
    sum = sqrt(sum1);
    rhs[*row] /= sum;
    val[*location] /= sum;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[*kmax - i] /= sum;
    }

    return 0;
} /* create_matrix_row_7pt_example2 */

/* ********************************************************************* */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example3(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1 ;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, z;
    static int k1;
    static double pi;
    static int icn;
    static double sum, sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D EXAMPLE3 operator on an n x n xn cube.   

    Example 3:   

    Uxx + Uyy + Uzz +100*X*Ux - Y*Uy + Z*Uz + 100*(X+Y+Z)*U/(X*Y*Z) = F   

    Where F is calculated for:   
      U=dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z) */
/* ----------------------------------------------------------------------- */

/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    h = 1. / (double) (global_2.n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    sum1 = 0.;
    icn = 0;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % global_2.n != global_2.n - 1) {
	val[k] = -1. - x * 50. * h;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % global_2.n != 0) {
	val[k] = x * 50. * h - 1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n;
    if (*row / global_2.n % global_2.n != global_2.n - 1) {
	val[k] = y * .5 * h - 1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n;
    if (*row / global_2.n % global_2.n != 0) {
	val[k] = -1. - y * .5 * h;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n * global_2.n;
    if (bindx[k] < global_2.n * global_2.n * global_2.n) {
	val[k] = -1. - z * .5 * h;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n * global_2.n;
    if (bindx[k] >= 0) {
	val[k] = z * .5 * h - 1.;
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6. - (x + y + z) * 100. * h * h / (x * y * z);

    f = exp(x * y * z) * sin(pi * y) * sin(pi * z) * (y * y * (z 
	    * z) * sin(pi * x) + pi * 2 * y * z * cos(pi * x) - pi * 
	    pi * sin(pi * x)) + exp(x * y * z) * sin(pi * x) * sin(pi * 
	    z) * (x * x * (z * z) * sin(pi * y) + pi * 2 * x * 
	    z * cos(pi * y) - pi * pi * sin(pi * y)) + exp(x * y * z) 
	    * sin(pi * x) * sin(pi * y) * (x * x * (y * y) * sin(
	    pi * z) + pi * 2 * x * y * cos(pi * z) - pi * pi * sin(pi 
	    * z)) + x * 100. * exp(x * y * z) * sin(pi * y) * sin(pi * 
	    z) * (y * z * sin(pi * x) + pi * cos(pi * x)) - y * exp(x * y 
	    * z) * sin(pi * x) * sin(pi * z) * (x * z * sin(pi * y) + 
	    pi * cos(pi * y)) + z * exp(x * y * z) * sin(pi * x) * sin(pi 
	    * y) * (x * y * sin(pi * z) + pi * cos(pi * z)) + (x + y + 
	    z) * 100. * exp(x * y * z) * sin(pi * x) * sin(pi * y) * sin(
	    pi * z) / (x * y * z);
    rhs[*row] = -f * h * h;
/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum1 += d1 * d1;
    sum = sqrt(sum1);
    rhs[*row] /= sum;
    val[*location] /= sum;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[*kmax - i] /= sum;
    }

    return 0;
} /* create_matrix_row_7pt_example3 */

/* ********************************************************************* */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example4(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1;
    static double pi;
    static int icn;
    static double sum, sum1;

/*   Purpose: */
/*   Add one row to an MSR matrix corresponding to a 7pt discrete */
/*   approximation to the 3D EXAMPLE4 operator on an n x n xn cube. */

/*    Example 4: */

/*    Uxx + Uyy + Uzz -10**5*X**2*(Ux + Uy + Uz) = F */

/*    where F is calculated for: */
/*    U=dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z) */
/* ----------------------------------------------------------------------- */

/*  Parameters: */
/*     row          == global row number of the new row to be added. */
/*     location     == local row where diagonal of the new row will be stored. */
/*     val,bindx    == (see user's guide). On output, val[] and bindx[] */
/*                     are appended such that the new row has been added. */

    n = global_1.nx;
/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    h = 1. / (double) (n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    f = exp(x * y * z) * sin(pi * y) * sin(pi * z) * (y * y * (z 
	    * z) * sin(pi * x) + pi * 2. * y * z * cos(pi * x) - pi * 
	    pi * sin(pi * x)) + exp(x * y * z) * sin(pi * x) * sin(pi * 
	    z) * (x * x * (z * z) * sin(pi * y) + pi * 2. * x * 
	    z * cos(pi * y) - pi * pi * sin(pi * y)) + exp(x * y * z) 
	    * sin(pi * x) * sin(pi * y) * (x * x * (y * y) * sin(
	    pi * z) + pi * 2. * x * y * cos(pi * z) - pi * pi * sin(
	    pi * z)) - x * x * 1.e5 * exp(x * y * z) * (sin(pi * y)
	     * sin(pi * z) * (y * z * sin(pi * x) + pi * cos(pi * x)) + 
	    sin(pi * x) * sin(pi * z) * (x * z * sin(pi * y) + pi * cos(
	    pi * y)) + sin(pi * x) * sin(pi * y) * (x * y * sin(pi * z) + 
	    pi * cos(pi * z)));
    rhs[*row] = -f * h * h;

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = x * 5.e4 * x * h - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -1. - x * 5.e4 * x * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = x * 5.e4 * x * h - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -1. - x * 5.e4 * x * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = x * 5.e4 * x * h - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -1. - x * 5.e4 * x * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6.;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example4 */


/* ******************************************************************** */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example5(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1 ;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1;
    static double pi;
    static int icn;
    static double sum, sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D EXAMPLE5 operator on an n x n xn cube.   

    Example 5:   
      Uxx + Uyy + Uzz -1000*(1+X**2)*Ux + 100*(Uy + Uz) = F   

    Where F is calculated for:   
      U=dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z) */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added.  */ 

    n = global_1.nx;
/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    h = 1. / (double) (n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    f = exp(x * y * z) * sin(pi * y) * sin(pi * z) * (y * y * (z 
	    * z) * sin(pi * x) + pi * 2. * y * z * cos(pi * x) - pi * 
	    pi * sin(pi * x) - (x * x + 1.) * 1.e3 * pi * cos(pi * x) 
	    + x * 100. * (y + z) * sin(pi * x)) + exp(x * y * z) * sin(pi 
	    * x) * sin(pi * z) * (x * x * (z * z) * sin(pi * y) 
	    + pi * 2. * x * z * cos(pi * y) - pi * pi * sin(pi * y) + 
	    pi * 100. * cos(pi * y)) + exp(x * y * z) * sin(pi * x) * sin(
	    pi * y) * (x * x * (y * y) * sin(pi * z) + pi * 2. *
	     x * y * cos(pi * z) - pi * pi * sin(pi * z) + pi * 
	    100. * cos(pi * z)) - (x * x + 1.) * 1.e3 * y * z * 
	    exp(x * y * z) * sin(pi * x) * sin(pi * y) * sin(pi * z);

    rhs[*row] = -f * h * h;

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = (x * x + 1.) * 500. * h - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -1. - (x * x + 1.) * 500. * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = -1. - h * 50.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = h * 50. - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = -1. - h * 50.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = h * 50. - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6.;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example5 */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int create_matrix_row_7pt_example6(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1;
    static double pi;
    static int icn;
    static double sum, sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D EXAMPLE6 operator on an n x n xn cube.   

    Example 6:   

    Uxx + Uyy + Uzz -1000*((1.-2*X)*Ux+(1-2*Y)*Uy+(1-2*Z)*Uz) = F   

    Where F is calculated for:   
              U=dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z) */
/* -------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    h = 1. / (double) (n + 1);
    i1 = n;
    k1 = *row / (i1 * i1);
    i1 = n;
    j = (*row - k1 * (i1 * i1)) / n;
    i1 = n;
    i = *row - k1 * (i1 * i1) - j * n;
    x = (i + 1) * h;
    y = (j + 1) * h;
    z = (k1 + 1) * h;
    f = exp(x * y * z) * sin(pi * y) * sin(pi * z) * (y * y * (z 
	    * z) * sin(pi * x) + pi * 2. * y * z * cos(pi * x) - pi * 
	    pi * sin(pi * x) - pi * 1.e3 * (1. - x * 2.) * cos(pi * x)) + 
	    exp(x * y * z) * sin(pi * x) * sin(pi * z) * (x * x * (
	    z * z) * sin(pi * y) + pi * 2. * x * z * cos(pi * y) - 
	    pi * pi * sin(pi * y) - pi * 1.e3 * (1. - y * 2.) * cos(pi * y)
	    ) + exp(x * y * z) * sin(pi * x) * sin(pi * y) * (x * x * 
	    (y * y) * sin(pi * z) + pi * 2. * x * y * cos(pi * z) - 
	    pi * pi * sin(pi * z) - pi * 1.e3 * (1. - z * 2.) * cos(pi 
	    * z)) - exp(x * y * z) * 1.e3 * sin(pi * x) * sin(pi * y) * 
	    sin(pi * z) * ((1. - x * 2.) * y * z + (1. - y * 2.) * x * 
	    z + (1. - z * 2.) * x * y);

    rhs[*row] = -f * h * h;

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = (1. - x * 2.) * 500. * h - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -1. - (1. - x * 2.) * 500. * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = (1. - y * 2.) * 500. * h - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -1. - (1. - y * 2.) * 500. * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = (1. - z * 2.) * 500. * h - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -1. - (1. - z * 2.) * 500. * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6.;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example6 */

/* ******************************************************************** */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example7(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1;
    static double pi;
    static int icn;
    static double sum, sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D EXAMPLE7 operator on an n x n xn cube.   

    Example 7:   

    Uxx + Uyy + Uzz -1000*X**2*Ux + 1000*U = F   

    Where F is calculated for:    
    U=dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z) */
/* -------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    h = 1. / (double) (n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    f = exp(x * y * z) * sin(pi * y) * sin(pi * z) * (y * y * (z 
	    * z) * sin(pi * x) + pi * 2. * y * z * cos(pi * x) - pi * 
	    pi * sin(pi * x)) + exp(x * y * z) * sin(pi * x) * sin(pi * 
	    z) * (x * x * (z * z) * sin(pi * y) + pi * 2. * x * 
	    z * cos(pi * y) - pi * pi * sin(pi * y)) + exp(x * y * z) 
	    * sin(pi * x) * sin(pi * y) * (x * x * (y * y) * sin(
	    pi * z) + pi * 2. * x * y * cos(pi * z) - pi * pi * sin(
	    pi * z)) - x * x * 1.e3 * exp(x * y * z) * (sin(pi * y)
	     * sin(pi * z) * (y * z * sin(pi * x) + pi * cos(pi * x))) + 
	    exp(x * y * z) * 1.e3 * sin(pi * x) * sin(pi * y) * sin(pi * z)
	    ;
    rhs[*row] = -f * h * h;

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = x * 500. * x * h - 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -1. - x * 500. * x * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6. - h * 1.e3 * h;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example7 */

/* ********************************************************************* */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example8(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1, icn;
    static double sum, sum1, sumf;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D EXAMPLE8 operator on an n x n xn cube.   

    Example 8 (problem F3D taken from Saad's book page 90):   

     -Uxx - Uyy - Uzz + 10*dexp(X*Y)*Ux +10*dexp(-X*Y)*Uy   
       +10*U*[Y*dexp(X*Y) - X*dexp(-X*Y)]= F   

    Where F is calculated from:   
    F=A*U ;  U=1. on the inner grid point, U=0 on the boundaries */

/* ----------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
    h = 1. / (double) (global_1.nx + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    icn = 0;
    sum = 0.;
    sumf = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = h * 5. * exp(x * y) - 1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -1. - h * 5. * exp(x * y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = h * 5. * exp(-x * y) - 1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -1. - h * 5. * exp(-x * y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = h * 10. * h * (y * exp(x * y) - x * exp(-x * y)) + 6.;
    sumf += val[*location];
    f = sumf;
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example8 */

/* ******************************************************************** */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_5pt_example9(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nxloc, int *nyloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, xm, ym, xp, yp;
    static int icn;
    static double uij, sum, sum1, uimj, uijm, uipj, uijp;
    double con; 

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5pt discrete   
    approximation to the 2D  EXAMPLE9 operator on an n x n  square.   

    Example 9 (taken from Kelley's book (1995) pages 26-27):   
                  -D(a(x,y)Du) = F   

    Where:        a(x,y)= cos(x);   

    and F is calculated for:   
                    u=10.*x*y*(1-x)*(1-y)*dexp(x**4.5) */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added.  */ 
    con=4.5;
    h = 1. / (double) (global_2.n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    xm = x - h;
    xp = x + h;
    yp = y + h;
    ym = y - h;
    uimj = 0.;
    uij = x * 10. * y * (1. - x) * (1. - y) * exp(pow(x, con));
    uipj = xp * 10. * y * (1. - xp) * (1. - y) * exp(pow(xp, con));
    if (*nxst != 1 || *nxst == 1 && i != 0) {
	uimj = xm * 10. * y * (1. - xm) * (1. - y) * exp(pow(xm, con));
    }
    uijp = x * 10. * yp * (1. - x) * (1. - yp) * exp(pow(x, con));
    uijm = x * 10. * ym * (1. - x) * (1. - ym) * exp(pow(x, con));
    f = uij * (cos(x) * 6. + cos(xp) + cos(xm)) - uipj * (cos(xp) + cos(x)) - 
	    uimj * (cos(xm) + cos(x)) - uijp * 2. * cos(x) - uijm * 2. * cos(
	    x);
    rhs[*row] = f;

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % global_2.n != global_2.n - 1) {
	val[k] = -cos(xp) - cos(x);
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % global_2.n != 0) {
	val[k] = -cos(xm) - cos(x);
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n;
    if (*row / global_2.n % global_2.n != global_2.n - 1) {
	val[k] = cos(x) * -2.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n;
    if (*row / global_2.n % global_2.n != 0) {
	val[k] = cos(x) * -2.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = cos(x) * 6. + cos(xp) + cos(xm);

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example9 */

/* ********************************************************************* */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_5pt_example10(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nxloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double h;
    static int i, k;
    static double x, ui, xm, xp;
    static int icn;
    static double uim, uip, xmm, sum, xpp, sum1, uimm, uipp;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to   
    the 5pt (1D) EXAMPLE10 operator on an n interval.   

    Example 10 (taken from the Ph.D Dissertation, July 2000, page 89-90   
    example GK4.16):   

    U=-0.05*X**2 +0.05*X */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    h = 1. / (double) (global_2.n + 1);
    i = *row;
    x = (double) (i + 1) * h;
    xm = x - h;
    xmm = xm - h;
    xp = x + h;
    xpp = xp + h;
    uim = 0.;
    uip = 0.;
    uimm = 0.;
    uipp = 0.;
    d1 = x;
    ui = d1 * d1 * -.05 + x * .05;
    if (*row < global_2.n - 1) {
	d1 = xp;
	uip = d1 * d1 * -.05 + xp * .05;
    }
    if (*row < global_2.n - 2) {
	d1 = xpp;
	uipp = d1 * d1 * -.05 + xpp * .05;
    }
    if (*row > 0) {
	d1 = xm;
	uim = d1 * d1 * -.05 + xm * .05;
    }
    if (*row > 1) {
	d1 = xmm;
	uimm = d1 * d1 * -.05 + xmm * .05;
    }
    if (*row == 0 || *row == global_2.n - 1) {
	rhs[*row] = ui * 5. - uip * 4. - uim * 4. + uipp + uimm;
    } else {
	rhs[*row] = ui * 6. - uip * 4. - uim * 4. + uipp + uimm;
    }

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row != global_2.n - 1) {
	val[k] = -4.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row != 0) {
	val[k] = -4.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + 2;
    if (*row < global_2.n - 2) {
	val[k] = 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 2;
    if (*row > 1) {
	val[k] = 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    if (*row == 0 || *row == global_2.n - 1) {
	val[*location] = 5.;
    } else {
	val[*location] = 6.;
    }

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example10 */

/* ************************************************************************ */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_5pt_example11(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nxloc, int *nyloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y;
    static int icn;
    static double sum, sum1, sumf;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5pt discrete   
    approximation to the 2D  EXAMPLE11 operator on an n x n  square.   

    Example 11 (taken from Saad's book  pages 89-90):   

         -Uxx -Uyy +D[(10*(x+y)*U]/Dx + D[(10*(x-y)*U]/Dy = F   

    Where  F is calculated from: F= A*U   
    and:   
      u=1.  on the inner grid points;  u=0 on the boundaries. */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    h = 1. / (double) (global_2.n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;

    icn = 0;
    sum = 0.;
    sumf = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % global_2.n != global_2.n - 1) {
	val[k] = h * 5. * (x + y) - 1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % global_2.n != 0) {
	val[k] = -1. - h * 5. * (x + y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n;
    if (*row / global_2.n % global_2.n != global_2.n - 1) {
	val[k] = h * 5. * (x - y) - 1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n;
    if (*row / global_2.n % global_2.n != 0) {
	val[k] = -1. - h * 5. * (x - y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 4.;
    sumf += val[*location];
    f = sumf;
    rhs[*row] = f;
/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example11 */

/* ************************************************************************* */
/* ************************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example12(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1;
    static double pi;
    static int icn;
    static double sum, sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D Poisson equation  on an n x n xn cube.   

    Example 12:   

    Uxx + Uyy + Uzz  = F   

    Where F is calculated for:   
         U=dexp(X*Y*Z)*dsin(PI*X)*dsin(PI*Y)*dsin(PI*Z) */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

/*  pi = 3.14159265359793; */
    pi = atan(1.) * 4.;
    n = global_1.nx;
    h = 1. / (double) (n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    f = exp(x * y * z) * sin(pi * y) * sin(pi * z) * (y * y * (z * z)
	    * sin(pi * x) + pi * 2. * y * z * cos(pi * x) - pi * pi 
	    * sin(pi * x)) + exp(x * y * z) * sin(pi * x) * sin(pi * z) 
	    * (x * x * (z * z) * sin(pi * y) + pi * 2. * x * 
	    z * cos(pi * y) - pi * pi * sin(pi * y)) + exp(x * y * z) 
	    * sin(pi * x) * sin(pi * y) * (x * x * (y * y) * sin(
	    pi * z) + pi * 2. * x * y * cos(pi * z) - pi * pi * sin(pi * z));

    rhs[*row] = -f * h * h;

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6.;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example12 */

/* ************************************************************************ */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_7pt_example13(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, z;
    static int k1, icn;
    static double sum, sum1, sumf;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D  EXAMPLE13 operator on an n x n x n  cube.   

    Example 13: (from Saad's book  pages 89-90 extended to 3D):   

        -Uxx -Uyy -Uzz +D[(10*(x+y)*U]/Dx + D[(10*(x-y)*U]/Dy = F   

    Where F is calculated from: F= A*U;   
    and:    
     u=1.  on the inner grid points;  u=0 on the boundaries. */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    h = 1. / (double) (global_2.n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    icn = 0;
    sum = 0.;
    sumf = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % global_2.n != global_2.n - 1) {
	val[k] = h * 5. * (x + y) - 1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % global_2.n != 0) {
	val[k] = -1. - h * 5. * (x + y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n;
    if (*row / global_2.n % global_2.n != global_2.n - 1) {
	val[k] = h * 5. * (x - y) - 1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n;
    if (*row / global_2.n % global_2.n != 0) {
	val[k] = -1. - h * 5. * (x - y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n * global_2.n;
    if (bindx[k] < global_2.n * global_2.n * global_2.n) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n * global_2.n;
    if (bindx[k] >= 0) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = 6.;
    sumf += val[*location];
    f = sumf;
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example13 */

/* ******************************************************************** */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example14(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1, icn;
    static double sum, sum1, sumf;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D EXAMPLE8 operator on an n x n xn cube.   

    Example 14 = example 8 modified; 1000 instead of 10   
    Example 8 (problem F3D taken from Saad's book page 90):   

     -Uxx - Uyy - Uzz + 1000*dexp(X*Y)*Ux +1000*dexp(-X*Y)*Uy   
       +1000*U*[Y*dexp(X*Y) - X*dexp(-X*Y)]= F   

    Where F is calculated from:   
    F=A*U ;  U=1. on the inner grid point, U=0 on the boundaries */
/* ----------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
    h = 1. / (double) (global_1.nx + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    z = (*nzst + k1) * h;

    icn = 0;
    sum = 0.;
    sumf = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = h * 500. * exp(x * y) - 1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -1. - h * 500. * exp(x * y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = h * 500. * exp(-x * y) - 1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -1. - h * 500. * exp(-x * y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = h * 1.e3 * h * (y * exp(x * y) - x * exp(-x * y)) + 6.;
    sumf += val[*location];
    f = sumf;
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example14 */

/* ********************************************************************* */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_5pt_example15(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nxloc, int *nyloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, v1, v2, v3, v4;
    static float xm, ym, xp, yp, acm, bcm, acp, bcp;
    static int icn;
    static double sum, sum1, sumf;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5pt discrete   
    approximation to the 2D EXAMPLE_15 operator on an n x n  square.   

    Example 15=problem F2DB from Saad's book "Iterative Methods   
      for Sparse Linear Systems" (page 90).   

    Domain: [0 1]*[0 1];   

     -(a*Ux)x - (b*Uy)y + (10*(x+y)*U)x +(10*(x-y)*U)y =R   

    Where:   
      for:  0.25< x,y < 0.75   a=b=1000.   
      else:                    a=b=1.   

    F is calculated from:   
    F=A*U ;  U=1. on all inner grid points.   

    With Dirichlet boundary condition: U=0 on all boundaries. */
/* -------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
    h = 1. / (double) (global_1.nx + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    xp = x + h * .5;
    xm = x - h * .5;
    yp = y + h * .5;
    ym = y - h * .5;
    v1 = 1.;
    v2 = 1.e3;
/*  v3= 1. */
/*  v4= 1000 */
    v3 = 1.e-4;
    v4 = 1.;
/*  v3= 10. */
/*  v3= 0.1 */
/*  v4= 10000. */
    acp = v1;
    acm = v1;
    bcp = v3;
    bcm = v3;
    if (xp > .25f && xp < .75f && (y > .25f && y < .75f)) {
	acp = v2;
    }
    if (xm > .25f && xm < .75f && (y > .25f && y < .75f)) {
	acm = v2;
    }
    if (yp > .25f && yp < .75f && (x > .25f && x < .75f)) {
	bcp = v4;
    }
    if (ym > .25f && ym < .75f && (x > .25f && x < .75f)) {
	bcm = v4;
    }

    icn = 0;
    sum = 0.;
    sumf = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = -acp + h * 5. * (x + y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -acm - h * 5. * (x + y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = -bcp + h * 5. * (x - y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -bcm - h * 5. * (x - y);
	d1 = val[k];
	sum += d1 * d1;
	sumf += val[k];
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = acp + acm + bcp + bcm;
    sumf += val[*location];
    f = sumf;
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example15 */

/* ********************************************************************* */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_5pt_example16(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nxloc, int *nyloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y;
    static double xm, ym, xp, yp, acm, bcm, acp, bcp;
    static int icn;
    static double sum, sum1, cons;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5pt discrete   
    approximation to the 2D EXAMPLE_16 operator on an nxn  square.   

    Example 16 = example #4 of Van Der Vorst (1992) Bi-CGSTAB paper.   

    Domain: [0 1]*[0 1];   

              -(a*Ux)x - (a*Uy)y  + b*Ux  =F   
    where:   
      a is discontinuous; see specification in the paper.   
      b = 2*exp(2*(x**2+y**2))   

    and F = 0 evrywhere except   
        F = 100 in the square [0.45 0.55]*[0.45 0.55]   

    With Dirichlet boundary conditions: u=0. on y=1. and u=1. on y=0 on   
    all other boundaries. */
/* -------------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
    h = 1. / (double) (global_1.nx + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;
    f = 0.;
    if (x > .45f && x < .55f && (y > .45f && y < .55f)) {
	d1 = h;
	f = d1 * d1 * 100.;
    }
    xp = x + h * .5;
    xm = x - h * .5;
    yp = y + h * .5;
    ym = y - h * .5;
    acp = 100.;
    acm = 100.;
    bcp = 100.;
    bcm = 100.;
/*    inner square */
    if (xp > .3f && xp < .7f && (y > .3f && y < .7f)) 
	acp = 1.e4;
    if (xm > .3f && xm < .7f && (y > .3f && y < .7f)) 
	acm = 1.e4;
    if (yp > .3f && yp < .7f && (x > .3f && x < .7f)) 
	bcp = 1.e4;
    if (ym > .3f && ym < .7f && (x > .3f && x < .7f)) 
	bcm = 1.e4;
/*  (inner) bounding strip */
/*  lower side */
    if (xp > .2f && xp < .8f && (y > .2f && y < .3f)) acp = 1.e-5;
    if (xm > .2f && xm < .8f && (y > .2f && y < .3f)) acm = 1.e-5;
    if (yp > .2f && yp < .3f && (x > .2f && x < .8f)) bcp = 1.e-5;
    if (ym > .2f && ym < .3f && (x > .2f && x < .8f)) bcm = 1.e-5;
/*  upper side (built of 2 strips) */
    if ((xp > .2f && xp < .49f || xp > .51f && xp < .8f) && (y > .7f && y < 
	    .8f))  acp = 1.e-5;
    if ((xm > .2f && xm < .49f || xm > .51f && xm < .8f) && (y > .7f && y < 
	    .8f))  acm = 1.e-5;
    if (yp > .7f && yp < .8f && (x > .2f && x < .49f || x > .51f && x < .8f)) 
	     bcp = 1.e-5;
    if (ym > .7f && ym < .8f && (x > .2f && x < .49f || x > .51f && x < .8f)) 
	     bcm = 1.e-5;
/*  2 side strips */
    if ((xp > .2f && xp < .3f || xp > .7f && xp < .8f) && (y > .3f && y < .7f))
        acp = 1.e-5;
    if ((xm > .2f && xm < .3f || xm > .7f && xm < .8f) && (y > .7f && y < .8f))
        acm = 1.e-5;
    if (yp > .3f && yp < .7f && (x > .2f && x < .3f || x > .7f && x < .8f)) 
	bcp = 1.e-5;
    if (ym > .3f && ym < .7f && (x > .2f && x < .3f || x > .7f && x < .8f)) 
	bcm = 1.e-5;

    icn = 0;
    sum = 0.;
/*  sumf = 0.;*/
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = -acp + h * exp((x * x + y * y) * 2.);
	d1 = val[k];
	sum += d1 * d1;
/*      sumf = sumf + val(k);*/
	++icn;
	++k;
    }
    if (*row % n == n - 1) {
	cons = -acp + h * exp((x * x + y * y) * 2.);
	f -= cons ;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -acm - h * exp((x * x + y * y) * 2.);
	d1 = val[k];
	sum += d1 * d1;
/*      sumf = sumf + val(k);*/
	++icn;
	++k;
    }
    if (*row % n == 0) {
	cons = -acm - h * exp((x * x + y * y) * 2.);
	f -= cons ;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = -bcp;
	d1 = val[k];
	sum += d1 * d1;
/*      sumf = sumf + val(k);*/
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -bcm;
	d1 = val[k];
	sum += d1 * d1;
/*      sumf = sumf + val(k);*/
	++icn;
	++k;
    }
    if (*row / n % n == 0) {
	cons = -bcm;
	f -= cons ;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = acp + acm + bcp + bcm;
/*  sumf = sumf +  val(location) */
/*  F = sumf */
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example16 */

/* ********************************************************************* */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_5pt_example17(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nxloc, int *nyloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y;
    static double v1, v2, dc, ec, xm, ym, xp, yp, acm, bcm, acp, bcp;
    static int icn;
    static double sum, sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5pt discrete   
    approximation to the 2D EXAMPLE_17 operator on an (n+2)x(n+2) square.   

    Example 17 = example #2 of Van Der Vorst (1992) Bi_CGSTAB paper.   

    Domain: [0 1]*[0 1];   

          -(d*Ux)x - (d*Uy)y  = 1   
    where:   
      for:   0.10 =< x,y <= 0.90   d=1000.   
      else:                        d=1.   

    With Dirichlet Boundary condition: u=0 on y=0 and   
    Neumann b.c. Du/Dn=0 on all other boundaries.  */ 
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
    h = 1. / (double) (global_1.nx + 1);
/*  F = H**2*(1.) */
    f = 1.;
/*  I = mod(location,NXloc) */
/*  J = mod(location/NXloc,NYloc) */
/*  X = (NXst + I)*H */
/*  Y = (NYst + J)*H */
    i = *row % (n + 2);
    j = *row / (n + 2) % (n + 2);
    x = i * h;
    y = j * h;

/*  Set eq. for the top boundary grid points (Neumann boundary condition). */

    if (*row > (n + 2) * (n + 1) - 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - (n + 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	f = 0.;
	goto L100;
    }

/*  Set eq. for the bottom boundary grid points (Dirichlet boundary condition). */

    if (*row <= n + 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	f = 0.;
	goto L100;
    }
/*  Set eq. for the left hand-side boundary grid points (Neumann boundary condition) */
    if (i == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + 1;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	f = 0.;
	goto L100;
    }

/*  Set eq. for the right hand-side boundary grid points (Neumann 
    boundary condition). */

    if (i == n + 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - 1;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	f = 0.;
	goto L100;
    }

/*  Set eq. for the inner grid points */

    xp = x + h * .5;
    xm = x - h * .5;
    yp = y + h * .5;
    ym = y - h * .5;

    v1 = 1.;
/*  v1=1000. */
/*  v2=1000. */
    v2 = 1.e4;
/*  v2=1. */
    d1 = h;
    acp = v1 / (d1 * d1);
    d1 = h;
    acm = v1 / (d1 * d1);
    d1 = h;
    bcp = v1 / (d1 * d1);
    d1 = h;
    bcm = v1 / (d1 * d1);
    if (xp >= .1f && xp <= .9f && (y >= .1f && y <= .9f)) {
	d1 = h;
	acp = v2 / (d1 * d1);
    }
    if (xm >= .1f && xm <= .9f && (y >= .1f && y <= .9f)) {
	d1 = h;
	acm = v2 / (d1 * d1);
    }
    if (yp >= .1f && yp <= .9f && (x >= .1f && x <= .9f)) {
	d1 = h;
	bcp = v2 / (d1 * d1);
    }
    if (ym >= .1f && ym <= .9f && (x >= .1f && x <= .9f)) {
	d1 = h;
	bcm = v2 / (d1 * d1);
    }
/*  dc=0. */
/*  ec=0. */
    dc = 100.;
    ec = 100.;

    icn = 0;
    sum = 0.;

    k = bindx[*location];
    bindx[k] = *row + 1;
    val[k] = -acp + dc * .5 / h;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 1;
    val[k] = -acm - dc * .5 / h;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (n + 2);
    val[k] = -bcp + ec * .5 / h;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (n + 2);
    val[k] = -bcm - ec * .5 / h;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = acp + acm + bcp + bcm;
L100:
    rhs[*row] = f;

/*  return 0;*/
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example17 */

/* ********************************************************************* */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example18(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1;
    static double v1, v2, dc, ec, fc;
    static int bk;
    static double xe, xm, ym, xp, yp, zp, xs, zm, acm, bcm, ccm, acp, bcp,
	     ccp;
    static int icn;
    static double sum1;

/*  Purpose: */
/*  Add one row to an MSR matrix corresponding to a 7pt discrete */
/*  approximation to the 3D EXAMPLE_18 operator on an n x n x n  cube. */

/*  Example 18: */
/*    The problem is taken from  I. G. Graham and M. J. Hagger paper, */
/*    SIAM J. Sci. Comp., 1999, pp. 2041-2066, "Unstructured additive */
/*    Schwarz-CG method for elliptic problems with highly discontinuous */
/*    coefficients". Our version enables also the addition of convectioni */
/*    terms. */

/*  Domain: [0 1]*[0 1]*[0 1]; */

/*   -(a*Ux)x - (b*Uy)y - (c*Uz)z + d*Ux + e*Uy + f*Uz + g*U =F */

/*  where: */
/*    d=e=f=100;  g=0; */

/*   and: */
/*    for: 0.333333 < x,y,z < 0.666666 a=b=c (=1.E+04) */
/*    else:                            a=b=c (=1.E+00) */

/*   F=0. */
/*   With Dirichlet boundary conditions: */
/*   U=1. on Y=0  and U=0 on all other planes */

/* ----------------------------------------------------------------------- */

/*  Parameters: */
/*     row          == global row number of the new row to be added. */
/*     location     == local row where diagonal of the new row will be stored. */
/*     val,bindx    == (see user's guide). On output, val[] and bindx[] */
/*                     are appended such that the new row has been added. */

/*  xs = 1./3.;
    xe = 2./3.;  */
    xs = 1.f/3.f;
    xe = 2.f/3.f;  

/*  v1= 1.e-01; */
/*  v2= 1.e-05; */
    v1 = 1.;
/*  v2= 1.e+06; */
    v2 = 1.e4;

    n = global_1.nx;
    h = 1. / (double) (global_1.nx + 1);
    i = *row % n;
    j = *row / n % n;
    k1 = *row / (n * n) % n;
    x = (i + 1) * h;
    y = (j + 1) * h;
    z = (k1 + 1) * h;
    xp = x + h * .5;
    xm = x - h * .5;
    yp = y + h * .5;
    ym = y - h * .5;
    zp = z + h * .5;
    zm = z - h * .5;

/*  Define the diffusion coefficients: */

    acp = v1;
    acm = v1;
    bcp = v1;
    bcm = v1;
    ccp = v1;
    ccm = v1;
    if (xp > xs && xp < xe && (y > xs && y < xe) && (z > xs && z < xe)) 
	acp = v2;
    if (xm > xs && xm < xe && (y > xs && y < xe) && (z > xs && z < xe)) 
	acm = v2;
    if (yp > xs && yp < xe && (x > xs && x < xe) && (z > xs && z < xe)) 
	bcp = v2;
    if (ym > xs && ym < xe && (x > xs && x < xe) && (z > xs && z < xe)) 
	bcm = v2;
    if (zp > xs && zp < xe && (x > xs && x < xe) && (y > xs && y < xe)) 
	ccp = v2;
    if (zm > xs && zm < xe && (x > xs && x < xe) && (y > xs && y < xe)) 
	ccm = v2;

/*  Define the convection coefficients: */

/*  dc=0. */
/*  ec=0. */
/*  fc=0. */
    dc = 200.;
    ec = 200.;
    fc = 200.;
/*  dc=100. */
/*  ec=100. */
/*  fc=100. */
/*  dc=10. */
/*  ec=10. */
/*  fc=10. */

    f =  0.;

/*  On the bottom  (U=1.0): */
    bk = *row - n * n;
    if (bk < 0 ) f = ccm + fc * .5 * h;

    icn = 0;
    sum1 = 0.;
/*  sumf=0.; */
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = -acp + dc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -acm - dc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = -bcp + ec * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -bcm - ec * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }
    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = -ccp + fc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -ccm - fc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }
    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = acp + acm + bcp + bcm + ccp + ccm;
/*  sumf= sumf + val[*location]; */
/*  f = sumf; */
    rhs[*row] = f;
/*  return 0; */

/*  Normalize the equations: */

    d1 = val[*location];
    sum1 += d1 * d1;
    sum1 = sqrt(sum1);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example18 */


/* ******************************************************************** */
/* ********************************************************************* */

/* Subroutine */ int create_matrix_row_7pt_example19(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1;
    static double v1, v2, dc, ec, fc;
    static int bk;
    static double xm, xp, xs, acm, bcm, ccm, acp, bcp, ccp;
    static int icn;
    static double sum1;

/*  Purpose: */
/*  Add one row to an MSR matrix corresponding to a 7pt discrete */
/*  approximation to the 3D EXAMPLE_19 operator on an n x n x n  cube. */

/*  Example 19:   
    The two domains model problem of L. Gerardo-Giorda, P. A. Tallec   
    and F. Nataf 2004 paper on advection-diffusion with discontinuous   
    coefficients.   

    Domain:  [0 1]*[0 1]*[0 1];   

     -(v(x)*Ux)x - (v(x)*Uy)y - (v(x)*Uz)z + d*Ux + e*Uy + f*Uz + g*U =F   

    Where:   
      g=1.;   
      d=e=f=prescribed ;   
    and:   
      for:   0. < x <= 0.5  v(x)=v1   
      for:   0.5< x < 1.0   v(x)=v2   

    F=0.   
    With Dirichlet boundary conditions:   
    U=1. on Y=0  and U=0 on all other boundaries. */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
    h = 1. / (double) (global_1.nx + 1);
    i = *row % n;
    j = *row / n % n;
    k1 = *row / (n * n) % n;
    x = (i + 1) * h;
    y = (j + 1) * h;
    z = (k1 + 1) * h;
    xp = x + h * .5;
    xm = x - h * .5;
    xs = .5;

/*  define the diffusion coefficients: */
    v1 = .1;
/*  v2= 1.e-07;*/
    v2 = 1e-5;
/*  v1= 1.e-05;*/
/*  v2= 1.e-01;*/

    acp = v1;
    acm = v1;
    bcp = v1;
    bcm = v1;
    ccp = v1;
    ccm = v1;
    if (xp > xs) {
	acp = v2;
    }
    if (xm > xs) {
	acm = v2;
    }
    if (x > xs) {
	bcp = v2;
	bcm = v2;
	ccp = v2;
	ccm = v2;
    }

/*   define the advection coefficients: */

/*  dc=1.;*/
/*  ec=0.;*/
/*  fc=0.;*/
    dc = 0.;
    ec = 1.;
    fc = 0.;

/*  dc=0.;*/
/*  ec=1.;*/
/*  fc=1.;*/

/*  dc=1.;*/
/*  ec=3.;*/
/*  fc=5.;*/

    f = 0.;
/*   on the bottom  (U=1.0): */
    bk = *row - n * n;
    if (bk < 0 ) f = ccm + fc * .5 * h;

    icn = 0;
    sum1 = 0.;
/*  sumf=0.;*/
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
        val[k] = -acp + dc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -acm - dc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = -bcp + ec * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -bcm - ec * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }
    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = -ccp + fc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -ccm - fc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }
    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = acp + acm + bcp + bcm + ccp + ccm + h * h;
/*  sumf= sumf + val[*location]; */
/*  f = sumf; */
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum1 += d1 * d1;
    sum1 = sqrt(sum1);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[*kmax - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example19 */

/* ************************************************************************ */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_5pt_example20(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nxloc, int *nyloc, 
	double *kn)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, pi;
    static int icn;
    static double sum, sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5pt discrete   
    approximation to the 2D EXAMPLE20 operator on an n x n  square.   

    Example20  taken from Erlangga, Vuik, Oosterlee, 2004 paper, in   
    applied Numerical Mathematics. Example #1 for Helemholtz eq..   

    Domain: [0 1]*[0 1]   

      Uxx + Uyy + K**2*U = (K**2 -5*PI**2)*sin(PI*X)*sin(2*PI*Y)   

    Dirichlet boundary conditions:  u=0 on the boundaries.   

    The analytic solution is:   
       u=sin(PI*X)*sin(2*PI*Y)   on the inner grid points. */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    pi = atan(1.) * 4.;
    h = 1. / (double) (global_2.n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    x = (*nxst + i) * h;
    y = (*nyst + j) * h;

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % global_2.n != global_2.n - 1) {
	val[k] = 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % global_2.n != 0) {
	val[k] = 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + global_2.n;
    if (*row / global_2.n % global_2.n != global_2.n - 1) {
	val[k] = 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - global_2.n;
    if (*row / global_2.n % global_2.n != 0) {
	val[k] = 1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
    }

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = *kn * *kn * h * h - 4.;
    f = (*kn * *kn - pi * 5. * pi) * sin(pi * x) * sin(pi * 2. * y);
    rhs[*row] = f * h * h;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[*kmax - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example20 */

/* ********************************************************************* */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_5pt_example21_2rows(int *row, 
	int *location, double *val, int *bindx, int *kmax, 
	double *rhs, int *nxst, int *nyst, int *nxloc, 
	int *nyloc, double *kn1, double *kn2, double *kn3)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y;
    static int ii;
    static double kn;
    static int icn;
    static double eps, sum, sum1;
    static int loca, rowa;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5pt discrete   
    approximation to the 2D EXAMPLE_21 operator on an (n+2)x(n+2) square.   

    Example 21 = similar to example #2 of Y.A. Erlangga, C. Vuik C.W. Oosterlee   
    2004 paper, in Applied Numerical Mathematics.   

    Domain: [0 1]*[0 1];   

            Uxx + Uyy + K**2*U = F   

    Where: F=0   

    With Dirichlet boundary condition on the bottom: u=0 on y=0, except at   
    the center point:   
              U=1   for x=0.5; y=0   
              U=0   for x#0.5; y=0   

    and Neumann b.c. Du/Dn -i*K*U=0 on all other boundaries. */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    eps = 1.e-5;
    h = 1. / (double) (global_2.n + 1);
    i = *location / 2 % *nxloc;
    j = *location / 2 / *nxloc;
    x = (*nxst + i - 1) * h;
    y = (*nyst + j - 1) * h;

/*  Note:  The results of the following single and double precision variable definitions may differ. 
    Though, the single def. in fortran = the single idef. in c. */

    if (y <= (1.f/3.f) ) kn = *kn1;
    if (y > (1.f/3.f) && y <= (2.f/3.f) ) *kn2;
    if (y > (2.f/3.f) ) *kn3;

/*  if (y <= (1./3.) ) kn = *kn1;
    if (y > (1./3.) && y <= (2./3.) ) kn = *kn2;
    if (y > (2./3.) ) kn = *kn3; */ 

    f = 0.;
/*  For a perturbation at the bottom center: */
    d1=x - .5;
    if (*nyst + j - 1 == 0 && fabs(d1) <= h * .5 + eps)  f = 1.;

/*  Set the  boundary conditions equations   
    on the upper boundary:  du/dn -i*Kn*u=0:  
    where du/dn is the outward derivative. */

    if (j + *nyst - 1 == global_2.n + 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - (global_2.n + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = kn * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*     the lower boundary (Dirichlet b.c.): */

    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*  the left hand-side boundary (Neumann b.c.) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < 
	    global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	bindx[k] = *row + 1;
	val[k] = kn * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*     the right hand-side boundary  (Neumann b.c. ) */
    if (i + *nxst - 1 == global_2.n + 1 && (j + *nyst - 1 > 0 && j + *nyst 
	    - 1 < global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	bindx[k] = *row + 1;
	val[k] = kn * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }
/* --------    end of the b.c. equations setting -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 2;
    val[k] = 1.;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 2;
    val[k] = 1.;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 2 << 1);
    val[k] = 1.;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (global_2.n + 2 << 1);
    val[k] = 1.;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = kn * kn * h * h - 4.;
L100:
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
/*  printf("\n%7d %13.10f %13.10f", *row, rhs[*row], val[*location]); */
/*  for (ii=1; ii<=i1; ++ii) printf(" %13.10f", val[*kmax - ii]) ; */
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[*kmax - ii] /= sum1;
    }

/* ----------------------------------- */
/* ----------------------------------- */
/*  Store the imaginary equation coefficients (B; A) */

    loca = *location + 1;
    rowa = *row + 1;
    f = 0.;

/*  Set the  b.c. equations: du/dn -i*Kn*u=0: */

/*  The upper boundary: */
    if (j + *nyst - 1 == global_2.n + 1) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - (global_2.n + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The lower boundary (Dirichlet b.c.): */

    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	*kmax = k;
	bindx[loca + 1] = k;
	val[loca] = 1.;
	goto L200;
    }

/*  The left hand-side boundary (Neumann b.c.): */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < 
	    global_2.n + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The right hand-side boundary  (Neumann b.c.): */
    if (i + *nxst - 1 == global_2.n + 1 && (j + *nyst - 1 > 0 && j + *nyst 
	    - 1 < global_2.n + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * h;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }
/* --------    end of the b.c. equation settings -------- */


/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    icn = 0;
    sum = 0.;
    k = bindx[loca];
    bindx[k] = rowa + 2;
    val[k] = 1.;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - 2;
    val[k] = 1.;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (global_2.n + 2 << 1);
    val[k] = 1.;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (global_2.n + 2 << 1);
    val[k] = 1.;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    *kmax = k;
    bindx[loca + 1] = k;
    val[loca] = kn * kn * h * h - 4.;
L200:
    rhs[rowa] = f;

/*  return */
/*  Normalize the equations: */

    d1 = val[loca];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[rowa] /= sum1;
    val[loca] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[*kmax - ii] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example21_2rows */

/* ********************************************************************* */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_7pt_example22(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nzst, int *nxloc, 
	int *nyloc, int *nzloc)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f, h;
    static int i, j, k, n;
    static double x, y, z;
    static int k1;
    static double v1, v2, v3, v4, v5, v6, v7, v8, dc, ec, fc;
    static int  bk;
    static double pi, xm, ym, xp, yp, zm, xs, ys, zs, zp, acm, bcm, ccm, 
	    acp, bcp, ccp;
    static int icn;
    static double sum1;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 7pt discrete   
    approximation to the 3D EXAMPLE_22 operator on an n x n x n  cube.   

    Example 22:   
    The two domains model problem of L. Gerardo-Giorda, P. A. Tallec   
    and F. Nataf 2004 paper on advection-diffusion with discontinuous   
    coefficients.   

    Domain:  [-0.5 0.5]*[0.5 -0.5]*[0 1];   

     -(v(x)*Ux)x - (v(x)*Uy)y - (v(x)*Uz)z + d*Ux + e*Uy + f*Uz + g*U =F   

    Where:   
      g=1.;   
      d=e=f=prescribed ;   

    F=0.   
    With Dirichlet boundary conditions:   
    U=1. on Z=0  and U=0 on all other boundaries. */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    n = global_1.nx;
    pi = atan(1.) * 4.;
    h = 1. / (double) (n + 1);
    i = *row % n;
    j = *row / n % n;
    k1 = *row / (n * n) % n;
/*  i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    k1 = *location / (*nxloc * *nyloc) % *nzloc; */
    x = (i + 1) * h - .5;
    y = (j + 1) * h - .5;
    z = (k1 + 1) * h;
    xp = x + h * .5;
    xm = x - h * .5;
    yp = y + h * .5;
    ym = y - h * .5;
    zp = z + h * .5;
    zm = z - h * .5;
/*  xs = 0.;
    ys = 0.;
    zs = .5;  */ 
    xs = 0.f;
    ys = 0.f;
    zs = .5f;   

/*  Define the diffusion coefficients: */
    v1 = .1;
    v2 = 1.e-6;
    v3 = .01; 
    v4 = 1.e-6;
    v5 = 1.e-6;
    v6 = .001; 
    v7 = 1.e-6;
    v8 = 1.e-4;

    if (xp < xs && y < ys && z < zs)     acp = v1;
    if (xm < xs && y < ys && z < zs)     acm = v1;
    if (xp < xs && y >= ys && z < zs)    acp = v2;
    if (xm < xs && y >= ys && z < zs)    acm = v2;   
    if (xp >= xs && y >= ys && z < zs)   acp = v3;
    if (xm >= xs && y >= ys && z < zs)   acm = v3;
    if (xp >= xs && y < ys && z < zs)    acp = v4;
    if (xm >= xs && y < ys && z < zs)    acm = v4; 


    if (xp < xs && y < ys && z >= zs)     acp = v5;
    if (xm < xs && y < ys && z >= zs)     acm = v5;
    if (xp < xs && y >= ys && z >= zs)    acp = v6;
    if (xm < xs && y >= ys && z >= zs)    acm = v6;
    if (xp >= xs && y >= ys && z >= zs)   acp = v7;
    if (xm >= xs && y >= ys && z >= zs)   acm = v7;
    if (xp >= xs && y < ys && z >= zs)    acp = v8;
    if (xm >= xs && y < ys && z >= zs)    acm = v8; 


    if (yp < ys && x < xs && z < zs)     bcp = v1; 
    if (ym < ys && x < xs && z < zs)     bcm = v1;
    if (yp >= ys && x < xs && z < zs)    bcp = v2;
    if (ym >= ys && x < xs && z < zs)    bcm = v2;
    if (yp >= ys && x >= xs && z < zs)   bcp = v3;
    if (ym >= ys && x >= xs && z < zs)   bcm = v3;
    if (yp < ys && x >= xs && z < zs)    bcp = v4;
    if (ym < ys && x >= xs && z < zs)    bcm = v4; 


    if (yp < ys && x < xs && z >= zs)     bcp = v5;
    if (ym < ys && x < xs && z >= zs)     bcm = v5;
    if (yp >= ys && x < xs && z >= zs)    bcp = v6;
    if (ym >= ys && x < xs && z >= zs)    bcm = v6;
    if (yp >= ys && x >= xs && z >= zs)   bcp = v7;
    if (ym >= ys && x >= xs && z >= zs)   bcm = v7;
    if (yp < ys && x >= xs && z >= zs)    bcp = v8;
    if (ym < ys && x >= xs && z >= zs)    bcm = v8;   


    if (zp < zs && x < xs && y < ys)     ccp = v1;
    if (zm < zs && x < xs && y < ys)     ccm = v1;
    if (zp < zs && x < xs && y >= ys)    ccp = v2;
    if (zm < zs && x < xs && y >= ys)    ccm = v2;
    if (zp < zs && x >= xs && y >= ys)   ccp = v3;
    if (zm < zs && x >= xs && y >= ys)   ccm = v3;
    if (zp < zs && x >= xs && y < ys)    ccp = v4;
    if (zm < zs && x >= xs && y < ys)    ccm = v4; 


    if (zp >= zs && x < xs && y < ys)     ccp = v5;
    if (zm >= zs && x < xs && y < ys)     ccm = v5;
    if (zp >= zs && x < xs && y >= ys)    ccp = v6;
    if (zm >= zs && x < xs && y >= ys)    ccm = v6;
    if (zp >= zs && x >= xs && y >= ys)   ccp = v7;
    if (zm >= zs && x >= xs && y >= ys)   ccm = v7;
    if (zp >= zs && x >= xs && y < ys)    ccp = v8;
    if (zm >= zs && x >= xs && y < ys)    ccm = v8;   

/*  Define the advection coefficients: */

    dc = pi * -2. * y;
    ec = pi * 2. * x;
    fc = sin(pi * 2. * x);

    f = 0.;
/*  On the bottom  (U=1.0): */
    bk = *row - n * n; 
    if (bk < 0 )  f = ccm + fc * .5 * h ;

    icn = 0;
    sum1 = 0.;
/*  sumf=0.; */
    k = bindx[*location];
    bindx[k] = *row + 1;
    if (*row % n != n - 1) {
	val[k] = -acp + dc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - 1;
    if (*row % n != 0) {
	val[k] = -acm - dc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row + n;
    if (*row / n % n != n - 1) {
	val[k] = -bcp + ec * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n;
    if (*row / n % n != 0) {
	val[k] = -bcm - ec * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }
    bindx[k] = *row + n * n;
    if (bindx[k] < n * n * n) {
	val[k] = -ccp + fc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }

    bindx[k] = *row - n * n;
    if (bindx[k] >= 0) {
	val[k] = -ccm - fc * .5 * h;
/*      sumf= sumf + val[k]; */
	d1 = val[k];
	sum1 += d1 * d1;
	++icn;
	++k;
    }
    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = acp + acm + bcp + bcm + ccp + ccm + h * h;
/*  sumf= sumf + val[*location]; */
/*  f = sumf; */
    rhs[*row] = f;
/*  return */

/*  Normalize the equations: */
    d1 = val[*location];
    sum1 += d1 * d1;
    sum1 = sqrt(sum1);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
/*  printf("\n%7d %13.10f %13.10f", *row, rhs[*row], val[*location]);  */ 
/*  if(i == 40 &&  k1 == 0) {
    printf("\n   %s", "    ");
    printf("\n%7d %13.10f %13.10f %7d %7d %7d %13.10f %13.10f %13.10f ", 
        *row, rhs[*row], val[*location], i, j, k1, x, y, z); 
    printf("\n%7d %13.10f %13.10f %3d %3d %3d %13.10f %13.10f %13.10f %13.10f %13.10f %13.10f", 
        *row, rhs[*row], val[*location], i, j, k1,  xp,  xm,  yp,  ym,  zp,  zm); 
    printf("\n%7d %13.10f %13.10f %3d %3d %3d %13.10f %13.10f %13.10f %13.10f %13.10f %13.10f", 
        *row, rhs[*row], val[*location], i, j, k1, acp, acm, bcp, bcm, ccp, ccm); 
    } */
    i1 = icn;
    for (i = 1; i <= i1; ++i) {
	val[k - i] /= sum1;
    }

    return 0;
} /* create_matrix_row_7pt_example22 */

/* ************************************************************************ */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_9pt_example23_2rows(int *row, 
	int *location, double *val, int *bindx, int *kmax, 
	double *rhs, int *nxst, int *nyst, int *nxloc, 
	int *nyloc, double *kn1, double *kn2, double *kn3, 
	int *order)
{
    /* System generated locals */
    int i1;
    double d1, d2;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, a0, ac;
    static int ii;
    static double as, kn, del;
    static int icn;
    static double eps, sum, hkn2, sum1, beta;
    static int loca, rowa;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 9pt discrete   
    approximation to the 2D EXAMPLE_23 operator on an (n+3)x(n+4) mesh.   

    Example 23 = similar to example #2 of Y.A. Erlangga, C. Vuik and   
    C.W. Oosterlee 2004 paper, in Applied Numerical Mathematics.   
    The variable Order defines the scheme order: 2nd / 4th/ 6th order scheme.   
    The 4th/6th order finite difference schemes follow Singer and   
    Turkel 2006, Erlanga and Turkel  2011? and Harrari and Turkel 1995.   
    with gamma=14/5, and DEL= -1.;0. or. 2. for the 6th order scheme.   

    Domain: [0 1]*[0 1];   

            Uxx + Uyy + K**2*U = F   

    Where: F=0   

    With Dirichlet boundary condition on the bottom: u=0 on y=0, except at   
    the center point:   
             U=1   for x=0.5; y=0   
             U=0   for x#0.5; y=0   

    and Neumann b.c. Du/Dn -i*K*U=0 on all other boundaries. */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */
    eps = 1.e-5;
    h = 1./ (double) (global_2.n + 1);
    i = *location / 2 % *nxloc;
    j = *location / 2 / *nxloc;
    x = (*nxst + i - 1) * h;
    y = (*nyst + j - 1) * h;

    if (y <= (1.f/3.f)) kn = *kn1;
    if (y > (1.f/3.f) && y <= (2.f/3.f) ) kn = *kn2;
    if (y > (2.f/3.f) ) kn = *kn3;

    beta = -kn;
    f = 0.;
    d1=x - .5;
    if (*nyst + j - 1 == 0 && fabs(d1) <= h * .5 + eps) f = 1.;
/*  if (*nyst + j - 1 == 0 && fabs(d1) <= h * .5 + eps) f = h*h; */

/*  Set the  boundary conditions equations */
/*  On the upper boundary:  du/dn -i*Kn*u=0: */
/*  where du/dn is the outward derivative. */

    if (j + *nyst - 1 == global_2.n + 2) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - (global_2.n + 4 << 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row - (global_2.n + 4 << 1) + 1;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * -2. * h * (1. - d1 * d1 / 6. + d2 * d2 
		    / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * -2. * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[*location] = 1.;
	goto L100;
    }

/*  The lower boundary (Dirichlet b.c.) */

    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*  The left hand-side boundary (Neumann b.c.: du/dn -i*Kn*u=0) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 <= 
	    global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part */
	bindx[k] = *row + 3;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
/* Computing 4th power */
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * -2. * h * (1. - d1 * d1 / 6. + d2 * d2 
		    / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * -2. * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	if (*order == 2) {
	    goto L35;
	}
	bindx[k] = *row + (global_2.n + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part */
	bindx[k] = *row + (global_2.n + 4 << 1) + 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the real part */
	bindx[k] = *row - (global_2.n + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginay part */
	bindx[k] = *row - (global_2.n + 4 << 1) + 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L35:
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[*location] = 1.;
	goto L100;
    }

/*  The right hand-side boundary  (Neumann b.c.: du/dn -i*Kn*u=0 ) */
    if (i + *nxst - 1 == global_2.n + 3 && (j + *nyst - 1 > 0 && j + *nyst 
	    - 1 <= global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginay part */
	bindx[k] = *row - 1;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
/* Computing 4th power */
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * -2. * h * (1. - d1 * d1 / 6. + d2 * d2 
		    / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * -2. * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	if (*order == 2) {
	    goto L45;
	}
	bindx[k] = *row + (global_2.n + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part */
	bindx[k] = *row + (global_2.n + 4 << 1) + 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the real part */
	bindx[k] = *row - (global_2.n + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginay part */
	bindx[k] = *row - (global_2.n + 4 << 1) + 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L45:
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[*location] = 1.;
	goto L100;
    }
/* --------    end of the b.c. equations setting -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    hkn2 = kn * kn * h * h;
    if (*order == 2) {
	a0 = hkn2 - 4.;
	as = 1.;
	ac = 0.;
    }
    if (*order == 4) {
	a0 = hkn2 * 67. / 90. - 10./3.;
	as = hkn2 * 2. / 45. + 2./3.;
	ac = hkn2 * 7. / 360. + 1./6.;
    }
    if (*order == 6) {
/*      del= 2. */
	del = -1.;
/*      del= 0. */
	d1 = hkn2;
	a0 = hkn2 * 67. / 90. - 10./3. + d1 * d1 * (del - 3.)  / 180.;
	d1 = hkn2;
	as = hkn2 * 2. / 45. + 2./3. + d1 * d1 * (3. - del *  2.) / 720.;
	d1 = hkn2;
	ac = hkn2 * 7. / 360. + 1./6. + d1 * d1 * del / 720.;
    }
    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 4 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L50;
    }

    bindx[k] = *row + (global_2.n + 4 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 4 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L50:

    bindx[k] = *row - (global_2.n + 4 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L55;
    }
    bindx[k] = *row - (global_2.n + 4 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (global_2.n + 4 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L55:

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = a0;
L100:
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[k - ii] /= sum1;
    }

/* --------------- */
/* --------------- */
/*  Store the imaginary equation coefficients (B; A) */

    loca = *location + 1;
    rowa = *row + 1;
    f = 0.;

/*  Set the  b.c. equations: du/dn -i*Kn*u=0: */
/*  The upper boundary: */

    if (j + *nyst - 1 == global_2.n + 2) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - (global_2.n + 4 << 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - (global_2.n + 4 << 1) - 1;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
/* Computing 4th power */
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * 2. * h * (1. - d1 * d1 / 6. + d2 * 
		    d2 / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2 * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[loca] = 1.;
	goto L200;
    }

/*  The lower boundary (Dirichlet b.c.) */

    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	*kmax = k;
	bindx[loca + 1] = k;
	val[loca] = 1.;
	goto L200;
    }

/*  The left hand-side boundary (Neumann b.c.) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 <= 
	    global_2.n + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa + 1;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
/* Computing 4th power */
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * 2. * h * (1. - d1 * d1 / 6. + d2 * 
		    d2 / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2. * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	if (*order == 2) {
	    goto L65;
	}

	bindx[k] = rowa + (global_2.n + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	bindx[k] = rowa + (global_2.n + 4 << 1) - 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	bindx[k] = rowa - (global_2.n + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa - (global_2.n + 4 << 1) - 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

L65:
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The right hand-side boundary  (Neumann b.c. ) */
    if (i + *nxst - 1 == global_2.n + 3 && (j + *nyst - 1 > 0 && j + *nyst 
	    - 1 <= global_2.n + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 3;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
/* Computing 4th power */
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * 2. * h * (1. - d1 * d1 / 6. + d2 *  d2 / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2. * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	if (*order == 2) {
	    goto L75;
	}

	bindx[k] = rowa + (global_2.n + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	bindx[k] = rowa + (global_2.n + 4 << 1) - 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	bindx[k] = rowa - (global_2.n + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa - (global_2.n + 4 << 1) - 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

L75:
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[loca] = 1.;
	goto L200;
    }
/* --------    end of the b.c. equation setting -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    icn = 0;
    sum = 0.;
    k = bindx[loca];
    bindx[k] = rowa + 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (global_2.n + 4 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L150;
    }

    bindx[k] = rowa + (global_2.n + 4 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (global_2.n + 4 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L150:

    bindx[k] = rowa - (global_2.n + 4 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L155;
    }

    bindx[k] = rowa - (global_2.n + 4 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (global_2.n + 4 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L155:

    *kmax = k;
    bindx[loca + 1] = k;
    val[loca] = a0;
L200:
    rhs[rowa] = f;

/*  return */
/*  Normalize the equations: */
    d1 = val[loca];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[rowa] /= sum1;
    val[loca] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[k - ii] /= sum1;
    }

    return 0;
} /* create_matrix_row_9pt_example23_2rows */
/* ********************************************************************** */
/* ********************************************************************** */
/* Subroutine */ int create_matrix_row_5_9_pt_example24_2rows(int *row, 
	int *location, double *val, int *bindx, int *kmax, 
	double *rhs, int *nxst, int *nyst, int *nxloc, 
	int *nyloc, double *cxy, double *xmin, double *xmax, 
	double *ymin, double *ymax, int *nx, int *ny, int 
	*order, double *freq)
{
    /* System generated locals */
    int i1;
    double d1, d2;

    /* Builtin functions */
/*  double atan(double), sqrt(double); */

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double ac_im1_jm1, ac_ip1_jm1, ac_im1_jp1, ac_ip1_jp1, 
	    x, y, kn_im1_jm1, kn_im1_jp1, kn_ip1_jm1, kn_ip1_jp1, a0, 
	    ac;
    static int ii;
    static double as, pi, hx, hy, del;
    static int icn;
    static double  sum, sum1, beta;
    static int loca, rowa;
    static double cons;
    static double kn_i_j, as_im1_j, as_i_jm1, as_ip1_j, 
	    as_i_jp1, kn_im1_j, kn_i_jm1, kn_ip1_j, kn_i_jp1;


/*  Purpose: */
/*  Add one row to an MSR matrix corresponding to a 5 / 9 pt discrete */
/*  approximation to the 2D Marmousi operator on an (nx+4)x(ny+4) domain. */

/*  Example 24 = example #2 of O. Schenk.  2009 */
/*  paper, SIAM J. Scientific Computing. */

/*  Domain: [Xmin-hx  Xmax+hx]*[Ymin-hy Ymax+hy]; */

/*          Uxx + Uyy + K**2*U = F */

/*          K=2*PI*freq/C */

/*   Neumann b.c. Du/Dn -i*K*U=gi on all  boundaries. */

/*  and: */
/*            f=1   for x=0.5, y=0 */
/*            f=0    else */


/* ----------------------------------------------------------------------- */

/*  Parameters: */
/*     row          == global row number of the new row to be added. */
/*     location     == local row where diagonal of the new row will be stored. */
/*     val,bindx    == (see user's guide). On output, val[] and bindx[] */
/*                     are appended such that the new row has been added. */


    hx = (*xmax - *xmin) / (*nx - 1);
    hy = (*ymax - *ymin) / (*ny - 1);
    i = *location / 2 % *nxloc;
    j = *location / 2 / *nxloc;
    x = *xmin - hx + (*nxst + i - 1) * hx;
    y = *ymin - hy + (*nyst + j - 1) * hy;
/* here we assume:  H = Hx = Hy ; */
    h = hx;

    pi = atan(1.) * 4.;
    cons = pi * 2. * *freq;
    kn_i_j = cons / cxy[*row / 2];
/* note: on the grid points outside the domain cxy() is defined as that of the nearest */
/*       grid point of the domain for which the finite difference scheme is applied. */
    if (*order > 2 && i + *nxst - 1 > 0 && i + *nxst - 1 < *nx + 3 && j + 
	    *nyst - 1 > 0 && j + *nyst - 1 < *ny + 3) {
	kn_im1_j = cons / cxy[*row / 2 - 1];
	kn_ip1_j = cons / cxy[*row / 2 + 1];
	kn_ip1_jp1 = cons / cxy[*row / 2 + (*nx + 4) + 1];
	kn_i_jp1 = cons / cxy[*row / 2 + (*nx + 4)];
	kn_im1_jp1 = cons / cxy[*row / 2 + (*nx + 4) - 1];
	kn_ip1_jm1 = cons / cxy[*row / 2 - (*nx + 4) + 1];
	kn_i_jm1 = cons / cxy[*row / 2 - (*nx + 4)];
	kn_im1_jm1 = cons / cxy[*row / 2 - (*nx + 4) - 1];
    }
/*  store the real equation coefficients (ai;-bi) */
    f = 0.;
/*  for a perturbation at the bottom center: */
/*  after the addition of the two external grid point at: (Xmin-h;*) and (Xmax+h;*): */
    if (*nyst + j - 1 == 1 && *nxst + i - 1 == (*nx + 3) / 2) {
	f = 1.;
    }

/*  Set the  boundary conditions equations: */
/* ----------------------------------------- */
/*  on the upper boundary:  du/dn -i*Kn*u=0: */
/*  where du/dn is the outward derivative. */

    if (j + *nyst - 1 == *ny + 3) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
/*  Note: the b.c. equation is set on the point on the boundary (not */
/*  the additional external point) */
	bindx[k] = *row - (*nx + 4 << 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row - (*nx + 4 << 1) + 1;
	beta = -kn_i_j;
	if (*order == 4 || *order == 6) {
	    d1 = beta * hy;
	    d2 = beta * hy, d2 *= d2;
	    val[k] = beta * -2 * hy * (1. - d1 * d1 / 6. + d2 * d2 /
		     120.);
	}
	if (*order == 2) {
	    val[k] = beta * -2. * hy;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[*location] = 1.;
	goto L100;
    }

/*     the lower boundary (Neumann b.c.  ) */

    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + (*nx + 4 << 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + (*nx + 4 << 1) + 1;
	beta = -kn_i_j;
	if (*order == 4 || *order == 6) {
	    d1 = beta * hy;
	    d2 = beta * hy, d2 *= d2;
	    val[k] = beta * -2. * hy * (1. - d1 * d1 / 6. + d2 * 
		    d2 / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * -2. * hy;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[*location] = 1.;
	goto L100;
    }

/*     the left hand-side boundary (Neumann b.c.) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 3)
	    ) {
	icn = 0;
	sum = 0.f;
	k = bindx[*location];
	bindx[k] = *row + 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 3;
	beta = -kn_i_j;
	if (*order == 4 || *order == 6) {
	    d1 = beta * hx;
	    d2 = beta * hx, d2 *= d2;
	    val[k] = beta * -2. * hx * (1. - d1 * d1 / 6. + d2 * 
		    d2 / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * -2. * hx;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	if (*order == 2) {
	    goto L35;
	}
/*  set zero value for the point exactly above the current (i=0) grid point: */
/*  the real part */
	bindx[k] = *row + (*nx + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part */
	bindx[k] = *row + (*nx + 4 << 1) + 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  set zero value for the point exactly below the grid point */
/*  the real part */
	bindx[k] = *row - (*nx + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginay part */
	bindx[k] = *row - (*nx + 4 << 1) + 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L35:
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element (the grid point outside the computational domain): */
	val[*location] = 1.;
	goto L100;
    }

/*     the right hand-side boundary  (Neumann b.c. ) */
    if (i + *nxst - 1 == *nx + 3 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *
	    ny + 3)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row - 1;
	beta = -kn_i_j;
	if (*order == 4 || *order == 6) {
	    d1 = beta * hx;
	    d2 = beta * hx, d2 *= d2;
	    val[k] = beta * -2. * hx * (1. - d1 * d1 / 6. + d2 * 
		    d2 / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * -2. * hx;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	if (*order == 2) {
	    goto L45;
	}
/*  set zero value for the point exactly above the grid point */
/*  the real part */
	bindx[k] = *row + (*nx + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part */
	bindx[k] = *row + (*nx + 4 << 1) + 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  set zero value for the point exactly below the grid point */
/*  the real part */
	bindx[k] = *row - (*nx + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginay part */
	bindx[k] = *row - (*nx + 4 << 1) + 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L45:
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[*location] = 1.;
	goto L100;
    }
/* --------    end of the b.c. equations setting -------- */

/*   set the eq. for the inner grid points. */
/*   In the present example there are no imaginary terms */
/*   for the inner grid points */

/* here we assume:  H = Hx = Hy ; */
    if (*order == 2) {
	a0 = kn_i_j * kn_i_j * h * h - 4.;
	as = 1.;
	ac = 0.;
	as_i_jm1 = as;
	as_i_jp1 = as;
	as_im1_j = as;
	as_ip1_j = as;
	ac_im1_jm1 = ac;
	ac_ip1_jm1 = ac;
	ac_im1_jp1 = ac;
	ac_ip1_jp1 = ac;
    }
    if (*order == 4) {
/*        DEL= 2. */
/*        DEL=-1. */
/*        DEL= 2.8 !for this value the coefficiens are:  67./90.; 2./45. and 7./360. */
	del = 0.;
/*        A0  = -10.D0/3.D0 + HKN2*67.D0/90.D0 */
	d1 = kn_i_j * h;
	a0 = d1 * d1 * (del + 24.) / 36. - 10./3.;
/*  symmetric points: */
	d1 = kn_i_jm1 * h;
	as_i_jm1 = d1 * d1 * (6. - del) / 72. + 2./3.;
	d1 = kn_i_jp1 * h;
	as_i_jp1 = d1 * d1 * (6. - del) / 72. + 2./3.;
	d1 = kn_im1_j * h;
	as_im1_j = d1 * d1 * (6. - del) / 72. + 2./3.;
	d1 = kn_ip1_j * h;
	as_ip1_j = d1 * d1 * (6. - del) / 72. + 2./3.;
/*  corner points: */
	d1 = kn_im1_jm1 * h;
	ac_im1_jm1 = d1 * d1 * del / 144. + 1./6.;
	d1 = kn_ip1_jm1 * h;
	ac_ip1_jm1 = d1 * d1 * del / 144. + 1./6.;
	d1 = kn_im1_jp1 * h;
	ac_im1_jp1 = d1 * d1 * del / 144. + 1./6.;
	d1 = kn_ip1_jp1 * h;
	ac_ip1_jp1 = d1 * d1 * del / 144. + 1./6.;
    }
/* endif 4th order */
    if (*order == 6) {
/*        DEL= 2. */
	del = -1.;
/*        DEL= 0. */
	d1 = kn_i_j * h;
	d2 = kn_i_j * h, d2 *= d2;
	a0 = d1 * d1 * 67. / 90. - 10./3. + d2 * d2 * (
		del - 3.) / 180.;

/*  symmetric points: */
	d1 = kn_i_jm1 * h;
	d2 = kn_i_jm1 * h, d2 *= d2;
	as_i_jm1 = d1 * d1 * 2. / 45. + 2./3. + d2 * 
		d2 * (3. - del * 2.) / 720.;
	d1 = kn_i_jp1 * h;
	d2 = kn_i_jp1 * h, d2 *= d2;
	as_i_jp1 = d1 * d1 * 2. / 45. + 2./3. + d2 * 
		d2 * (3. - del * 2.) / 720.;
	d1 = kn_im1_j * h;
	d2 = kn_im1_j * h, d2 *= d2;
	as_im1_j = d1 * d1 * 2. / 45. + 2./3. + d2 * 
		d2 * (3. - del * 2.) / 720.;
	d1 = kn_ip1_j * h;
	d2 = kn_ip1_j * h, d2 *= d2;
	as_ip1_j = d1 * d1 * 2. / 45. + 2./3. + d2 * 
		d2 * (3. - del * 2.) / 720.;
/*   corner points: */
	d1 = kn_im1_jm1 * h;
	d2 = kn_im1_jm1 * h, d2 *= d2;
	ac_im1_jm1 = d1 * d1 * 7. / 360. + 1./6. + d2 * 
		d2 * del / 720.;
	d1 = kn_ip1_jm1 * h;
	d2 = kn_ip1_jm1 * h, d2 *= d2;
	ac_ip1_jm1 = d1 * d1 * 7. / 360. + 1./6. + d2 * 
		d2 * del / 720.;
	d1 = kn_im1_jp1 * h;
	d2 = kn_im1_jp1 * h, d2 *= d2;
	ac_im1_jp1 = d1 * d1 * 7. / 360. + 1./6. + d2 * 
		d2 * del / 720.;
	d1 = kn_ip1_jp1 * h;
	d2 = kn_ip1_jp1 * h, d2 *= d2;
	ac_ip1_jp1 = d1 * d1 * 7. / 360. + 1./6. + d2 * 
		d2 * del / 720.;
    }
/* endif 6th order */
    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 2;
    val[k] = as_ip1_j;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 2;
    val[k] = as_im1_j;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (*nx + 4 << 1);
    val[k] = as_i_jp1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (*nx + 4 << 1);
    val[k] = as_i_jm1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L50;
    }

    bindx[k] = *row + (*nx + 4 << 1) + 2;
    val[k] = ac_ip1_jp1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (*nx + 4 << 1) - 2;
    val[k] = ac_im1_jp1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (*nx + 4 << 1) + 2;
    val[k] = ac_ip1_jm1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (*nx + 4 << 1) - 2;
    val[k] = ac_im1_jm1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;
L50:

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = a0;
L100:
    rhs[*row] = f;

/*     return */
/*  normalize the equations: */
    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[*kmax - ii] /= sum1;
    }

/* --------------- */
/* --------------- */
/*    store the imaginary equation coefficients (bi; ai) */

    loca = *location + 1;
    rowa = *row + 1;
    f = 0.;

/*    Set the  b.c. equations: du/dn -i*Kn*u=0: */

/*     the upper boundary: */

    if (j + *nyst - 1 == *ny + 3) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - (*nx + 4 << 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - (*nx + 4 << 1) - 1;
	beta = -kn_i_j;
	if (*order == 4 || *order == 6) {
	    d1 = beta * hy;
	    d2 = beta * hy, d2 *= d2;
	    val[k] = beta * 2. * hy * (1. - d1 * d1 / 6. + d2 * d2 
		    / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2 * hy;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[loca] = 1.;
	goto L200;
    }

/*     the lower boundary: du/dn -i*Kn*u=0 */

    if (j + *nyst - 1 == 0) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + (*nx + 4 << 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa + (*nx + 4 << 1) - 1;
	beta = -kn_i_j;
	if (*order == 4 || *order == 6) {
	    d1 = beta * hy;
	    d2 = beta * hy, d2 *= d2;
	    val[k] = beta * 2. * hy * (1. - d1 * d1 / 6. + d2 * d2 
		    / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2. * hy;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
/* the additional point - outside the computational domain: */
	val[loca] = 1.;
	goto L200;
    }

/*     the left hand-side boundary (Neumann b.c.) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 3)
	    ) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
/*        bindx(k)  = rowA - 1 */
	bindx[k] = rowa + 1;
	beta = -kn_i_j;
	if (*order == 4 || *order == 6) {
	    d1 = beta * hx;
	    d2 = beta * hx, d2 *= d2;
	    val[k] = beta * 2. * hx * (1. - d1 * d1 / 6. + d2 * d2 
		    / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2. * hx;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	if (*order == 2) {
	    goto L65;
	}

	bindx[k] = rowa + (*nx + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa + (*nx + 4 << 1) - 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa - (*nx + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa - (*nx + 4 << 1) - 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L65:
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*     the right hand-side boundary  (Neumann b.c. ) */
    if (i + *nxst - 1 == *nx + 3 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *
	    ny + 3)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 3;
	beta = -kn_i_j;
	if (*order == 4 || *order == 6) {
	    d1 = beta * hx;
	    d2 = beta * hx, d2 *= d2;
	    val[k] = beta * 2. * hx * (1. - d1 * d1 / 6. + d2 * d2 
		    / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2. * hx;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	if (*order == 2) {
	    goto L75;
	}

	bindx[k] = rowa + (*nx + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa + (*nx + 4 << 1) - 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa - (*nx + 4 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa - (*nx + 4 << 1) - 1;
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L75:

	*kmax = k;
	bindx[loca + 1] = k;
/* the diagonal element: */
/* the additional point - outside the computational domain: */
	val[loca] = 1.;
	goto L200;
    }
/* --------    end of the b.c. equation settings -------- */


/*   set the eq. for the inner grid points. */
/*   In the present example there are no imaginary terms */
/*   for the inner grid points */

    icn = 0;
    sum = 0.;
    k = bindx[loca];
    bindx[k] = rowa + 2;
    val[k] = as_ip1_j;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - 2;
    val[k] = as_im1_j;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (*nx + 4 << 1);
    val[k] = as_i_jp1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (*nx + 4 << 1);
    val[k] = as_i_jm1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L150;
    }

    bindx[k] = rowa + (*nx + 4 << 1) + 2;
    val[k] = ac_ip1_jp1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (*nx + 4 << 1) - 2;
    val[k] = ac_im1_jp1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (*nx + 4 << 1) + 2;
    val[k] = ac_ip1_jm1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (*nx + 4 << 1) - 2;
    val[k] = ac_im1_jm1;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;
L150:
    *kmax = k;
    bindx[loca + 1] = k;
    val[loca] = a0;
L200:
    rhs[rowa] = f;

/*     return */
/*  normalize the equations: */

    d1 = val[loca];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[rowa] /= sum1;
    val[loca] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[*kmax - ii] /= sum1;
    }

    return 0;
} /* create_matrix_row_5_9_pt_example24_2rows */

/* ************************************************************************ */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_5pt_example24_2rows(int *row, 
	int *location, double *val, int *bindx, int *kmax, 
	double *rhs, int *nxst, int *nyst, int *nxloc, 
	int *nyloc, double *kn, double *xmin, double *xmax, 
	double *ymin, double *ymax, int *nx, int *ny)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f;
    static int i, j, k;
    static double x, y;
    static int ii;
    static double hx, hy;
    static int icn;
    static double  sum, sum1;
    static int loca, rowa;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5pt discrete   
    approximation to the 2D Marmousi operator on an (nx+2)x(ny+2) domain.   

    Example 24 = example #2 of O. Schenk.  2009   
    paper, SIAM J. Scientific Computing.   

    Domain: [Xmin-hx  Xmax+hx]*[Ymin-hy Ymax+hy];   

            Uxx + Uyy + K**2*U = F   

            K=2*PI*freq/C   

    Neumann b.c. Du/Dn -i*K*U=gi on all  boundaries.   
    and:   
              f=1   for x=0.5, y=0   
              f=0    else */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    hx = (*xmax - *xmin) / (*nx - 1);
    hy = (*ymax - *ymin) / (*ny - 1);
    i = *location / 2 % *nxloc;
    j = *location / 2 / *nxloc;
    x = *xmin - hx + (*nxst + i - 1) * hx;
    y = *ymin - hy + (*nyst + j - 1) * hy;

/*  Store the real equation coefficients (ai;-bi) */
    f = 0.;
    if (*nyst + j - 1 == 1 && *nxst + i - 1 == (*nx + 1) / 2) 
	f = 1. / (hx * hy);
/*    1   F = 1. */

/*  Set the  boundary conditions equations   
    on the upper boundary:  du/dn -i*Kn*u=0:   
    where du/dn is the outward derivative. */

    if (j + *nyst - 1 == *ny + 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - (*nx + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = *kn * hy;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  The lower boundary (Neumann b.c.  ) */
    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + (*nx + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = *kn * hy;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  The left hand-side boundary (Neumann b.c.) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 1)
	    ) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = *kn * hx;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  The right hand-side boundary (Neumann b.c. ) */
    if (i + *nxst - 1 == *nx + 1 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *
	    ny + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = *kn * hx;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }
/* --------    end of the b.c. equations setting -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 2;
    val[k] = 1. / (hx * hx);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 2;
    val[k] = 1. / (hx * hx);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (*nx + 2 << 1);
    val[k] = 1. / (hy * hy);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (*nx + 2 << 1);
    val[k] = 1. / (hy * hy);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = -2. / (hx * hx) - 2. / (hy * hy) + *kn * *kn;
L100:
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[*kmax - ii] /= sum1;
    }

/* ------------------------------ */
/* ------------------------------ */
/*  Store the imaginary equation coefficients (bi; ai) */

    loca = *location + 1;
    rowa = *row + 1;
    f = 0.;

/*  Set the  b.c. equations: du/dn -i*Kn*u=0: */

/*  The upper boundary: */
    if (j + *nyst - 1 == *ny + 1) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - (*nx + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -(*kn) * hy;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The lower boundary: du/dn -i*Kn*u=0 */

    if (j + *nyst - 1 == 0) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + (*nx + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -(*kn) * hy;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The left hand-side boundary (Neumann b.c.) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 1)
	    ) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -(*kn) * hx;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The right hand-side boundary  (Neumann b.c. ) */
    if (i + *nxst - 1 == *nx + 1 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *
	    ny + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -(*kn) * hx;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }
/* --------    end of the b.c. equation settings -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    icn = 0;
    sum = 0.;
    k = bindx[loca];
    bindx[k] = rowa + 2;
    val[k] = 1. / (hx * hx);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - 2;
    val[k] = 1. / (hx * hx);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (*nx + 2 << 1);
    val[k] = 1. / (hy * hy);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (*nx + 2 << 1);
    val[k] = 1. / (hy * hy);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    *kmax = k;
    bindx[loca + 1] = k;
    val[loca] = -2. / (hx * hx) - 2. / (hy * hy) + *kn * *kn;
L200:
    rhs[rowa] = f;

/*  return */
/*  Normalize the equations: */
    d1 = val[loca];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[rowa] /= sum1;
    val[loca] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[*kmax - ii] /= sum1;
    }

    return 0;
} /* create_matrix_row_5pt_example24_2rows */

/* ************************************************************************ */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_7pt_example25_2rows(int *row, 
	int *location, double *val, int *bindx, int *kmax, 
	double *rhs, int *nxst, int *nyst, int *nzst, int 
	*nxloc, int *nyloc, int *nzloc, double *kn1, double *
	kn2, double *kn3, double *xmin, double *xmax, double *
	ymin, double *ymax, double *zmin, double *zmax, int *
	nx, int *ny, int *nz)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static double f;
    static int i, j, k;
    static double x, y, z;
    static int ii, kk;
    static double kn, hx, hy, hz;
    static int icn;
    static double sum, sum1;
    static int loca;
/*  static double suma, sumb; */
    static int rowa;
    static double ztot;   
    double z13,z23;


/*  Purpose: */
/*  Add one row to an MSR matrix corresponding to a 7pt discrete */
/*  approximation to the 3D operator on an (nx+2)x(ny+2)x(nz+2) domain. */

/*  Example 25 is similar to example #3 of O. Schenk.  2009 */
/*  paper, in SIAM J. Scientific Computing, (3D, 3 layer case). */

/*  Domain: [Xmin-hx  Xmax+hx]*[Ymin-hy Ymax+hy]*[Zmin-hz Zmax+hz]; */

/*          Uxx + Uyy + Uzz + K**2*U = F */

/*          K=K1; K2; K3 for the 3 different layers. */

/*   Neumann b.c. Du/Dn -i*K*U = gi  on all boundaries. */

/*  and:      gi=0 */
/*            f=1   for x=0.5, y=0, z=0.5 */
/*            f=0    else */


/* ----------------------------------------------------------------------- */

/*  Parameters: */
/*     row          == global row number of the new row to be added. */
/*     location     == local row where diagonal of the new row will be stored. */
/*     val,bindx    == (see user's guide). On output, val[] and bindx[] */
/*                     are appended such that the new row has been added. */

    hx = (*xmax - *xmin) / (*nx - 1);
    hy = (*ymax - *ymin) / (*ny - 1);
    hz = (*zmax - *zmin) / (*nz - 1);
    i = *location / 2 % *nxloc;
    j = *location / 2 / *nxloc % *nyloc;
    kk = *location / 2 / (*nxloc * *nyloc) % *nzloc;
    x = *xmin - hx + (*nxst + i - 1) * hx;
    y = *ymin - hy + (*nyst + j - 1) * hy;
    z = *zmin - hz + (*nzst + kk - 1) * hz;
    ztot = *zmax - *zmin;   
    z13= (1./3.)*ztot ;
    z23= (2./3.)*ztot ;
    
    if (z < *zmin + z13) {
	kn = *kn1;
    }
    if (z >= *zmin + z13 && z < *zmin + z23 ) {
	kn = *kn2;
    }
    if (z >= *zmin + z23 ) {
	kn = *kn3;
    }

/*  Store the real equation coefficients (ai;-bi) */
    f = 0.;
/*     for a perturbation at the bottom center: */
    if (*nyst + j - 1 == 1 && *nxst + i - 1 == (*nx + 1) / 2 && *nzst + kk 
	    - 1 == (*nz + 1) / 2) {
	f = 1. / (hx * hy * hz);
    }
/*   1  (NZst+KK-1.eq.(NZ+1)/2)) F = +1. */

/*  Set the  boundary conditions equations */
/*  on the upper boundary:  du/dn -i*Kn*u=0: */
/*  where du/dn is the outward derivative. */

    if (j + *nyst - 1 == *ny + 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - (*nx + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = kn * hy;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  The lower boundary (Neumann b.c.): */

    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + (*nx + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = kn * hy;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  The left hand-side boundary (Neumann b.c.): */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 1)
	     && (kk + *nzst - 1 > 0 && kk + *nzst - 1 < *nz + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = kn * hx;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  The right hand-side boundary  (Neumann b.c.): */
    if (i + *nxst - 1 == *nx + 1 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *
	    ny + 1) && (kk + *nzst - 1 > 0 && kk + *nzst - 1 < *nz + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = kn * hx;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }
/*  The back side boundary (Neumann b.c.): */
    if (kk + *nzst - 1 == 0 && (i + *nxst - 1 >= 0 && i + *nxst - 1 <= *
	    nx + 1) && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row + (*nx + 2 << 1) * (*ny + 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = kn * hz;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  The front side boundary (Neumann b.c.): */
    if (kk + *nzst - 1 == *nz + 1 && (i + *nxst - 1 >= 0 && i + *nxst - 1 
	    <= *nx + 1) && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - (*nx + 2 << 1) * (*ny + 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = *row + 1;
	val[k] = kn * hz;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }
/* --------    end of the b.c. equations setting -------- */

/*  Set the eq. for the inner grid points. */
/*  In the present example there are no imaginary terms */
/*  for the inner grid points. */

    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 2;
    val[k] = 1. / (hx * hx);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 2;
    val[k] = 1. / (hx * hx);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (*nx + 2 << 1);
    val[k] = 1. / (hy * hy);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (*nx + 2 << 1);
    val[k] = 1. / (hy * hy);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (*nx + 2 << 1) * (*ny + 2);
    val[k] = 1. / (hz * hz);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (*nx + 2 << 1) * (*ny + 2);
    val[k] = 1. / (hz * hz);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = -2. / (hx * hx) - 2. / (hy * hy) - 2. / (hz * hz) + kn 
	    * kn;
L100:
    rhs[*row] = f;

/*  return 0;*/
/*  Normalize the equations: */
    d1 = val[*location];
    sum += d1 * d1;
/*  suma = sum; */
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[*kmax - ii] /= sum1;
    }
/*  printf("\n%7d %13.10f %13.10f", *row, rhs[*row], val[*location]);
    for (ii=1; ii<=i1; ++ii) printf(" %13.10f", val[*kmax - ii]) ; */

/* -------------------------------- */
/* -------------------------------- */
/*  Store the imaginary equation coefficients (bi; ai) */

    loca = *location + 1;
    rowa = *row + 1;
    f = 0.;

/*  Set the  b.c. equations: du/dn -i*Kn*u=0: */

/*  The upper boundary: */
    if (j + *nyst - 1 == *ny + 1) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - (*nx + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * hy;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The lower boundary: du/dn -i*Kn*u=0 */
    if (j + *nyst - 1 == 0) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + (*nx + 2 << 1);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * hy;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The left hand-side boundary (Neumann b.c.): */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 1)
	     && (kk + *nzst - 1 > 0 && kk + *nzst - 1 < *nz + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * hx;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The right hand-side boundary  (Neumann b.c.): */
    if (i + *nxst - 1 == *nx + 1 && (j + *nyst - 1 > 0 && j + *nyst - 1 < *
	    ny + 1) && (kk + *nzst - 1 > 0 && kk + *nzst - 1 < *nz + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - 2;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * hx;
/*        val(k) =  -Kn*H */
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }
/*  The back side boundary (Neumann b.c.): */
    if (kk + *nzst - 1 == 0 && (i + *nxst - 1 >= 0 && i + *nxst - 1 <= *
	    nx + 1) && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa + (*nx + 2 << 1) * (*ny + 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * hz;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }


/*  The front side boundary (Neumann b.c.): */
    if (kk + *nzst - 1 == *nz + 1 && (i + *nxst - 1 >= 0 && i + *nxst - 1 
	    <= *nx + 1) && (j + *nyst - 1 > 0 && j + *nyst - 1 < *ny + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - (*nx + 2 << 1) * (*ny + 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 1;
	val[k] = -kn * hz;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }
/* --------    end of the b.c. equation settings -------- */


/*  Set the eq. for the inner grid points. */
/*  In the present example there are no imaginary terms */
/*  for the inner grid points. */

    icn = 0;
    sum = 0.;
    k = bindx[loca];
    bindx[k] = rowa + 2;
    val[k] = 1. / (hx * hx);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - 2;
    val[k] = 1. / (hx * hx);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (*nx + 2 << 1);
    val[k] = 1. / (hy * hy);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (*nx + 2 << 1);
    val[k] = 1. / (hy * hy);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (*nx + 2 << 1) * (*ny + 2);
    val[k] = 1. / (hz * hz);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (*nx + 2 << 1) * (*ny + 2);
    val[k] = 1. / (hz * hz);
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    *kmax = k;
    bindx[loca + 1] = k;
    val[loca] = -2. / (hx * hx) - 2. / (hy * hy) - 2. / (hz * hz) + kn * kn;
L200:
    rhs[rowa] = f;

/*  return 0 ;*/
/*  Normalize the equations: */

    d1 = val[loca];
    sum += d1 * d1;
/*  sumb = sum; */
    sum1 = sqrt(sum);
    rhs[rowa] /= sum1;
    val[loca] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[*kmax - ii] /= sum1;
    }
/*  printf("\n%7d %13.10f %13.10f", rowa, rhs[rowa], val[loca]);
    for (ii=1; ii<=i1; ++ii) printf(" %13.10f", val[*kmax - ii]) ;  */
    return 0;
} /* create_matrix_row_7pt_example25_2rows */


/* ************************************************************************ */
/* ************************************************************************ */

/* Subroutine */ int create_matrix_row_9pt_example26_2rows(int *row, 
	int *location, double *val, int *bindx, int *kmax, 
	double *rhs, int *nxst, int *nyst, int *nxloc, 
	int *nyloc, double *kn1, double *kn2, double *kn3, 
	int *order)
{
    /* System generated locals */
    int i1;
    double d1, d2;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, a0, ac;
    static int ii;
    static double as, pi, kn, del;
    static int icn;
    static double eps, sum, hkn2, sum1, beta;
    static int loca, rowa;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5point/9point discrete   
    approximation to the 2D EXAMPLE_26 operator on an (n+3)*(n+2) rectangle.   
    Note: An external line of points is added to implement a high order scheme   
    of the boundary condition on the right side of the domain.   

    Example 26 = similar to the example  of Y.A. Erlangga,  and   
    E. Turkel 2010? paper.   
    We use a 4th/6th order finite difference scheme following Singer and   
    Turkel 2006, Erlanga and Turkel  2011? and Harrari and Turkel 1995.        
    The 6th order schem uses gamma=14/5 and DEl=-1. / 2./ 0.   

    Domain: [0 1]*[-0.5 0.5];   

            Uxx + Uyy + K**2*U = F   

    where: F=0   

    With Dirichlet boundary condition on the bottom: u=0 on y=-0.5; on   
    the top:  u=0 on y=0.5, and on the left side: u=cos(PI*y) on x=0.   
    On the right side, Neumann boundary condition: Du/Dx +i*beta*U=0 on x=1.   

    The analytic solution of this problem is:   
        u(x,y)=cos(PI*y)*exp(-i*beta*x)    where: PI**2+beta**2=K**2 */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. */

    pi = atan(1.) * 4.;
    h = 1. / (double) (global_2.n + 1);
    i = *location / 2 % *nxloc;
    j = *location / 2 / *nxloc;
    x = (*nxst + i - 1) * h;
    y = (*nyst + j - 1) * h - .5;

    if (y <= 1.f/3.f-0.5f ) kn = *kn1;
    if (y > 1.f/3.f-0.5f && y <= 2.f/3.f-0.5f ) kn = *kn2;
    if (y > 2.f/3.f-0.5f ) kn = *kn3;

    d1 = kn;
    d2 = pi;
    beta = sqrt(d1 * d1 - d2 * d2);
    f = 0.;
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < 
	    global_2.n + 1)) f = cos(pi * y);

/*  Set the  boundary conditions equations */

/*  On the upper boundary (Dirichlet b.c.): u=0 */
    if (j + *nyst - 1 == global_2.n + 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  On the lower boundary (Dirichlet b.c.): u=0 */
    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*  The left hand-side boundary (Dirichlet b.c.): u=cos(y) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < 
	    global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*  The right hand-side boundary  (Neumann b.c.): Du/Dx+i*beta*u=0 */
    if (i + *nxst - 1 == global_2.n + 2 && (j + *nyst - 1 > 0 && j + *nyst 
	    - 1 < global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginay part */
	bindx[k] = *row - 1;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
/* Computing 4th power */
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * -2. * h * (1. - d1 * d1 / 6. + d2 * d2 
		    / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * -2. * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[*location] = 1.;
	goto L100;
    }
/* --------    end of the b.c. equations setting -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points */

    hkn2 = kn * kn * h * h;
    if (*order == 2) {
	a0 = hkn2 - 4.;
	as = 1.;
	ac = 0.;
    }
    if (*order == 4) {
	a0 = hkn2 * 67. / 90. - 10./3.;
	as = hkn2 * 2. / 45. + 2./3.;
	ac = hkn2 * 7. / 360. + 1./6.;
    }
    if (*order == 6) {
/*      del= 2. */
	del = -1.;
/*      del= 0. */
	d1 = hkn2;
	a0 = hkn2 * 67. / 90. - 10./3. + d1 * d1 * (del - 3.) / 180.; 
	d1 = hkn2;
	as = hkn2 * 2. / 45. + 2./3. + d1 * d1 * (3. - del * 2.) / 720.;
	d1 = hkn2;
	ac = hkn2 * 7. / 360. + 1./6. + d1 * d1 * del / 720.;
    }
    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 3 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L50;
    }

    bindx[k] = *row + (global_2.n + 3 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 3 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L50:

    bindx[k] = *row - (global_2.n + 3 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L55;
    }
    bindx[k] = *row - (global_2.n + 3 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (global_2.n + 3 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L55:

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = a0;
L100:
    rhs[*row] = f;

/*  return */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[k - ii] /= sum1;
    }

/* --------------- */
/* --------------- */
/*  Store the imaginary equation coefficients (B; A) */

    loca = *location + 1;
    rowa = *row + 1;
    f = 0.;

/*  Set the  b.c. equations: */

/*  The upper boundary u=0: */
    if (j + *nyst - 1 == global_2.n + 1) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The lower boundary (Dirichlet b.c.) u=0: */
    if (j + *nyst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	*kmax = k;
	bindx[loca + 1] = k;
	val[loca] = 1.;
	goto L200;
    }

/*  The left hand-side boundary (Dirichlet b.c.): u=cos(PI*y) */
    if (i + *nxst - 1 == 0 && (j + *nyst - 1 > 0 && j + *nyst - 1 < 
	    global_2.n + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The right hand-side boundary  (Neumann b.c. ) */
    if (i + *nxst - 1 == global_2.n + 2 && (j + *nyst - 1 > 0 && j + *nyst 
	    - 1 < global_2.n + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - 4;
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - 3;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
/* Computing 4th power */
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * 2. * h * (1. - d1 * d1 / 6. + d2 * d2 / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2 * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[loca] = 1.;
	goto L200;
    }
/* --------    end of the b.c. equation setting -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    icn = 0;
    sum = 0.;
    k = bindx[loca];
    bindx[k] = rowa + 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (global_2.n + 3 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L150;
    }

    bindx[k] = rowa + (global_2.n + 3 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (global_2.n + 3 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L150:

    bindx[k] = rowa - (global_2.n + 3 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L155;
    }

    bindx[k] = rowa - (global_2.n + 3 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (global_2.n + 3 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L155:

    *kmax = k;
    bindx[loca + 1] = k;
    val[loca] = a0;
L200:
    rhs[rowa] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[loca];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[rowa] /= sum1;
    val[loca] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[k - ii] /= sum1;
    }

    return 0;
} /* create_matrix_row_9pt_example26_2rows */

/* ************************************************************************ */
/* ******************************************************************** */

/* Subroutine */ int example26_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr1, double *sumerr2, 
	double *kn1, double *kn2, double *kn3, int *order, 
        int *ipr)
{
    /* System generated locals */
    int i1, i2;
    double d1, d2;

    /* Local variables */
    static double h;
    static int i, j;
    static double x, y;
    static int ii;
    static double kn, pi, beta, u_imag, u_real;

/*  Purpose: calculate the error for a 2D Helmholtz problem.   
    Example26 is similar to the example of Erlangga and E. Turkel,   
    in the 2011? paper.   
    The equations are solved using a 2nd/4th/6th order difference scheme.   

    Domain: [0 1]*[-1/2 1/2]   

          Uxx + Uyy + K**2*U = f   

    Where: f=0   
    With Dirichlet boundary conditions:   
    On the bottom: U=0 on y=-1/2,  On the top:  U=0 on y=1/2;   
    On the left side: U=cos(PI*y) on x=0.   
    and Neumann b.c. on the right side: DU/Dx +i*beta*U=0 on x=1.,   
    Where: beta**2+PI**2=K**2   

    The analytic solution is: U=cos(PI*y)*exp(-i*beta*x)   
                               =cos(PI*y)*[cos(-beta*x)+ i*sin(-beta*x)] */
/* ----------------------------------------------------------------------- */

    if(*ipr == 0) goto L80;
    printf("\n  %s", "    ");
    printf("\n  %s", "    ");
    printf("\n  %s","  Results for example 26:");
    printf("\n  %s","  Similar to the example in Erlangga and Turkel, 2011 paper.");
    printf("\n  %s", "    ");
    if (*order == 2) {
       printf("\n  %s","  Using a 2nd order finite difference scheme.     ");
    }
    if (*order == 4) {
       printf("\n  %s","  Using a 4th order finite difference scheme.     ");
    }
    if (*order == 6) {
       printf("\n  %s","  Using a 6th order finite difference scheme.     ");
    }
    printf("\n  %s", "    ");
    printf("\n  %s","  Domain: [0 1]*[-1/2 1/2]  ");
    printf("\n  %s", "    ");
    printf("\n  %s","  Uxx+Uyy+K**2*U= f  ");
    printf("\n  %s", "    ");
    printf("\n  %s","  Where:  f=0.           ");
    printf("\n  %s", "    ");
    printf("\n  %s","  Dirichlet boundary conditions:      ");
    printf("\n  %s","     U=0 on y=-1/2 ; U=0 on  y=1/2; U=cos(PI*y) on  x=0");
    printf("\n  %s","  Neumann b.c.:  DU/Dx +i*beta*U=0 on x=1.;           ");
    printf("\n  %s","  Where: beta**2+PI**2= K**2 "                         );
    printf("\n  %s","   " );
    printf("\n  %s","   The analytic solution is: U=cos(PI*y)*exp(-i*beta*x) " );
    printf("\n  %s","------------------------------------------------------------------" );
    printf("\n  %s","   ");
    printf("\n  %s","   ");
L80:

    pi = atan(1.) * 4.;
    h = 1. / (double) (*nx + 1);
    *sumerr1 = 0.;
    *sumerr2 = 0.;
    ii = 0;
    i1 = *ny + 1;
    for (j = 0; j <= i1; ++j) {
	y = (double) j * h - .5;

	if (y <= (1./3.)-0.5 ) kn = *kn1;
	if (y > (1./3.)-0.5 && y <= (2./3.)-0.5 ) kn = *kn2;
	if (y > (2./3.)-0.5 ) kn = *kn3;
	d1 = kn;
	d2 = pi;
	beta = sqrt(d1 * d1 - d2 * d2);
	i2 = *nx + 2;
	for (i = 0; i <= i2; ++i) {
	    x = (double) i * h;
	    u_real = cos(pi * y) * cos(-beta * x);
	    u_imag = cos(pi * y) * sin(-beta * x);
	    d1 = u_real;
	    *sumerr1 += d1 * d1;
	    d1 = u_imag;
	    *sumerr2 += d1 * d1;
	    err[ii] = xnew[ii] - u_real;
	    ++ii;
	    err[ii] = xnew[ii] - u_imag;
	    ++ii;
	}
    }
    *sumerr1 = sqrt(*sumerr1);
    *sumerr2 = sqrt(*sumerr2);
    if(*ipr == 1)
    printf("\n  SUMERR1= %f    SUMERR2= %f", *sumerr1, *sumerr2);
    return 0;
} /* example26_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int create_matrix_row_9pt_example27_2rows(int *row, 
	int *location, double *val, int *bindx, int *kmax, 
	double *rhs, int *nxst, int *nyst, int *nxloc, 
	int *nyloc, double *kn1, double *kn2, double *kn3, 
	int *order)
{
    /* System generated locals */
    int i1;
    double d1, d2;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, a0, ac;
    static int ii;
    static double as, pi, kn, del;
    static int icn;
    static double sum, hkn2, sum1, beta;
    static int loca, rowa;


/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5point/9point discrete   
    approximation to the 2D EXAMPLE_27 operator on an (n+2)*(n+3) rectangle.   
    Note: An external line of points is added to implement a high order scheme   
    of the boundary condition on the right side of the domain.   

    Example 27 = similar to the example  of Y.A. Erlangga,  and   
    E. Turkel 2011? paper.   
    The variable Order defines the scheme order: 2nd / 4th/ 6th order scheme.   
    The 4th/6th order finite difference schemes follow Singer and   
    Turkel 2006, Erlanga and Turkel  2011? and Harrari and Turkel 1995.   
    The 6th order schem uses gamma=14/5 and DEl=-1. / 2./ 0.   

    Domain: [-0.5 0.5]*[0 1];   

           Uxx + Uyy + K**2*U = F   

    Where: F=0   

    With Dirichlet boundary condition on the left side: u=0 on x=-0.5; on   
    the right side:  u=0 on x=0.5, and on the bottom: u=cos(PI*x) on y=0.   
    On the top Neumann boundary condition: Du/Dy +i*beta*U=0 on y=1.   

    The analytic solution of this problem is:   
           u(x,y)=cos(PI*x)*exp(-i*beta*y)    where: PI**2+beta**2=K**2 */
/* ----------------------------------------------------------------------- */
/*  Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added. 
       Order        == scheme order */

    pi = atan(1.) * 4.;
    h = 1. / (double) (global_2.n + 1);
    i = *location / 2 % *nxloc;
    j = *location / 2 / *nxloc;
    x = (*nxst + i - 1) * h - .5;
    y = (*nyst + j - 1) * h;

    if (x <= 1.f/3.f-0.5f) kn = *kn1;
    if (x > 1.f/3.f-0.5 && y <= 2.f/3.f-0.5f) kn = *kn2;
    if (x > 2.f/3.f-0.5f ) kn = *kn3;

    d1 = kn;
    d2 = pi;
    beta = sqrt(d1 * d1 - d2 * d2);

    f = 0.;
    if (j + *nyst - 1 == 0 && (i + *nxst - 1 > 0 && i + *nxst - 1 < 
	    global_2.n + 1)) f = cos(pi * x);

/*  Set the  boundary conditions equations */

/*  On the right side boundary (Dirichlet b.c.): u=0 */

    if (i + *nxst - 1 == global_2.n + 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	if (*order == 2) {
	    goto L35;
	}
	if ((*order == 4 || *order == 6) && (j + *nyst - 1 == 0 || j + *nyst 
		- 1 == global_2.n + 2)) {
	    goto L35;
	}
	bindx[k] = *row + (global_2.n + 2 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = *row - (global_2.n + 2 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L35:
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  On the left side boundary (Dirichlet b.c.): u=0 */

    if (i + *nxst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	if (*order == 2) {
	    goto L45;
	}
	if ((*order == 4 || *order == 6) && (j + *nyst - 1 == 0 || j + *nyst 
		- 1 == global_2.n + 2)) {
	    goto L45;
	}
	bindx[k] = *row + (global_2.n + 2 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = *row - (global_2.n + 2 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L45:
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*  The  bottom boundary (Dirichlet b.c.): u=cos(Pi*x) */
    if (j + *nyst - 1 == 0 && (i + *nxst - 1 > 0 && i + *nxst - 1 < 
	    global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*  The top boundary  (Neumann b.c.): Du/Dy+i*beta*u=0 */
    if (j + *nyst - 1 == global_2.n + 2 && (i + *nxst - 1 > 0 && i + *
	    nxst - 1 < global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	bindx[k] = *row - (global_2.n + 2 << 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginay part */
	bindx[k] = *row - (global_2.n + 2 << 1) + 1;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
/* Computing 4th power */
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * -2. * h * (1. - d1 * d1 / 6. + d2 * d2 / 120.); 
	}
	if (*order == 2) {
	    val[k] = beta * -2. * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[*location] = 1.;
	goto L100;
    }
/* --------    end of the b.c. equations setting -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    hkn2 = kn * kn * h * h;
    if (*order == 2) {
	a0 = hkn2 - 4.;
	as = 1.;
	ac = 0.;
    }
    if (*order == 4) {
	a0 = hkn2 * 67. / 90. - 10./3.;
	as = hkn2 * 2. / 45. + 2./3.;
	ac = hkn2 * 7. / 360. + 1./6.;
    }
    if (*order == 6) {
/*      del= 2. */
	del = -1.;
/*      del= 0. */
	d1 = hkn2;
	a0 = hkn2 * 67. / 90. - 10./3. + d1 * d1 * (del - 3.) / 180.;
	d1 = hkn2;
	as = hkn2 * 2. / 45. + 2./3. + d1 * d1 * (3. - del *  2.) / 720.;
	d1 = hkn2;
	ac = hkn2 * 7. / 360. + 1./6. + d1 * d1 * del / 720.;
    }
    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 2 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L50;
    }

    bindx[k] = *row + (global_2.n + 2 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 2 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L50:

    bindx[k] = *row - (global_2.n + 2 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L55;
    }
    bindx[k] = *row - (global_2.n + 2 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (global_2.n + 2 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L55:

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = a0;
L100:
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[k - ii] /= sum1;
    }

/* ------------------------------------ */
/* ------------------------------------ */
/*  Store the imaginary equation coefficients (B; A) */

    loca = *location + 1;
    rowa = *row + 1;
    f = 0.;

/*  Set the  b.c. equations: */

/*  The right side boundary u=0: */

    if (i + *nxst - 1 == global_2.n + 1) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	if (*order == 2) {
	    goto L56;
	}
	if ((*order == 4 || *order == 6) && (j + *nyst - 1 == 0 || j + *nyst 
		- 1 == global_2.n + 2)) {
	    goto L56;
	}
	bindx[k] = rowa + (global_2.n + 2 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa - (global_2.n + 2 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L56:
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The left side boundary (Dirichlet b.c.) u=0: */

    if (i + *nxst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	if (*order == 2) {
	    goto L58;
	}
	if ((*order == 4 || *order == 6) && (j + *nyst - 1 == 0 || j + *nyst 
		- 1 == global_2.n + 2)) {
	    goto L58;
	}
	bindx[k] = rowa + (global_2.n + 2 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = rowa - (global_2.n + 2 << 1);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L58:
	*kmax = k;
	bindx[loca + 1] = k;
	val[loca] = 1.;
	goto L200;
    }

/*  The bottom boundary (Dirichlet b.c.): u=cos(Pi*x) */
    if (j + *nyst - 1 == 0 && (i + *nxst - 1 > 0 && i + *nxst - 1 < 
	    global_2.n + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
	val[loca] = 1.;
	goto L200;
    }

/*  The top boundary  (Neumann b.c. ): */
    if (j + *nyst - 1 == global_2.n + 2 && (i + *nxst - 1 > 0 && i + *
	    nxst - 1 < global_2.n + 1)) {
/*  the real part: */
	icn = 0;
	sum = 0.;
	k = bindx[loca];
	bindx[k] = rowa - (global_2.n + 2 << 2);
	val[k] = -1.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
/*  the imaginary part: */
	bindx[k] = rowa - (global_2.n + 2 << 1) - 1;
	if (*order == 4 || *order == 6) {
	    d1 = beta * h;
	    d2 = beta * h, d2 *= d2;
	    val[k] = beta * 2. * h * (1. - d1 * d1 / 6. + d2 * d2 / 120.);
	}
	if (*order == 2) {
	    val[k] = beta * 2. * h;
	}
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
	*kmax = k;
	bindx[loca + 1] = k;
/*  the diagonal element: */
/*  the additional point - outside the computational domain: */
	val[loca] = 1.;
	goto L200;
    }
/* --------    end of the b.c. equation setting -------- */


/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points */

    icn = 0;
    sum = 0.;
    k = bindx[loca];
    bindx[k] = rowa + 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - 2;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (global_2.n + 2 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L150;
    }

    bindx[k] = rowa + (global_2.n + 2 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa + (global_2.n + 2 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L150:

    bindx[k] = rowa - (global_2.n + 2 << 1);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L155;
    }

    bindx[k] = rowa - (global_2.n + 2 << 1) + 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = rowa - (global_2.n + 2 << 1) - 2;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L155:

    *kmax = k;
    bindx[loca + 1] = k;
    val[loca] = a0;
L200:
    rhs[rowa] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[loca];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[rowa] /= sum1;
    val[loca] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[k - ii] /= sum1;
    }

    return 0;
} /* create_matrix_row_9pt_example27_2rows */

/* ************************************************************************ */
/* ******************************************************************** */

/* Subroutine */ int example27_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr1, double *sumerr2, 
	double *kn1, double *kn2, double *kn3, int *order,
         int *ipr )
{
    /* System generated locals */
    int i1, i2;
    double d1, d2;

    /* Local variables */
    static double h;
    static int i, j;
    static double x, y;
    static int ii;
    static double kn, pi, beta, u_imag, u_real;

/*  Purpose: calculate the error for a 2D Helmholtz problem.   
    Example27 is  example of Erlangga and Turkel,   
    2011 paper.   
    The equations are solved using a 2nd/4th/6th order difference scheme.   

    Domain: [-1/2 1/2]*[0 1]   

           Uxx + Uyy + K**2*U = f   

    Where: f=0   

    With Dirichlet boundary conditions:   
    On the left side: u=0 on x=-1/2,  On the right side:  u=0 on x=1/2;   
    On the bottom   : u=cos(PI*x) on y=0.   
    Neumann b.c. on the top: Du/Dy +i*beta*U=0 on y=1.,   
    where: beta**2+PI**2=K**2   

    The analytic solution is: u=cos(PI*x)*exp(-i*beta*y)   
                               =cos(PI*x)*[cos(-beta*y)+ i*sin(-beta*y)] */
/* ----------------------------------------------------------------------- */

    if(*ipr == 0) goto L80;
    printf("\n  %s", "    ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Results for example 27:");
    printf("\n  %s","   Similar to the example in Erlangga and Turkel, 2011 paper.");
    printf("\n  %s", "    ");
    if (*order == 2) {
    printf("\n  %s","   Using a 2nd order finite difference scheme.     ");
    }
    if (*order == 4) {
    printf("\n  %s","   Using a 4th order finite difference scheme.     ");
    }
    if (*order == 6) {
    printf("\n  %s","   Using a 6th order finite difference scheme.     ");
    }
    printf("\n  %s", "    ");
    printf("\n  %s","   Domain: [-1/2 1/2]*[0 1]  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Uxx+Uyy+K**2*U= f  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Where:  f=0.           ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Dirichlet boundary conditions:      ");
    printf("\n  %s","      U=0 on x=-1/2 ; U=0 on  x=1/2; U=cos(PI*x) on  y=0");
    printf("\n  %s","   Neumann b.c.: DU/Dy +i*beta*U=0 on y=1.;           ");
    printf("\n  %s","   Where: beta**2+PI**2= K**2 "                         );
    printf("\n  %s","   " );                                                
    printf("\n  %s","   The analytic solution is: U=cos(PI*x)*exp(-i*beta*y) "    );
    printf("\n  %s","---------------------------------------------------------------" );
    printf("\n  %s","   ");                                                
    printf("\n  %s","   ");                                                
L80:
    pi = atan(1.) * 4.;
    h = 1. / (double) (*nx + 1);
    *sumerr1 = 0.;
    *sumerr2 = 0.;
    ii = 0;
    i1 = *ny + 2;
    for (j = 0; j <= i1; ++j) {
	y = (double) j * h;

	if (y <= 1./3.) kn = *kn1;
	if (y > 1./3. && y <= 2./3.) kn = *kn2;
	if (y > 2./3.) kn = *kn3;
	 
	d1 = kn;
	d2 = pi;
	beta = sqrt(d1 * d1 - d2 * d2);
	i2 = *nx + 1;
	for (i = 0; i <= i2; ++i) {
	    x = (double) i * h - .5;
	    u_real = cos(pi * x) * cos(-beta * y);
	    u_imag = cos(pi * x) * sin(-beta * y);
	    d1 = u_real;
	    *sumerr1 += d1 * d1;
	    d1 = u_imag;
	    *sumerr2 += d1 * d1;
	    err[ii] = xnew[ii] - u_real;
	    ++ii;
	    err[ii] = xnew[ii] - u_imag;
	    ++ii;
	}
    }
    *sumerr1 = sqrt(*sumerr1);
    *sumerr2 = sqrt(*sumerr2);
    if(*ipr == 1)
    printf("\n  SUMERR1= %f    SUMERR2= %f", *sumerr1, *sumerr2);
    return 0;
} /* example27_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int create_matrix_row_9pt_example28(int *row, int *
	location, double *val, int *bindx, int *kmax, double *
	rhs, int *nxst, int *nyst, int *nxloc, int *nyloc, 
	double *kn, int *order)
{
    /* System generated locals */
    int i1;
    double d1, d2;

    /* Local variables */
    static double f, h;
    static int i, j, k;
    static double x, y, a0, ac;
    static int ii;
    static double as, pi, del;
    static int icn;
    static double sum, hkn2, sum1, beta;

/*  Purpose:   
    Add one row to an MSR matrix corresponding to a 5point/9point discrete   
    approximation to the 2D EXAMPLE_28 operator on an (n+2)*(n+2) rectangle.   
    Note: An external line of points is added to implement a high order scheme   
    of the boundary condition on the right side of the domain.   

    Example 28 - is our example   

    The variable Order defines the scheme order: 2nd / 4th/ 6th order scheme.   
    The 4th/6th order finite difference schemes follow Singer and   
    Turkel 2006, Erlanga and Turkel  2011? and Harrari and Turkel 1995.        
    The 6th order scheme uses gamma=14/5, and DEl=-1. / 2./ 0.   

    Domain: [0 1]*[0 1];   

            Uxx + Uyy + K**2*U = F   

    Where: F=0   
    Case A:   
    With Dirichlet boundary condition on the left side: u=0 on X=0.; on   
    The right side:  u=0 on x=1.0, and on the bottom: u=0. on y=0.   
    On the top: u=sin(PI*X)*sin(beta) on y=1.   

    The analytic solution is:   
          u(x,y)= sin(PI*X)*sin(beta*Y)  with: PI**2+beta**2=K**2 */
/* ----------------------------------------------------------------------- */
/*   Parameters:   
       row          == global row number of the new row to be added.   
       location     == local row where diagonal of the new row will be stored.   
       val,bindx    == (see user's guide). On output, val[] and bindx[]   
                       are appended such that the new row has been added.   
       Order        == scheme order */


    pi = atan(1.) * 4.;
    h = 1. / (double) (global_2.n + 1);
    i = *location % *nxloc;
    j = *location / *nxloc % *nyloc;
    x = (*nxst + i - 1) * h;
    y = (*nyst + j - 1) * h;
    d1 = *kn;
    d2 = pi;
    beta = sqrt(d1 * d1 - d2 * d2);

    f = 0.;
/*  CASE A: */
/* --------- */
/*   if (J+NYst-1.eq.0.and.(I+NXst-1.gt.0.and.I+NXst-1.lt.n+1)) */
/*  1    F=0.D0 */
/*   if (J+NYst-1.eq.n+1.and.(I+NXst-1.gt.0.and.I+NXst-1.lt.n+1)) */
/*  1    F=dsin(PI*X)*dsin(beta) */
/*   if (I+NXst-1.eq.n+1) F=0.D0 */
/*   if (I+NXst-1.eq.0) F=0.D0 */
/*  CASE B: */
/* --------- */
/*   const=Kn/dsqrt(2.D0) */
/*   if (J+NYst-1.eq.0.and.(I+NXst-1.gt.0.and.I+NXst-1.lt.n+1)) */
/*  1    F=dsin(const*X) */
/*   if (J+NYst-1.eq.n+1.and.(I+NXst-1.gt.0.and.I+NXst-1.lt.n+1)) */
/*  1    F=dsin(const*X)*dcos(const) */
/*   if (I+NXst-1.eq.n+1) F=0.D0 */
/*   if (I+NXst-1.eq.0) F=dsin(const)*dcos(const*Y) */
/*  CASE C: */
/* --------- */
    if (j + *nyst - 1 == 0 && (i + *nxst - 1 > 0 && i + *nxst - 1 < 
	    global_2.n + 1)) f = sin(pi * x);

    if (j + *nyst - 1 == global_2.n + 1 && (i + *nxst - 1 > 0 && i + *
	    nxst - 1 < global_2.n + 1)) f = sin(pi * x) * cos(beta);

    if (i + *nxst - 1 == global_2.n + 1) f = 0.;
    if (i + *nxst - 1 == 0) f = 0.;
/* --------- */

/*  Set the  boundary conditions equations */

/*  On the right side boundary (Dirichlet b.c.): u=0 */

    if (i + *nxst - 1 == global_2.n + 1) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	if (*order == 2) {
	    goto L35;
	}
	if ((*order == 4 || *order == 6) && (j + *nyst - 1 == 0 || j + *nyst 
		- 1 == global_2.n + 1)) {
	    goto L35;
	}
	bindx[k] = *row + (global_2.n + 2);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = *row - (global_2.n + 2);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L35:
	*kmax = k;
	bindx[*location + 1] = k;
/*  the diagonal element: */
	val[*location] = 1.;
	goto L100;
    }

/*  On the left side boundary (Dirichlet b.c.): u=0 */

    if (i + *nxst - 1 == 0) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];

	if (*order == 2) {
	    goto L45;
	}
	if ((*order == 4 || *order == 6) && (j + *nyst - 1 == 0 || j + *nyst 
		- 1 == global_2.n + 1)) {
	    goto L45;
	}
	bindx[k] = *row + (global_2.n + 2);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;

	bindx[k] = *row - (global_2.n + 2);
	val[k] = 0.;
	d1 = val[k];
	sum += d1 * d1;
	++icn;
	++k;
L45:
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*  The bottom boundary (Dirichlet b.c.): */
    if (j + *nyst - 1 == 0 && (i + *nxst - 1 > 0 && i + *nxst - 1 < 
	    global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }

/*  The top boundary: */
    if (j + *nyst - 1 == global_2.n + 1 && (i + *nxst - 1 > 0 && i + *
	    nxst - 1 < global_2.n + 1)) {
	icn = 0;
	sum = 0.;
	k = bindx[*location];
	*kmax = k;
	bindx[*location + 1] = k;
	val[*location] = 1.;
	goto L100;
    }
/* --------    end of the b.c. equations setting -------- */

/*  Set the eq. for the inner grid points.   
    In the present example there are no imaginary terms   
    for the inner grid points. */

    hkn2 = *kn * *kn * h * h;
    if (*order == 2) {
	a0 = hkn2 - 4.;
	as = 1.;
	ac = 0.;
    }
    if (*order == 4) {
	a0 = hkn2 * 67. / 90. - 10./3.;
	as = hkn2 * 2. / 45. + 2./3. ;
	ac = hkn2 * 7. / 360. + 1./6.;
    }
    if (*order == 6) {
/*      del= 2. */
	del = -1.;
/*      del= 0. */
	d1 = hkn2;
	a0 = hkn2 * 67. / 90. - 10./3.  + d1 * d1 * (del - 3.) 
		/ 180.;
	d1 = hkn2;
	as = hkn2 * 2. / 45. + 2./3.  + d1 * d1 * (3. - del * 
		2.) / 720.;
	d1 = hkn2;
	ac = hkn2 * 7. / 360. + 1./6.  + d1 * d1 * del / 720.;
    }
    icn = 0;
    sum = 0.;
    k = bindx[*location];
    bindx[k] = *row + 1;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - 1;
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 2);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L50;
    }

    bindx[k] = *row + (global_2.n + 2) + 1;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row + (global_2.n + 2) - 1;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L50:

    bindx[k] = *row - (global_2.n + 2);
    val[k] = as;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    if (*order == 2) {
	goto L55;
    }
    bindx[k] = *row - (global_2.n + 2) + 1;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

    bindx[k] = *row - (global_2.n + 2) - 1;
    val[k] = ac;
    d1 = val[k];
    sum += d1 * d1;
    ++icn;
    ++k;

L55:

    *kmax = k;
    bindx[*location + 1] = k;
    val[*location] = a0;
L100:
    rhs[*row] = f;

/*  return 0; */
/*  Normalize the equations: */

    d1 = val[*location];
    sum += d1 * d1;
    sum1 = sqrt(sum);
    rhs[*row] /= sum1;
    val[*location] /= sum1;
    i1 = icn;
    for (ii = 1; ii <= i1; ++ii) {
	val[k - ii] /= sum1;
    }

    return 0;
} /* create_matrix_row_9pt_example28 */

/* ************************************************************************ */
/* ******************************************************************** */

/* Subroutine */ int example28_err(double *xnew, double *err, 
	int *nx, int *ny, double *sumerr, double *kn, int 
	*order, int *ipr)
{
    /* System generated locals */
    int i1, i2;
    double d1, d2;

    /* Local variables */
    static double h;
    static int i, j;
    static double u, x, y;
    static int ii;
    static double pi, beta;

/*  Purpose: calculate the error for a 2D Helmholtz problem.   
    Example28 - is our example   
    The equations are solved using a 2nd/4th/6th order difference scheme.   

    Domain: [0  1]*[0 1]   

          Uxx + Uyy + K**2*U = f   

    Where: f=0   

    Case A:   
    With Dirichlet boundary conditions: On the left side: U=0 on X=0.; On   
    the right side:  U=0 on x=1.0;  On the bottom: U=0. on y=0.   
    On the top: U=sin(PI*X)*sin(beta) on y=1.   

    The analytic solution is:   
            U(x,y)= sin(PI*X)*sin(beta*Y)  with: PI**2+beta**2=K**2 */
/* ---------------------------------------------------------------------- */

    if(*ipr == 0) goto L80;
    printf("\n  %s", "    ");
    printf("\n  %s", "    ");
/*  printf("\n  %s","   Results for example 28, case A:  "); */
/*  printf("\n  %s","   Results for example 28, case B:  "); */
    printf("\n  %s","   Results for example 28, case C:  ");
    printf("\n  %s", "    ");
    if (*order == 2) {
       printf("\n  %s","   Using a 2nd order finite difference scheme.     ");
    }
    if (*order == 4) {
       printf("\n  %s","   Using a 4th order finite difference scheme.     ");
    }
    if (*order == 6) {
       printf("\n  %s","   Using a 6th order finite difference scheme.     ");
    }
    printf("\n  %s", "    ");
    printf("\n  %s","   Domain: [-1/2 1/2]*[0 1]  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Uxx+Uyy+K**2*U= f  ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Where:  f=0.           ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Dirichlet boundary condition on all boundaries:      ");
    printf("\n  %s", "    ");
    printf("\n  %s","   Case C: Analytic solution: U=sin(PI*x)*cos(beta*y)");
    printf("\n  %s","   Where: PI**2+beta**2=K**2                             ");
    printf("\n  %s","-------------------------------------------------------" );
    printf("\n  %s", "    ");
L80:
    
    pi = atan(1.) * 4.;
    h = 1. / (double) (*nx + 1);
    d1 = *kn;
    d2 = pi;
    beta = sqrt(d1 * d1 - d2 * d2);
    *sumerr = 0.;
    ii = 0;
    i1 = *ny + 1;
    for (j = 0; j <= i1; ++j) {
	y = (double) j * h;
	i2 = *nx + 1;
	for (i = 0; i <= i2; ++i) {
	    x = (double) i * h;
/*  case A: */
/* --------- */
/*          U=dsin(PI*X)*dsin(beta*Y) */
/*  case B: */
/* --------- */
/*          const=Kn/dsqrt(2.D0) */
/*          U=dsin(const*X)*dcos(const*Y) */
/*  case C: */
/* --------- */
            u = sin(pi * x) * cos(beta * y);
	    d1 = u;
	    *sumerr += d1 * d1;
	    err[ii] = xnew[ii] - u;
/*          write(6,100) ii,x,y,u, xnew(ii),err(ii),sumerr,beta */
	    ++ii;
	}
    }
    *sumerr = sqrt(*sumerr);
    if(*ipr == 1)
    printf("\n  SUMERR= %f", *sumerr);
    return 0;
} /* example28_err */

/* ******************************************************************** */
/* ******************************************************************** */

/* Subroutine */ int transf_3d(double *q, double *ecoef, int *
	update_inv, int *n_update, int *comm_cart, int *me, 
	int *rank_source_x, int *rank_source_y, int *
	rank_source_z, int *rank_dest_x, int *rank_dest_y, 
	int *rank_dest_z, int *index_s_xm, int *index_s_xp, 
	int *index_s_ym, int *index_s_yp, int *index_s_zm, 
	int *index_s_zp, int *index_r_xm, int *index_r_xp, 
	int *index_r_ym, int *index_r_yp, int *index_r_zm, 
	int *index_r_zp, int *ixp, int *ixm, int *iyp, 
	int *iym, int *izp, int *izm, int *ndim)
{
    /* System generated locals */
    int i1;

    /* Local variables */

    static int tagx_min, tagy_min, tagz_min, i, tagx_plus, 
	    tagy_plus, tagz_plus;
    MPI_Status status;

    /* Parameter adjustments */
    --index_r_zp;
    --index_r_zm;
    --index_r_yp;
    --index_r_ym;
    --index_r_xp;
    --index_r_xm;
    --index_s_zp;
    --index_s_zm;
    --index_s_yp;
    --index_s_ym;
    --index_s_xp;
    --index_s_xm;   

    /* Function Body */
    tagx_min = 10;
    tagx_plus = 20;
    tagy_min = 30;
    tagy_plus = 40;
    tagz_min = 50;
    tagz_plus = 60;

/*   Define the vectors to be transfered to the neighboring processors */
    i1 = *ixm;
    for (i = 1; i <= i1; ++i) {
	buf_s_xm[i - 1] = q[index_s_xm[i] - 1];
    }

    i1 = *ixp;
    for (i = 1; i <= i1; ++i) {
	buf_s_xp[i - 1] = q[index_s_xp[i] - 1];
    }

    if (*ndim == 1) {
	goto L10;
    }
    i1 = *iym;
    for (i = 1; i <= i1; ++i) {
	buf_s_ym[i - 1] = q[index_s_ym[i] - 1];
    }

    i1 = *iyp;
    for (i = 1; i <= i1; ++i) {
	buf_s_yp[i - 1] = q[index_s_yp[i] - 1];
    }

    if (*ndim == 2) {
	goto L10;
    }
    i1 = *izm;
    for (i = 1; i <= i1; ++i) {
	buf_s_zm[i - 1] = q[index_s_zm[i] - 1];
    }

    i1 = *izp;
    for (i = 1; i <= i1; ++i) {
	buf_s_zp[i - 1] = q[index_s_zp[i] - 1];
    }
L10:

/* send the correction to the neighbors */

    if (*rank_source_x > -1) {
	MPI_Sendrecv(buf_s_xm,*ixm,MPI_DOUBLE_PRECISION,*rank_source_x,tagx_min, 
		buf_r_xm,*ixm,MPI_DOUBLE_PRECISION,*rank_source_x,tagx_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_x > -1) {
	MPI_Sendrecv(buf_s_xp,*ixp,MPI_DOUBLE_PRECISION,*rank_dest_x,tagx_plus, 
		buf_r_xp,*ixp,MPI_DOUBLE_PRECISION,*rank_dest_x,tagx_min, 
		*comm_cart,&status );
    }

    if (*ndim == 1) {
	goto L20;
    }
    if (*rank_source_y > -1) {
	MPI_Sendrecv(buf_s_ym,*iym,MPI_DOUBLE_PRECISION,*rank_source_y,tagy_min, 
		buf_r_ym,*iym,MPI_DOUBLE_PRECISION,*rank_source_y,tagy_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_y > -1) {
	MPI_Sendrecv(buf_s_yp,*iyp,MPI_DOUBLE_PRECISION,*rank_dest_y,tagy_plus, 
		buf_r_yp,*iyp,MPI_DOUBLE_PRECISION,*rank_dest_y,tagy_min, 
		*comm_cart,&status );
    }

    if (*ndim == 2) {
	goto L20;
    }
    if (*rank_source_z > -1) {
	MPI_Sendrecv(buf_s_zm,*izm,MPI_DOUBLE_PRECISION,*rank_source_z,tagz_min, 
		buf_r_zm,*izm,MPI_DOUBLE_PRECISION,*rank_source_z,tagz_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_z > -1) {
	MPI_Sendrecv(buf_s_zp,*izp,MPI_DOUBLE_PRECISION,*rank_dest_z,tagz_plus, 
		buf_r_zp,*izp,MPI_DOUBLE_PRECISION,*rank_dest_z,tagz_min, 
		*comm_cart,&status );
    }
L20:

/*  Update qnewl() of the processor's inner points due to the   
    neighboring processors */

    i1 = *ixm;
    for (i = 1; i <= i1; ++i) {
	q[index_r_xm[i] - 1] += buf_r_xm[i - 1];
    }

    i1 = *ixp;
    for (i = 1; i <= i1; ++i) {
	q[index_r_xp[i] - 1] += buf_r_xp[i - 1];
    }

    if (*ndim == 1) {
	goto L30;
    }
    i1 = *iym;
    for (i = 1; i <= i1; ++i) {
	q[index_r_ym[i] - 1] += buf_r_ym[i - 1];
    }

    i1 = *iyp;
    for (i = 1; i <= i1; ++i) {
	q[index_r_yp[i] - 1] += buf_r_yp[i - 1];
    }

    if (*ndim == 2) {
	goto L30;
    }
    i1 = *izm;
    for (i = 1; i <= i1; ++i) {
	q[index_r_zm[i] - 1] += buf_r_zm[i - 1];
    }

    i1 = *izp;
    for (i = 1; i <= i1; ++i) {
	q[index_r_zp[i] - 1] += buf_r_zp[i - 1];
    }
L30:

/*  Get  the new values of the inner boundary points: */

    i1 = *n_update - 1;
    for (i = 0; i <= i1; ++i) {
	q[update_inv[i]] *= ecoef[update_inv[i]];
    }

    i1 = *ixm;
    for (i = 1; i <= i1; ++i) {
	buf_s_xm[i - 1] = q[index_r_xm[i] - 1];
    }
    i1 = *ixp;
    for (i = 1; i <= i1; ++i) {
	buf_s_xp[i - 1] = q[index_r_xp[i] - 1];
    }
    if (*ndim == 1) {
	goto L40;
    }
    i1 = *iym;
    for (i = 1; i <= i1; ++i) {
	buf_s_ym[i - 1] = q[index_r_ym[i] - 1];
    }
    i1 = *iyp;
    for (i = 1; i <= i1; ++i) {
	buf_s_yp[i - 1] = q[index_r_yp[i] - 1];
    }
    if (*ndim == 2) {
	goto L40;
    }
    i1 = *izm;
    for (i = 1; i <= i1; ++i) {
	buf_s_zm[i - 1] = q[index_r_zm[i] - 1];
    }
    i1 = *izp;
    for (i = 1; i <= i1; ++i) {
	buf_s_zp[i - 1] = q[index_r_zp[i] - 1];
    }
L40:

/* Send the new values of XNEW of the inner points to the neighbors */

    if (*rank_source_x > -1) {
	MPI_Sendrecv(buf_s_xm,*ixm,MPI_DOUBLE_PRECISION,*rank_source_x,tagx_min, 
		buf_r_xm,*ixm,MPI_DOUBLE_PRECISION,*rank_source_x,tagx_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_x > -1) {
	MPI_Sendrecv(buf_s_xp,*ixp,MPI_DOUBLE_PRECISION,*rank_dest_x,tagx_plus, 
		buf_r_xp,*ixp,MPI_DOUBLE_PRECISION,*rank_dest_x,tagx_min, 
		*comm_cart,&status );
    }

    if (*ndim == 1) {
	goto L50;
    }
    if (*rank_source_y > -1) {
	MPI_Sendrecv(buf_s_ym,*iym,MPI_DOUBLE_PRECISION,*rank_source_y,tagy_min, 
		buf_r_ym,*iym,MPI_DOUBLE_PRECISION,*rank_source_y,tagy_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_y > -1) {
	MPI_Sendrecv(buf_s_yp,*iyp,MPI_DOUBLE_PRECISION,*rank_dest_y,tagy_plus, 
		buf_r_yp,*iyp,MPI_DOUBLE_PRECISION,*rank_dest_y,tagy_min, 
		*comm_cart,&status );
    }

    if (*ndim == 2) {
	goto L50;
    }
    if (*rank_source_z > -1) {
	MPI_Sendrecv(buf_s_zm,*izm,MPI_DOUBLE_PRECISION,*rank_source_z,tagz_min, 
		buf_r_zm,*izm,MPI_DOUBLE_PRECISION,*rank_source_z,tagz_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_z > -1) {
	MPI_Sendrecv(buf_s_zp,*izp,MPI_DOUBLE_PRECISION,*rank_dest_z,tagz_plus, 
		buf_r_zp,*izp,MPI_DOUBLE_PRECISION,*rank_dest_z,tagz_min, 
		*comm_cart,&status );
    }
L50:

/*  Update the external points */

    i1 = *ixm;
    for (i = 1; i <= i1; ++i) {
	q[index_s_xm[i] - 1] = buf_r_xm[i - 1];
    }
    i1 = *ixp;
    for (i = 1; i <= i1; ++i) {
	q[index_s_xp[i] - 1] = buf_r_xp[i - 1];
    }
    if (*ndim == 1) {
	goto L60;
    }
    i1 = *iym;
    for (i = 1; i <= i1; ++i) {
	q[index_s_ym[i] - 1] = buf_r_ym[i - 1];
    }
    i1 = *iyp;
    for (i = 1; i <= i1; ++i) {
	q[index_s_yp[i] - 1] = buf_r_yp[i - 1];
    }
    if (*ndim == 2) {
	goto L60;
    }
    i1 = *izm;
    for (i = 1; i <= i1; ++i) {
	q[index_s_zm[i] - 1] = buf_r_zm[i - 1];
    }
    i1 = *izp;
    for (i = 1; i <= i1; ++i) {
	q[index_s_zp[i] - 1] = buf_r_zp[i - 1];
    }
L60:
    return 0;
} /* transf_3d */

/* ******************************************************************** */
/* ************************************************************************ */

/* Subroutine */ int transf_2d(double *q, double *ecoef, int *
	update_inv, int *n_update, int *comm_cart, int *me, 
	int *rank_source_x, int *rank_source_y, int *
	rank_dest_x, int *rank_dest_y, int *index_s_xm, int 
	*index_s_xp, int *index_s_ym, int *index_s_yp, int *
	index_r_xm, int *index_r_xp, int *index_r_ym, int *
	index_r_yp, int *ixp, int *ixm, int *iyp, int *iym, 
	int *ndim)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    static int tagx_min, tagy_min, i, tagx_plus, tagy_plus;
    MPI_Status status;

    /* Parameter adjustments */
    --index_r_yp;
    --index_r_ym;
    --index_r_xp;
    --index_r_xm;
    --index_s_yp;
    --index_s_ym;
    --index_s_xp;
    --index_s_xm;

    /* Function Body */
    tagx_min = 10;
    tagx_plus = 20;
    tagy_min = 30;
    tagy_plus = 40;

/*  Define the vectors to be transfered to the neighboring processors */
    i1 = *ixm;
    for (i = 1; i <= i1; ++i) {
	buf_s_xm[i - 1] = q[index_s_xm[i] - 1];
    }

    i1 = *ixp;
    for (i = 1; i <= i1; ++i) {
	buf_s_xp[i - 1] = q[index_s_xp[i] - 1];
    }

    i1 = *iym;
    for (i = 1; i <= i1; ++i) {
	buf_s_ym[i - 1] = q[index_s_ym[i] - 1];
    }

    i1 = *iyp;
    for (i = 1; i <= i1; ++i) {
	buf_s_yp[i - 1] = q[index_s_yp[i] - 1];
    }

/* Send the correction to the neighbors */

    if (*rank_source_x > -1) {
	MPI_Sendrecv(buf_s_xm,*ixm,MPI_DOUBLE_PRECISION,*rank_source_x,tagx_min, 
		buf_r_xm,*ixm,MPI_DOUBLE_PRECISION,*rank_source_x,tagx_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_x > -1) {
	MPI_Sendrecv(buf_s_xp,*ixp,MPI_DOUBLE_PRECISION,*rank_dest_x,tagx_plus, 
		buf_r_xp,*ixp,MPI_DOUBLE_PRECISION,*rank_dest_x,tagx_min, 
		*comm_cart,&status );
    }

    if (*rank_source_y > -1) {
	MPI_Sendrecv(buf_s_ym,*iym,MPI_DOUBLE_PRECISION,*rank_source_y,tagy_min, 
		buf_r_ym,*iym,MPI_DOUBLE_PRECISION,*rank_source_y,tagy_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_y > -1) {
	MPI_Sendrecv(buf_s_yp,*iyp,MPI_DOUBLE_PRECISION,*rank_dest_y,tagy_plus, 
		buf_r_yp,*iyp,MPI_DOUBLE_PRECISION,*rank_dest_y,tagy_min, 
		*comm_cart,&status );
    }

/*  Update qnewl() of the processor's inner points due to the 
    neighboring processors */

    i1 = *ixm;
    for (i = 1; i <= i1; ++i) {
	q[index_r_xm[i] - 1] += buf_r_xm[i - 1];
    }

    i1 = *ixp;
    for (i = 1; i <= i1; ++i) {
	q[index_r_xp[i] - 1] += buf_r_xp[i - 1];
    }

    i1 = *iym;
    for (i = 1; i <= i1; ++i) {
	q[index_r_ym[i] - 1] += buf_r_ym[i - 1];
    }

    i1 = *iyp;
    for (i = 1; i <= i1; ++i) {
	q[index_r_yp[i] - 1] += buf_r_yp[i - 1];
    }

/*  Get  the new values of the inner boundary points: */

    i1 = *n_update - 1;
    for (i = 0; i <= i1; ++i) {
	q[update_inv[i]] *= ecoef[update_inv[i]];
    }

    i1 = *ixm;
    for (i = 1; i <= i1; ++i) {
	buf_s_xm[i - 1] = q[index_r_xm[i] - 1];
    }
    i1 = *ixp;
    for (i = 1; i <= i1; ++i) {
	buf_s_xp[i - 1] = q[index_r_xp[i] - 1];
    }
    i1 = *iym;
    for (i = 1; i <= i1; ++i) {
	buf_s_ym[i - 1] = q[index_r_ym[i] - 1];
    }
    i1 = *iyp;
    for (i = 1; i <= i1; ++i) {
	buf_s_yp[i - 1] = q[index_r_yp[i] - 1];
    }

    if (*rank_source_x > -1) {
	MPI_Sendrecv(buf_s_xm,*ixm,MPI_DOUBLE_PRECISION,*rank_source_x,tagx_min, 
		buf_r_xm,*ixm,MPI_DOUBLE_PRECISION,*rank_source_x,tagx_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_x > -1) {
	MPI_Sendrecv(buf_s_xp,*ixp,MPI_DOUBLE_PRECISION,*rank_dest_x,tagx_plus, 
		buf_r_xp,*ixp,MPI_DOUBLE_PRECISION,*rank_dest_x,tagx_min, 
		*comm_cart,&status );
    }

    if (*rank_source_y > -1) {
	MPI_Sendrecv(buf_s_ym,*iym,MPI_DOUBLE_PRECISION,*rank_source_y,tagy_min, 
		buf_r_ym,*iym,MPI_DOUBLE_PRECISION,*rank_source_y,tagy_plus, 
		*comm_cart,&status );
    }

    if (*rank_dest_y > -1) {
	MPI_Sendrecv(buf_s_yp,*iyp,MPI_DOUBLE_PRECISION,*rank_dest_y,tagy_plus, 
		buf_r_yp,*iyp,MPI_DOUBLE_PRECISION,*rank_dest_y,tagy_min, 
		*comm_cart,&status );
    }

/*  Update the external points */

    i1 = *ixm;
    for (i = 1; i <= i1; ++i) {
	q[index_s_xm[i] - 1] = buf_r_xm[i - 1];
    }
    i1 = *ixp;
    for (i = 1; i <= i1; ++i) {
	q[index_s_xp[i] - 1] = buf_r_xp[i - 1];
    }
    i1 = *iym;
    for (i = 1; i <= i1; ++i) {
	q[index_s_ym[i] - 1] = buf_r_ym[i - 1];
    }
    i1 = *iyp;
    for (i = 1; i <= i1; ++i) {
	q[index_s_yp[i] - 1] = buf_r_yp[i - 1];
    }
    return 0;
} /* transf_2d */

/* ******************************************************************** */
/* ************************************************************************ */
