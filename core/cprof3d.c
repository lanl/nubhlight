#include "cprof3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#define TRUE  1
#define FALSE 0

#define MYH5CHECK(ierr) \
    if(ierr < 0) { \
        fprintf(stderr, "HDF5 Error. %s:%d\n", __FILE__, __LINE__); \
        fflush(stderr); \
        exit(1); \
    }

#define SIGN(X) ((X<0)?(-1):(1))

typedef struct {
    double origin[3];
    double delta[3];
    int n[3];
} Reflevel;

typedef struct {
    Reflevel * levels;
    int nlevels;
} Grid;

typedef struct {
    Grid const * grid;
    double ** data;
} DataSet;

typedef struct {
    hid_t file_id;
    Grid * grid;
} DataFile;

void cprof3d_cmap_default(
        double x_in, double y_in, double z_in,
        double * x_out, double * y_out, double * z_out) {
    *x_out = x_in;
    *y_out = y_in;
    *z_out = z_in;
}
void cprof3d_cmap_reflecting_xy(
        double x_in, double y_in, double z_in,
        double * x_out, double * y_out, double * z_out) {
    *x_out = x_in;
    *y_out = y_in;
    *z_out = fabs(z_in);
}

double cprof3d_transf_default(
        double var, double x, double y, double z) {
    (void)x;
    (void)y;
    (void)z;
    return var;
}
double cprof3d_transf_reflecting_xy(
        double var, double x, double y, double z) {
    (void)x;
    (void)y;
    return SIGN(z)*var;
}

void * cprof3d_open_file(
        char const * fname) {
    char gname[BUFSIZ];
    herr_t ierr;
    hsize_t dims[3];

    DataFile * self = malloc(sizeof(DataFile));
    self->file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
        MYH5CHECK(self->file_id);

    self->grid = malloc(sizeof(Grid));
    self->grid->nlevels = 0;
    for(int il = 0; ; ++il) {
        snprintf(gname, BUFSIZ, "reflevel=%d", il);
        htri_t group_exists = H5Lexists(self->file_id, gname, H5P_DEFAULT);
        if(group_exists == TRUE) {
            self->grid->nlevels++;
        }
        else {
            break;
        }
    }
    if(0 == self->grid->nlevels) {
        fprintf(stderr, "%s is not a valid profile file!", fname);
        fflush(stderr);
        exit(1);
    }

    self->grid->levels = malloc(self->grid->nlevels*sizeof(Reflevel));
    for(int il = 0; il < self->grid->nlevels; ++il) {
        snprintf(gname, BUFSIZ, "reflevel=%d", il);

        ierr = H5LTget_attribute_double(self->file_id, gname, "delta",
                &self->grid->levels[il].delta[0]);
            MYH5CHECK(ierr);

        double extent[6];
        ierr = H5LTget_attribute_double(
                self->file_id, gname, "extent", &extent[0]);
            MYH5CHECK(ierr);

        hid_t group_id = H5Gopen(self->file_id, gname, H5P_DEFAULT);
            MYH5CHECK(group_id);
        ierr = H5LTget_dataset_info(group_id, "rho", &dims[0], NULL, NULL);
            MYH5CHECK(ierr);
        H5Gclose(group_id);

        for(int d = 0; d < 3; ++d) {
            self->grid->levels[il].origin[d] = extent[2*d];
            self->grid->levels[il].n[d] = dims[d];
        }

        /* Assume that we are dealing with cell centered data from Carpet */
        if(0 == il) {
            for(int d = 0; d < 3; ++d) {
                self->grid->levels[il].origin[d] +=
                    self->grid->levels[il].delta[d];
            }
        }
        else {
            for(int d = 0; d < 3; ++d) {
                self->grid->levels[il].origin[d] +=
                    0.5*self->grid->levels[il].delta[d];
            }
        }
    }

    return (void *)self;
}

void cprof3d_close_file(
        void * dfile) {
    DataFile * self = (DataFile *)dfile;
    H5Fclose(self->file_id);
    free(self->grid->levels);
    free(self->grid);
    free(self);
}

void * cprof3d_read_dset(
        void const * _dfile,
        char const * vname) {
    DataFile * dfile = (DataFile *)_dfile;

    DataSet * dset = malloc(sizeof(DataSet));
    dset->data = malloc(dfile->grid->nlevels*sizeof(double *));
    dset->grid = dfile->grid;

    char gname[BUFSIZ];
    herr_t ierr;
    hid_t group_id;
    hsize_t dims[3];
    for(int il = 0; il < dfile->grid->nlevels; ++il) {
        snprintf(gname, BUFSIZ, "reflevel=%d", il);
        group_id = H5Gopen(dfile->file_id, gname, H5P_DEFAULT);
            MYH5CHECK(group_id);

        ierr = H5LTget_dataset_info(group_id, vname, &dims[0], NULL, NULL);
            MYH5CHECK(ierr);
        dset->data[il] = malloc(dims[0]*dims[1]*dims[2]*sizeof(double));
        ierr = H5LTread_dataset_double(group_id, vname, dset->data[il]);
            MYH5CHECK(ierr);

        H5Gclose(group_id);
    }

    return (void *)dset;
}

void cprof3d_del_dset(
        void * dset) {
    DataSet * self = (DataSet *)dset;
    for(int i = 0; i < self->grid->nlevels; ++i) {
        free(self->data[i]);
    }
    free(self->data);
    free(self);
}

static inline int index3d(
        Reflevel const * rlevel,
        int i, int j, int k) {
    return k + rlevel->n[2]*(j + rlevel->n[1]*i);
}

static inline double grid_point(
        Reflevel const * rlevel,
        int const d, int const idx) {
    return rlevel->origin[d] + rlevel->delta[d]*idx;
}

static bool interp_on_level(
        Reflevel const * rlevel,
        double const * idata,
        cprof3d_cmap cmap,
        cprof3d_transf transf,
        int const np,
        double const * xp,
        double const * yp,
        double const * zp,
        double * odata,
        bool * interpolated) {
    int const stride[3] = {
        index3d(rlevel, 1, 0, 0),
        index3d(rlevel, 0, 1, 0),
        index3d(rlevel, 0, 0, 1),
    };
#pragma omp parallel for
    for(int ip = 0; ip < np; ++ip) {
        if(interpolated[ip]) {
            continue;
        }
        double rp[3];
        cmap(xp[ip], yp[ip], zp[ip], &rp[0], &rp[1], &rp[2]);

        interpolated[ip] = true;
        for(int d = 0; d < 3; ++d) {
            interpolated[ip] = interpolated[ip] &&
                (rp[d] >= rlevel->origin[d]);
            interpolated[ip] = interpolated[ip] &&
                (rp[d] <= grid_point(rlevel, d, rlevel->n[d]-1));
        }
        if(!interpolated[ip]) {
            continue;
        }

        int idx[3]; 
        double weight[3][2];
        for(int d = 0; d < 3; ++d) {
            idx[d] = (rp[d] - rlevel->origin[d])/rlevel->delta[d];
            weight[d][1] = (rp[d] - grid_point(rlevel, d, idx[d])) /
                rlevel->delta[d];
            weight[d][0] = (1.0 - weight[d][1]);
        }

        int const ijk = index3d(rlevel, idx[0], idx[1], idx[2]);
        double sum = 0;
        for(int di = 0; di < 2; ++di)
        for(int dj = 0; dj < 2; ++dj)
        for(int dk = 0; dk < 2; ++dk) {
            sum += weight[0][di]*weight[1][dj]*weight[2][dk]*
                idata[ijk + di*stride[0] + dj*stride[1] + dk*stride[2]];
        }
        odata[ip] = transf(sum, xp[ip], yp[ip], zp[ip]);
    }
    bool done = true;
    for(int ip = 0; ip < np; ++ip) {
        done = done && interpolated[ip];
    }
    return done;
}

bool cprof3d_interp(
        void const * dset,
        cprof3d_cmap cmap,
        cprof3d_transf transf,
        int const np,
        double const * xp,
        double const * yp,
        double const * zp,
        double * data) {
    DataSet const * self = (DataSet const *)dset;

    bool * interpolated = malloc(np*sizeof(bool));
    memset(interpolated, 0, np*sizeof(bool));

    bool done = false;
    for(int il = self->grid->nlevels - 1; il >= 0; --il) {
        done = interp_on_level(&self->grid->levels[il], self->data[il],
                cmap, transf, np, xp, yp, zp, data, interpolated);
        if(done) {
            break;
        }
    }

    free(interpolated);

    return done;
}
