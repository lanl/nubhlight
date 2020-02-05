/******************************************************************************
 *                                                                            *
 * xdmf_output.c                                                              *
 *                                                                            *
 * Write the metadata for xdmf for Visit                                      *
 *                                                                            *
 ******************************************************************************/

/*
 * For information on xdmf see:
 * http://www.xdmf.org/index.php/XDMF_Model_and_Format
 * and
 * https://www.visitusers.org/index.php?title=Using_XDMF_to_read_HDF5
 * And this was a particularly useful example:
 * https://stackoverflow.com/questions/36718593/describing-5-dimensional-hdf5-matrix-with-xdmf
 *
 * Note that visit supports vectors and tensors... and in principle xdmf does
 * too however discussions online indicate that combining the two is dodgy, so I
 * treat everything as a scalar using hyperslabs ~JMM
 */

#include "decs.h"

static int  iX1_max;
static char dname[STRLEN], gname[STRLEN], name[STRLEN], fname[STRLEN];

void geom_meta(FILE *fp);
void coord_meta(FILE *fp, int indx);
void prim_meta(FILE *fp, const char *vnams[NVAR], int indx);
void scalar_meta(
    FILE *fp, const char *name, const char *sourcename, int precision);
void vec_meta(FILE *fp, const char *name);
void tensor_meta(
    FILE *fp, const char *name, const char *sourcename, int precision);

void vec_component(FILE *fp, const char *name, int indx);
void tensor_component(FILE *fp, const char *name, const char *sourcename,
    int precision, int mu, int nu);

void write_xml_file(int dump_id, double t, const char *vnams[NVAR]) {
// I know I don't need to do this every time... but it doesn't cost anything
// ~JMM
#if METRIC == MKS
  iX1_max = bl_i_of_r(Rout_vis);
#else
  iX1_max = N1TOT;
#endif // METRIC

  sprintf(dname, "../dump_%08d.h5", dump_id);
  sprintf(gname, "../grid.h5");
  sprintf(fname, "dump_%08d.xmf", dump_id);
  strcpy(name, xmfdir);
  strcat(name, fname);
  int full_dump = (dump_id % DTf == 0);

  fprintf(stdout, "XMF %s\n", name);

  FILE *fp = fopen(name, "w");

  // header
  fprintf(fp, "<?xml version=\"1.0\" ?>\n");
  fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fp, "<Xdmf Version=\"3.0\">\n");
  fprintf(fp, "  <Domain>\n");
  fprintf(fp, "    <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
  fprintf(fp, "      <Time Value=\"%16.14e\"/>\n", t);
  fprintf(fp,
      "      <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d "
      "%d\"/>\n",
      iX1_max + 1, N2TOT + 1, N3TOT + 1);

  geom_meta(fp); // Geometry
  fprintf(fp, "\n");

  // Jacobians
  fprintf(fp, "      <!-- JACOBIANS -->\n");
  fprintf(fp, "      <!-- contravariant -->\n");
  tensor_meta(fp, "Lambda_h2cart_con", gname, 64);
  fprintf(fp, "      <!-- covariant -->\n");
  tensor_meta(fp, "Lambda_h2cart_cov", gname, 64);
  fprintf(fp, "\n");

  // Metric
  fprintf(fp, "      <!-- METRIC -->\n");
  fprintf(fp, "      <!-- contravariant -->\n");
  tensor_meta(fp, "gcon", gname, 64);
  fprintf(fp, "      <!-- covariant -->\n");
  tensor_meta(fp, "gcov", gname, 64);
  fprintf(fp, "      <!-- determinant -->\n");
  scalar_meta(fp, "gdet", gname, 64);
  fprintf(fp, "      <!-- lapse -->\n");
  scalar_meta(fp, "alpha", gname, 64);
  fprintf(fp, "\n");

  // Variables
  fprintf(fp, "      <!-- PRIMITIVES -->\n");
  PLOOP prim_meta(fp, vnams, ip);
  fprintf(fp, "\n");

  if (full_dump) {
    fprintf(fp, "      <!-- DERIVED VARS -->\n");
    scalar_meta(fp, "divb", dname, 32);
    fprintf(fp, "      <!-- jcon -->\n");
    vec_meta(fp, "jcon");
#if OUTPUT_EOSVARS
    {
      scalar_meta(fp, "PRESS", dname, 32);
      scalar_meta(fp, "ENT", dname, 32);
      scalar_meta(fp, "TEMP", dname, 32);
      scalar_meta(fp, "CS2", dname, 32);
    }
#endif
#if ELECTRONS
    {
      scalar_meta(fp, "Qvisc", dname, 32);
#if RADIATION
      { scalar_meta(fp, "Qcoul", dname, 32); }
#endif // RADIATION
    }
#endif // ELECTRONS
#if RADIATION
    {
      scalar_meta(fp, "nph", dname, 32);
      fprintf(fp, "      <!-- Jrad -->\n");
      vec_meta(fp, "Jrad");
      fprintf(fp, "      <!-- Rmunu -->\n");
      tensor_meta(fp, "Rmunu", dname, 32);
    }
#endif // RADIATION
  }

  // footer
  fprintf(fp, "    </Grid>\n");
  fprintf(fp, "  </Domain>\n");
  fprintf(fp, "</Xdmf>\n");

  fclose(fp);
}

void geom_meta(FILE *fp) {
  fprintf(fp, "      <!-- GRID DEFINITION -->\n");
  fprintf(fp, "      <Geometry GeometryType=\"X_Y_Z\">\n");
  for (int d = 1; d < NDIM; d++)
    coord_meta(fp, d);
  fprintf(fp, "      </Geometry>\n");
}

void vec_meta(FILE *fp, const char *name) {
  DLOOP1 vec_component(fp, name, mu);
}

void tensor_meta(
    FILE *fp, const char *name, const char *sourcename, int precision) {
  DLOOP2 tensor_component(fp, name, sourcename, precision, mu, nu);
}

void coord_meta(FILE *fp, int indx) {
  fprintf(fp,
      "        <DataItem ItemType=\"Hyperslab\" Dimensions=\"%d %d %d\" "
      "Type=\"Hyperslab\">\n",
      iX1_max + 1, N2TOT + 1, N3TOT + 1);
  fprintf(fp, "          <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
  fprintf(fp, "            0 0 0 %d\n", indx);
  fprintf(fp, "            1 1 1 1\n");
  fprintf(fp, "            %d %d %d 1\n", iX1_max + 1, N2TOT + 1, N3TOT + 1);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp,
      "          <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" "
      "Precision=\"64\" Format=\"HDF\">\n",
      N1TOT + 1, N2TOT + 1, N3TOT + 1, NDIM);
  fprintf(fp, "            %s:/XFcart\n", gname);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp, "        </DataItem>\n");
}

void prim_meta(FILE *fp, const char *vnams[NVAR], int indx) {
  fprintf(fp,
      "      <Attribute Name=\"%s\" AttributeType=\"Scalar\" "
      "Center=\"Cell\">\n",
      vnams[indx]);
  fprintf(fp,
      "        <DataItem ItemType=\"Hyperslab\" Dimensions=\"%d %d %d\" "
      "Type=\"Hyperslab\">\n",
      iX1_max, N2TOT, N3TOT);
  fprintf(fp, "          <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
  fprintf(fp, "            0 0 0 %d\n", indx);
  fprintf(fp, "            1 1 1 1\n");
  fprintf(fp, "            %d %d %d 1\n", iX1_max, N2TOT, N3TOT);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp,
      "          <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" "
      "Precision=\"32\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT, NVAR);
  fprintf(fp, "            %s:/P\n", dname);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp, "        </DataItem>\n");
  fprintf(fp, "      </Attribute>\n");
}

void scalar_meta(
    FILE *fp, const char *name, const char *sourcename, int precision) {
  fprintf(fp,
      "      <Attribute Name=\"%s\" AttributeType=\"Scalar\" "
      "Center=\"Cell\">\n",
      name);
  fprintf(fp,
      "        <DataItem ItemType=\"Hyperslab\" Dimensions=\"%d %d %d\" "
      "Type=\"Hyperslab\">\n",
      iX1_max, N2TOT, N3TOT);
  fprintf(fp, "          <DataItem Dimensions=\"3 3\" Format=\"XML\">\n");
  fprintf(fp, "            0 0 0\n");
  fprintf(fp, "            1 1 1\n");
  fprintf(fp, "            %d %d %d\n", iX1_max, N2TOT, N3TOT);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp,
      "          <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" "
      "Precision=\"%d\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT, precision);
  fprintf(fp, "            %s:/%s\n", sourcename, name);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp, "        </DataItem>\n");
  fprintf(fp, "      </Attribute>\n");
}

void vec_component(FILE *fp, const char *name, int indx) {
  fprintf(fp,
      "      <Attribute Name=\"%s_%d\" AttributeType=\"Scalar\" "
      "Center=\"Cell\">\n",
      name, indx);
  fprintf(fp,
      "        <DataItem ItemType=\"Hyperslab\" Dimensions=\"%d %d %d\" "
      "Type=\"Hyperslab\">\n",
      iX1_max, N2TOT, N3TOT);
  fprintf(fp, "          <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
  fprintf(fp, "            0 0 0 %d\n", indx);
  fprintf(fp, "            1 1 1 1 \n");
  fprintf(fp, "            %d %d %d 1\n", iX1_max, N2TOT, N3TOT);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp,
      "          <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" "
      "Precision=\"32\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT, NDIM);
  fprintf(fp, "            %s:/%s\n", fname, name);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp, "        </DataItem>\n");
  fprintf(fp, "      </Attribute>\n");
}

void tensor_component(FILE *fp, const char *name, const char *sourcename,
    int precision, int mu, int nu) {
  fprintf(fp,
      "      <Attribute Name=\"%s_%d%d\" AttributeType=\"Scalar\" "
      "Center=\"Cell\">\n",
      name, mu, nu);
  fprintf(fp,
      "        <DataItem ItemType=\"Hyperslab\" Dimensions=\"%d %d %d\" "
      "Type=\"Hyperslab\">\n",
      iX1_max, N2TOT, N3TOT);
  fprintf(fp, "          <DataItem Dimensions=\"3 5\" Format=\"XML\">\n");
  fprintf(fp, "            0 0 0 %d %d\n", mu, nu);
  fprintf(fp, "            1 1 1 1 1\n");
  fprintf(fp, "            %d %d %d 1 1\n", iX1_max, N2TOT, N3TOT);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp,
      "          <DataItem Dimensions=\"%d %d %d %d %d\" NumberType=\"Float\" "
      "Precision=\"%d\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT, NDIM, NDIM, precision);
  fprintf(fp, "            %s:/%s\n", sourcename, name);
  fprintf(fp, "          </DataItem>\n");
  fprintf(fp, "        </DataItem>\n");
  fprintf(fp, "      </Attribute>\n");
}
