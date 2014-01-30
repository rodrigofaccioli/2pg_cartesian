#include <malloc.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define TAM_BLOCO 400
#define MAX(A,B) (((A)>=(B)) ? (A) : (B))
#define MIN(A,B) (((A)<=(B)) ? (A) : (B))
#define MAX_COMMAND 700
#define MAX_PATH_FILE_NAME 220
#define MAX_FILE_NAME 40
#define MAX_PATH 200
#define MAX_LINE_FILE 500
#define MAX_INT_NUMBER_RANDOM 999999
/* This macro calculates the size of a array - idea from GROMACS source-code*/
#define asize(a) (sizeof(a)/sizeof((a)[0]))
#define PI 3.14159265
#define MAX_STR_DIHEDRAL_ANGLE_TYPE 20
#define MAX_RANDOM_STRING 50
#define PREFIX_PDB_FILE_NAME "PROT_"
#define PREFIX_PDB_FILE_NAME_EA "PROT_IND_"
#define PREFIX_PDB_FILE_NAME_EA_FINAL "PROT_IND_FINAL_"
#define PREFIX_PDB_FILE_NAME_AUX "PROT_AUX_"

/* These defines are used for distance crowding*/
#define MIN_DIST -999999999
#define MAX_DIST 999999999