#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <stddef.h>
#include "CRC32.cpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define  OUTPUT_LINE_MAX  1<<7
#define  BUF_SIZE  1<<20

struct buf_result {
	char *buf;
	int count;
};
typedef struct buf_result BUF_RESULT;

void print_estimation_br_nodisplay(int collasion_mat[], double *jac_sim_mat, int b, int r);
double *jac_sim(unsigned int (*)[493]); int up_tri_pos(int i, int j);
int next_uniq(int p, unsigned int * sh);
int collasion_check(int c_h, int (*min_hashs)[c_h], int fst, int sed, int b, int r);
void *collasion_mats_gen(int c_br);
void *shingle(char (*bufs)[500]);
int q_compare(const void *p1, const void *p2);
void *min_hash(unsigned int (*shingles)[493], int c_h);
void *hash_params_gen(int c_h);
void *txt_gen(void);
int *change_point_gen(void);
char **buf_gen(void);
char char_gen(int *char_tab, int origin);
int *char_tab_gen(void);
extern unsigned int CRC32(void *pData, size_t iLen);
void buf_write(char *line, BUF_RESULT *, FILE *fp);
void buf_flush(BUF_RESULT *, FILE *);

void print_estimation_br(int collasion_mat[], double *jac_sim_mat, int b, int r);
int print_estimation(int (*collasion_mats)[500 * 999], double *jac_sim_mat, int (*br)[2], int c_br);
void *mh_lsh(unsigned int (*shingles)[493], int (*br)[2], int c_br);
void lsh(int c_h, int (*min_hashs)[c_h], int (*collasion_mats)[500 * 999], int (*br)[2], int c_br);
