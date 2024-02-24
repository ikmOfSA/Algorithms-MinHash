#include "min_hash.h"

int main()
{
	unsigned sec = (unsigned)time(NULL);				//seeding random number generator
	srand(sec);

	clock_t t_start;						//used to monitor function running time
	double t_txt, t_shingle, t_jac_sim, t_mh_lsh, t_estimation, t_results;

	int br[6][2] = {{2, 15},					//b * r, the LSH parameters
		    {3, 10},
		    {5, 6},
		    {6, 5},
		    {10, 3},
		    {15, 2}};
	int c_br = sizeof(br)/sizeof(br[0]);

	t_start = clock();
	char (*bufs)[500];
	bufs = txt_gen();						//generate txts
	t_txt = (double)(clock() - t_start) / CLOCKS_PER_SEC;
	
	t_start = clock();
	unsigned int (*shingles)[493];
	shingles = shingle(bufs);					//shingle txts
	t_shingle = (double)(clock() - t_start) / CLOCKS_PER_SEC;

	t_start = clock();
	double *jac_sim_mat = jac_sim(shingles);			//calculate jaccard similarity between two shingles
	t_jac_sim = (double)(clock() - t_start) / CLOCKS_PER_SEC;

	t_start = clock();
	int (*collasion_mats)[500 * 999] = mh_lsh(shingles, br, c_br);	//1000 rounds of minhash-LSH
	t_mh_lsh = (double)(clock() - t_start) / CLOCKS_PER_SEC;

	t_start = clock();
	print_estimation(collasion_mats, jac_sim_mat, br, c_br);	//print results, write to file at the same time.
	t_results = (double)(clock() - t_start) / CLOCKS_PER_SEC;

	printf("Txts generation time = %fs\n", t_txt);
	printf("Shingling time = %fs\n", t_shingle);
	printf("Jaccard similarity calculating time = %fs\n", t_jac_sim);
	printf("Minhash-LSH time = %fs\n", t_mh_lsh);
	printf("Results output time = %fs\n", t_results);
	printf("----------\n");
	printf("Done.\n");
}

/*do 1000 * minhash-LSH cycle*/
void *mh_lsh(unsigned int (*shingles)[493], int (*br)[2], int c_br)
{
	int c_h = br[0][0] * br[0][1];

	int (*collasion_mats)[500 * 999] = collasion_mats_gen(c_br);	//upper triangular matrix to hold the collasion count of minhash/LSH
	int (*min_hashs)[c_h];

	printf("Calculating min-hash and LSH:\n");

	for (int i = 0; i < 1000; i++) {	//1000 * minhash-LSH rounds

		min_hashs = min_hash(shingles, c_h);	//do minhash
		lsh(c_h, min_hashs, collasion_mats, br, c_br);	//do LSH, put collasion in collasion matrix
		free(min_hashs);

		if (i % 10 == 0)
			printf("\t%d rounds finished\n", i);
	}
	printf("\t1000 rounds minhash-LSH finished\n");
	printf("----------\n");

	return (void *)collasion_mats;
}

/*generate txts*/
void *txt_gen()
{
	int *change_point = change_point_gen();
	char (*bufs)[500] = malloc(sizeof(char[1000][500]));
	int *char_tab = char_tab_gen();
	int c;

	struct stat st = {0};
	if (stat("txts", &st) == -1) {
		mkdir("txts", 0700);
	}

	for (int buf_pos = 0, buf_i; buf_pos < 500; buf_pos++) {

		c = char_gen(char_tab, 0);
		for(buf_i = 0; buf_i < 1000; buf_i++) {
			if(change_point[buf_i] == buf_pos)
				c = char_gen(char_tab, c);
			bufs[buf_i][buf_pos] = c;
		}
	}


	
	char * f_name = malloc(sizeof(char) * 20);	//write txts into files
	FILE *fp;
	for (int i = 0; i < 1000; i++) {
		sprintf(f_name, "txts/%d.txt", i+1);

		fp = fopen(f_name, "w");
		fwrite(bufs[i], sizeof(char), 500, fp);
		fclose(fp);
	}

	printf("Txt generating finished\n");
	return (void *)bufs;
}

/*shingling txts*/
void *shingle(char (*bufs)[500])
{
	unsigned int (*shingles)[493] = malloc(sizeof(unsigned int[1000][493]));

	for (int buf_i = 0, buf_pos; buf_i < 1000; buf_i++) {
		for (buf_pos = 0; buf_pos < 493; buf_pos++) {
			shingles[buf_i][buf_pos] = CRC32(&bufs[buf_i][buf_pos], 8);
		}

		qsort(shingles[buf_i], 493, sizeof(unsigned int), q_compare);	//sort the shingle txts, for use in Jaccard similarity calculation
	}

	printf("Shingling finished\n");
	return (void *)shingles;
}

/*calculate Jaccard similarity*/
double *jac_sim(unsigned int (*shingles)[493])
{	
	printf("Calculating jaccard similarity:\n");
	double jac_sim;
	double *jac_sim_mat = malloc(sizeof(double) * 500 * 999);
	unsigned int *sh1, *sh2;
	int p1, p2;
	int u, n;
	int pos;

	for (int i = 0, j; i < 999; i++)
		for (j = i + 1; j < 1000; j++) {
			sh1 = shingles[i];
			sh2 = shingles[j];

			p1 = 0;
			p2 = 0;

			u = 1;
			n = 0;

			/*Intuition: check two shingle vectors in one row, if two shingles are different, 
			  then one more different element found, union + 1, then more to check the next element; 
			  or if they are the same, then one more intersection found, intersection + 1.*/
			while (1) {

				if (sh1[p1] < sh2[p2]) {
					u++;
					p1 = next_uniq(p1, sh1);
					if (p1 == 0)
						break;
				}
				else if (sh1[p1] > sh2[p2]) {
					u++;
					p2 = next_uniq(p2, sh2);
					if (p2 == 0)
						break;
				}
				else if (sh1[p1] == sh2[p2]) {
					n++;
					p2 = next_uniq(p2, sh2);
					if (p2 ==0)
						break;
				}
			}

			if (p1 == 0) {
				while (1) {
					p2 = next_uniq(p2, sh2);
					if (p2 != 0)
						u++;
					else
						break;
				}
			}
			else if (p2 == 0) {
				while (1) {
					p1 = next_uniq(p1, sh1);
					if (p1 != 0)
						u++;
					else
						break;
				}
			}

			jac_sim = (double)n / u;
			pos = up_tri_pos(i, j);
			jac_sim_mat[pos] = jac_sim;

			if (pos % 100000 == 0)
				printf("\t%d%%\n", pos / 5000);
		}

	printf ("\tProcessing completes.\n");
	return jac_sim_mat;
}

/*hold the collasion count between two shingles*/
void *collasion_mats_gen(int c_br)
{
	int d = 999 * 500;
	int (*collasion_mats)[d] = malloc(c_br * sizeof(int) * d);
	
	for (int i = 0, j; i < c_br; i++)
		for (j = 0; j < d; j++)
			collasion_mats[i][j] = 0;	

	return (void *)collasion_mats;
}

/*do min-hash*/
void *min_hash(unsigned int (*shingles)[493], int c_h)
{
	int (*params)[2] = hash_params_gen(c_h);
	
	long p = 4294967311L;	//prime greater then 2^32.
	long n = 4294967296L;	//N, size of the hash table.
	unsigned int h;

	unsigned int (*min_hashs)[c_h] = malloc(sizeof(unsigned int[1000][c_h]));
	for (int i = 0, j; i < 1000; i++)
		for (j = 0; j < c_h; j++)
			min_hashs[i][j] = UINT_MAX;

	for (int i_buf = 0, i_pos, i_hash; i_buf < 1000; i_buf++) {
		for (i_pos = 0; i_pos < 493; i_pos++) {
			for (i_hash = 0; i_hash < c_h; i_hash++) {

				h = (long)(params[i_hash][0] * shingles[i_buf][i_pos] + params[i_hash][1]) % p;
				
				if (h < min_hashs[i_buf][i_hash])
					min_hashs[i_buf][i_hash] = h;
			}
		}
	}	

	return (void *)min_hashs;
}

/*do LSH*/
void lsh(int c_h, int (*min_hashs)[c_h], int (*collasion_mats)[500 * 999], int (*br)[2], int c_br)
{
	int b, r;
	int mat_pos;

	for (int i = 0; i < c_br; i++) {
		b = br[i][0];
		r = br[i][1];
		
		for (int fst = 0, sed; fst < 999; fst++)
			for (sed = fst + 1; sed < 1000; sed++) {
				if (collasion_check(c_h, min_hashs, fst, sed, b, r)) {
					mat_pos = up_tri_pos(fst, sed);
					collasion_mats[i][mat_pos]++;
				}
			}
	}
}

/*print retults and write to file in "results/.cvs"*/
int print_estimation(int (*collasion_mats)[500 * 999], double *jac_sim_mat, int (*br)[2], int c_br)
{
	int b, r;
	char display = 0;

	printf("Prepared to write results into files in \"results/*.cvs\", do you also want them to display in the monitor? (y, n):\n");
	scanf("%c", &display);

	if (display == 'y' || display == 'Y') {

		for (int i_br = 0; i_br < c_br; i_br++) {

			b = br[i_br][0];
			r = br[i_br][1];
			print_estimation_br(collasion_mats[i_br], jac_sim_mat, b, r);
		}
	}
	else {

		printf("Writing files...\n");
		printf("----------\n");
		for (int i_br = 0; i_br < c_br; i_br++) {

			b = br[i_br][0];
			r = br[i_br][1];
			print_estimation_br_nodisplay(collasion_mats[i_br], jac_sim_mat, b, r);
		}
	}
}

void print_estimation_br(int collasion_mat[], double *jac_sim_mat, int b, int r)
{
	struct stat st = {0};
	if (stat("results", &st) == -1) {
		mkdir("results", 0700);
	}

	char * f_name = malloc(sizeof(char) * 20);
	sprintf(f_name, "results/b=%d, r=%d.csv", b, r);
	FILE *fp = fopen(f_name, "w");

	BUF_RESULT *buf_result = malloc(sizeof(BUF_RESULT));
	buf_result->buf = malloc(sizeof(char) * BUF_SIZE);
	buf_result->count = 0;

	double pr;
	int pos;
	double prediction;
	double jac_sim;
	char *output_line = malloc(sizeof(char) * LINE_MAX);

	for (int i = 0, j; i < 999; i++)
		for (j = i + 1; j < 1000; j++) {

			pos = up_tri_pos(i, j);

			pr = collasion_mat[pos] / 1000.0;
			jac_sim = jac_sim_mat[pos];
			prediction = 1 - pow(1 - pow(jac_sim, r), b); 

			sprintf(output_line, "%d,%d,%8.6f,%8.6f,%8.6f\n", b, r, pr, jac_sim, prediction);
			buf_write(output_line, buf_result, fp);

			printf("(%d, %d)\b\b\tprobability = %8.4f%%\tsimilarity = %8.4f%%\tprediction = %8.4f%%\n", i, j, pr * 100, jac_sim * 100, prediction * 100);
		}	

	buf_flush(buf_result, fp);
	fclose(fp);
	free(buf_result->buf);
	free(buf_result);
	
	printf("----------\n");
	printf("Above is for b = %d, r = %d\n", b, r);
	printf("----------\n");
}

void print_estimation_br_nodisplay(int collasion_mat[], double *jac_sim_mat, int b, int r)
{
	struct stat st = {0};
	if (stat("results", &st) == -1) {
		mkdir("results", 0700);
	}

	char * f_name = malloc(sizeof(char) * 20);
	sprintf(f_name, "results/b=%d, r=%d.csv", b, r);
	FILE *fp = fopen(f_name, "w");

	BUF_RESULT *buf_result = malloc(sizeof(BUF_RESULT));
	buf_result->buf = malloc(sizeof(char) * BUF_SIZE);
	buf_result->count = 0;

	double pr;
	int pos;
	double prediction;
	double jac_sim;
	char *output_line = malloc(sizeof(char) * LINE_MAX);

	for (int i = 0, j; i < 999; i++)
		for (j = i + 1; j < 1000; j++) {

			pos = up_tri_pos(i, j);

			pr = collasion_mat[pos] / 1000.0;
			jac_sim = jac_sim_mat[pos];
			prediction = 1 - pow(1 - pow(jac_sim, r), b); 

			sprintf(output_line, "%d,%d,%8.6f,%8.6f,%8.6f\n", b, r, pr, jac_sim, prediction);
			buf_write(output_line, buf_result, fp);
		}	

	buf_flush(buf_result, fp);
	fclose(fp);
	free(buf_result->buf);
	free(buf_result);
}

/*return upper triangular matrix position by (row, column) index*/
int up_tri_pos(int i, int j)
{
	int pos = 0;
	int tmp;

	if (i > j) {
		
		tmp = i;
		i = j;
		j = tmp;
	}

	if (i == j)
		return 0;
	
	pos = (1999 - i) * i / 2 + j - i - 1;
	
	return pos;
}

/*return next uniq shingle position for ordered shingle matrix*/
int next_uniq(int p, unsigned int * sh)
{
	p++;

	while(1) {
		if (p == 493) {
			p = 0;
			break;
		}
		else if (sh[p] == sh[p - 1]) {
			p++;
		}
		else
			break;
	}

	return p;
}

/*use buffer to write to file to reduce I/O cost*/
void buf_write(char *line, BUF_RESULT *buf_ptr, FILE *fp)
{
	char *top = buf_ptr->buf + buf_ptr->count;
	buf_ptr->count += sprintf(top, "%s", line);

	if (buf_ptr->count + OUTPUT_LINE_MAX > BUF_SIZE) {
		fwrite(buf_ptr->buf, sizeof(char), buf_ptr->count, fp);
		buf_ptr->count = 0;
	}
}

/*flush remaining content in the buffer to file*/
void buf_flush(BUF_RESULT *buf_ptr, FILE *fp)
{
	fwrite(buf_ptr->buf, sizeof(char), buf_ptr->count, fp);
	buf_ptr->count = 0;
}

/*Check if two min-hash sigs have collasion*/
int collasion_check(int c_h, int (*min_hashs)[c_h], int fst, int sed, int b, int r)
{
	int collasion;
	int front;	//Mark the begining of the examined band

	for (int i = 0, j; i < b; i++) {

		collasion = 1;
		front = i * r;
		
		for (j = 0; j < r; j++) {
			int c1 = min_hashs[fst][front + j];
			int c2 = min_hashs[sed][front + j];

			if (c1 != c2) {
				collasion = 0;
				break;
			}
		}

		if (collasion == 1)
			break;
	}

	return collasion;
}

/*used by qsort()*/
int q_compare(const void *p1, const void *p2)
{
	long fst = *(unsigned int *)p1;
	long sed = *(unsigned int *)p2;
	int compare;

	if (fst > sed)
		compare = 1;
	else if (fst < sed)
		compare = -1;
	else
		compare = 0;

	return compare;
}

/*generate hash function parameters*/
void *hash_params_gen(int c_h)
{
	int c = c_h * 2;
	int *hash_params = malloc(sizeof(int) * c);

	for (int i = 0; i < c; i++)
		hash_params[i] = rand();

	return (void *)hash_params;
}

/*generate difference points between txts*/
int *change_point_gen()
{
	int r;
	int *change_point = malloc(sizeof(int) * 1000);

	change_point[0] = INT_MAX;
	for (int i = 1; i < 1000; i++) {
		r = rand();
		r = (int)floor(r / (RAND_MAX + 1.0) * 500);
		change_point[i] = r;
	}

	return change_point;
}

/*generate random alphabet character, case sensitive*/
char char_gen(int *char_tab, int origin)
{
	int r;
	int c = origin;

	while (c == origin) {

		r = rand();
		r = (int)floor(r / (RAND_MAX + 1.0) * 52);
		c = char_tab[r];
	}

	return (char)c;
}

/*generate alphabet character table, used by char_gen()*/
int *char_tab_gen()
{
	int *char_tab = malloc(sizeof(int) * 52);

	for (int i = 0, j = 65; i < 52; i++, j++) {
		while(j > 'Z' && j < 'a')
		       j++;
		char_tab[i] = j;	
	}
	
	return char_tab;
}
