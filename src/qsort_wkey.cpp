// Syungkwon Ra : Modified from sanos source index

//
// qsort.c
//
// Quick sort
//
// Copyright (C) 2002 Michael Ringgaard. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 
// 1. Redistributions of source code must retain the above copyright 
//    notice, this list of conditions and the following disclaimer.  
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.  
// 3. Neither the name of the project nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission. 
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
// SUCH DAMAGE.
// 


static void shortsort(char *lo, char *hi, unsigned int width, unsigned int* lo_key, unsigned int* hi_key, int (*comp)(const void *, const void *));
static void swap(char *p, char *q, unsigned int width, unsigned int *p_key, unsigned int *q_key);

/* num = size of base array */
/* width = sizeof(base[0]) */

void qsort_wkey(void *base, unsigned int num, unsigned int width, unsigned int *key, int (*comp)(const void *, const void *))
{

#define CUTOFF 8

	char *lo, *hi;
	char *mid;
	char *loguy, *higuy;
	unsigned int size;
	char *lostk[30], *histk[30];
	unsigned int *lostk_key[30], *histk_key[30];
	int stkptr;
	unsigned int *lo_key, *hi_key, *mid_key;
	unsigned int *loguy_key, *higuy_key;

	if (num < 2 || width == 0) return;
	stkptr = 0;

	lo = (char *) base;
	hi = (char *) base + width * (num - 1);

	lo_key = key;
	hi_key = key + (num - 1);

recurse:
	size = ((unsigned int)(hi - lo)) / width + 1;

	if (size <= CUTOFF) {
		shortsort(lo, hi, width, lo_key, hi_key, comp);
	}
	else {
		mid = lo + (size / 2) * width;
		mid_key = lo_key + (size / 2);

		swap(mid, lo, width, mid_key, lo_key);

		loguy = lo;
		higuy = hi + width;

		loguy_key = lo_key;
		higuy_key = hi_key + 1;

		for (;;) {
			do { 
				loguy += width; 
				loguy_key++;
			} while (loguy <= hi && comp(loguy, lo) <= 0);

			do { 
				higuy -= width; 
				higuy_key--;
			} while (higuy > lo && comp(higuy, lo) >= 0);

			if (higuy < loguy) break;

			swap(loguy, higuy, width, loguy_key, higuy_key);
		}

		swap(lo, higuy, width, lo_key, higuy_key);

		if (higuy - 1 - lo >= hi - loguy) {
			if (lo + width < higuy) {
				lostk[stkptr] = lo;
				histk[stkptr] = higuy - width;

				lostk_key[stkptr] = lo_key;
				histk_key[stkptr] = higuy_key - 1;

				++stkptr;
			}

			if (loguy < hi) {
				lo = loguy;
				lo_key = loguy_key;
				goto recurse;
			}
		}
		else {
			if (loguy < hi) {
				lostk[stkptr] = loguy;
				histk[stkptr] = hi;

				lostk_key[stkptr] = loguy_key;
				histk_key[stkptr] = hi_key;

				++stkptr;
			}

			if (lo + width < higuy) {
				hi = higuy - width;
				hi_key = higuy_key - 1;
				goto recurse;
			}
		}
	}

	--stkptr;
	if (stkptr >= 0) {
		lo = lostk[stkptr];
		hi = histk[stkptr];
		lo_key = lostk_key[stkptr];
		hi_key = histk_key[stkptr];

		goto recurse;
	}
	else
		return;
}

static void shortsort(char *lo, char *hi, unsigned int width, unsigned int *lo_key, unsigned int *hi_key, int (*comp)(const void *, const void *))
{
	char *p, *max;
	unsigned int *max_key;

	while (hi > lo) {
		max = lo;
		max_key = lo_key;
		for (p = lo+width; p <= hi; p += width) {
			if (comp(p, max) > 0) { 
				max = p;
				max_key = lo_key + ((unsigned int)(p-lo)) / width;
			}
		}

		swap(max, hi, width, max_key, hi_key);

		hi -= width;
		hi_key--;
	}
}

static void swap(char *a, char *b, unsigned int width, unsigned int *a_key, unsigned int *b_key)
{
	char tmp;
	unsigned int tmp_key;

	if (a != b) {
		while (width--) {
			tmp = *a;
			*a++ = *b;
			*b++ = tmp;
		}
		tmp_key = *a_key;
		*a_key = *b_key;
		*b_key = tmp_key;
	}
}

int
dqsort_cmp_decrease(const void *a, const void *b)
{
	double *data_a, *data_b;
	data_a = (double *)a;
	data_b = (double *)b;

	if (*data_a > *data_b) return 1;
	else if (*data_a < *data_b) return -1;
	else return 0;
}

int
dqsort_cmp_increase(const void *a, const void *b)
{
	double *data_a, *data_b;
	data_a = (double *)a;
	data_b = (double *)b;

	if (*data_a < *data_b) return 1;
	else if (*data_a > *data_b) return -1;
	else return 0;
}

/*

int _tmain(int argc, _TCHAR* argv[])
{
#define SIZE 53
	unsigned int size = SIZE;
	double data1[SIZE];
	double data2[SIZE];
	int i;
	unsigned int key[SIZE];

	for (i = 0; i < SIZE; i++) {
		key[i] = i;
	}

	double tmp;
	for (i = 0; i < SIZE; i++) {
		tmp = ((double)rand()) / ((double)RAND_MAX);
		data1[i] = tmp;
		data2[i] = tmp;
	}

	printf("data : ");
	for(i = 0; i < SIZE; i++) {
		printf("%d %lf\n", i, data1[i]);
	}

	qsort_wkey(data1, SIZE, sizeof(double), key, bigger);
	qsort(data2, SIZE, sizeof(double), bigger);

	printf("data : ");
	for(i = 0; i < SIZE; i++) {
		printf("%d %lf\n", i, data1[i]);
	}
	printf("\n");

	printf("key : ");
	for( int i = 0; i < SIZE; i++) {
		printf("%d ", key[i]);
	}
	printf("\n");

	printf("sorted data : ");
	for(i = 0; i < SIZE; i++) {
		if (data1[i] != data2[i]) printf("-");
	}
	printf("\n");

	return 0;
}

*/
