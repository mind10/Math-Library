#ifndef QSORT_WKEY
#define QSORT_WKEY

int dqsort_cmp_decrease(const void *a, const void *b);
int dqsort_cmp_increase(const void *a, const void *b);

/* num = size of base array */
/* width = sizeof(base[0]) */
/* the size of key = num */

void qsort_wkey(void *base, unsigned int num, unsigned int width, unsigned int *key, int (*comp)(const void *, const void *));

#endif /* QSORT_WKEY */

