/* Copyright 2020 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
#include <err.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pcg_variants.h"

// count bits set in a 64 bit int
#define popcount64(x) __builtin_popcountll(x)

static inline size_t
rank(size_t n, uint64_t *v, uint64_t x)
{
	size_t i;
	for (i=0; i<=n-2; i++) {
		if ((v[i] & x) == x)
			return i;
	}
	// shouldn't get here
	errx(1, "x has no ancestor in v");
}

/*
 * RNNI distance between trees T and R.
 * Collienne & Gavryushkin (2020) https://arxiv.org/abs/2007.12307
 */
uint64_t
t64_rnni_distance(size_t n, uint64_t * T, uint64_t * R)
{
	size_t i, r;
	uint64_t d = 0;
	uint64_t *T1;

	assert(n > 1 && n <= 64);

	if ((T1 = malloc(sizeof(*T1) * (n-1))) == NULL) {
		err(1, "malloc");
	}
	memcpy(T1, T, sizeof(*T1) * (n-1));

	for (i=n-2; i<n-1; i--) {
		r = rank(n, T1, R[i]);
		while (r > i) {
			uint64_t u, v;
			v = T1[r];
			u = T1[r-1];
			if (v & u) {
				uint64_t w, x, y;
				// u is a child of v
				// now find children of u
				if (popcount64(u) == 2) {
					// both children are leaves; make x the
					// child with the right-most bit set
					x = u & -u;
				} else {
					size_t j;
					// descend through T1 until we find a
					// child, x, of u
					for (j=r-1; j<=r-1; j--) {
						x = T1[j];
						if ((u & x) == x)
							break;
					}
					if (j > r-1) {
						// shouldn't get here
						errx(1, "x has no ancestor in u");
					}
				}
				y = u - x;
				w = v - u;
				if ((R[i] & (x + w)) == x + w) {
					T1[r-1] = x + w;
				} else {
					T1[r-1] = y + w;
				}
			} else {
				// swap u and v
				T1[r] = u;
				T1[r-1] = v;
			}
			d++;
			r--;
		}
	}

	free(T1);

	return d;
}

/*
 * Generate a random tree with n leaves.
 * Each node in the tree is represented as a 64 bit int, where the binary bits
 * in the integer correspond to the node's decendent leaves.
 */
void
t64_random(pcg32_random_t * rng, size_t n, uint64_t * v)
{
	size_t i, j;
	size_t a, b;
	uint64_t *u;

	assert(n > 1 && n <= 64);

	if ((u = malloc(sizeof(*u) * n)) == NULL) {
		err(1, "malloc");
	}

	// leaves of the tree
	for (i=0; i<n; i++) {
		u[i] = 1ull << i;
	}

	for (j=n, i=0; i<n-1; i++) {
		// pick two nodes, a and b, to be joined
		a = pcg32_boundedrand_r(rng, j);
		b = pcg32_boundedrand_r(rng, --j);

		v[i] = u[a];
		// put last entry into a's position so that u[a] != u[b]
		u[a] = u[j];
		// combine a and b into new node
		v[i] |= u[b];
		// replace b's position with the new node
		u[b] = v[i];
	}

	free(u);
}

void
t64_print(size_t n, uint64_t * v)
{
	size_t i, j;

	assert(n > 1 && n <= 64);

	for (i=n-2; i<=n-2; i--) {
		for (j=0; j<n; j++) {
			printf("%d", (uint8_t)(v[i] >> j & 0x1));
		}
		printf("\n");
	}
}

int
main()
{
	const size_t n = 50; // number of leaves
	const size_t nreps = 100000;
	pcg32_random_t rng;
	uint64_t v[n-1], u[n-1];
	uint64_t d = 0;
	size_t i;

	pcg32_srandom_r(&rng, 12u, 34u);
	for (i=0; i<nreps; i++) {
		t64_random(&rng, n, v);
		t64_random(&rng, n, u);
		d += t64_rnni_distance(n, v, u);
	}
	//t64_print(n, v);
	printf("%lf\n", (double)d/nreps);

	return 0;
}
