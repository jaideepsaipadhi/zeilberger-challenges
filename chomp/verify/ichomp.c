/* ichomp.c - independent Chomp verifier, written from scratch (Oct 2026)
 * without reference to chomp.c / chomp2.c / chomp3.c.
 *
 * Position: row lengths l[0] >= l[1] >= ... >= l[a-1], 0 <= l[i] <= b, l[0] >= 1.
 * Cell (0,0) is poisoned.  The position consisting of the poison alone is P.
 * Bite (r,c) with c < l[r], (r,c) != (0,0):  l[k] = min(l[k], c) for k >= r.
 *
 * Rank: lexicographic order of (l[0],...,l[a-1]).  A bite never increases any
 * row, so a child is lexicographically smaller and its rank is smaller: one
 * increasing sweep classifies everything.  rank(l) = sum_i C(a-i-1+l[i], a-i)
 * (count of sequences that agree up to row i-1 and are smaller at row i).
 *
 * The ranks are recomputed from scratch for every child (O(a) per child), with
 * no incremental tricks, so this shares no shortcut with the production engine.
 *
 * Usage:  ./ichomp A B            solve the A x B box, report every bar inside
 *         ./ichomp A B a' b'      ... and print the winning bites of a' x b'
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int A, B;
static uint64_t binom[80][80];
static uint8_t *Pbit;            /* bit set => P-position */

static inline uint64_t rank_of(const int *l) {
    uint64_t r = 0;
    for (int i = 0; i < A; i++) {
        int k = A - i - 1;
        r += binom[k + l[i]][k + 1];
    }
    return r;
}
static inline int isP(uint64_t r) { return (Pbit[r >> 3] >> (r & 7)) & 1; }
static inline void setP(uint64_t r) { Pbit[r >> 3] |= (uint8_t)(1u << (r & 7)); }

/* returns 1 if position l (already ranked smaller children done) is P */
static int classify(const int *l) {
    int ch[64];
    if (l[0] == 1 && (A == 1 || l[1] == 0)) return 1;     /* poison only */
    for (int r = 0; r < A && l[r] > 0; r++) {
        for (int c = 0; c < l[r]; c++) {
            if (r == 0 && c == 0) continue;
            memcpy(ch, l, sizeof(int) * A);
            for (int k = r; k < A; k++) if (ch[k] > c) ch[k] = c;
            if (isP(rank_of(ch))) return 0;              /* move to P => N */
        }
    }
    return 1;
}

static void bites(int a2, int b2) {
    int l[64], ch[64], n = 0;
    for (int i = 0; i < A; i++) l[i] = (i < a2) ? b2 : 0;
    printf("%d x %d winning bites (row,col, 1-based):", a2, b2);
    for (int r = 0; r < a2; r++)
        for (int c = 0; c < b2; c++) {
            if (r == 0 && c == 0) continue;
            memcpy(ch, l, sizeof(int) * A);
            for (int k = r; k < A; k++) if (ch[k] > c) ch[k] = c;
            if (isP(rank_of(ch))) { printf(" (%d,%d)", r + 1, c + 1); n++; }
        }
    printf("   total %d\n", n);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s A B [a b]\n", argv[0]); return 1; }
    A = atoi(argv[1]); B = atoi(argv[2]);
    for (int n = 0; n < 80; n++) { binom[n][0] = 1; for (int k = 1; k <= n; k++) binom[n][k] = binom[n-1][k-1] + (k <= n-1 ? binom[n-1][k] : 0); }
    uint64_t total = binom[A + B][A];
    fprintf(stderr, "box %dx%d: %llu positions, %.2f GB\n", A, B,
            (unsigned long long)total, total / 8.0 / 1e9);
    Pbit = calloc(total / 8 + 1, 1);
    if (!Pbit) { fprintf(stderr, "out of memory\n"); return 1; }

    /* enumerate in lex order: odometer over non-increasing sequences */
    int l[64]; for (int i = 0; i < A; i++) l[i] = 0;
    uint64_t r = 0, nP = 0;
    for (;;) {
        if (l[0] > 0) {
            if (rank_of(l) != r) { fprintf(stderr, "rank mismatch at %llu\n", (unsigned long long)r); return 2; }
            if (classify(l)) { setP(r); nP++; }
        }
        r++;
        /* next sequence in lex order */
        int i = A - 1;
        while (i >= 0 && l[i] == (i == 0 ? B : l[i-1])) i--;
        if (i < 0) break;
        l[i]++;
        for (int k = i + 1; k < A; k++) l[k] = 0;
        if ((r & ((1ULL << 30) - 1)) == 0) fprintf(stderr, "  %.1f%%\n", 100.0 * r / total);
    }
    if (r != total) { fprintf(stderr, "count mismatch %llu vs %llu\n", (unsigned long long)r, (unsigned long long)total); return 2; }
    fprintf(stderr, "done: %llu P-positions\n", (unsigned long long)nP);

    /* census: number of winning first moves for every bar inside the box */
    if (argc >= 5) { bites(atoi(argv[3]), atoi(argv[4])); return 0; }
    for (int a2 = 2; a2 <= A; a2++)
        for (int b2 = 2; b2 <= B; b2++) {
            int lb[64], ch[64], n = 0;
            for (int i = 0; i < A; i++) lb[i] = (i < a2) ? b2 : 0;
            for (int rr = 0; rr < a2; rr++)
                for (int c = 0; c < b2; c++) {
                    if (rr == 0 && c == 0) continue;
                    memcpy(ch, lb, sizeof(int) * A);
                    for (int k = rr; k < A; k++) if (ch[k] > c) ch[k] = c;
                    n += isP(rank_of(ch));
                }
            if (n != 1) printf("%dx%d: %d winning moves\n", a2, b2, n);
        }
    return 0;
}
