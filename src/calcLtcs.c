///////////////////////////////////////////////////////////////////////////////
// Author: Matthias Gerstl 
// Email: matthias.gerstl@acib.at 
// Company: Austrian Centre of Industrial Biotechnology (ACIB) 
// Web: http://www.acib.at Copyright
// (C) 2015 Published unter GNU Public License V3
///////////////////////////////////////////////////////////////////////////////
//Basic Permissions.
// 
// All rights granted under this License are granted for the term of copyright
// on the Program, and are irrevocable provided the stated conditions are met.
// This License explicitly affirms your unlimited permission to run the
// unmodified Program. The output from running a covered work is covered by
// this License only if the output, given its content, constitutes a covered
// work. This License acknowledges your rights of fair use or other equivalent,
// as provided by copyright law.
// 
// You may make, run and propagate covered works that you do not convey,
// without conditions so long as your license otherwise remains in force. You
// may convey covered works to others for the sole purpose of having them make
// modifications exclusively for you, or provide you with facilities for
// running those works, provided that you comply with the terms of this License
// in conveying all material for which you do not control copyright. Those thus
// making or running the covered works for you must do so exclusively on your
// behalf, under your direction and control, on terms that prohibit them from
// making any copies of your copyrighted material outside their relationship
// with you.
// 
// Disclaimer of Warranty.
// 
// THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE
// LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR
// OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND,
// EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
// ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.
// SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY
// SERVICING, REPAIR OR CORRECTION.
// 
// Limitation of Liability.
// 
// IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL
// ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE
// PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
// GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE
// OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA
// OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
// PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
// EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGES.
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#include "generalFunctions.c"
#include "efmMethods.c"
#include "bitmakros.h"

#define BITSIZE        CHAR_BIT

#define MAX_ARGS       11
#define ARG_INPUT      0
#define ARG_LTCS_OUT   1
#define ARG_SFILE      2
#define ARG_RFILE      3
#define ARG_RVFILE     4
#define ARG_LOOPS      5
#define ARG_ZERO       6 
#define ARG_THREADS    7 
#define ARG_FULL_OUT   8
#define ARG_CSV_OUT    9
#define ARG_ANALYSIS   10

#define ERROR_ARGS     1
#define ERROR_ZERO_NR  2
#define ERROR_FILE     3
#define ERROR_RAM      4
#define ERROR_EFM      5

struct thread_args
{
    int max_threads;
    int thread_id;
    int bitarray_size;
    unsigned long ltcs_count;
    unsigned long efm_count;
    unsigned int reaction;
    char** ltcs;
    char** matrix;
    unsigned long* new_ltcs_count;
    char** new_ltcs;
};

/*
 * return actual time as string
 */
char* getTime()
{
    static char time_string[20];
    time_t now = time (0);
    strftime (time_string, 100, "%Y-%m-%d %H:%M:%S", localtime (&now));
    return time_string;
}

/**
 * calculates needed number of chars for bit support of ltcs
 */
unsigned long getBitsize(unsigned long count)
{
    unsigned long bitsize = count / BITSIZE;
    int modulo  = count % BITSIZE;
    if (modulo > 0)
    {
        bitsize++;
    }
    return bitsize;
}

/**
 * find subsets of LTCS and set memory of those LTCS free
 */
void filterLtcs(char** ltcs, char** notAnLtcs, unsigned long ltcs_count,
        unsigned long efm_count)
{
    unsigned long bitarray_size = getBitsize(ltcs_count);
    char* m_notAnLtcs = calloc(1, bitarray_size);
    unsigned long li, lj, lk;
    for (li = 0; li < (ltcs_count - 1); li++) {
        if (!BITTEST(m_notAnLtcs, li))
        {
            for (lj = (li+1); lj < ltcs_count; lj++) {
                if (!BITTEST(m_notAnLtcs, lj))
                {
                    int a = 0;
                    int b = 0;
                    for (lk = 0; lk < efm_count; lk++) {
                        if (BITTEST(ltcs[li],lk) && !BITTEST(ltcs[lj],lk))
                        {
                            a = 1;
                        }
                        else if (!BITTEST(ltcs[li],lk) && BITTEST(ltcs[lj],lk))
                        {
                            b = 1;
                        }
                        if (a > 0 && b > 0)
                        {
                            lk = efm_count;
                        }
                    }
                    // ltcs 2 is a subset of ltcs 1
                    if (a > 0 && b == 0)
                    {
                        BITSET(m_notAnLtcs, lj);
                        free(ltcs[lj]);
                    }
                    // ltcs 1 is a subset of ltcs 2
                    else if (a == 0 && b > 0)
                    {
                        BITSET(m_notAnLtcs, li);
                        free(ltcs[li]);
                        lj = ltcs_count;
                    }
                }
            }
        }
    }
    *notAnLtcs = m_notAnLtcs;
}

/*
 * calculate the cardinality of a LTCS
 */
unsigned long getLtcsCardinality(char* ltcs, unsigned long efm_count)
{
    unsigned long res = 0;
    unsigned long li;
    for (li = 0; li < efm_count; li++) {
        if (BITTEST(ltcs, li))
        {
            res++;
        }
    }
    return(res);
}

void performAnalysis(double*** result, char** ltcs, char** mat, char**
        reaction, unsigned long ltcs_count, unsigned long efm_count, unsigned
        int rx_count)
{
    unsigned long li, lj, lk;
    double** m_result = calloc(rx_count, sizeof(double*));
    for (li = 0; li < rx_count; li++) {
        m_result[li] = calloc(ltcs_count, sizeof(double));
    }
    unsigned long bitarray_size = getBitsize(rx_count * 2);
    unsigned long* counts = calloc(rx_count, sizeof(unsigned long));
    char* n_vect = calloc(1, bitarray_size);
    for (li = 0; li < ltcs_count; li++) {
        for (lj = 0; lj < bitarray_size; lj++) {
            BITCLEAR(counts, lj);
        }
        for (lj = 0; lj < rx_count; lj++) {
            counts[lj] = 0;
        }
        double c = (double) getLtcsCardinality(ltcs[li], efm_count);
        for (lj = 0; lj < efm_count; lj++) {
            if (BITTEST(ltcs[li], lj))
            {
                for (lk = 0; lk < rx_count; lk++) {
                    if (BITTEST(mat[lj], 2*lk))
                    {
                        counts[lk]++;
                    }
                    else if (BITTEST(mat[lj], 2*lk+1))
                    {
                        BITSET(n_vect, lk);
                        counts[lk]++;
                    }
                }
            }
        }
        for (lj = 0; lj < rx_count; lj++) {
            double val = (double)counts[lj]*100/c;
            if (BITTEST(n_vect, lj))
            {
                val *= -1;
            }
            m_result[lj][li] = val;
        }
    }
    free(n_vect);
    free(counts);
    *result = m_result;
}


/*
 * print efm matrix for debugging
 */
void printMatrix(char** matrix, unsigned long efm_count, unsigned int rx_count)
{
    unsigned long i;
    unsigned int j;
    unsigned int bits = rx_count * 2;
    for (i = 0; i < efm_count; i++)
    {
        for (j = 0; j < bits; j+=2)
        {
            if (BITTEST(matrix[i],j))
            {
                printf("  1");
                if (BITTEST(matrix[i], j+1))
                {
                    quitError("ERROR: fwd and rev is active in EFM\n", ERROR_EFM);
                }
            }
            else if (BITTEST(matrix[i],j+1))
            {
                printf(" -1");
            }
            else
            {
                printf("  0");
            }
        }
        printf("\n");
    }
}

/*
 * print a bit vector for debugging
 */
void printBitset (char* bs, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        if (BITTEST(bs, i))
        {
            printf("%d ", i);
        }
    }
    printf("\n");
}

/*
 * return number of reversible reactions
 */
unsigned int getReversibleReactionCount (char* reversible_reactions, unsigned int rx_count)
{
    unsigned int rev_rx_count = 0;
    unsigned int i;
    for (i = 0; i < rx_count; i++) 
    {
        if (BITTEST(reversible_reactions, i))
        {
            rev_rx_count++;
        }
    }
    return rev_rx_count;
}

/*
 * stores an EFM into a bitset if it is not an internal loop
 * every reaction is mapped to two bits
 * first bit is if reaction has a positive flux
 * second bit is set if reaction has a negative flux
 * none is set if reaction has no flux
 */
int loadEfm (char* matrix, char* line, char* reversible_reactions, char*
        rxs_fwd, char* rxs_rev, int rx_count, char* exchange_reaction, int
        checkLoops, double threshold)
{
    double n_thres = -1 * threshold;
    char *ptr;
    ptr = strtok(line, "\n\t ");
    int i = 0;
    int j = 0;
    while(ptr != NULL) 
    {
        double x = atof(ptr);
        if (checkLoops > 0 && BITTEST(exchange_reaction, i))
        {
            // check if EFM is an internal loop
            if (x >= threshold || x <= n_thres)
            {
                checkLoops = 0;
            }
        }
        if (BITTEST(reversible_reactions,i))
        {
            // positive flux
            if (x >= threshold)
            {
                BITSET(matrix, 2*j);
                BITSET(rxs_fwd,j);
            }
            // negative flux
            else if (x <= n_thres)
            {
                BITSET(matrix, (2*j)+1);
                BITSET(rxs_rev, j);
            }
            j++;
        }
        i++;
        ptr = strtok(NULL, "\n\t ");
    }
    if (checkLoops > 0)
    {
        // EFM is a loop - free memory
        free(matrix);
        matrix=NULL;
        return(0);
    }
    else
    {
        return(1);
    }
}

/*
 * increase memory size of efm matrix and loop vector
 */
void reallocInitialMatrix (char*** matrix, char** loops, unsigned long index)
{
    int alloc_size = 10000;
    unsigned long new_size = 0;
    if (index % alloc_size == 0)
    {
        new_size = index + alloc_size;
        *matrix = (char**) realloc(*matrix, new_size * sizeof(char*));
        unsigned long bitarray_size = getBitsize(new_size);
        *loops = (char*) realloc(*loops, bitarray_size);
        if (NULL == matrix || NULL == loops) 
        {
            quitError("Not enough free memory in reallocInitialMatrix\n", ERROR_RAM);
        }
    }
}

/*
 * allocate memory for an EFM
 */
void allocateEfm (char** vector, unsigned int reactions)
{
    unsigned long bitarray_size = getBitsize(2*reactions);
    *vector = calloc(1, bitarray_size);
    if (NULL == vector)
    {
        quitError("Not enough free memory\n", ERROR_RAM);
    }
}

/*
 * clean matrix
 * keep only those reactions that have a positive and a negative flux in any of
 * the EFMs
 */
void getCleanedMatrix(char **initial_mat, unsigned long efms, unsigned int
        init_rx_count, char ***cleaned_mat, char *rxs_both, unsigned int
        result_rx_count, char* loops) {
    unsigned int i;
    unsigned int j;
    unsigned long ul;
    char **m_cleaned_mat = calloc(efms, sizeof(char*));
    unsigned long bitarray_size = getBitsize(2 * result_rx_count);
    for (ul = 0; ul < efms; ul++) 
    {
        if (!BITTEST(loops, ul))
        {
            m_cleaned_mat[ul] = calloc(1, bitarray_size);
            if (NULL == m_cleaned_mat[ul])
            {
                quitError("Not enough free memory in getCleanedMatrix\n", ERROR_RAM);
            }
            j = 0;
            for (i = 0; i < init_rx_count; i++) 
            {
                if (BITTEST(rxs_both,i))
                {
                    if (BITTEST(initial_mat[ul], 2*i))
                    {
                        BITSET(m_cleaned_mat[ul], j);
                    }
                    j++;
                    if (BITTEST(initial_mat[ul], 2*i+1))
                    {
                        BITSET(m_cleaned_mat[ul], j);
                    }
                    j++;
                }
            }
        }
    }
    *cleaned_mat = m_cleaned_mat;
}

/*
 * read file containing reversibilty of reactions in form of:
 * 1 0 0 1 1
 * where 1 describes a reversible and 0 a non reversible reaction
 */
void readReversibleFile(char *filename, char *reversible_reactions, unsigned
        int rx_count)
{
    char* line = NULL;
    size_t len = 0;
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        quitError("Error in opening file\n", ERROR_FILE);
    }
    while ( getline(&line, &len, file) != -1)
    {
        char *ptr;
        ptr = strtok(line, "\n\t ");
        int i = 0;
        while(ptr != NULL) 
        {
            int x = atoi(ptr);
            if (x == 1)
            {
                BITSET(reversible_reactions, i);
            }
            else if (x == 0)
            {
                BITCLEAR(reversible_reactions, i);
            }
            else
            {
                quitError("not a correct reversibility file\n", ERROR_FILE);
            }
            i++;
            ptr = strtok(NULL, "\n\t ");
        }
        if (i != rx_count)
        {
            quitError("not a correct reversibility file\n", ERROR_FILE);
        }
    }
    fclose(file);
    free(line);
}

void readReactions(char *filename, char*** rx_names, unsigned int rx_count)
{
    char* line = NULL;
    size_t len = 0;
    FILE *file = fopen(filename, "r");
    char** m_rx_names = calloc(rx_count, sizeof(char*));
    if (!file)
    {
        quitError("Error in opening file\n", ERROR_FILE);
    }
    if ( getline(&line, &len, file) != -1)
    {
        char *ptr;
        ptr = strtok(line, "\"\n\t ");
        int i = 0;
        while(ptr != NULL) 
        {
            if (i >= rx_count)
            {
                quitError("not correct reaction file: more reactions than in EFM file\n", ERROR_FILE);
            }
            size_t len = strlen(ptr);
            m_rx_names[i] = calloc(len + 1, sizeof(char));
            strcpy(m_rx_names[i], ptr);
            ptr = strtok(NULL, "\"\n\t ");
            i++;
        }
        if (i != rx_count)
        {
            printf("Number of reactions in EFM file: %d\n", rx_count);
            printf("Number of reactions in reaction file: %d\n", i);
            quitError("not a correct reaction file\n", ERROR_FILE);
        }
    }
    fclose(file);
    free(line);
    *rx_names = m_rx_names;
}

/*
 * read stoichiometric matrix file and searches exchange reactions
 * exchange reactions are defined as those reactions that have only positive
 * or only negative values
 */
void defineExchangeReactions(int *cols, char** exchange_reaction, double
        threshold, char* filename)
{
    double n_thres = -1 * threshold;
    FILE *rx_file = fopen(filename, "r");
    if (!rx_file)
    {
        quitError("Error in opening file\n", ERROR_FILE);
    }
    int rx_count = getRxCount(rx_file);
    fclose(rx_file);
    char* line = NULL;
    size_t len = 0;
    unsigned long bitarray_size = getBitsize(rx_count);
    char* m_exchange_reaction = calloc(1, bitarray_size);
    char* m_pos = calloc(1, bitarray_size);
    char* m_neg = calloc(1, bitarray_size);
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        quitError("Error in opening file\n", ERROR_FILE);
    }
    while ( getline(&line, &len, file) != -1)
    {
        char *ptr;
        ptr = strtok(line, "\n\t ");
        int i = 0;
        while(ptr != NULL) 
        {
            double val = atof(ptr);
            // positive value
            if (val > threshold)
            {
                BITSET(m_pos, i);
            }
            // negative value
            else if (val < n_thres)
            {
                BITSET(m_neg,i);
            }
            i++;
            ptr = strtok(NULL, "\n\t ");
        }
    }
    free(line);
    line = NULL;
    fclose(file);
    int i;
    for (i = 0; i < rx_count; i++) {
        // check if reaction have positive and negative values
        if ( (BITTEST(m_pos,i) && !BITTEST(m_neg,i)) || (!BITTEST(m_pos,i) && BITTEST(m_neg,i)) )
        {
            BITSET(m_exchange_reaction, i);
        }
    }
    free(m_pos);
    free(m_neg);
    *cols = rx_count;
    *exchange_reaction = m_exchange_reaction;
}

/*
 * read EFM file
 * store EFMs in bit vectors
 * clean matrix by removing reactions that have a flux in only one direction
 * store internal loops in EFMs
 */
void readInitialMatrix(int rx_count, unsigned long* efm_count, char*
        reversible_reactions, unsigned int* rev_rx_count, char***
        initial_mat, char*** full_mat, int keep_full_mat, char*
        exchange_reaction, char** loops, int checkLoops, double threshold, char
        *filename)
{
    printf("%s loading EFMs\n", getTime());
    char** m_initial_mat = NULL;
    char* m_loops = NULL;
    char* line = NULL;
    size_t len = 0;
    unsigned long mat_ix = 0;

    unsigned int m_rev_rx_count =
        getReversibleReactionCount(reversible_reactions, rx_count);
    unsigned long rxs_bitarray_size = getBitsize(m_rev_rx_count);
    char* rxs_fwd  = calloc(1, rxs_bitarray_size);
    char* rxs_rev  = calloc(1, rxs_bitarray_size);
    char* rxs_both = calloc(1, rxs_bitarray_size);

    // read EFMs
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        quitError("Error in opening file\n", ERROR_FILE);
    }
    while ( getline(&line, &len, file) != -1)
    {
        // allocate memory and load EFM
        reallocInitialMatrix(&m_initial_mat, &m_loops, mat_ix);
        allocateEfm(&m_initial_mat[mat_ix], m_rev_rx_count);
        int loaded = loadEfm(m_initial_mat[mat_ix], line, reversible_reactions,
                rxs_fwd, rxs_rev, rx_count, exchange_reaction, checkLoops,
                threshold);
        // check if EFM is an internal loop
        if (loaded > 0)
        {
            BITCLEAR(m_loops, mat_ix);
        }
        else
        {
            BITSET(m_loops, mat_ix);
        }
        mat_ix++;
    }
    free(line);
    line = NULL;
    fclose(file);
    unsigned long bitarray_size = getBitsize(mat_ix);
    m_loops = (char*) realloc(m_loops, bitarray_size);

    // define reactions that have a flux in both directions
    unsigned int i = 0;
    unsigned int result_rx_count = 0;
    for (i = 0; i < m_rev_rx_count; i++) 
    {
        if (BITTEST(rxs_fwd,i) && BITTEST(rxs_rev,i))
        {
            BITSET(rxs_both,i);
            result_rx_count++;
        }
    }
    free(rxs_fwd);
    free(rxs_rev);
    rxs_fwd = NULL;
    rxs_rev = NULL;

    // clean matrix by removing reactions that are active in only one direction
    printf("%s optimizating EFM matrix\n", getTime());
    getCleanedMatrix(m_initial_mat, mat_ix, m_rev_rx_count, initial_mat,
            rxs_both, result_rx_count, m_loops);
    if (keep_full_mat < 1)
    {
        unsigned long ul;
        for (ul = 0; ul < mat_ix; ul++) 
        {
            if (!BITTEST(m_loops,ul))
            {
                free(m_initial_mat[ul]);
            }
        }
        free(m_initial_mat);
        m_initial_mat = NULL;
    }
    free(rxs_both);
    rxs_both = NULL;

    // set return pointer
    *efm_count = mat_ix;
    *rev_rx_count = result_rx_count;
    *loops = m_loops;
    *full_mat = m_initial_mat;
}

/*
 * thread function to find new ltcs
 * splits an ltcs into two ltcs based on positive and negative flux values in
 * EFMs at given reaction index
 */
void* findLtcsThread(void *pointer_thread_args)
{
    // unpack given arguments
    struct thread_args* thread_args = (struct thread_args*) pointer_thread_args;
    int thread_id                 = thread_args->thread_id;
    int max_threads               = thread_args->max_threads;
    int bitarray_size             = thread_args->bitarray_size;
    unsigned long ltcs_count      = thread_args->ltcs_count;
    unsigned long efm_count       = thread_args->efm_count;
    unsigned int reaction         = thread_args->reaction;
    char** ltcs                   = thread_args->ltcs;
    char** matrix                 = thread_args->matrix;
    unsigned long* t_ltcs_count   = thread_args->new_ltcs_count;
    char** t_ltcs                 = thread_args->new_ltcs;

    t_ltcs_count[0] = 0;
    if (NULL == t_ltcs)
    {
        quitError("Not enough free memory in findLtcs", ERROR_RAM);
    }
    unsigned long ui, uj;
    for (ui = 0; ui < ltcs_count; ui++) 
    {
        // check if thread is responsible
        if (ui%max_threads == thread_id)
        {
            int pos = 0;
            int neg = 0;
            char *p_vect = (char *)calloc(1, bitarray_size);
            char *n_vect = (char *)calloc(1, bitarray_size);
            for (uj = 0; uj < efm_count; uj++) 
            {
                // if efm is in ltcs
                if (BITTEST(ltcs[ui],uj))
                {
                    // efms with positive flux to new ltcs 1
                    if (BITTEST(matrix[uj],2*reaction))
                    {
                        BITSET(p_vect, uj);
                        pos = 1;
                    } 
                    // efms with negative flux to new ltcs 2
                    else if (BITTEST(matrix[uj], 2*reaction+1))
                    {
                        BITSET(n_vect, uj);
                        neg = 1;
                    } 
                    // efms with zero flux to both new ltcs
                    else
                    {
                        BITSET(p_vect, uj);
                        BITSET(n_vect, uj);
                    }
                }
            }
            // create new ltcs 1 if EFMs with positive fluxes were found, else
            // set memory free
            if (pos > 0)
            {
                t_ltcs[t_ltcs_count[0]] = p_vect;
                t_ltcs_count[0]++;
            }
            else
            {
                free(p_vect);
            }
            // create new ltcs 2 if EFMs with negative fluxes were found, else
            // set memory free
            if (neg > 0)
            {
                t_ltcs[t_ltcs_count[0]] = n_vect;
                t_ltcs_count[0]++;
            }
            else
            {
                free(n_vect);
            }
        }
    }
    return((void*)NULL);
}

/*
 * find ltcs
 * for each reaction start threads to find new ltcs
 * merge ltcs found by threads to new set of ltcs until all reactions are
 * processed
 */
void findLtcs(unsigned int rx_count, unsigned long efm_count, char **mat, char*
        loops, char ***ltcs, unsigned long *ltcs_count, int max_threads)
{
    unsigned int i;
    unsigned long ui;
    unsigned long bitarray_size = getBitsize(efm_count);
    // initialize ltcs with 1 ltcs containing all EFMs that are not internal
    // loops
    unsigned long m_ltcs_count = 1;
    char **m_ltcs = (char**) calloc(1, sizeof(char*));
    m_ltcs[0] = calloc(1, bitarray_size);
    for (ui = 0; ui < efm_count; ui++) {
        if (!BITTEST(loops, ui))
        {
            BITSET(m_ltcs[0], ui);
        }
    }
    if (NULL == m_ltcs)
    {
        quitError("Not enough free memory in findLtcs", ERROR_RAM);
    }

    // process each single reaction
    for (i = 0; i < rx_count; i++) 
    {
        printf("%s ", getTime());
        printf("iteration %d/%d: ", i+1, rx_count);

        // prepare threads
        int num_threads = (m_ltcs_count > max_threads) ? max_threads : m_ltcs_count;
        pthread_t thread[num_threads];
        struct thread_args thread_args[num_threads];
        int ti;
        for (ti = 0; ti < num_threads; ti++) {
            thread_args[ti].thread_id = ti;
            thread_args[ti].max_threads = num_threads;
            thread_args[ti].bitarray_size = bitarray_size;
            thread_args[ti].ltcs_count = m_ltcs_count;
            thread_args[ti].efm_count = efm_count;
            thread_args[ti].reaction = i;
            thread_args[ti].ltcs = m_ltcs;
            thread_args[ti].matrix = mat;
            thread_args[ti].new_ltcs_count = calloc(1, sizeof(unsigned long));
            unsigned long size = m_ltcs_count * 2 / num_threads + num_threads;
            char** new_ltcs = (char**) calloc(size, sizeof(char*));
            thread_args[ti].new_ltcs = new_ltcs;
        }

        // start threads
        for (ti = 0; ti < num_threads; ti++) {
            pthread_create(&thread[ti], NULL, findLtcsThread,
                    (void*)&thread_args[ti]);
        }

        // join threads
        unsigned long new_ltcs_count = 0;
        for (ti = 0; ti < num_threads; ti++) {
            pthread_join(thread[ti], NULL);
            new_ltcs_count += *thread_args[ti].new_ltcs_count;
        }

        // prepare memory for new ltcs sets
        for (ui = 0; ui < m_ltcs_count; ui++) {
            free(m_ltcs[ui]);
        }
        free(m_ltcs);
        m_ltcs = NULL;
        m_ltcs = calloc(new_ltcs_count, sizeof(char*));
        unsigned long index = 0;
        unsigned long tj;
        // merge ltcs found by threads to new set
        for (ti = 0; ti < num_threads; ti++) {
            for (tj = 0; tj < *thread_args[ti].new_ltcs_count; tj++) {
                m_ltcs[index] = thread_args[ti].new_ltcs[tj];
                index++;
            }
            free(thread_args[ti].new_ltcs);
            free(thread_args[ti].new_ltcs_count);
        }
        m_ltcs_count = new_ltcs_count;
        printf("%lu ltcs\n", m_ltcs_count);
    }
    *ltcs = m_ltcs;
    *ltcs_count = m_ltcs_count;
}


int main (int argc, char *argv[])
{
    printf("Start: %s\n", getTime());

    //================================================== 
    // define arguments and usage
    char *optv[MAX_ARGS] = { "-i", "-o", "-s", "-r", "-v", "-l", "-z", "-t", "-f", "-c", "-a" };
    char *optd[MAX_ARGS] = {"efm file  (tab separated like:  0.4\t0\t-0.24)",
                            "output file [default: ltcs.out]",
                            "stoichiometric matrix file [optional, needed to find internal loops]",
                            "reaction file [optional, but necessary if option -a is set]",
                            "reversibility file [optional, but saves memory!]",
                            "loops output [only if stoichiometric matrix is given, default: loops.out]",
                            "zero threshold [default: 1e-10]",
                            "number of threads [default: 1]",
                            "full output [yes/no; default: yes] if set to no only summary is printed and ltcs are not saved",
                            "print ltcs in csv format [yes/no; default: yes] if set to no ltcs output will be e.g. 10011 instead of 1,0,0,1,1",
                            "analysis output file - needs option -r"};
    char *optr[MAX_ARGS];
    char *description = "Calculate largest thermodynamically consistent sets of "
        "EFMs\nbased only on the reversibility of the reactions";
    char *usg = "\n   calcLtcs -i efms.txt -o ltcs.out\n   calcLtcs -i efms.txt -o ltcs.out -s sfile -v rvfile -z 1e-6 -t 8 -a yes -r rfile";
    // end define arguments and usage
    //================================================== 

    //================================================== 
    // read arguments
    readArgs(argc, argv, MAX_ARGS, optv, optr);
    if ( !optr[ARG_INPUT] )
    {
        usage(description, usg, MAX_ARGS, optv, optd);
        quitError("Missing argument\n", ERROR_ARGS);
    }
    char* ltcsout = optr[ARG_LTCS_OUT] ? optr[ARG_LTCS_OUT] : "ltcs.out";
    char* arg_rvfile = optr[ARG_RVFILE] ? optr[ARG_RVFILE] : "not available";
    char* arg_sfile = optr[ARG_SFILE] ? optr[ARG_SFILE] : "not available";
    char* arg_rfile = optr[ARG_RFILE] ? optr[ARG_RFILE] : "not available";
    char* loopout = optr[ARG_SFILE] ? (optr[ARG_LOOPS] ? optr[ARG_LOOPS] : "loops.out") : "not available";
    double threshold = optr[ARG_ZERO] ? atof(optr[ARG_ZERO]) : 1e-10;
    int threads = optr[ARG_THREADS] ? atoi(optr[ARG_THREADS]) : 1;
    int full_out = optr[ARG_FULL_OUT] ? (!strcmp(optr[ARG_FULL_OUT], "no") ? 0 : 1) : 1;
    int csv_out = optr[ARG_CSV_OUT] ? (!strcmp(optr[ARG_CSV_OUT], "no") ? 0 : 1) : 1;
    int arg_analysis = optr[ARG_ANALYSIS] ? (optr[ARG_RFILE] ? 1 : 0) : 0;
    char* arg_analysisfile = optr[ARG_ANALYSIS] ? optr[ARG_ANALYSIS] : "not available";
    // end read arguments
    //================================================== 

    //================================================== 
    // print arguments summary
    printf("\n");
    printf("Input:            %s\n", optr[ARG_INPUT]);
    printf("Output:           %s\n", ltcsout);
    printf("sfile:            %s\n", arg_sfile);
    printf("rvfile:           %s\n", arg_rvfile);
    printf("rfile:            %s\n", arg_rfile);
    printf("Zero threshold:   %.2e\n", threshold);
    printf("Threads:          %d\n", threads);
    printf("Loops output:     %s\n", full_out > 0 ? loopout : "no");
    printf("Full output:      %s\n", full_out > 0 ? "yes" : "no");
    if (full_out > 0)
    {
        printf("csv output:       %s\n", csv_out > 0 ? "yes (1,0,1,1)" : "no (1011)");
    }
    printf("Perform analysis: %s\n", arg_analysis > 0 ? "yes" : "no");
    if (arg_analysis > 0)
    {
        printf("analysis file:    %s\n", arg_analysisfile);
    }
    printf("\n");
    // end print arguments summary
    //================================================== 

    //================================================== 
    // check files
    FILE *file = fopen(optr[ARG_INPUT], "r");
    FILE *fileout = NULL;
    FILE *fileloops = NULL;
    FILE *fileanalysis = NULL;
    if (!file)
    {
        quitError("Error in opening input file\n", ERROR_FILE);
    }
    if (full_out > 0)
    {
        fileout = fopen(ltcsout, "w");
        if (!fileout)
        {
            quitError("Error in opening output file\n", ERROR_FILE);
        }
        if (optr[ARG_SFILE])
        {
            fileloops = fopen(loopout, "w");
            if (!fileloops)
            {
                quitError("Error in opening loop outputfile\n", ERROR_FILE);
            }
        }
    }
    if (arg_analysis > 0)
    {
        fileanalysis = fopen(arg_analysisfile, "w");
        if (!fileanalysis)
        {
            quitError("Error in opening analysis file\n", ERROR_FILE);
        }
    }
    // end check files
    //================================================== 
    
    //================================================== 
    // read number of reactions
    int rx_count = getRxCount(file);
    fclose(file);
    if (rx_count < 1)
    {
        quitError("Error in EFM file format; number of reactions < 1\n",
                ERROR_FILE);
    }
    // end read number of reactions
    //================================================== 

    //================================================== 
    // define exchange reactions
    char* exchange_reaction = NULL;
    int s_cols = 0;
    if (optr[ARG_SFILE])
    {
        defineExchangeReactions(&s_cols, &exchange_reaction, threshold,
                optr[ARG_SFILE]);
        if (s_cols != rx_count)
        {
            quitError("Error in EFM of stoichiometric file format. Number of reactions is not equal\n", ERROR_FILE);
        }
    }
    // end define exchange reactions
    //================================================== 

    //================================================== 
    // define reversible reactions
    unsigned int   rev_rx_count  = 0;
    unsigned long rev_rx_bitsize = getBitsize(rx_count);
    char* reversible_reactions = calloc(1, rev_rx_bitsize);
    int i;
    for (i = 0; i < rx_count; i++) 
    {
        BITSET(reversible_reactions, i);
    }
    if (optr[ARG_RVFILE])
    {
        readReversibleFile(optr[ARG_RVFILE], reversible_reactions, rx_count);
    }
    else
    {
        for (i = 0; i < rx_count; i++) 
        {
            BITSET(reversible_reactions, i);
        }
    }
    // end define reversible reactions
    //================================================== 
    
    //================================================== 
    // read reaction names
    char** reaction_names = NULL;
    if (arg_analysis > 0)
    {
        readReactions(optr[ARG_RFILE], &reaction_names, rx_count);
    }
    //
    //================================================== 

    //================================================== 
    // read efm matrix
    char* loops = NULL;
    char** initial_mat = NULL;
    char** full_mat = NULL;
    unsigned long efm_count = 0;
    int checkLoops = optr[ARG_SFILE] ? 1 : 0;
    readInitialMatrix(rx_count, &efm_count, reversible_reactions,
            &rev_rx_count, &initial_mat, &full_mat, arg_analysis,
            exchange_reaction, &loops, checkLoops, 1e-8, optr[ARG_INPUT]);
    // end read efm matrix
    //================================================== 

    //================================================== 
    // find ltcs
    unsigned long ltcs_count = 0;
    char** ltcs = NULL;
    findLtcs(rev_rx_count, efm_count, initial_mat, loops, &ltcs, &ltcs_count, threads);
    // end find ltcs
    //================================================== 
    
    //================================================== 
    // filter ltcs to remove subsets
    printf("%s filter LTCS to remove subsets\n", getTime());
    char* notAnLtcs = NULL;
    filterLtcs(ltcs, &notAnLtcs, ltcs_count, efm_count);
    // end filter ltcs to remove subsets
    //================================================== 

    //================================================== 
    // print ltcs to file
    if (full_out > 0)
    {
        printf("%s save ltcs\n", getTime());
    }
    unsigned long li;
    unsigned long ul;
    unsigned long result_ltcs_count = 0;
    int makeSep = 0;
    for (li = 0; li < ltcs_count; li++) {
        if (!BITTEST(notAnLtcs, li))
        {
            makeSep = 0;
            result_ltcs_count++;
            if (full_out > 0)
            {
                for (ul = 0; ul < efm_count; ul++) {
                    if (makeSep > 0 && csv_out > 0)
                    {
                        fprintf(fileout, ",");
                    }
                    makeSep = 1;
                    if (BITTEST(ltcs[li],ul))
                    {
                        fprintf(fileout, "1");
                    }
                    else
                    {
                        fprintf(fileout, "0");
                    }
                }
                fprintf(fileout, "\n");
            }
        }
    }
    if (full_out > 0)
    {
        fclose(fileout);
    }
    // end print ltcs to file
    //================================================== 

    //================================================== 
    // perform analysis
    if (arg_analysis > 0)
    {
        printf("%s perform and save analysis of LTCS\n", getTime());
        double** analysis_values = NULL;
        performAnalysis(&analysis_values, ltcs, full_mat, reaction_names,
                ltcs_count, efm_count, rx_count);
        unsigned long lj;
        for (li = 0; li < rx_count; li++) {
            fprintf(fileanalysis, "%s", reaction_names[li]);
            for (lj = 0; lj < ltcs_count; lj++) {
                fprintf(fileanalysis, ",%.2f", analysis_values[li][lj]);
            }
            fprintf(fileanalysis, "\n");
        }
        fclose(fileanalysis);

        for (li = 0; li < rx_count; li++) {
            free(analysis_values[li]);
            free(reaction_names[li]);
        }
        free(analysis_values);
        free(reaction_names);
    }
    // end perform analysis
    //================================================== 
    
    //================================================== 
    // set initial_mat memory free
    for (li=0; li<efm_count; li++)
    {
        if (!BITTEST(loops, li))
        {
            free(initial_mat[li]);
        }
    }
    free(initial_mat);
    initial_mat = NULL;
    if (arg_analysis > 0)
    {
        for (li=0; li<efm_count; li++)
        {
            if (!BITTEST(loops, li))
            {
                free(full_mat[li]);
            }
        }
        free(full_mat);
        full_mat = NULL;
    }
    // end set initial_mat memory free
    //================================================== 

    //================================================== 
    // count and print loops to file
    unsigned long loop_count = 0;
    if (optr[ARG_SFILE])
    {
        if (full_out > 0)
        {
            printf("%s save internal loops\n", getTime());
            for (ul = 0; ul < efm_count; ul++) {
                if (ul > 0 && csv_out > 0)
                {
                    fprintf(fileloops, ",");
                }
                if (BITTEST(loops, ul))
                {
                    fprintf(fileloops, "1");
                    loop_count++;
                }
                else
                {
                    fprintf(fileloops, "0");
                }
            }
            fclose(fileloops);
        }
        else
        {
            for (ul = 0; ul < efm_count; ul++) {
                if (BITTEST(loops, ul))
                {
                    loop_count++;
                }
            }
        }
    }
    // end count and print loops to file
    //================================================== 

    //================================================== 
    // print summary
    printf("\n");
    printf("Nr of EFMS:           %lu\n", efm_count);
    printf("Nr of LTCS:           %lu\n", result_ltcs_count);
    if (optr[ARG_SFILE])
    {
        printf("Nr of internal loops: %lu\n", loop_count);
    }
    printf("\nSizes of LTCS:\n");
    for (li = 0; li < ltcs_count; li++) {
        if (!BITTEST(notAnLtcs, li))
        {
            if (li > 0)
            {
                printf(",");
            }
            printf("%lu", getLtcsCardinality(ltcs[li], efm_count));
        }
    }
    printf("\n\n");
    printf("End: %s\n", getTime());
    // end print summary
    //================================================== 

    //================================================== 
    // set memory free
    free(loops);
    free(exchange_reaction);
    loops = NULL;
    for (li = 0; li < ltcs_count; li++) {
        if (!BITTEST(notAnLtcs, li))
        {
            free(ltcs[li]);
        }
    }
    free(ltcs);
    free(notAnLtcs);
    ltcs = NULL;
    free(reversible_reactions);
    reversible_reactions = NULL;
    // end set memory free
    //================================================== 

    return EXIT_SUCCESS;
}
