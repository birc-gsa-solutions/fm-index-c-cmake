#include "fasta.h"
#include "fastq.h"
#include "sam.h"
#include <cstr.h>
#include <stdio.h>
#include <string.h>

static const char *suffix = "fmindex";
static char *gen_preproc_fname(const char *fasta_fname)
{
    // length of fasta file name, '.', the suffix, + '\0'
    char *buf = cstr_malloc(strlen(fasta_fname) + 1 + strlen(suffix) + 1);
    sprintf(buf, "%s.%s", fasta_fname, suffix);
    return buf;
}

static void write_name(FILE *f, const char *x)
{
    size_t len = strlen(x);
    fwrite(&len, sizeof len, 1, f);
    fwrite(x, len, 1, f);
}

static char *read_name(FILE *f)
{
    size_t len;
    fread(&len, sizeof len, 1, f);
    char *buf = cstr_malloc(len + 1);
    fread(buf, len, 1, f);
    buf[len] = '\0';
    return buf;
}

// Check if we are at the end of file before we attempt
// to read more.
static inline bool more_data(FILE *f)
{
    int c = fgetc(f);
    if (c == EOF)
        return false;
    ungetc(c, f);
    return true;
}

static void preprocess_genome(const char *genome_fname)
{
    struct fasta_records *genome = load_fasta_records(genome_fname);
    if (!genome)
    {
        perror("Could not open fasta file");
        exit(EXIT_FAILURE);
    }

    char *preproc_fname = gen_preproc_fname(genome_fname);
    FILE *f = fopen(preproc_fname, "w");
    free(preproc_fname);

    if (!f)
    {
        perror("Could not open preprocessing file");
        exit(EXIT_FAILURE);
    }
    struct fm_index *fm_index = 0;
    for (struct fasta_record *fa_rec = fasta_records(genome);
         fa_rec;
         fa_rec = fa_rec->next)
    {
        write_name(f, fa_rec->name);
        cstr_bwt_preproc *tables = cstr_bwt_preprocess(fa_rec->seq);
        cstr_write_bwt_tables(f, tables);
        cstr_free_bwt_preproc(tables);
    }

    fclose(f);
}

struct fm_index
{
    const char *chr_name;     // We don't own this
    cstr_bwt_preproc *tables; // We do own this (and must free it)
    struct fm_index *next;    // Link to next chromosome.
};

static struct fm_index *read_preprocessed_genome(const char *genome_fname)
{
    char *preproc_fname = gen_preproc_fname(genome_fname);
    FILE *f = fopen(preproc_fname, "r");
    free(preproc_fname);

    if (!f)
    {
        perror("Could not open preprocessing file");
        exit(EXIT_FAILURE);
    }

    struct fm_index *fm_index = 0;
    while (more_data(f))
    {
        struct fm_index *new_fm_index = cstr_malloc(sizeof *new_fm_index);
        new_fm_index->next = fm_index;
        fm_index = new_fm_index;

        fm_index->chr_name = read_name(f);
        fm_index->tables = cstr_read_bwt_tables(f);
    }
    return fm_index;
}

static void free_fm_index(struct fm_index *fm_index)
{
    struct fm_index *next;
    for (; fm_index; fm_index = next)
    {
        next = fm_index->next;
        cstr_free_bwt_preproc(fm_index->tables);
        free(fm_index);
    }
}

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr,
                "Usage: %s -p genome\n       %s genome reads\n",
                argv[0], argv[0]);
        return 1;
    }

    if (strcmp("-p", argv[1]) == 0)
    {
        // preprocessing
        const char *genome_fname = argv[2];
        fprintf(stderr, "Preprocessing genome '%s'\n", genome_fname);
        char *preproc_fname = gen_preproc_fname(genome_fname);
        preprocess_genome(genome_fname);
        free(preproc_fname);
    }
    else
    {
        // mapping
        const char *genome_fname = argv[1];
        const char *reads_fname = argv[2];
        fprintf(stderr, "Mapping in genome '%s' for reads in '%s'\n", genome_fname, reads_fname);

        fprintf(stderr, "Reading preprocesed genome\n");
        struct fm_index *chromosomes = read_preprocessed_genome(genome_fname);

        fprintf(stderr, "Processing genome.\n");
        struct fastq_iter fqiter;
        struct fastq_record fqrec;
        static const size_t CIGAR_BUF_SIZE = 2048;
        char cigarbuf[CIGAR_BUF_SIZE];

        FILE *fq = fopen(reads_fname, "r");
        if (!fq)
        {
            abort(); // I can always implement better checking another day...
        }
        init_fastq_iter(&fqiter, fq);

        while (next_fastq_record(&fqiter, &fqrec))
        {
            snprintf(cigarbuf, CIGAR_BUF_SIZE, "%lldM", fqrec.seq.len);

            for (struct fm_index *fm_index = chromosomes; fm_index; fm_index = fm_index->next)
            {
                cstr_exact_matcher *m = cstr_fmindex_search(fm_index->tables, fqrec.seq);
                for (long long i = cstr_exact_next_match(m); i != -1; i = cstr_exact_next_match(m))
                {
                    print_sam_line(stdout, (const char *)fqrec.name.buf,
                                   fm_index->chr_name, i, cigarbuf,
                                   (const char *)fqrec.seq.buf);
                }
                cstr_free_exact_matcher(m);
            }
        }

        fclose(fq);
        free_fm_index(chromosomes);
    }
    return 0;
}

/*

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "%s genome reads", argv[0]);
        return 1;
    }
    const char *genome_fname = argv[1];
    const char *reads_fname = argv[2];

    // Initialise global alphabet.
    // NB! This assumes that we know the alphabet already and that it is "acgt" + $.
    // If we don't know it, we need to obtain it from the chromosomes.
    cstr_const_sslice letters = CSTR_SLICE_STRING0((const char *)"acgt");
    cstr_init_alphabet(&ALPHA, letters);

    struct fasta_records *genome = load_fasta_records(genome_fname);
    struct sa *suf_arrays = 0;
    for (struct fasta_record *fa_rec = fasta_records(genome); fa_rec; fa_rec = fa_rec->next)
    {
        suf_arrays = process_fasta_rec(fa_rec, suf_arrays);
    }

    struct fastq_iter fqiter;
    struct fastq_record fqrec;
    static const size_t CIGAR_BUF_SIZE = 2048;
    char cigarbuf[CIGAR_BUF_SIZE];

    FILE *fq = fopen(reads_fname, "r");
    if (!fq)
    {
        abort(); // I can always implement better checking another day...
    }
    init_fastq_iter(&fqiter, fq);

    while (next_fastq_record(&fqiter, &fqrec))
    {
        snprintf(cigarbuf, CIGAR_BUF_SIZE, "%lldM", fqrec.seq.len);

        for (struct sa *sa = suf_arrays; sa; sa = sa->next)
        {
            cstr_exact_matcher *m = cstr_sa_bsearch(*sa->sa, sa->x, fqrec.seq);
            for (long long i = cstr_exact_next_match(m); i != -1; i = cstr_exact_next_match(m))
            {
                print_sam_line(stdout, (const char *)fqrec.name.buf,
                               sa->chr_name, i, cigarbuf,
                               (const char *)fqrec.seq.buf);
            }
            cstr_free_exact_matcher(m);
        }
    }

    fclose(fq);
    free_sa(suf_arrays);
    free_fasta_records(genome);

    return 0;
}
*/