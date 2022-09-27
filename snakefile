"""
Right now the results are different every time the pipeline is run. 
This is because gargammel fragSim creates different DNA fragments each 
time, it is not possible to set a seed for this step. Therefore, when
publishing the results, the fragments should be made publicly available.

AdapterRemoval uses a seed when choosing between bases with equal Phred
scores when merging the reads. The other trimming methods don't specify
what happens in such a case.

TO DO
- add gzip
- delete unneccessary variables
- Check if the merged output file really contains only merged reads by
  comparing with leonardos files (header of the reads):
  bbmerge, adna-trim 

Unclear Phred Score offset of adna-trim, bbmerge
"""


### Set global variables


SEED = 2718

# fragment variables
LENGTHS = [100]  # list(range(1, 250+1)) + [1000]
NUMFRAGS = 2  # 1000000

# paths for output directories
OUTDIR = "/net/node07/home/projects/DNA_reconstruct/output/"
OUTDIR_BEN = OUTDIR + "benchmarks/"
OUTDIR_SIM = OUTDIR + "simulations/"
OUTDIR_REC = OUTDIR + "reconstructions/"
OUTDIR_EVA = OUTDIR + "evaluation/"

# paths to tools
FRAGSIM = "/home/ctools/gargammel/src/fragSim"
ADPTSIM = "/home/ctools/gargammel/src/adptSim"
ART = "/home/ctools/gargammel/art_src_MountRainier_Linux/art_illumina"

LEEHOM = "/home/ctools/leeHom-1.2.15/src/leeHom"
ADPTREM = "/home/ctools/adapterremoval-2.3.2/build/AdapterRemoval"
CLIPMERGE = "/home/ctools/ClipAndMerge-1.7.8/build/libs/ClipAndMerge-1.7.8.jar"
SEQTK = "/home/ctools/seqtk-1.3/seqtk"
ADNA = "/home/ctools/adna/adna-trim"
BBMERGE = "/home/ctools/bbmap_38_91/bbmerge.sh"
FASTP = "/home/ctools/fastp/fastp"
SEQPREP = "/home/ctools/SeqPrep-1.3.2/SeqPrep"

# path to genome and adapter sequences
IN_FILE = "/home/databases/genomes/Homo_sapiens/CHM13_T2T/CHM13_T2T.fa"
ADPT_FILE = "/home/projects/DNA_reconstruct/gabrieldir/adapters.fa"

# adapter sequences
ADPT1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG"
ADPT2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT"


### Run all


rule all:
    input:
        expand(
            OUTDIR_EVA + "leeHom/gen_n{n}_l{l}_leehom.csv",
            n=NUMFRAGS,
            l=LENGTHS,
        ),


### Simulation


rule simulate_fragments:
    """
    gargammel fragSim:
    simulation of ancient DNA fragments being retrieved at random from 
    the genome.
    """
    input:
        IN_FILE,
    output:
        OUTDIR_SIM + "gen_n{n}_l{l}.fa",
    benchmark:
        OUTDIR_BEN + "simulate_fragments/gen_n{n}_l{l}.tsv"
    wildcard_constraints:
        n="\d+",
        l="\d+",
    shell:
        "{FRAGSIM} -n {wildcards.n} -l {wildcards.l} {input} > {output}"


rule add_adapters:
    """
    gargammel adptSim:
    adding of adapters to create raw Illumina reads (without errors and
    quality scores)

    Output reads as ART (unzipped fasta) with wrap-around for paired-end
    mode.
    """
    input:
        OUTDIR_SIM + "gen_n{n}_l{l}.fa",
    output:
        OUTDIR_SIM + "gen_n{n}_l{l}_adpt.fa",
    benchmark:
        OUTDIR_BEN + "add_adapters/gen_n{n}_l{l}_adpt.tsv"
    wildcard_constraints:
        n="\d+",
        l="\d+",
    shell:
        # -l        : Desired read length
        # -artp     : Output reads as ART with wrap-around
        ("{ADPTSIM}" " -l 125" " -artp" " {output} {input}")


rule simulate_reads:
    """
    add sequencing errors and corresponding quality scores
    """
    input:
        OUTDIR_SIM + "gen_n{n}_l{l}_adpt.fa",
    output:
        OUTDIR_SIM + "gen_n{n}_l{l}_s1.fq",
        OUTDIR_SIM + "gen_n{n}_l{l}_s2.fq",
    params:
        out_prefix=OUTDIR_SIM + "gen_n{n}_l{l}_s",
    benchmark:
        OUTDIR_BEN + "simulate_reads/gen_n{n}_l{l}_reads.tsv"
    wildcard_constraints:
        n="\d+",
        l="\d+",
    shell:
        # --insRate     : insertion rate
        # -dr           : deletion rate
        # --seqSys      : HS25: Illumina HiSeq 2500 (125bp, 150bp)
        # --len         : read length
        # --rcount      : number of read pairs to be generated per sequence
        # --paired      : paired-end read simulation
        # --amplicon    : amplicon sequencing simulation
        # --noALN       : do not output ALN alignment file
        # --quiet       : turn off end of run summary
        # --rndSeed     : the seed for random number generator, use a
        #                 fixed seed to generate two identical datasets
        #                 from different runs
        (
            "{ART}"
            " --insRate 0 --insRate2 0"
            " -dr 0 -dr2 0"
            " --seqSys HS25"
            " --len 125"
            " --rcount 1"
            " --paired"
            " --amplicon"
            " --noALN"
            " --quiet"
            " --rndSeed {SEED}"
            " -i {input} -o {params.out_prefix}"
        )


### Reconstruction aka trimming


rule leeHom:
    """Reconstruction using leeHom"""
    input:
        s1=OUTDIR_SIM + "gen_n{n}_l{l}_s1.fq",
        s2=OUTDIR_SIM + "gen_n{n}_l{l}_s2.fq",
    output:
        OUTDIR_REC + "leeHom/gen_n{n}_l{l}_lh.fq.gz",
        OUTDIR_REC + "leeHom/gen_n{n}_l{l}_lh_r1.fq.gz",
        OUTDIR_REC + "leeHom/gen_n{n}_l{l}_lh_r2.fq.gz",
        OUTDIR_REC + "leeHom/gen_n{n}_l{l}_lh.fail.fq.gz",
        OUTDIR_REC + "leeHom/gen_n{n}_l{l}_lh_r1.fail.fq.gz",
        OUTDIR_REC + "leeHom/gen_n{n}_l{l}_lh_r2.fail.fq.gz",
    wildcard_constraints:
        n="\d+",
        l="\d+",
    params:
        out_prefix=OUTDIR_REC + "leeHom/gen_n{n}_l{l}_lh",
    benchmark:
        OUTDIR_BEN + "leeHom/gen_n{n}_l{l}_lh.tsv"
    shell:
        # (default : PHRED33)
        # -t [threads]            Use multiple cores (default : 1)
        (
            "{LEEHOM}"
            " --adapterFirstRead {ADPT1}"
            " --adapterSecondRead {ADPT2}"
            " --ancientdna"
            " -fq1 {input.s1}"
            " -fq2 {input.s2}"
            " -fqo {params.out_prefix}"
        )


rule AdapterRemoval:
    """Reconstruction using AdapterRemoval"""
    input:
        s1=OUTDIR_SIM + "gen_n{n}_l{l}_s1.fq",
        s2=OUTDIR_SIM + "gen_n{n}_l{l}_s2.fq",
    output:
        OUTDIR_REC + "AdapterRemoval/gen_n{n}_l{l}_ar.pair1.truncated",
        OUTDIR_REC + "AdapterRemoval/gen_n{n}_l{l}_ar.pair2.truncated",
        OUTDIR_REC + "AdapterRemoval/gen_n{n}_l{l}_ar.singleton.truncated",
        OUTDIR_REC + "AdapterRemoval/gen_n{n}_l{l}_ar.discarded",
        OUTDIR_REC + "AdapterRemoval/gen_n{n}_l{l}_ar.collapsed.truncated",
        OUTDIR_REC + "AdapterRemoval/gen_n{n}_l{l}_ar.collapsed",
    wildcard_constraints:
        n="\d+",
        l="\d+",
    params:
        out_prefix=OUTDIR_REC + "AdapterRemoval/gen_n{n}_l{l}_ar",
    benchmark:
        OUTDIR_BEN + "AdapterRemoval/gen_n{n}_l{l}_ar.tsv"
    shell:
        # --seed SEED
        #   Sets the RNG seed used when choosing between bases with 
        #   equal Phred scores.
        # Phred is by default 33
        # number of threades is by default 1
        (
            "{ADPTREM}"
            #" --gzip"
            " --collapse"
            " --minlength 1"
            " --adapter1 {ADPT1}"
            " --adapter2 {ADPT2}"
            " --file1 {input.s1}"
            " --file2 {input.s2}"
            " --basename {params.out_prefix}"
            " --seed {SEED}"
        )


rule ClipAndMerge:
    """Reconstruction using ClipAndMerge"""
    input:
        s1=OUTDIR_SIM + "gen_n{n}_l{l}_s1.fq",
        s2=OUTDIR_SIM + "gen_n{n}_l{l}_s2.fq",
    output:
        OUTDIR_REC + "ClipAndMerge/gen_n{n}_l{l}_cm.fq",
    wildcard_constraints:
        n="\d+",
        l="\d+",
    benchmark:
        OUTDIR_BEN + "ClipAndMerge/gen_n{n}_l{l}_cm.tsv"
    shell:
        # Phred Score offset default: 33
        # -rm_no_partner : Remove reads with no pairing partner after
        #                  adapter clipping. (default: false)
        (
            "java -jar {CLIPMERGE}"
            " -in1 {input.s1}"
            " -in2 {input.s2}"
            " -f {ADPT1}"
            " -r {ADPT2}"
            " -o {output}"
            " -l 1"
            " -rm_no_partner"
        )


rule seqtk_adna_trim:
    """Reconstruction using seqtk and adna-trim"""
    input:
        s1=OUTDIR_SIM + "gen_n{n}_l{l}_s1.fq",
        s2=OUTDIR_SIM + "gen_n{n}_l{l}_s2.fq",
    output:
        OUTDIR_REC + "seqtk_adna-trim/gen_n{n}_l{l}_at.fq",
    wildcard_constraints:
        n="\d+",
        l="\d+",
    #params:
    #    out_pe=OUTDIR_REC + "seqtk_adna_trim/gen_n{n}_l{l}_at",
    benchmark:
        OUTDIR_BEN + "seqtk_adna-trim/gen_n{n}_l{l}_at.tsv",
    shell:
        # seqtk mergepe: interleave two paired-end FASTA/Q files
        # adna-trim:
        # -l INT       min read/fragment length to output [default: 30]
        # -t INT       number of threads [default: 2]
        # -p STR       output PE reads to STR.R[12].fq.gz [default: discard pe]
        (
            "{SEQTK} mergepe"
            " {input.s1} {input.s2} |"
            " {ADNA}"
            " -l 1"
            " -t 1"
            #" -p {params.out_pe}"
            " -"
            " > {output}"
        )


rule bbmerge:
    """Reconstruction using BBMerge"""
    input:
        s1=OUTDIR_SIM + "gen_n{n}_l{l}_s1.fq",
        s2=OUTDIR_SIM + "gen_n{n}_l{l}_s2.fq",
    output:
        m=OUTDIR_REC + "bbmerge/gen_n{n}_l{l}_bb.fq",
        #u1=OUTDIR_REC + "bbmerge/gen_n{n}_l{l}_bb.R1.fq"
        #u2=OUTDIR_REC + "bbmerge/gen_n{n}_l{l}_bb.R2.fq"
    wildcard_constraints:
        n="\d+",
        l="\d+",
    params:
    benchmark:
        OUTDIR_BEN + "bbmerge/gen_n{n}_l{l}_bb.tsv",
    shell:
        # t=1           : Set threads to 1
        # outu=<file>   : File for unmerged reads.
        (
            "{BBMERGE}"
            " in1={input.s1}"
            " in2={input.s2}"
            " out={output.m}"
            #" outu1={output.u1}"
            #" outu2={output.u2}"
            " adapter={ADPT_FILE}"
            " t=1"
            " mininsert=1"
            " mininsert0=1"
        )


rule fastp:
    """Reconstruction using fastp"""
    input:
        s1=OUTDIR_SIM + "gen_n{n}_l{l}_s1.fq",
        s2=OUTDIR_SIM + "gen_n{n}_l{l}_s2.fq",
    output:
        m=OUTDIR_REC + "fastp/gen_n{n}_l{l}_fp.fq",
        #u1=OUTDIR_REC + "fastp/gen_n{n}_l{l}_fp.R1.fq",
        #u2=OUTDIR_REC + "fastp/gen_n{n}_l{l}_fp.R2.fq",
        #json=OUTDIR_REC + "fastp/gen_n{n}_l{l}_fp_report.json",
        #html=OUTDIR_REC + "fastp/gen_n{n}_l{l}_fp_report.html",
    wildcard_constraints:
        n="\d+",
        l="\d+",
    params:
        out_prefix=OUTDIR_REC + "fastp/gen_n{n}_l{l}",
    benchmark:
        OUTDIR_BEN + "fastp/gen_n{n}_l{l}_fp.tsv",
    shell:
        # takes phread33 as input
        (
            "{FASTP}"
            " --merge "
            " --in1 {input.s1}"
            " --in2 {input.s2}"
            " --adapter_sequence {ADPT1}"
            " --adapter_sequence_r2 {ADPT2}"
            " --merged_out {output.m}"
            " --disable_length_filtering"
            " --length_required 1"
            #" --out1 {output.u1}"
            #" --out2 {output.u2}"
            " --json /dev/null"
            " --html /dev/null"
        )
    

rule SeqPrep:
    """Reconstruction using SeqPrep"""
    input:
        s1=OUTDIR_SIM + "gen_n{n}_l{l}_s1.fq",
        s2=OUTDIR_SIM + "gen_n{n}_l{l}_s2.fq",
    output:
        m=OUTDIR_REC + "SeqPrep/gen_n{n}_l{l}_sp.fq.gz",
        #u1=OUTDIR_REC + "SeqPrep/gen_n{n}_l{l}_sp.R1fq.gz",
        #u2=OUTDIR_REC + "SeqPrep/gen_n{n}_l{l}_sp.R2fq.gz",
    wildcard_constraints:
        n="\d+",
        l="\d+",
    benchmark:
        OUTDIR_BEN + "SeqPrep/gen_n{n}_l{l}_sp.tsv"
    shell:
        # The output is always gziped compressed.
        # -f <first read input fastq filename>
        # -r <second read input fastq filename>
        # -s <perform merging and output the merged reads to this file>
        # -1 <first read output fastq filename>
        # -2 <second read output fastq filename>
        # -L <Minimum length of a trimmed/merged read; default = 30>
        # -A <forward read primer/adapter sequence to trim>
        # -B <reverse read primer/adapter sequence to trim>
        (
            "{SEQPREP}"
            " -f {input.s1}"
            " -r {input.s2}"
            " -s {output.m}"
            " -1 /dev/null"
            " -2 /dev/null"
            " -L 1"
            " -A {ADPT1}"
            " -B {ADPT2}"
        )


### Evaluate trimming performance


rule evaluate:
    input:
        orig=OUTDIR_SIM + "gen_n{n}_l{l}.fa",
        lh=OUTDIR_REC + "leeHom/gen_n{n}_l{l}_lh.fq.gz",
        ar=OUTDIR_REC + "AdapterRemoval/gen_n{n}_l{l}_ar.collapsed",
        cm=OUTDIR_REC + "ClipAndMerge/gen_n{n}_l{l}_cm.fq",
        at=OUTDIR_REC + "seqtk_adna-trim/gen_n{n}_l{l}_at.fq",
        bb=OUTDIR_REC + "bbmerge/gen_n{n}_l{l}_bb.fq",
        fp=OUTDIR_REC + "fastp/gen_n{n}_l{l}_fp.fq",
        sp=OUTDIR_REC + "SeqPrep/gen_n{n}_l{l}_sp.fq.gz",
    output:
        lh_out=OUTDIR_EVA + "leeHom/gen_n{n}_l{l}_leehom.csv",
    run:
        shell("touch {output.lh_out}")

