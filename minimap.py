#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import argparse
from collections import OrderedDict

from nextmeta.common import check_path, check_paths, mkdir, read_tsv, read_files, get_version
from dagflow import DAG, Task, ParallelTask, do_dag


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "1131978210qq.com"
__version__ = "v2.2.0"


MINIMAP_BIN = "/export/personal/software/software/minimap2/v2.17/"
SAMTOOLS_BIN = "/export/personal/software/software/samtools/v1.10/bin/"
SEQUENCER = OrderedDict([
    ("RSII", {"minimap2": "-ax map-pb"}
    ),
    ("Sequel", {"minimap2": "-ax map-pb"}
    ),
    ("GridION", {"minimap2": "-ax map-ont"}
    ),
    ("PromethION", {"minimap2": "-ax map-ont"}
    ),
    ("illumina", {"minimap2": "-ax sr"}
    ),
    ("mgi", {"minimap2": "-ax sr"}
    ),
])

SOFTWARE_VERSION = {
    "minimap2": {
        "GETVER": "%s/minimap2 --version" % MINIMAP_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.17",
    },
    "samtools": {
        "GETVER": "%s/samtools --version" % SAMTOOLS_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.10",
    },
}


def read_fasta(file):

    if file.endswith(".gz"):
        fa = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fa = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''
    for line in fa:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            seq = seq.split('\n')
            if len(seq)==2:
                yield seq[0], seq[1]
            seq = ''
            line = line.strip(">").split()[0]
            seq += "%s\n" % line
            continue
        seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield seq[0], seq[1]
    fa.close()


def read_fastq(file):

    if file.endswith(".gz"):
        fp = gzip.open(file, 'r')
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith("@") and (len(seq)==0 or len(seq)>=5):
            seq = []
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith("@") and len(seq)==4:
            yield seq[0], seq[1]
            seq = []
            seq.append(line.strip("@").split()[0])
            continue
        seq.append(line)

    if len(seq)==4:
        yield seq[0], seq[1]
    fp.close()


def gmk2pb(string):

    string = string.lower().strip()

    if string.endswith('g') or string.endswith('gb'):
        base = string.split('g')[0]
        base = float(base)*1e+09
    elif string.endswith('m') or string.endswith('mb'):
        base = string.split('m')[0]
        base = float(base)*1e+06
    elif string.endswith('k') or string.endswith('kb'):
        base = string.split('k')[0]
        base = float(base)*1e+03
    else:
        base = float(string)

    return int(base)


def allot_base(file, name, base, output_dir):

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    n = 0
    lable = 1
    base = gmk2pb(base)
    fname = "%s/%s.%s.fa" % (output_dir, name, str(lable))
    output = open(fname, 'w')
    file_list = []

    for seqid, seq in fh:
        if n>=base:
            output.close()
            file_list.append(check_path(fname))
            n = len(seq)
            lable += 1
            fname = "%s/%s.%s.fa" % (output_dir, name, str(lable))
            output = open(fname, 'w')
            output.write(">%s\n%s\n" % (seqid, seq))
        else:
            output.write(">%s\n%s\n" % (seqid, seq))
            n += len(seq)

    output.close()
    file_list.append(check_path(fname))

    return file_list


def create_minimap_task(genomes, name, read1, read2, sequencer, thread, job_type, job_option, work_dir="", out_dir="", split=""):

    option = OrderedDict()
    option["minimap2"] = {
        "version": get_version(SOFTWARE_VERSION["minimap2"]),
        "option": SEQUENCER[sequencer]["minimap2"]
    }
    option["samtools"] = {
        "version": get_version(SOFTWARE_VERSION["samtools"]),
        "option": "samtools view |samtools sort"
    }

    prefix = [os.path.basename(i) for i in genomes]
    id = "minimap"

    if sequencer in ['llumina', 'mgi']:
        reads = "%s %s" % (read1, read2)
    else:
        reads = read1

    temps = []
    for i in range(len(genomes)):
        if split=="split":
            temps.append("--split-prefix temp%s" % i)
        else:
            temps.append("")

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option=job_option,
        script="""
export PATH={samtools}:{minimap}:$PATH
samtools faidx {{genomes}}
minimap2 {{temps}} -t {thread} {x} {{genomes}} {reads} \\
|samtools view --threads {thread} -bS -T {{genomes}}.fai \\
|samtools sort --threads {thread} -m {memory}G -o {{prefix}}.sort.bam
""".format(minimap=MINIMAP_BIN,
           samtools=SAMTOOLS_BIN,
           x=SEQUENCER[sequencer]["minimap2"],
           reads=reads,
           thread=thread,
           memory=thread*4),
        temps=temps,
        genomes=genomes,
        prefix=prefix
    )

    join_task = Task(
        id="merge_minimap",
        work_dir=work_dir,
        type=job_type,
        option=job_option,
        script="""
export PATH={samtools}:$PATH
samtools merge -f -c --threads {thread} {name}.sorted.bam ./minimap/*.sort.bam
samtools index {name}.sorted.bam
""".format(samtools=SAMTOOLS_BIN,
           name=name,
           thread=thread)
    )
    join_task.set_upstream(*tasks)

    return tasks, join_task, option, os.path.join(work_dir, "%s.sorted.bam" % name)


def run_minimap(genome, name, read1, read2, sequencer, thread, job_type, work_dir="", out_dir="", split=""):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    genome = check_path(genome)
    read1 = check_path(read1)

    if sequencer in ['llumina', 'mgi']:
        read2 = check_path(read2)

    dag = DAG("minimap")
    work_dict = {
        "split": "00_data",
        "minimap": "01_minimap",
    }

    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))

    genomes = allot_base(file=genome, name=name, base='10M', output_dir=work_dict["split"])

    if job_type == "local":
        job_option = ""
    else:
        job_option = "-l vf=%sG -pe smp %s" % (str(thread*6) , thread)
    tasks, join_task, option, bam =  create_minimap_task(
        genomes=genomes,
        name=name,
        read1=read1,
        read2=read2,
        sequencer=sequencer,
        thread=thread,
        job_type=job_type,
        job_option=job_option,
        work_dir=work_dict["minimap"],
        out_dir=out_dir,
        split=split)

    dag.add_task(*tasks)
    dag.add_task(join_task)

    return dag, option, bam


def minimap(genome, name, read1, read2, sequencer, thread, job_type, concurrent, refresh, work_dir, out_dir, split=""):

    dag, option, bam = run_minimap(genome, name, read1, read2, sequencer, thread, job_type, work_dir, out_dir, split)
    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    with open(os.path.join(out_dir, "minimap2.json"), "w") as fh:
        json.dump(option, fh, indent=2)


def minimap_hlep_args(parser):

    parser.add_argument("-r1", "--read1", metavar='FILE', nargs='+', type=str, required=True,
        help="Input R1 data, if it is three-generation sequencing, enter the reads sequence.")
    parser.add_argument("-r2", "--read2", metavar='FILE', nargs='+', type=str, default='',
        help="Input R2 data, not input if it is three-generation sequencing.")
    parser.add_argument("-n", "--name", metavar="STR", type=str, default="out",
        help="Input sample name.")
    parser.add_argument("-g", "--genome", metavar='FILE', type=str, required=True,
        help="Input genome file.")
    parser.add_argument("-seq", "--sequencer", choices=["PromethION", "GridION", "RSII", "Sequel", "illumina", "mgi"], default="PromethION",
        help="Select sequencing platform, default=PromethION")
    parser.add_argument("-s", "--split", choices=["split", ""], default="",
        help="Cache output, default=")
    parser.add_argument("-t", "--thread", metavar='INT', type=int, default=4,
        help="Set the running thread, default=4")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", default=".",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", default=".",
        help="Output directory (default: current directory)")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    minimap.py Compare reads to genome.

attention:
    minimap.py --read1 ngs.r1.fq --read2 ngs.r2.fq --genome genome.fa
    minimap.py --read1 tgs.read.fq --genome genome.fa

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = minimap_hlep_args(parser).parse_args()

    minimap(args.genome, args.name, args.read1, args.read2, args.sequencer, args.thread, args.job_type, args.concurrent, args.refresh, args.work_dir, args.out_dir, args.split)


if __name__ == "__main__":

    main()
