#!/usr/bin/env python3

import sys, getopt, re, subprocess, os, argparse, shutil, pickle, fnmatch, yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML._paml import PamlError

debug_mode = False
def debug(msg):
    if debug_mode:
        print(msg, file=sys.stderr)


parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True, nargs="+", 
                   help = "Path to yaml file with config and task " +
                   "definitions.")
parser.add_argument("--verbose", default=False, help="Be verbose", action='store_true')

#parser.add_argument("--tree", required=True, nargs="+", help="Treefile to supply to CODEML (one or more)")
#parser.add_argument("--codeml-ctl", required=True, nargs="+", help="CodeML control file (one or more files)")
#parser.add_argument("--homology-file", required=True, help="File with homology definitions")
#parser.add_argument("--macse", required=True, help="Path to MACSE jar")
#parser.add_argument("--data", default="./codeml_data",
#                    help="Data spool dir for results from subtasks (defaults to ./codeml_data)")
#parser.add_argument("--chunk-size", default="10", help="Array job chunk size",
#                    type=int)
parser.add_argument("--outfile", help="Store output to this file, otherwise, print to STDOUT")

mut_group=parser.add_mutually_exclusive_group(required=True)
mut_group.add_argument('-c','--calculate', help="Run calculations for given chunk",
                        nargs='?', type=int)

mut_group.add_argument('-n', '--get-tasks', help="Get number of chunks",
                       action="store_true")
mut_group.add_argument('-e', '--extract',
                       help="List of datapaths to extract (see for value tree). Wildcards are allowed",
                       nargs="*")
mut_group.add_argument("-l", "--list-paths", help="List available result paths",
                       action="store_true")


args=parser.parse_args()


if args.verbose:
    debug_mode = True

debug("Verbose mode enabled")
    
tmp = os.getenv("TMPDIR")
if tmp == None:
    tmp = "."

debug("Using %s for temporary files" % tmp)

config = yaml.load( open(args.config[0]).read(), yaml.Loader )

# Read homology file
homologies = [x.rstrip() for x in open( config['common']['homology_file'] )]
num_tasks = int(len(homologies) / config['common'].get('chunk_size', 5)) + 1

# Report number of tasks
if args.get_tasks:
    debug("Calculating number of tasks")
    print(num_tasks)
    sys.exit(0)

elif args.calculate:
    debug("Doing homologies %d:%d" % ((args.calculate - 1) * config['common'].get('chunk_size', 5), 
                                      args.calculate * config['common'].get('chunk_size', 5)))
    # Run actual calculations
    debug("Opening sequence files")
    seq_files = [ SeqIO.index(x,"fasta") for x in config['common']['fasta_files'] ]

 
    
    # Regular expression to parse sequence description
    rex = re.compile(r'''(?x) (\w+) : ( ' (?: \\. | [^'] )* ' | " (?: \\. | [^"] )* "|[^'"\s]+)''')

    paml_results = list()
    n=1  
    for line in homologies[(args.calculate - 1) * config['common'].get('chunk_size', 5) : 
                           args.calculate * config['common'].get('chunk_size', 5)]:
        # Get sequences
        seqs = [ x[0][x[1]] for x in zip(seq_files,line.split()) ]

        # Extract metadata from sequence descriptions
        tokenised_descs = [{k:v.strip("'\"") for k, v in re.findall(rex, y.description)} for y in seqs]

        cds = [ x[0][x[1]:] for x in zip(seqs, [int(t["offset"]) for t in tokenised_descs])]
            
        # Get protein sequences
        prots = [ c.seq.translate(to_stop=True) for c in cds ]

        # Truncate cds before stop codon
        cds = [ c[0][:c[1]*3] for c in zip(cds, [len(p) for p in prots]) ]

        # Base name for alignment files
        fnbase = tmp + "/" + "_".join([c.id for c in cds])
        
        # Write CDS to fasta
        SeqIO.write(cds,open(fnbase + ".fasta","w"),"fasta")

        # Run MACSE
        macse_args=['java', '-Xms1G', '-Xmx2G', 
                    '-jar', config['common']['macse'], '-prog', 'alignSequences',
                    '-seq', fnbase + ".fasta"]
        if debug_mode:
            macse = subprocess.Popen(macse_args)
        else:
            macse = subprocess.Popen(macse_args, stdout=subprocess.PIPE)

        macse.wait()

        # Parse alignment
        aln = AlignIO.read(open(fnbase + "_NT.fasta"), "fasta")

        # Rename sequences in alignment for CODEML and convert "!" output by
        # macse to "-". Also rid the sequence of internal stops and flag
        # alignment as possible pseudogene
        pseudo=False
        cnt = 1;
        for a in aln:
            a.id="seq" + str(cnt)
            a.seq = Seq(str(a.seq).replace("!","-"))
            tseq = Seq(str(a.seq).replace("-","N"))
            pr = str(tseq.translate())
            list_seq = list(str(a.seq))
            
            for i in range(len(pr)):
                if pr[i] == "*":
                    print("!!!!!! Internal STOP after alignment, possible pseudogene %d\n" % i)
                    pseudo = True
                    list_seq[ i*3 : i*3+3 ] = list("NNN")
            if pseudo:
                a.seq = Seq("".join(list_seq))
            cnt += 1
        try:
            # Write as phylip
            AlignIO.write(aln, fnbase + ".phylip", "phylip-sequential")
        except ValueError:
            paml_results.append(dict( paml_results = dict(),
                                  seq_ids      = [s.id for s in seqs],
                                  seq_lens     = [len(s) for s in cds],
                                  aln_len      = aln.get_alignment_length(),
                                  pseudo       = pseudo ))
            continue
        
        # Run the different combinations of tree and codeml.ctl
        p_res=dict()
        for job_name, job in config['jobs'].items():
            # Copy treefile
            shutil.copy( job['tree'], fnbase + ".tree")
    
            # run Codeml
            cml = codeml.Codeml(alignment = fnbase + ".phylip", tree = fnbase + ".tree", out_file = fnbase + ".out", working_dir = tmp + "/scratch")
            cml.read_ctl_file( job['ctl'] )
            try:
                cml_res=cml.run(verbose=debug_mode)
            
                def to_path(path, d):
                    res = list()
                    for k,v in d.items():
                        if isinstance(v, dict):
                            res += to_path(path + "/" + str(k), v)
                        else:
                            res.append((path + "/" + str(k), v))
                    return res
                 
                for p, cml in to_path("", cml_res):
                    p_res[job_name + str(p)] = cml
            
            except PamlError:
                debug("CodeML failed")
            except ValueError:
                debug("CodeML failed")
            except OSError:
                debug("CodeML failed")

        paml_results.append(dict( paml_results = p_res,
                                  seq_ids      = [s.id for s in seqs],
                                  seq_lens     = [len(s) for s in cds],
                                  aln_len      = aln.get_alignment_length(),
                                  pseudo       = pseudo ))
    pickle.dump(paml_results, open(config['common']['data_dir'] + "/" + str(args.calculate) + ".pickle", "wb"))

else:
    # Extract or list paths. Both need to read the whole dataset to find available paths/data
    paml_results = list()
    for i in range(1, num_tasks + 1):
        try:	
            r = pickle.load(open(config['common']['data_dir'] + "/" + str(i) + ".pickle", "rb"))
            paml_results += r
        except:
            debug("File not found: " + str(i) + ".pickle")

    paths = set()
    for ri in range(len(paml_results)):
        r=paml_results[ri]
        for p in list(r["paml_results"].keys()): # Copy before iterate
            # Expand array-like results
            if type(r["paml_results"][p]) == list:
                for pair in enumerate(r["paml_results"][p]):
                    tp=p + "/" + str(pair[0])
                    paths |= set([ tp ])
                    paml_results[ri]["paml_results"][tp]=pair[1]
                    
            paths |= set([p])

    paths = list(paths)
    paths.sort()

    paths_l = [ p.split("/") for p in paths ]

    i = 0
    to_discard = set()
    try:
        while True:
            col = set([ p[i] for p in paths_l ])
            if len(col) == 1:
                to_discard |= { i }
            i += 1
    except IndexError:
        pass

    new_paths = [ "/".join( [ v for k, v in zip(range(len(p)),p) \
                              if k not in to_discard ]) for p in paths_l ]
    o2n = dict(zip(paths, new_paths))

    if args.list_paths:
        for p in new_paths:
            print(p)

    else:
        # Filter paml_results according to filter rules
        if len(args.extract):
            to_extract = args.extract
        else:
            to_extract = [ j + '/' + p for j in config['jobs'].keys() 
                                       for p in config['jobs'][j]['export']]

        display_keys = set()
        for p in paml_results:
            p["paml_results"] = {o2n[k]:v for k,v in p["paml_results"].items() if
                                 len([x for x in to_extract
                                      if fnmatch.fnmatch(o2n[k],x)])}
            display_keys |= set(p["paml_results"].keys())
        display_keys = list(display_keys)
        display_keys.sort()

        num_seqs = len(paml_results[0]["seq_lens"])

        
        if args.outfile:
            out = open(args.outfile, "w")
        else:
            out = sys.stdout
            
        # Print header
        print("%s\t%s\t%s\t%s\t%s" % (
            "\t".join(["seq" + str(x + 1) + "id" for x in range(num_seqs)]),
            "\t".join(["seq" + str(x + 1) + "len" for x in range(num_seqs)]),
            "aln_len",
            "pseudo",
            "\t".join(display_keys)),
            file=out)
            

        for p in paml_results:
            print("%s\t%s\t%d\t%s\t%s" % (
                "\t".join(p["seq_ids"]),
                "\t".join([str(x) for x in p["seq_lens"]]),
                p["aln_len"],
                ("YES" if p["pseudo"] else "NO"),
                "\t".join([str(p["paml_results"][x]) if x in p["paml_results"] else "NA" for x in display_keys])),
                file=out)
#        'aln_len': 2079, 'paml_results': {}, 'seq_lens': [2527, 4755, 2711], 'seq_ids': ['BN884_T00000004_1', 'BN885_T00007199_1', 'SCLET00007519_1']
        
        
