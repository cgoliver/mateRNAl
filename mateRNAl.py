###############################################################################
#Copyright (c) 2017, Carlos Oliver, Vladimir Reinharz,                        #
#& Jerome Waldispuhl                                                          #
#All rights reserved.                                                         #
#                                                                             #
#Redistribution and use in source and binary forms, with or without           #
#modification, are permitted provided that the following conditions are met:  #
#* Redistributions of source code must retain the above copyright             #
#notice, this list of conditions and the following disclaimer.                #
#* Redistributions in binary form must reproduce the above copyright          #
#notice, this list of conditions and the following disclaimer in the          #
#documentation and/or other materials provided with the distribution.         #
#* Neither the name of the <organization> nor the                             #
#names of its contributors may be used to endorse or promote products         #
#derived from this software without specific prior written permission.        #
#                                                                             #
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  #
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    #
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   #
#ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY       #
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   #
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  #
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT   #
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS#
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE                  #
###############################################################################

#author: Carlos Oliver

import begin
import sys
import multiprocessing
import os
import copy
from collections import namedtuple
import re
import subprocess
import os
import tempfile
import random
import logging
import numpy as np
import multiprocessing as mp
import shutil
import math
import RNAstats

BASES = ["A", "U", "C", "G"]
R = 0.0019872041
T = 310.15

CSV_HEADER = "generation,sequence,structure,energy,probability,gc,mutations,fitness,id,parent\n"
#sequence class
class RNA:
    def __init__(self, seq, struc, eng, prob, muts, id, parent, ancestor):
        self.sequence = seq
        self.structure = struc
        self.energy = eng
        self.probability = prob
        self.gc = gc_content(seq)
        self.id = id
        self.mutations = muts
        self.parent = parent
        self.ancestor = ancestor
        pass
    def __str__(self):
        return ",".join(map(str, [self.sequence, self.structure, self.energy, self.probability\
               , self.gc, self.mutations, self.fitness, self.id, self.parent]))

    def set_fitness(self, fit="energy", beta=1,target=None):
        if fit == "energy":
            #careful. energy is negative so no need to multiply by -1
            #self.fitness = math.exp((-1 * beta * self.energy) / (R * T))
            self.fitness = math.exp((-1 * beta * self.energy) / (R * T))
        if fit == "target":
            self.fitness = math.exp((-1 * beta * bp_dist(self.structure,\
            target)) / len(self.structure))
        pass

def hamming(s1, s2):
    diffs = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            diffs = diffs + 1
    return diffs

def bp_dist(structure, target):
    bps_struc = ss_to_bp(structure)
    bps_target = ss_to_bp(target)

    return len(bps_struc.symmetric_difference(bps_target))

def ss_to_bp(structure):
    bps = set()
    l = []
    for i , x in enumerate(structure):
        if x == "(":
            l.append(i)
        elif x == ")":
            bps.add((l.pop(), i))
    return bps

#fold method
def fold(sequences):
    with tempfile.TemporaryDirectory() as tmpdir:
        home_dir = os.getcwd()
        #go to tempdir
        os.chdir(tmpdir)
        #write sequences to input file
        with open("seq.in", "w+") as seqfile:
            for i, s in enumerate(sequences):
                seqfile.write("> %s \n" % (i))
                seqfile.write(s + "\n")

        #make output file and call RNAfold
        with open("seq.in", "r") as seqin, open("rnafold.out", "w") as rnaout:
            p = subprocess.run(["RNAfold", "-d", "0", "--noPS", "-p"],\
                stdin=seqin, stdout=rnaout)
        rnas = []
        rnatup = namedtuple("rnatup", "sequence structure energy probability")
        with open("rnafold.out", "r") as rnaout:
            lines = rnaout.readlines()
            for i, l in enumerate(lines):
                if ">" in l:
                    sequence = lines[i+1].strip()
                    struc_info = lines[i+2].split()
                    structure = struc_info[0]
                    energy = float(re.findall("-?\d+.\d+",  struc_info[-1])[0])
                    p = lines[i+5].split()[6].strip(";")
                    if 'e' in p:
                        n = p.split('e')
                        prob = math.pow(float(n[0]), float(n[1]))
                    else:
                        prob = float(p)
                    rnas.append(rnatup(sequence=sequence, structure=structure,\
                        energy=energy, probability=prob))
        os.chdir(home_dir)

    return rnas

def gc_content(s):
    return len([n for n in s if n == "G" or n == "C"]) / float(len(s))
def populate(size, length, gc, fit="energy", target=None):
    pop = []
    while len(pop) < size:
        seq = None
        while True:
            seq = "".join([np.random.choice(BASES, p=[(1-gc)/2, (1- gc)/2,\
                gc/2, gc/2]) for _ in range(length)])
            if gc - 0.1 <= gc_content(seq) <= gc + 0.1:
                break
        pop.append(seq)
    folded = fold(pop)
    rnas = [RNA(f.sequence, f.structure, f.energy, f.probability, 0, i, i,\
                f.sequence) for i, f in enumerate(folded)]
    for r in rnas:
        r.set_fitness(fit=fit, target=target)
    return rnas

def mutate(rna, mutation_rate):
    mutations = 0
    new_seq = ""
    for s in rna.sequence:
        r = random.random()
        if r < mutation_rate:
            mutations = mutations + 1
            new_bases = [b for b in BASES if b != s]
            new_seq = new_seq + np.random.choice(new_bases)
        else:
            new_seq = new_seq + s
    return (new_seq, mutations)

def write_pops(pops, dest, verbose=False):
    with open(dest, "w+") as d:
        d.write(CSV_HEADER)
        if verbose:
            for gen,p in enumerate(pops):
                for rna in p:
                    d.write("{0},{1}\n".format(gen, rna))
        else:
            for gen, p in enumerate(pops):
                gen_uniques = {}
                for rna in p:
                    gen_uniques[rna.sequence] = rna
                for seq in gen_uniques:
                    d.write("{0},{1}\n".format(gen,gen_uniques[seq]))


def pop_stat(pop, stat):
    stat_list = np.mean([getattr(s, stat) for s in pop])
    return (np.mean(stat_list), np.std(stat_list))

# select
def select(pop, mutation_rate, gc, fit="energy", target=None, density=False,\
        K=10, fixpop=True):
    #normalize(pop, density=density, K=K)
    fitnesses = [r.fitness for r in pop]
    tot = np.sum(fitnesses)
    parents = []
    if fixpop:
        parents = np.random.choice(pop, p=[rna.fitness/tot for rna in pop], replace=True,\
            size=len(pop))
    else:
        strucs = {}
        families = {}
        #get phenotype counts
        for rna in pop:
            hp = RNAstats.loop_counter(rna.structure)['hairpin']
            families.setdefault(hp, 0)
            families[hp] += 1
            strucs.setdefault(rna.structure, 0)
            strucs[rna.structure] += 1

        print(families)
        #add offspring for each parent to parents list
        for rna in pop:
            P = strucs[rna.structure]
            offspring = rna.fitness * P  * (1 - (P/K))
            for _ in range(int(offspring)):
                parents.append(rna)
    next_gen = []
    re_fold = []
    for i, p in enumerate(parents):
        child_seq = None
        mutations = 0
        while True:
            #mutate the sequence
            child_seq, mutations = mutate(p, mutation_rate)
            #compute mutated GC
            child_gc = gc_content(child_seq)
            #this loop breaks only once the sequence gc is in the gc range
            if gc - 0.1 <= child_gc <= gc + 0.1:
                break
        child_obj = RNA(child_seq, p.structure, p.energy, p.probability,\
                int(hamming(p.sequence, p.ancestor)), i, p.id, p.ancestor)
        if mutations:
            re_fold.append((child_obj, i))
        next_gen.append(child_obj)
    #re-fold rnas that had mutations in their sequences
    re_folded = fold([s[0].sequence for s in re_fold])
    for i, r in enumerate(re_folded):
        update_seq = next_gen[re_fold[i][1]]

        update_seq.structure = r.structure
        update_seq.sequence = r.sequence
        update_seq.energy = r.energy
        update_seq.probability = r.probability
    #set the fitness of each individual in the next generation
    for c in next_gen:
        c.set_fitness(fit=fit, target=target)
    return next_gen
def normalize(pop, density=False, K=10):
    raw_fit = [rna.fitness for rna in pop]
    tot = np.sum(raw_fit)
    for rna in pop:
        #rna.fitness = rna.fitness / tot
        rna.fitness = 1 / (1 + math.exp(rna.energy))
    if density == False:
        return
    else:
        #count number of structures per family
        families = {0: 0, 1: 0, 2: 0, 3: 0}
        strucs = {}
        for rna in pop:
            ss = RNAstats.loop_counter(rna.structure)
            hp = ss['hairpin']
            rna.hp = hp
            families[hp] += 1
            strucs.setdefault(rna.structure, 0)
            strucs[rna.structure] += 1
        print(families)
        print(len(strucs))
        for rna in pop:
            fit_before = rna.fitness
            rna.fitness = max(0.0, rna.fitness * (1 - (strucs[rna.structure] / K)))
            #print("before: %s, after: %s, delta: %s, hp: %s" \
               # % (fit_before, rna.fitness, fit_before - rna.fitness, rna.hp))
            normalize(pop, density=False)
    pass
#evolve
def evolve(args):
    generations, size, length, fit, gc, mutation_rate, dest,\
            target, verbose, density, K, fixpop = args
    current_pop = populate(size, length, gc, fit=fit, target=target)
    pops = [current_pop]
    for g in range(generations):
        #density normalize current_pop
        next_pop = select(current_pop, mutation_rate, gc, fit=fit,\
                target=target, density=density, K=K, fixpop=fixpop)
        current_pop = next_pop
        pops.append(current_pop)
    write_pops(pops, dest,verbose=verbose)
    pass

def dest_format(dest, i):
    head, tail = os.path.split(dest)
    prefix, suffix = tail.split(".")
    prefix = prefix + "_" + str(i)
    new_tail = prefix + "." + suffix
    return os.path.join(head, new_tail)
@begin.start
@begin.convert(generations=int, size=int, length=int, fit=str, gc=float,\
        mutation_rate=float, runs=int, procs=int, target=str, dest=str,\
        density=bool, K=int, fixpop=bool)
def start(generations=20, size=10, length=50, fit='energy', gc=0.5,\
        mutation_rate=0.1, runs=1, procs=1, beta=1, target=None, dest="maternal.csv",\
        verbose=False, density=False, K=10, fixpop=True):

    if fit =="target":
        length=len(target)
    todo = ((generations, size, length, fit, gc, mutation_rate,\
        dest_format(dest, i), target, verbose, density, K, fixpop) for i in range(runs))
    if runs > 1:
        with multiprocessing.Pool(procs) as pool:
            pool.map(evolve, todo)
    else:
        evolve(next(todo))
