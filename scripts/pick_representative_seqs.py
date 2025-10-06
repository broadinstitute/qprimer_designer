#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import numpy as np
import time
from Bio import SeqIO


# In[ ]:


"""Classes and methods for applying locality-sensitive hashing.

This code is modified from CATCH (https://github.com/broadinstitute/catch/blob/master/catch/utils/lsh.py),
which is released under the following license:
MIT License

Copyright (c) 2018 Hayden Metsky, Broad Institute, Inc., and Massachusetts
Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from collections import defaultdict
import hashlib
import heapq
import logging
import math
import random

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


class HammingDistanceFamily:
    """An LSH family that works with Hamming distance by sampling bases."""

    def __init__(self, dim):
        self.dim = dim

    def make_h(self):
        """Construct a random hash function for this family.

        Returns:
            hash function
        """
        i = random.randint(0, self.dim - 1)
        def h(x):
            assert len(x) == self.dim
            return x[i]
        return h

    def P1(self, dist):
        """Calculate lower bound on probability of collision for nearby sequences.

        Args:
            dist: Hamming distance; suppose two sequences are within this
                distance of each other

        Returns:
            lower bound on probability that two sequences (e.g., probes) hash
            to the same value if they are within dist of each other
        """
        return 1.0 - float(dist) / float(self.dim)


class HammingWithGUPairsDistanceFamily:
    """An LSH family that works with a pseudo-Hamming distance, which
    tolerates G-U base pairs, by sampling and transforming bases.

    This family can be used with near-neighbor queries if the goal
    is to tolerate G-U pairing between guide and target sequence,
    assuming there are subsequent checks using a distance function
    more conservative (and accurate) than the relaxed one used here.
    However, it should be used with caution with clustering, as the
    distance function used here is relaxed (see below for explanation).

    This family hashes sequences such that similar sequences hash
    to the same value, while tolerating G-U differences. This is
    meant to accommodate G-U base pairing:
    An RNA guide with U can bind to target RNA with G, and vice-versa.
    Note that a guide sequence that ends up being synthesized to RNA
    is the reverse-complement of the guide sequence constructed here, so
    that it will hybridize to the target (here, the guide sequences
    are designed to match the target). If the RNA guide were to have U and
    the RNA target were to have G, then the guide sequence here would be
    A and the target would be G. If the RNA guide were to have G
    and the RNA target were to have U, then the guide sequence here were
    would be C and the target would be T. Thus, to allow G-U pairing,
    we count a base X in the guide sequence as matching a base Y in
    the target if either of the following is true:
      - X == 'A' and Y == 'G'
      - X == 'C' and Y == 'T'

    This is tricky here because the distance function should have
    symmetry -- i.e., the distance between x and y should be equal
    to the distance between y and x. But, given how we defined the
    distance above (in terms of how to count which bases match), this
    would be violated: the distance would depend on whether x is
    a guide sequence or a target sequence (and whether y is a guide
    sequence or target). To solve this, the family of hash functions
    defined here is meant to work with the following distance metric:
      The distance between sequence x and sequence y is the number
      of mismatched bases between x and y after transforming all
      'A' bases in each sequence to 'G' and all 'C' bases in each
      sequence to 'T'. This is equivalent to Hamming distance after
      the 'A'->'G' and 'C'->'T' transformations.

    For example, we will treat the sequences 'AA' and 'GG' as identical
    and ensure they hash to the same value. Likewise, 'CC' and 'TT'
    must hash to the same value. That is, the distance between them
    should be treated as 0 (regardless of which represents a guide
    sequence and which a target).

    This is a more relaxed distance function than what we might need.
    That is, it may report two sequences x and y as being closer to
    each other than we ought to be considering them to be. For example,
    if x represents a guide sequence and x='GG' and y represents a
    target sequence and y='AA', then the distance between x and y
    will be 0 even though we would not want to think of x as binding
    to y. However, this is OK for the purposes of LSH for near-neighbor
    lookups between guide and target sequence, assuming the queries
    have subsequent checks using a more conservative and accurate
    (albeit asymmetric) distance function; the hash family will
    simply lead to more false positive near neighbor results than
    otherwise, but these will still subsequently be pruned by the
    query function.

    The value of P1 is identical to the value defined above for
    (regular) Hamming distance. To see this, let two sequences x
    and y be within distance d of each other, where d is an 'accurate'
    asymmetric distance that accounts for whether x is a guide
    sequence or target (and likewise for y). Let d' be a corresponding
    upper bound on the distance between x and y according to the
    distance function defined above (i.e., after 'A'->'G' and
    'C'->'T' transformations of x and y). Note that d' <= d. The
    probability that x and y hash to different values is <= d'/dim.
    So the probability of collision for x and y, P1(d'), is
    P1(d') >= 1.0 - d'/dim. Since d' <= d, P1(d') >= P1(d). Hence,
    even if d is used as a substitute for d' (as d is the value
    by which we ultimately want to compare x and y), then P1(d)
    is still a valid lower bound on the probability of collision.
    """

    def __init__(self, dim):
        self.dim = dim

    def make_h(self):
        """Construct a random hash function for this family.

        Returns:
            hash function
        """
        i = random.randint(0, self.dim - 1)
        def h(x):
            assert len(x) == self.dim
            b = x[i]
            # Perform 'A'->'G'and 'C'->'T' transformations, as
            # described above; ignore ambiguity codes (for now) and
            # treat them as mismatches
            if b == 'A':
                b = 'G'
            if b == 'C':
                b = 'T'
            return b
        return h

    def P1(self, dist):
        """Calculate lower bound on probability of collision for nearby sequences.

        Args:
            dist: suppose two sequences are within this distance of each other.
                The distance can be either a Hamming distance that does not
                account for G-U pairing, or a Hamming distance after 'A'->'G'
                and 'C'->'T' transformations (see above for explanation of
                why the value returned will be a valid lower bound regardless)

        Returns:
            lower bound on probability that two sequences (e.g., guides) hash
            to the same value if they are within dist of each other
        """
        return 1.0 - float(dist) / float(self.dim)


class MinHashFamily:
    """An LSH family that works by taking the minimum N elements of a
    permutation of k-mers in a string/sequence (MinHash).

    See (Broder et al. 1997) and (Andoni and Indyk 2008) for details.

    The signature produced here also has similarity to Mash (Ondov et al. 2016).
    """

    def __init__(self, kmer_size, N=1, use_fast_str_hash=False):
        """
        Args:
            kmer_size: length of each k-mer to hash
            N: represent the signature of a sequence using hash values of
                the N k-mers in the sequence that have the smallest hash
                values
            use_fast_str_hash: if True, use a faster (~10x faster) inner
                 hash function for strings (k-mers) that is not deterministic
                 and varies across Python processes; this may not be
                 appropriate if the family is shared across processes and
                 different processes expect the same hash
        """
        self.kmer_size = kmer_size
        self.N = N
        self.use_fast_str_hash = use_fast_str_hash

    def make_h(self):
        """Construct a random hash function for this family.

        Here, we treat a sequence as being a set of k-mers. We calculate
        a hash value for each k-mer and the hash function on the sequence
        returns the N smallest of these in sorted order.

        Returns:
            hash function
        """
        # First construct a random hash function for a k-mer that
        # is a universal hash function (effectively a "permutation"
        # of the k-mer)
        # Constrain all values to be in [0, 2^31 - 1] to have a bound
        # on the output of the universal hash function; this upper bound
        # is nice because it is also a prime, so we can simply work
        # modulo (2^31 - 1)
        p = 2**31 - 1
        # Let the random hash function be:
        #   (a*x + b) mod p
        # for random integers a, b (a in [1, p] and b in [0, p])
        a = random.randint(1, p)
        b = random.randint(0, p)
        if self.use_fast_str_hash:
            def kmer_hash(x):
                # Hash a k-mer x with the random hash function
                # hash(..) uses a random seed in Python 3.3+, so its output
                #   varies across Python processes; thus, this may not be
                #   suitable if used across processes
                x_hash = abs(hash(x))
                return (a * x_hash + b) % p
        else:
            def kmer_hash(x):
                # Hash a k-mer x with the random hash function
                # hashlib.md5(..) gives a deterministic hash value of the k-mer
                #   but is ~10x slower than hash(..)
                x_hash = int(hashlib.md5(x.encode('utf-8')).hexdigest(), 16)
                return (a * x_hash + b) % p

        def h(s):
            # For a string/sequence s, have the MinHash function be the minimum
            # N hashes, in sorted order, over all the k-mers in it
            assert self.kmer_size <= len(s)
            if self.kmer_size >= len(s) / 2:
                logger.warning(("The k-mer size %d is large (> (1/2)x) "
                    "compared to the size of a sequence to hash (%d), which "
                    "might make it difficult for MinHash to find similar "
                    "sequence"), self.kmer_size, len(s))
            num_kmers = len(s) - self.kmer_size + 1
            if num_kmers < self.N:
                logger.warning(("The number of k-mers (%d) in a given "
                    "sequence is too small to produce a signature of "
                    "size %d; the MinHash family might provide unreliable "
                    "distances against the sequence.") % (num_kmers, self.N))
            def kmer_hashes():
                num_yielded = 0
                # Keep yielding k-mers until at least self.N
                #   have been yielded (helpful when s is short)
                while num_yielded < self.N:
                    for i in range(num_kmers):
                        kmer = s[i:(i + self.kmer_size)]
                        num_yielded += 1
                        yield kmer_hash(kmer)
            if self.N == 1:
                # Speed the special case where self.N == 1
                return tuple([min(kmer_hashes())])
            else:
                return tuple(sorted(heapq.nsmallest(self.N, kmer_hashes())))
        return h

    def P1(self, dist):
        """Calculate lower bound on probability of collision for nearby sequences.

        Args:
            dist: Jaccard distance (1 minus Jaccard similarity); suppose
                two sequences are within this distance of each other. The
                Jaccard similarity can be thought of as the overlap in k-mers
                between the two sequences

        Returns:
            lower bound on probability that two sequences (e.g., probes) hash
            to the same value if they are within dist of each other
        """
        # With MinHash, the collision probability is the Jaccard similarity
        return 1.0 - dist

    def estimate_jaccard_dist(self, hA, hB):
        """Estimate Jaccard distance between two MinHash signatures.

        Args:
            hA: signature output by h(A) for a hash function h given by
                make_h() and a sequence A (specifically, an ordered list of
                hash values of k-mers in A)
            hB: signature output by h(B) for the same hash function used
                to generate hA, and a sequence B (specifically, an ordered
                list of hash values of k-mers in B)

        Returns:
            estimate of Jaccard distance between signatures A and B
        """
        # Let X = h(hA \union hB). This is equivalent to h(A \union B)
        # where A and B, here, represent the set of k-mers in their
        # respective sequences. X is a random sample of the hashes of
        # the k-mers in (A \union B), and (Y = X \intersect hA \intersect hB)
        # are the hashes in both X and (A \intersect B). So |Y|/|X| is
        # an estimate of the Jaccard similarity. Since X has self.N
        # elements (it is the output of h(.)), this is |Y|/self.N

        # Iterate over hA and hB simultaneously, to count the number of
        # hash values in common between them, among the first self.N of
        # them
        # This is possible because hA and hB are in sorted order
        hA_i, hB_i = 0, 0
        intersect_count = 0
        union_count = 0
        while hA_i < len(hA) and hB_i < len(hB):
            if union_count == self.N:
                # We found all the hash values in X, as defined above
                break
            elif hA[hA_i] < hB[hB_i]:
                hA_i += 1
                union_count += 1
            elif hA[hA_i] > hB[hB_i]:
                hB_i += 1
                union_count += 1
            else:
                # We found a new unique hash value in common:
                # hA[hA_i] == hB[hB_i]
                intersect_count += 1
                union_count += 1
                hA_i += 1
                hB_i += 1

        # It should be that union_count == self.N, except perhaps in
        # some edge cases (e.g., if a hash value is repeated in a signature)
        # And |Y| == intersect_count, where Y is defined above
        similarity = float(intersect_count) / union_count
        return 1.0 - similarity


class HashConcatenation:
    """Concatenated hash functions (AND constructions)."""

    def __init__(self, family, k, join_as_str=False):
        """
        Args:
            family: hash family object; must have a make_h() function
            k: number of hash functions to concatenate
            join_as_str: if True, concatenate the output of the k hash
                functions into a string before returning the concatenated
                result; if False (default), simply return a tuple of the
                k outputs
        """
        self.family = family
        self.k = k
        self.join_as_str = join_as_str
        self.hs = [family.make_h() for _ in range(k)]

    def g(self, x):
        """Evaluate random hash functions and concatenate the result.

        Args:
            x: point (e.g., probe) on which to evaluate hash functions

        Returns:
            concatenation of the result of the self.k random hash functions
            evaluated at x
        """
        if self.join_as_str:
            return ''.join(h(x) for h in self.hs)
        else:
            return tuple([h(x) for h in self.hs])


class NearNeighborLookup:
    """Support for approximate near neighbor lookups.

    This implements the R-near neighbor reporting problem described in
    Andoni and Indyk 2008.
    """

    def __init__(self, family, k, dist_thres, dist_fn, reporting_prob,
            hash_idx=None, join_concat_as_str=False):
        """
        This selects a number of hash tables (defined as L in the above
        reference) according to the strategy it outlines: we want any
        neighbor (within dist_thres) of a query to be reported with
        probability at least reporting_prob; thus, the number of
        tables should be [log_{1 - (P1)^k} (1 - reporting_prob)]. In
        the above reference, delta is 1.0 - reporting_prob.

        Args:
            family: object giving family of hash functions
            k: number of hash functions from family to concatenate
            dist_thres: consider any two objects within this threshold
                of each other to be neighbors
            dist_fn: function f(a, b) that calculates the distance between
                a and b, to compare against dist_thres
            reporting_prob: report any neighbor of a query with
                probability at least equal to this
            hash_idx: if set, the inserted points are tuples and should
                be key'd on the hash_idx'd index; e.g., (A, B, C) might
                be a point and if hash_idx is 0, it is hashed only based on A,
                B and C simply store additional information along with A,
                and queries are based on distance to A
            join_concat_as_str: if True, have concatenated hash functions
                return a string rather than tuple (this can be more
                efficient, but works only if each hash function from the
                family returns a string)
        """
        self.family = family
        self.k = k
        self.dist_thres = dist_thres
        self.dist_fn = dist_fn
        self.hash_idx = hash_idx

        P1 = self.family.P1(dist_thres)
        if P1 == 1.0:
            # dist_thres might be 0, and any number of hash tables can
            # satisfy the reporting probability
            self.num_tables = 1
        else:
            self.num_tables = math.log(1.0 - reporting_prob, 1.0 - math.pow(P1, k))
            self.num_tables = int(math.ceil(self.num_tables))

        # Setup self.num_tables hash tables, each with a corresponding
        # function for hashing into it (the functions are concatenations
        # of k hash functions from the given family)
        self.hashtables = []
        self.hashtables_masked = []
        self.hashtables_g = []
        for j in range(self.num_tables):
            g = HashConcatenation(self.family, self.k,
                join_as_str=join_concat_as_str)
            self.hashtables += [defaultdict(set)]
            self.hashtables_masked += [defaultdict(set)]
            self.hashtables_g += [g]

    def add(self, pts):
        """Insert given points into each of the hash tables.

        Args:
            pts: collection of points (e.g., probes) to add to the hash
                tables
        """
        for j in range(self.num_tables):
            ht = self.hashtables[j]
            g = self.hashtables_g[j].g
            for p in pts:
                p_key = p[self.hash_idx] if self.hash_idx is not None else p
                ht[g(p_key)].add(p)

    def mask(self, mask_idx, mask_val):
        """Mask points from the hash tables that meet given criteria.

        The points stored must be tuples. This moves the points that meet
        the given criteria (according to mask_idx and mask_val) into separate
        hash tables, so that they can be returned to the main hash tables
        when they should be unmasked.

        Args:
            mask_idx: mask according to this index in a tuple specifying
                a point
            mask_val: mask all points p where p[mask_idx] == mask_val
        """
        for j in range(self.num_tables):
            ht = self.hashtables[j]
            ht_masked = self.hashtables_masked[j]
            keys_to_del = set()
            for k in ht.keys():
                # Find all points p in bucket with key k to mask
                p_to_mask = set()
                for p in ht[k]:
                    if p[mask_idx] == mask_val:
                        p_to_mask.add(p)

                for p in p_to_mask:
                    # Delete p from bucket with key k
                    ht[k].remove(p)
                    # Add p to the mask hash table
                    ht_masked[k].add(p)

                if len(ht[k]) == 0:
                    keys_to_del.add(k)
            # Delete empty buckets in ht
            for k in keys_to_del:
                del ht[k]

    def unmask_all(self):
        """Unmask all points that have been masked.

        This moves all points from the mask hash tables into the main
        hash tables.
        """
        for j in range(self.num_tables):
            ht = self.hashtables[j]
            ht_masked = self.hashtables_masked[j]
            for k in ht_masked.keys():
                for p in ht_masked[k]:
                    ht[k].add(p)

            # Reset this mask hash table
            self.hashtables_masked[j] = defaultdict(set)

    def query(self, q):
        """Find neighbors of a query point.

        Args:
            q: query point (e.g., probe); if self.hash_idx is set and
                the inserted points are tuples, q should only be the
                key of what to search for (i.e., distance is measured
                between q and p[self.hash_idx] for a stored point p)

        Returns:
            collection of stored points that are within self.dist_thres of
            q; all returned points are within this distance, but the
            returned points might not include all that are
        """
        neighbors = set()
        for j in range(self.num_tables):
            ht = self.hashtables[j]
            g = self.hashtables_g[j].g
            for p in ht[g(q)]:
                p_key = p[self.hash_idx] if self.hash_idx is not None else p
                if self.dist_fn(q, p_key) <= self.dist_thres:
                    neighbors.add(p)
        return neighbors


# In[ ]:


"""Functions for clustering sequences before input.

This includes computing a distance matrix using MinHash, and
clustering that matrix.
"""

from collections import defaultdict
from collections import OrderedDict
import logging
import operator

import numpy as np
from scipy.cluster import hierarchy

# instead of importing everything with ADAPT, only use lsh.py
# import lsh 
# change this when using as a .py file

from pathlib import Path as _Path
from collections import OrderedDict as _OrderedDict

def _load_fasta_to_dict(_fasta_path):
    """Minimal FASTA reader -> OrderedDict[header] = sequence (uppercased)."""
    if not isinstance(_fasta_path, (str, bytes, _Path)):
        raise TypeError("fasta_path must be a path-like string")
    path = str(_fasta_path)
    seqs = _OrderedDict()
    header = None
    chunks = []
    with open(path, 'r') as _fh:
        for _line in _fh:
            if _line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(chunks).replace(' ', '').replace('\t', '').upper()
                header = _line[1:].strip()
                chunks = []
            else:
                s = _line.strip()
                if s:
                    chunks.append(s)
    if header is not None:
        seqs[header] = ''.join(chunks).replace(' ', '').replace('\t', '').upper()
    return seqs


__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def make_signatures_with_minhash(family, seqs):
    """Construct a signature using MinHash for each sequence.

    Args:
        family: lsh.MinHashFamily object
        seqs: dict mapping sequence header to sequences

    Returns:
        dict mapping sequence header to signature
    """
    # Construct a single hash function; use the same for all sequences
    h = family.make_h()

    signatures = {}
    for name, seq in seqs.items():
        signatures[name] = h(seq)
    return signatures


def create_condensed_dist_matrix(n, dist_fn):
    """Construct a 1d condensed distance matrix for scipy.

    Args:
        n: number of elements whose pairwise distances to store in the
            matrix
        dist_fn: function such that dist_fn(i, j) gives the distance
            between i and j, for all i<j<n

    Returns:
        condensed 1d distance matrix for input to scipy functions
    """
    def idx(i, j):
        # Compute index in 1d vector for pair (i, j)
        return int((-1 * i*i)/2 + i*n - 3*i/2 + j - 1)

    dist_matrix_len = int(n*(n-1)/2)
    dist_matrix = np.zeros(dist_matrix_len)

    for j in range(n):
        for i in range(j):
            dist_matrix[idx(i, j)] = dist_fn(i, j)

    return dist_matrix


def cluster_from_dist_matrix(dist_matrix, threshold, num_clusters=None):
    """Use scipy to cluster a distance matrix.

    Args:
        dist_matrix: distance matrix, represented in scipy's 1d condensed form
        threshold: maximum inter-cluster distance to merge clusters (higher
            results in fewer clusters)
        num_clusters: if set, cluster such that this is the maximum number
            of clusters; threshold is ignored and must be None

    Returns:
        list c such that c[i] is a collection of all the observations
        (whose pairwise distances are indexed in dist) in the i'th
        cluster, in sorted order by descending cluster size
    """
    if threshold is None and num_clusters is None:
        raise ValueError(("threshold or num_clusters must be set"))
    elif threshold is not None and num_clusters is not None:
        raise ValueError(("Only one of threshold or num_clusters can be "
            "set"))

    linkage = hierarchy.linkage(dist_matrix, method='average')
    if threshold:
        clusters = hierarchy.fcluster(linkage, threshold, criterion='distance')
    elif num_clusters:
        clusters = hierarchy.fcluster(linkage, num_clusters,
            criterion='maxclust')

    # clusters are numbered starting at 1, but base the count on
    # first_clust_num just in case this changes
    first_clust_num = min(clusters)
    num_clusters_found = max(clusters) + 1 - first_clust_num
    elements_in_cluster = defaultdict(list)
    for i, clust_num in enumerate(clusters):
        elements_in_cluster[clust_num].append(i)
    cluster_sizes = {c: len(elements_in_cluster[c])
                     for c in range(first_clust_num,
                                    num_clusters_found + first_clust_num)}

    elements_in_cluster_sorted = []
    for clust_num, _ in sorted(cluster_sizes.items(),
            key=operator.itemgetter(1), reverse=True):
        elements_in_cluster_sorted += [elements_in_cluster[clust_num]]
    return elements_in_cluster_sorted


def cluster_with_minhash_signatures(seqs, k=12, N=100, threshold=0.1,
        num_clusters=None, return_dist_matrix_and_indices=False):
    """Cluster sequences based on their MinHash signatures.

    Args:
        seqs: dict mapping sequence header to sequences
        k: k-mer size to use for k-mer hashes (smaller is likely more
            sensitive for divergent genomes, but may lead to false positives
            in determining which genomes are close)
        N: number of hash values to use in a signature (higher is slower for
            clustering, but likely more sensitive for divergent genomes)
        threshold: maximum inter-cluster distance to merge clusters, in
            average nucleotide dissimilarity (1-ANI, where ANI is
            average nucleotide identity); higher results in fewer
            clusters
        num_clusters: if set, cluster such that this is the maximum number
            of clusters; threshold is ignored and must be None
        return_dist_matrix_and_indices: if set, return the pairwise distance
            matrix and sequences indices in each cluster; when used,
            seqs should be an OrderedDict

    Returns:
        if return_dist_matrix_and_indices:
            tuple (dm, c) such that dm is the pairwise distance matrix and
            c is a list such that c[i] gives a collection of sequence
            indices (corresponding to indices in dm) in the same cluster
            (note dm is a 1d condensed matrix in scipy's form); clusters in
            c are sorted in descending order of size
        else:
            list c such that c[i] gives a collection of sequence headers
            in the same cluster, and the clusters in c are sorted
            in descending order of size
    """
    # try FASTA first
    if isinstance(seqs, (str, bytes, _Path)):
        seqs = _load_fasta_to_dict(seqs)
        
    if len(seqs) == 1:
        if return_dist_matrix_and_indices:
            return None, [[0]]
        else:
            return [[list(seqs.keys())[0]]]

#     family = lsh.MinHashFamily(k, N=N) 
    # this is a temporary fix because it's defined in our notebook
    # but when running it in a .py file, comment this lower one out and use the lsh.MinHashFamily 
    family = MinHashFamily(k, N=N)
    signatures_map = make_signatures_with_minhash(family, seqs)

    # Map each sequence header to an index (0-based), and get
    # the signature for the corresponding index
    num_seqs = len(seqs)
    seq_headers = []
    signatures = []
    for name, seq in seqs.items():
        seq_headers += [name]
        signatures += [signatures_map[name]]

    # Eq. 4 of the Mash paper (Ondov et al. 2016) shows that the
    # Mash distance, which is shown to be closely related to 1-ANI, is:
    #  D = (-1/k) * ln(2*j/(1+j))
    # where j is a Jaccard similarity. Solving for j:
    #  j = 1/(2*exp(k*D) - 1)
    # So, for a desired distance D in terms of 1-ANI, the corresponding
    # Jaccard distance is:
    #  1.0 - 1/(2*exp(k*D) - 1)
    # We can use this to calculate a clustering threshold in terms of
    # Jaccard distance
    if threshold is not None:
        jaccard_dist_threshold = 1.0 - 1.0/(2.0*np.exp(k*threshold) - 1)
    else:
        # Ignore inter-cluster distance; use num_clusters instead
        jaccard_dist_threshold = None

    def jaccard_dist(i, j):
        # Return estimated Jaccard dist between signatures at
        # index i and index j
        return family.estimate_jaccard_dist(
            signatures[i], signatures[j])

    dist_matrix = create_condensed_dist_matrix(num_seqs, jaccard_dist)
    clusters = cluster_from_dist_matrix(dist_matrix,
        jaccard_dist_threshold, num_clusters=num_clusters)

    if return_dist_matrix_and_indices:
        return dist_matrix, clusters
    else:
        seqs_in_cluster = []
        for cluster_idxs in clusters:
            seqs_in_cluster += [[seq_headers[i] for i in cluster_idxs]]
        return seqs_in_cluster


def find_representative_sequences(seqs, k=12, N=100, threshold=0.1,
        num_clusters=None, frac_to_cover=1.0):
    """Find a set of representative sequences.

    This clusters seqs, and then determines a medoid for each cluster.
    It returns the medoids.

    This will not return representative sequences with ambiguity or NNNs.

    Args:
        seqs, k, N, threshold, num_clusters: see cluster_with_minhash_signatures()
        frac_to_cover: return medoids from clusters that collectively
            account for at least this fraction of all sequences; this
            allows ignoring representative sequences for outlier
            clusters

    Returns:
        tuple (set of sequence headers representing cluster medoids,
        fraction of all sequences contained in cluster)
    """
    # Accept a FASTA filepath in addition to a dict of sequences
    if isinstance(seqs, (str, bytes, _Path)):
        seqs = _load_fasta_to_dict(seqs)
    seqs = OrderedDict(seqs)
    dist_matrix, clusters = cluster_with_minhash_signatures(
            seqs, k=k, N=N, threshold=threshold, num_clusters=num_clusters,
            return_dist_matrix_and_indices=True)

    seqs_items = list(seqs.items())
    n = len(seqs)
    def idx(i, j):
        # Compute index in 1d vector for pair (i, j)
        if i > j:
            i, j = j, i
        return int((-1 * i*i)/2 + i*n - 3*i/2 + j - 1)

    rep_seqs = []
    rep_seqs_frac = []
    num_seqs_accounted_for = 0
    for cluster_idxs in clusters:
        # Stop if we have already accounted for frac_to_cover of the
        # sequences
        # Note that clusters should be sorted in descending order of
        # size, so any clusters after this one will be the same size
        # or smaller
        if float(num_seqs_accounted_for) / len(seqs) >= frac_to_cover:
            break

        # Find the medoid of this cluster
        # Simply look over all pairs in the cluster (there are faster
        # algorithms, though not linear)
        curr_medoid = None
        curr_medoid_dist_total = None
        for i in cluster_idxs:
            # Only allow i to be the medoid if it does not have ambiguity
            seq = seqs_items[i][1]
            if sum(seq.count(b) for b in ('A','C','G','T')) != len(seq):
                # Has ambiguity or NNNs; skip
                continue

            # Compute the total distance to all other sequences in this
            # cluster, and check if this is the medoid
            dist_total = 0
            for j in cluster_idxs:
                if i == j:
                    continue
                dist_total += dist_matrix[idx(i,j)]
            if curr_medoid is None or dist_total < curr_medoid_dist_total:
                curr_medoid = i
                curr_medoid_dist_total = dist_total
        if curr_medoid is not None:
            rep_seqs += [curr_medoid]
            rep_seqs_frac += [float(len(cluster_idxs)) / len(seqs)]
            num_seqs_accounted_for += len(cluster_idxs)
        else:
            # All sequences have ambiguity or NNNs; raise a warning and
            # skip this cluster
            logger.warning(("Cannot find medoid for cluster of size %d "
                "because all sequences have ambiguity or NNNs; skipping "
                "this cluster"),
                len(cluster_idxs))

    return ([seqs_items[i][0] for i in rep_seqs], rep_seqs_frac)


# In[ ]:


def parse_params(paramFile):
    params = {}
    for l in open(paramFile):
        if '=' in l:
            name, value = l.split('=')
            try:
                params[name.strip()] = int(value.strip())
            except ValueError:
                params[name.strip()] = value.strip()
    window = params['DESIGN_WINDOW']
    return window


# In[ ]:


def main():
    parser = argparse.ArgumentParser(prog='python -u generate_primers.py', 
                                     description='Generate primer candidates')
    parser.add_argument('--in', dest='msa', required=True, help='FASTA file of target sequences')
    parser.add_argument('--out', dest='representative', required=True, help='FASTA file of primers')
    parser.add_argument('--params', dest='param_file', required=True, help='Parameters file')
    parser.add_argument('--name', dest='name', required=True, help='Pathogen name')
    args = parser.parse_args()
    window = parse_params(args.param_file)
    alns = [ (s.id, str(s.seq)) for s in SeqIO.parse(args.msa, 'fasta') ]
    
    print(f'Picking representatives out of {len(alns)} sequences in {args.msa}...')
    startTime = time.time()
    
    alnlen = len(alns[0][1])
    gapcnts = {}
    for st in range(0, alnlen, window//2):
        if st + window > alnlen:
            st = alnlen - window
            if st in gapcnts:
                break
        gapcnts[st] = [ aln[st:st+window].count('-') for sid, aln in alns ]
    
    bestSt = min(gapcnts.keys(), key=lambda x:(np.median(gapcnts[x]), np.average(gapcnts[x])))
    trimFa = args.representative + '.tmp'
    with open(trimFa, 'wt') as out:
        for sid, aln in alns:
            trimseq = aln[bestSt:].replace('-','')[:window].upper()
            out.write(f'>{sid}\n{trimseq}\n')
    
    trimseqs = { s.id:str(s.seq) for s in SeqIO.parse(trimFa, 'fasta') }
    with open(args.representative, 'wt') as out:
        sids, fracs = find_representative_sequences(trimFa, frac_to_cover=.99)
        for sid, frac in zip(sids, fracs):
            out.write(f'>{sid}\tcov:%.3f\n{trimseqs[sid]}\n' % frac)
            
    runtime = (time.time() - startTime)
    print(f'Wrote {len(sids)} sequences to {args.representative} (%.1f sec).' % runtime)
    
    
if __name__ == '__main__':
    main()

