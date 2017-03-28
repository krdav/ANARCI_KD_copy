#    ANARCI - Antibody Numbering and Antigen Receptor ClassIfication
#    Copyright (C) 2016 Oxford Protein Informatics Group (OPIG)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.#
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
ANARCI - Antigen Receptor Numbering And ClassIfication

Oxford Protein Informatics Group (OPIG). 2015

ANARCI performs alignments of sequences to databases of Hidden Markov Models (HMMs).
Those that align with a significant score are classified by species and chain type.
They are then numbered with a scheme of the user's choosing. 

Currently implemented schemes: 
    IMGT
    Chothia (IGs only)
    Kabat (IGs only)
    Martin / Enhanced Chothia (IGs only)

Currently recognisable species (chains):
    Human (heavy, kappa, lambda, alpha, beta)
    Mouse (heavy, kappa, lambda, alpha, beta)
    Rat (heavy, kappa, lambda)
    Rabbit (heavy, kappa, lambda)
    Pig (heavy, kappa, lambda)
    Rhesus Monkey (heavy, kappa)
    
Notes:
 o Currently, mouse and rat are not reliably distinguishable (use assign_germline to
   refine the species annotation.
 o The AHo scheme is implemented heuristically.
 o Fragments of antibodies can be recognised. Play with the bit_score_threshold in
   _parse_hmmer_query if you want to investigate this feature.
 o Really long CDR3s will not be identified. although the v-gene segment will be.
 o Unusual N and C terminal residues are often not considered as part of the domain by
   HMMER. Thus, they are not included in the numbering.


'''

import os
import sys
import tempfile
from textwrap import wrap
from subprocess import Popen, PIPE
from itertools import groupby

# Import the HMMER parser from the distributed version of Biopython.
from .Bio.SearchIO.HmmerIO import Hmmer3TextParser as HMMERParser

# Import from the schemes submodule
from schemes import *
from germlines import all_germlines
all_species = all_germlines['V']['H'].keys()

amino_acids = sorted(list("QWERTYIPASDFGHKLCVNM"))
set_amino_acids = set(amino_acids)
anarci_path  = os.path.split(__file__)[0]

scheme_short_to_long = { "m":"martin", "c":"chothia", "k":"kabat","imgt":"imgt", "kabat":"kabat", "chothia":"chothia", "martin":"martin", "i":"imgt", "a":"aho","aho":"aho","wolfguy":"wolfguy", "w":"wolfguy"}

scheme_names = scheme_short_to_long.keys() 
chain_type_to_class = {"H":"H", "K":"L", "L":"L", "A":"A", "B":"B", "G":"G", "D":"D"}

HMM_path =  os.path.join( anarci_path, "dat", "HMMs" )

all_reference_states = range( 1, 129 ) # These are the IMGT reference states (matches)

class HMMscanError(Exception):
    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        super(HMMscanError, self).__init__(message)

## Utility functions ##
def read_fasta(filename):
    """
    Read a sequence file and parse as description, string 
    """
    return [ r for r in fasta_iter(filename) ]

def fasta_iter(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    https://www.biostars.org/p/710/
    """
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


def write_fasta(sequences, f):
    """
    Write a list of sequences to file. 
    
    should be a list of name, sequence tuples
    
    f should be an open file
    """
    for name, sequence in sequences:
        print >> f, ">%s"%name
        print >> f, '\n'.join(['\n'.join(wrap(block, width=80)) for block in sequence.splitlines()])
    
    
def validate_sequence(sequence):
    """
    Check whether a sequence is a protein sequence or if someone has submitted something nasty.
    """
    assert len(sequence) < 10000, "Sequence too long."
    assert not (set( sequence.upper() ) - set_amino_acids ), "Unknown amino acid letter found in sequence: %s"% ", ".join(list((set( sequence.upper() ) - set_amino_acids )))    
    return True

def validate_numbering((numbering, start, end), name_seq=[]):
    """
    Wrapper to do some basic validation of the numbering.
    
    Further validation could be done but at the moment we just check that the numbering indices are incremental (they should be)
    """

    # Validation.
    # 1. The indices of the numbering increase
    # 2. A contiguous segment of the sequence has been numbered. # Commented out now - should be added if further developed for sanity check
#    name, seq = name_seq
#    last = -1
#    nseq=""
#    for (index, _), a in numbering:
#        assert index >= last, "Numbering was found to decrease along the sequence %s. Something is wrong...debugging time"%name
#        last = index
#        nseq +=a.replace("-","")
#    assert nseq in seq, "The algorithm did not number a contiguous segment for sequence %s ... debugging time"%name
    return numbering, start, end

def anarci_output(numbered, sequences, alignment_details, outfile, sequence_id=None, domain_id=None):
    """
    Output to open file

    If sequence_id is specified as an integer then only this sequence will be printed. 
    Otherwise all sequences will be printed.

    If domain_id is specified as an integer then only this domain will be printed. 
    Otherwise all domains will be printed.

    If domain_id is specified then sequence_id must also be specified. 
    """       
    assert (sequence_id is not None) or (sequence_id is None and domain_id is None), "If domain_id is specified, sequence_id must also be specified."

    for i in xrange(len(numbered)):
        if sequence_id is None:
            print >> outfile, "# %s"%sequences[i][0] # print the name
        if numbered[i] is not None:
            if sequence_id is not None:
                if i != sequence_id: continue
            print >> outfile, "# ANARCI numbered"
            for j in xrange( len(numbered[i])): # Iterate over domains
                if domain_id is not None:
                    if j != domain_id: continue
                print >> outfile, "# Domain %d of %d"%(j+1, len(numbered[i]) )
                print >> outfile, "# Most significant HMM hit"
                print >> outfile, "#|species|chain_type|e-value|score|seqstart_index|seqend_index|"
                alignment_details[i][j]["evalue"] = str( alignment_details[i][j]["evalue"] )
                print >> outfile, "#|%s|%s|%s|%.1f|%d|%d|"%tuple( [alignment_details[i][j][field] for field in 
                                                                     ["species","chain_type","evalue","bitscore"]] 
                                                                   +[ numbered[i][j][1], numbered[i][j][2]] )
                
                if 'germlines' in alignment_details[i][j]:
                    print >> outfile, '# Most sequence-identical germlines'
                    print >> outfile, '#|species|v_gene|v_identity|j_gene|j_identity|'
                    (species, vgene), vid =alignment_details[i][j]['germlines']['v_gene']
                    if vgene is None:
                        vgene, vid = 'unknown', 0
                    (_,jgene), jid =alignment_details[i][j]['germlines']['j_gene']
                    if jgene is None:
                        jgene, jid = 'unknown', 0
                    print >> outfile, '#|%s|%s|%.2f|%s|%.2f|'%(species, vgene, vid, jgene, jid )	
                chain_type = chain_type_to_class[  alignment_details[i][j]["chain_type"] ]
                print >> outfile, "# Scheme = %s"%alignment_details[i][j]["scheme"]
                for (index, insertion), aa in numbered[i][j][0]:
                    print >> outfile, chain_type, ("%d"%index).ljust(5), insertion, aa
        print >> outfile, "//"


## Parsing and recognising domain hits from hmmscan ##

def _domains_are_same(dom1, dom2):
    """
    Check to see if the domains are overlapping.
    @param dom1: 
    @param dom2: 

    @return: True or False  
    """
    dom1, dom2 = sorted( [dom1, dom2], key=lambda x: x.query_start  )
    if dom2.query_start >= dom1.query_end:
        return False
    return True


def _parse_hmmer_query(query, bit_score_threshold=80):
    """
    
    @param query: hmmer query object from Biopython
    @param bit_score_threshold: the threshold for which to consider a hit a hit. 
    
    The function will identify multiple domains if they have been found and provide the details for the best alignment for each domain.
    This allows the ability to identify single chain fvs and engineered antibody sequences as well as the capability in the future for identifying constant domains. 

    """
    hit_table = [ ['id', 'description', 'evalue', 'bitscore', 'bias', 
                    'query_start', 'query_end' ] ]

    # Find the best hit for each domain in the sequence.

    top_descriptions, domains,state_vectors = [], [], []

    if query.hsps: # We have some hits
        for hsp in sorted(query.hsps, key=lambda x: x.evalue): # Iterate over the matches of the domains in order of their e-value (most significant first)
            new=True
            if hsp.bitscore >= bit_score_threshold: # Only look at those with hits that are over the threshold bit-score.
                for i in xrange( len(domains) ): # Check to see if we already have seen the domain
                    if _domains_are_same( domains[i], hsp ):
                        new = False
                        break      
                hit_table.append( [ hsp.hit_id, hsp.hit_description, hsp.evalue, hsp.bitscore, hsp.bias, hsp.query_start, hsp.query_end] )
                if new: # It is a new domain and this is the best hit. Add it for further processing.
                    domains.append( hsp )
                    top_descriptions.append(  dict( zip(hit_table[0], hit_table[-1]) ) ) # Add the last added to the descriptions list. 

        # Reorder the domains according to the order they appear in the sequence.         
        ordering = sorted( range(len(domains)), key=lambda x: domains[x].query_start)
        domains = [ domains[_] for _ in ordering ]
        top_descriptions = [ top_descriptions[_] for _ in ordering ]         

    for i in xrange(len(domains)): # If any significant hits were identified parse and align them to the reference state.
        domains[i].order = i
        species, chain = top_descriptions[i]["id"].split("_")
        state_vectors.append( _hmm_alignment_to_states(domains[i]) ) # Alignment to the reference states.
        top_descriptions[i][ "species"] = species # Reparse
        top_descriptions[i][ "chain_type"] = chain
        top_descriptions[i][ "query_start"] = state_vectors[-1][0][-1] # Make sure the query_start agree if it was changed

    return hit_table, state_vectors, top_descriptions


def _hmm_alignment_to_states(hsp):
    """
    Take a hit hsp and turn the alignment into a state vector with sequence indices
    """

    # Extract the strings for the reference states and the posterior probability strings     
    reference_string = hsp.aln_annotation["RF"]
    state_string = hsp.aln_annotation["PP"]
    
    assert len(reference_string) == len(state_string), "Aligned reference and state strings had different lengths. Don't know how to handle"
    
    # Extract the start an end points of the hmm states and the sequence
    # These are python indices i.e list[ start:end ] and therefore start will be one less than in the text file
    _hmm_start = hsp.hit_start
    _hmm_end = hsp.hit_end
     
    _seq_start = hsp.query_start
    _seq_end = hsp.query_end

    # Handle cases where there are n terminal modifications.
    # In most cases the user is going to want these included in the numbered domain even though they are not 'antibody like' and 
    # not matched to the germline. Only allow up to a maximum of 5 unmatched states at the start of the domain
    # Adds a bug here if there is a very short linker between a scfv domains with a modified n-term second domain
    # Thus this is only done for the first identified domain ( hence order attribute on hsp )
    if hsp.order == 0 and _hmm_start and _hmm_start < 5: 
        n_extend = _hmm_start 
        if _hmm_start > _seq_start:
            n_extend = min( _seq_start , _hmm_start - _seq_start )
        state_string = '8'*n_extend + state_string  
        reference_string = 'x'*n_extend + reference_string
        _seq_start = _seq_start - n_extend
        _hmm_start = _hmm_start - n_extend


    # Generate lists for the states and the sequence indices that are included in this alignment
    hmm_states = all_reference_states[ _hmm_start : _hmm_end ] 
    sequence_indices = range(_seq_start,  _seq_end)
    h, s = 0, 0 # initialise the current index in the hmm and the sequence
    
    state_vector = []
    # iterate over the state string (or the reference string)
    for i in xrange( len(state_string) ):
        if reference_string[i] == "x": # match state
            state_type = "m"
        else: # insert state
            state_type = "i"
        
        if state_string[i] == ".": # overloading if deleted relative to reference. delete_state
            state_type = "d"
            sequence_index = None
        else:
            sequence_index = sequence_indices[s]    
        # Store the alignment as the state identifier (uncorrected IMGT annotation) and the index of the sequence
        
        state_vector.append(  ((hmm_states[h], state_type),  sequence_index )  )        

        # Updates to the indices         
        if state_type == "m":
            h+=1
            s+=1
        elif state_type == "i":
            s+=1
        else: # delete state
            h+=1
        
    return state_vector


def parse_hmmer_output(filedescriptor="", bit_score_threshold=80):
    """
    Parse the output of HMMscan and return top alignment and the score table for each input sequence.
    """
    results  = []
    if type(filedescriptor) is str:
        openfile = open
    elif type(filedescriptor) is int:
        openfile = os.fdopen
    
    with openfile(filedescriptor) as inputfile:
        p = HMMERParser( inputfile )
        for query in p:
            results.append(_parse_hmmer_query(query,bit_score_threshold=bit_score_threshold ))
    return results


def run_hmmer(sequence_list,hmm_database="ALL",hmmerpath="", ncpu=None, bit_score_threshold=80):
    """
    Run the sequences in sequence list against a precompiled hmm_database.

    Those sequence that have a significant hit with a bit score over a threshold will
    be recognised and an alignment given. The alignment will be used to number the 
    sequence.

    @param sequence_list: a list of (name, sequence) tuples. Both are strings
    @param hmm_database: The hmm database to use. Currently, all hmms are in the ALL database.
                         The code to develop new models is in build_pipeline in the git repo.
    @param hmmerpath: The path to hmmer binaries if not in the path
    @param ncpu: The number of cpu's to allow hmmer to use.
    """

    # Check that hmm_database is available
    
    assert hmm_database in ["ALL"], "Unknown HMM database %s"%hmm_database    
    HMM = os.path.join( HMM_path, "%s.hmm"%hmm_database )


    # Create a fasta file for all the sequences. Label them with their sequence index
    # This will go to a temp file
    fasta_filehandle, fasta_filename =  tempfile.mkstemp( ".fasta", text=True )
    with os.fdopen(fasta_filehandle,'w') as outfile:
        write_fasta(sequence_list, outfile)

    output_filehandle, output_filename =  tempfile.mkstemp( ".txt", text=True )

    # Run hmmer as a subprocess
    if hmmerpath:
        hmmscan = os.path.join(hmmerpath,"hmmscan")
    else:
        hmmscan = "hmmscan"
    try:
        if ncpu is None:
            command = [ hmmscan, "-o", output_filename, HMM,  fasta_filename]
        else:
            command = [ hmmscan, "-o", output_filename, "--cpu", str(ncpu), HMM,  fasta_filename]
        process = Popen( command, stdout=PIPE, stderr=PIPE  )
        _, pr_stderr = process.communicate()

        if pr_stderr:
            _f = os.fdopen(output_filehandle) # This is to remove the filedescriptor from the os. I have had problems with it before.
            _f.close()
            
            raise HMMscanError(pr_stderr)
        results = parse_hmmer_output(output_filehandle, bit_score_threshold=bit_score_threshold)
        
    finally:
        # clear up
        os.remove(fasta_filename)
        os.remove(output_filename)
        
    return results


def number_sequence_from_alignment(state_vector, sequence, scheme="imgt", chain_type=None):
    """
    Given you have an alignment. Give back the numbering
    
    @param state_vector: List of states from the hmm. Effectively these are imgt columns but CDR3 has not been redone. 
    @param sequence: The original sequence string or list.
    @param scheme: The numbering scheme to apply
    @param chain_type: The type of chain to apply numbering for. Some schemes do not require this (IMGT). Others (e.g. Chothia/Wolfguy) do.
    
    @return: A list of numbering identifier / amino acids tuples over the domain that has been numbered. The indices of the start (inclusive) and end point (exclusive) in the sequence for the numbering 
    """
    scheme=scheme.lower()
    if scheme == "imgt":
        return number_imgt(state_vector, sequence)
    elif scheme == "chothia":
        if chain_type == "H":
            return number_chothia_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_chothia_light(state_vector, sequence)
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))
    elif scheme == "kabat":
        if chain_type == "H":
            return number_kabat_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_kabat_light(state_vector, sequence)
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))
    elif scheme == "martin":
        if chain_type == "H":
            return number_martin_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_martin_light(state_vector, sequence)
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))
    elif scheme == "aho":
        return number_aho(state_vector, sequence, chain_type) # requires the chain type to heuristically put the CDR1 gap in position.
    elif scheme == "wolfguy":
        if chain_type == "H":
            return number_wolfguy_heavy( state_vector, sequence )     
        elif chain_type in "KL":
            return number_wolfguy_light( state_vector, sequence )  
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))      
    else:
        raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))


def get_identity( state_sequence, germline_sequence ):
    """
    Get the partially matched sequence identity between two aligned sequences. 
    Partial in the sense that gaps can be in the state_sequence.
    """
    # Ensure that the sequences are the expected length
    assert len( state_sequence) == len(germline_sequence ) == 128
    n, m = 0, 0
    for i in xrange( 128 ):
        if germline_sequence[i] == "-":continue
        if state_sequence[i].upper() == germline_sequence[i]: m+=1
        n+=1

    if not n:
        return 0    
    return float(m)/n
    

def run_germline_assignment(state_vector, sequence, chain_type, allowed_species=None ):
    """
    Find the closest sequence identity match.
    """
    genes={'v_gene': [None,None],
           'j_gene': [None,None],
         }


    # Extract the positions that correspond to match (germline) states. 
    state_dict = dict( ((i, 'm'),None) for i in xrange(1,129))
    state_dict.update(dict(state_vector))
    state_sequence = "".join([ sequence[state_dict[(i, 'm')]] if state_dict[(i,'m')] is not None else "-" for i in xrange(1,129) ])

    # Iterate over the v-germline sequences of the chain type of interest.
    # The maximum sequence identity is used to assign the germline 
    if chain_type in all_germlines["V"]:
        if allowed_species is not None:
            assert all( [ sp in all_germlines['V'][chain_type] for sp in allowed_species ] ), 'Unknown species specified. Choose from %s'%','.join( all_species )
        else:
            allowed_species = all_species
        seq_ids = {}
        for species in allowed_species:
            if species not in all_germlines["V"][ chain_type ]: continue # Previously bug.
            for gene, germline_sequence in all_germlines["V"][ chain_type ][ species ].iteritems():
                seq_ids[ (species, gene) ] = get_identity( state_sequence , germline_sequence )
        genes['v_gene' ][0] = max( seq_ids, key=lambda x: seq_ids[x] )
        genes['v_gene' ][1] = seq_ids[ genes['v_gene' ][0] ]
        
        # Use the assigned species for the v-gene for the j-gene. 
        # This assumption may affect exotically engineered abs but in general is fair.
        species = genes['v_gene' ][0][0]       
        if chain_type in all_germlines["J"]:
            if species in all_germlines["J"][chain_type]:
                seq_ids = {}
                for gene, germline_sequence in all_germlines["J"][ chain_type ][ species ].iteritems():
                    seq_ids[ (species, gene) ] = get_identity( state_sequence , germline_sequence )
                genes['j_gene' ][0] = max( seq_ids, key=lambda x: seq_ids[x] )
                genes['j_gene' ][1] = seq_ids[ genes['j_gene' ][0] ]
     
    return genes



##################################
# High level numbering functions #
##################################

# Main function for ANARCI
def anarci(sequences, scheme="imgt", database="ALL", output=False, outfile=None, allow=set(["H","K","L","A","B","G","D"]), hmmerpath="", ncpu=None, assign_germline=False, allowed_species=None, bit_score_threshold=80):
    """
    The main function for anarci. This can be used to identify domains, number them and classify their receptor type (and species).
    Multiple sequences can be handed to the function at once. 
    For a more basic interface to the algorithm see the "number" function                
    
    @param sequences: A list of name and amino acid sequence string tuples. 
                      e.g. [ ("seq1","EVQLQQSGAEVVRSG ..."),
                             ("seq2","DIVMTQSQKFMSTSV ...") ]
    @param scheme:    The numbering scheme that should be applied. Choose from imgt, chothia, kabat or martin
    @param database:  The HMMER database that should be used. Normally not changed unless a custom db is created.
    @param output:    Boolean flag to say whether the result should be output.
    @param outfile:   The name of the file to output to. If output is True and outfile is None then output is printed
                      to stdout.
    @param allow:     A set containing the chain types that should be recognised. If chothia, kabat or martin is used
                      as the scheme, anarci will ignore tcr chains.
    @param hmmerpath: The path to hmmscan. If left unspecified then the PATH will be searched. 
    @param ncpu:      The number of cpu's that hmmer should be allowed to use. If not specified then the hmmscan 
                      default is used.

    @return: Three lists. Numbered, Alignment_details and Hit_tables.
             Each list is in the same order as the input sequences list.
             A description of each entry in the three lists is as followed.
               o Numbered: will be None if no domain was found for that sequence or a list of domains with their 
                           numbering, start and finish indices.
               o Alignment_details: will be None if no domain was found for that sequence or a dictionary for each
                           domain identified containing the details of the alignment (chain type, e-value, species etc).
               o Hit_tables: None if no domain was found for that sequence or a nested list for each domain containing
                           the hit table from hmmscan.
    
    """
    
    # Validate the input scheme
    try:
        scheme = scheme_short_to_long[scheme.lower()]
    except KeyError:
        raise AssertionError, "Unrecognised or unimplemented scheme: %s"%scheme        
        
    # Perform the alignments of the sequences to the hmm database
    alignments = run_hmmer( sequences, hmm_database=database, hmmerpath=hmmerpath, ncpu=ncpu, bit_score_threshold=bit_score_threshold )    
    
    # Iteration over the sequence alignments performing the desired numbering 
    numbered = []
    alignment_details = []
    hit_tables = []
    for i in xrange(len(sequences)):
        # Unpack
        hit_table, state_vectors, detailss = alignments[i] # We may have multiple domains per sequence (e.g. single chain fvs). 

        # Iterate over all the domains in the sequence that have been recognised (typcially only 1 with the current hmms available)
        hit_numbered, hit_details = [], []
        for di in xrange( len( state_vectors ) ):
            state_vector = state_vectors[di]
            details      = detailss[di]
            details["scheme"]=scheme
            details["query_name"]=sequences[i][0]            
            if state_vector and details["chain_type"] in allow: # only number things that are allowed. We still keep the alignment details and hit_table
                try:
                    # Do the numbering and validate (for development purposes)
                    hit_numbered.append( validate_numbering(number_sequence_from_alignment(state_vector, sequences[i][1], scheme=scheme, chain_type=details["chain_type"]), sequences[i] ) )
                    if assign_germline:
                        details["germlines"] = run_germline_assignment( state_vector, sequences[i][1], details["chain_type"], allowed_species=allowed_species)
                    hit_details.append( details )
                except AssertionError, e: # Handle errors. Those I have implemented should be assertion.
                    print >> sys.stderr, str(e)
                    raise e # Validation went wrong. Error message will go to stderr. Want this to be fatal during development.
                except Exception, e:
                    print >> sys.stderr, "Error: Something really went wrong that has not been handled"
                    print >> sys.stderr, str(e)
                    raise e
                
        if hit_numbered: 
            numbered.append( hit_numbered )
            alignment_details.append( hit_details )
        else: 
            numbered.append( None )
            alignment_details.append( None )
        hit_tables.append(hit_table)
    

    if output: 
        outto, close=sys.stdout, False
        if outfile:
            outto, close = open(outfile,'w'), True
        anarci_output(numbered, sequences, alignment_details, outto)
        if close:
            outto.close()

    return numbered, alignment_details, hit_tables
                

# Wrapper function for ANARCI for back compatibility with sabdab annotate functions. 
def number(sequence, scheme="imgt", database="ALL", allow=set(["H","K","L","A","B","G","D"])):
    """
    Given a sequence string, use anarci to number it using the scheme of choice.
    Only the first domain will be recognised and numbered

    @param sequence: An amino acid sequence string
    @param scheme: The numbering scheme that should be applied. Choose from imgt, chothia, kabat or martin
    @param database: The HMMER database that should be used. Normally not changed unless a custom db is created.
    @param allow: A set containing the chain types that should be recognised. If chothia, kabat or martin is used
                  as the scheme, anarci will ignore tcr chains.

    @return: If the sequence can be numbered, a list containing the numbering and sequence; and the chain type. 
             Otherwise both are False.
    """
    
    try:
        validate_sequence( sequence )  
        scheme = scheme_short_to_long[scheme.lower()]
    except KeyError:
        raise AssertionError, "Unrecognised to unimplemented scheme: %s"%scheme
    
    if len(sequence) < 70: # Length check. ANARCI can number fragments of chains well. Encourage full domain numbering. 
        return False, False
   
    try:
        numbered, alignment_details, _ = anarci( [("sequence_0", sequence)], scheme=scheme, database=database, output=False, allow=allow )
    except AssertionError: # Catch where the user has tried to number a TCR with an antibody scheme
        return False, False
    

    # We return the numbering list and the chain type where kappa and lambda chains are both "L" for light
    if numbered[0]:
        return numbered[0][0][0], chain_type_to_class[alignment_details[0][0]["chain_type"]]
    else:
        return False, False

if __name__ == "__main__":
    # Test and example useage of the anarci function. 
    sequences = [ ("12e8:H","EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAP"),
                  ("12e8:L","DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASV"),
                  ("scfv:A","DIQMTQSPSSLSASVGDRVTITCRTSGNIHNYLTWYQQKPGKAPQLLIYNAKTLADGVPSRFSGSGSGTQFTLTISSLQPEDFANYYCQHFWSLPFTFGQGTKVEIKRTGGGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFDFSRYDMSWVRQAPGKRLEWVAYISSGGGSTYFPDTVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCARQNKKLTWFDYWGQGTLVTVSSHHHHHH"),
                  ("lysozyme:A","KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL")]

    results = anarci(sequences, scheme="imgt", output=True)
    numbering, alignment_details, hit_tables = results

    expect_one_VH_domain_numbering, expect_one_VL_domain_numbering, expect_VH_then_VL_numbering, expect_None = numbering
    assert  len(expect_one_VH_domain_numbering) == 1
    assert  len(expect_one_VL_domain_numbering) == 1
    assert  len(expect_VH_then_VL_numbering)    == 2
    assert  expect_None                         == None




