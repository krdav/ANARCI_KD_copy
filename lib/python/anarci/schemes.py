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
Module containing functions to convert hmm alignment to a numbering scheme. 

Currently implemented

For IG's
IMGT
Chothia
Kabat
Martin (Extended Chothia)
Aho 
Wolfguy

For TR's
IMGT
(Aho)

---------------------------------------------------------------------------------------------------------------------
Functions are written to a template:

There are 128 match states in the HMMs (these are the IMGT states). The alignment to these states must be converted to
correspond to the scheme of choice. 

We define:
  - a state string consisting of 'X' and 'I' where:
    X  means that for the state there is an equivalent position in the numbering scheme.
    I  means that for the state there is not an equivalent position in the numbering scheme. It should therefore be 
       considered as an insertion in the scheme.
       
  - a region string consisting of characters (integers in the currently implemented schemes). Each character 
corresponds to a contiguous region. Therefore each state can be assigned a region according to the scheme. 

  - a mapping between region characters and region indices as a dictionary. e.g. the first region character maps
to 0, second to 1 ...

  - a dictionary containing the difference between state number (imgt) and scheme number at the *beginning* of 
each region using the region indices as keys and the difference as values.

  - the number of regions defined

  - a list for which delete states should not be included in the numbering (typically those for the cdrs). This
will allow the length of the region to be the number of residues found instead of the number of possible states plus
insertions.

 
This all goes into the _number_regions function along with the sequence and the state_vector (the alignment from the
HMM).

_number regions will then divide the aligned part of the sequence into as many regions as defined above. Within each 
region it will give a numbering according to the input parameters. A list of lists will be returned containing the 
numbered sequence for each region.

Some of the regions will not be numbered correctly according to the scheme. For example the insertions for the CDRs
will not necessarily be on the correct residue. For each different scheme these regions are then modified (see code
for implementation) 

Finally the full numbered sequence is compiled and returned to the calling function.
---------------------------------------------------------------------------------------------------------------------

Other schemes can be implemented following the template above. 


'''

# Alphabet used for insertion (last (-1th) is a blank space for no insertion)
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ "

# Blosum62 matrix. Used in some annotation methods to recognise pre-defined motifs
blosum62 = {('B', 'N'): 3, ('W', 'L'): -2, ('G', 'G'): 6, ('X', 'S'): 0, ('X', 'D'): -1, ('K', 'G'): -2, ('S', 'E'): 0, ('X', 'M'): -1, ('Y', 'E'): -2, ('W', 'R'): -3, ('I', 'R'): -3, ('X', 'Z'): -1, ('H', 'E'): 0, ('V', 'M'): 1, ('N', 'R'): 0, ('I', 'D'): -3, ('F', 'D'): -3, ('W', 'C'): -2, ('N', 'A'): -2, ('W', 'Q'): -2, ('L', 'Q'): -2, ('S', 'N'): 1, ('Z', 'K'): 1, ('V', 'N'): -3, ('Q', 'N'): 0, ('M', 'K'): -1, ('V', 'H'): -3, ('G', 'E'): -2, ('S', 'L'): -2, ('P', 'R'): -2, ('D', 'A'): -2, ('S', 'C'): -1, ('E', 'D'): 2, ('Y', 'G'): -3, ('W', 'P'): -4, ('X', 'X'): -1, ('Z', 'L'): -3, ('Q', 'A'): -1, ('V', 'Y'): -1, ('W', 'A'): -3, ('G', 'D'): -1, ('X', 'P'): -2, ('K', 'D'): -1, ('T', 'N'): 0, ('Y', 'F'): 3, ('W', 'W'): 11, ('Z', 'M'): -1, ('L', 'D'): -4, ('M', 'R'): -1, ('Y', 'K'): -2, ('F', 'E'): -3, ('M', 'E'): -2, ('S', 'S'): 4, ('X', 'C'): -2, ('Y', 'L'): -1, ('H', 'R'): 0, ('P', 'P'): 7, ('K', 'C'): -3, ('S', 'A'): 1, ('P', 'I'): -3, ('Q', 'Q'): 5, ('L', 'I'): 2, ('P', 'F'): -4, ('B', 'A'): -2, ('Z', 'N'): 0, ('M', 'Q'): 0, ('V', 'I'): 3, ('Q', 'C'): -3, ('I', 'H'): -3, ('Z', 'D'): 1, ('Z', 'P'): -1, ('Y', 'W'): 2, ('T', 'G'): -2, ('B', 'P'): -2, ('P', 'A'): -1, ('C', 'D'): -3, ('Y', 'H'): 2, ('X', 'V'): -1, ('B', 'B'): 4, ('Z', 'F'): -3, ('M', 'L'): 2, ('F', 'G'): -3, ('S', 'M'): -1, ('M', 'G'): -3, ('Z', 'Q'): 3, ('S', 'Q'): 0, ('X', 'A'): 0, ('V', 'T'): 0, ('W', 'F'): 1, ('S', 'H'): -1, ('X', 'N'): -1, ('B', 'Q'): 0, ('K', 'A'): -1, ('I', 'Q'): -3, ('X', 'W'): -2, ('N', 'N'): 6, ('W', 'T'): -2, ('P', 'D'): -1, ('B', 'C'): -3, ('I', 'C'): -1, ('V', 'K'): -2, ('X', 'Y'): -1, ('K', 'R'): 2, ('Z', 'R'): 0, ('W', 'E'): -3, ('T', 'E'): -1, ('B', 'R'): -1, ('L', 'R'): -2, ('Q', 'R'): 1, ('X', 'F'): -1, ('T', 'S'): 1, ('B', 'D'): 4, ('Z', 'A'): -1, ('M', 'N'): -2, ('V', 'D'): -3, ('F', 'A'): -2, ('X', 'E'): -1, ('F', 'H'): -1, ('M', 'A'): -1, ('K', 'Q'): 1, ('Z', 'S'): 0, ('X', 'G'): -1, ('V', 'V'): 4, ('W', 'D'): -4, ('X', 'H'): -1, ('S', 'F'): -2, ('X', 'L'): -1, ('B', 'S'): 0, ('S', 'G'): 0, ('P', 'M'): -2, ('Y', 'M'): -1, ('H', 'D'): -1, ('B', 'E'): 1, ('Z', 'B'): 1, ('I', 'E'): -3, ('V', 'E'): -2, ('X', 'T'): 0, ('X', 'R'): -1, ('R', 'R'): 5, ('Z', 'T'): -1, ('Y', 'D'): -3, ('V', 'W'): -3, ('F', 'L'): 0, ('T', 'C'): -1, ('X', 'Q'): -1, ('B', 'T'): -1, ('K', 'N'): 0, ('T', 'H'): -2, ('Y', 'I'): -1, ('F', 'Q'): -3, ('T', 'I'): -1, ('T', 'Q'): -1, ('P', 'L'): -3, ('R', 'A'): -1, ('B', 'F'): -3, ('Z', 'C'): -3, ('M', 'H'): -2, ('V', 'F'): -1, ('F', 'C'): -2, ('L', 'L'): 4, ('M', 'C'): -1, ('C', 'R'): -3, ('D', 'D'): 6, ('E', 'R'): 0, ('V', 'P'): -2, ('S', 'D'): 0, ('E', 'E'): 5, ('W', 'G'): -2, ('P', 'C'): -3, ('F', 'R'): -3, ('B', 'G'): -1, ('C', 'C'): 9, ('I', 'G'): -4, ('V', 'G'): -3, ('W', 'K'): -3, ('G', 'N'): 0, ('I', 'N'): -3, ('Z', 'V'): -2, ('A', 'A'): 4, ('V', 'Q'): -2, ('F', 'K'): -3, ('T', 'A'): 0, ('B', 'V'): -3, ('K', 'L'): -2, ('L', 'N'): -3, ('Y', 'N'): -2, ('F', 'F'): 6, ('L', 'G'): -4, ('B', 'H'): 0, ('Z', 'E'): 4, ('Q', 'D'): 0, ('X', 'B'): -1, ('Z', 'W'): -3, ('S', 'K'): 0, ('X', 'K'): -1, ('V', 'R'): -3, ('K', 'E'): 1, ('I', 'A'): -1, ('P', 'H'): -2, ('B', 'W'): -4, ('K', 'K'): 5, ('H', 'C'): -3, ('E', 'N'): 0, ('Y', 'Q'): -1, ('H', 'H'): 8, ('B', 'I'): -3, ('C', 'A'): 0, ('I', 'I'): 4, ('V', 'A'): 0, ('W', 'I'): -3, ('T', 'F'): -2, ('V', 'S'): -2, ('T', 'T'): 5, ('F', 'M'): 0, ('L', 'E'): -3, ('M', 'M'): 5, ('Z', 'G'): -2, ('D', 'R'): -2, ('M', 'D'): -3, ('W', 'H'): -2, ('G', 'C'): -3, ('S', 'R'): -1, ('S', 'I'): -2, ('P', 'Q'): -1, ('Y', 'A'): -2, ('X', 'I'): -1, ('E', 'A'): -1, ('B', 'Y'): -3, ('K', 'I'): -3, ('H', 'A'): -2, ('P', 'G'): -2, ('F', 'N'): -3, ('H', 'N'): 1, ('B', 'K'): 0, ('V', 'C'): -1, ('T', 'L'): -1, ('P', 'K'): -1, ('W', 'S'): -3, ('T', 'D'): -1, ('T', 'M'): -1, ('P', 'N'): -2, ('K', 'H'): -1, ('T', 'R'): -1, ('Y', 'R'): -2, ('L', 'C'): -1, ('B', 'L'): -4, ('Z', 'Y'): -2, ('W', 'N'): -4, ('G', 'A'): 0, ('S', 'P'): -1, ('E', 'Q'): 2, ('C', 'N'): -3, ('H', 'Q'): 0, ('D', 'N'): 1, ('Y', 'C'): -2, ('L', 'H'): -3, ('E', 'C'): -4, ('Z', 'H'): 0, ('H', 'G'): -2, ('P', 'E'): -1, ('Y', 'S'): -2, ('G', 'R'): -2, ('B', 'M'): -3, ('Z', 'Z'): 4, ('W', 'M'): -1, ('Y', 'T'): -2, ('Y', 'P'): -3, ('Y', 'Y'): 7, ('T', 'K'): -1, ('Z', 'I'): -3, ('T', 'P'): -1, ('V', 'L'): 1, ('F', 'I'): 0, ('G', 'Q'): -2, ('L', 'A'): -1, ('M', 'I'): 1}

# General function to give annotations for regions that have direct mappings onto the hmm alignment (imgt states)

def _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions):
    """
    General function to number a sequence and divide it into different regions  
    
    @param sequence: The sequence string
    @param state_vector: The list of states from the aligned hmm
    @param state_string: A string of states for the scheme relative to IMGT (this is X for a direct equivalence, I if needs to be treated as insertion)
    @param region_string: A string of characters that indicate which hmm states are in each regions for this scheme (i.e. how should the sequence be divided up)
    @param region_index_dict: A dictionary converting the characters in region string to an index of the regions. 
    @param rels: The difference of the numbering integer at the *start* of each region
    @param n_regions: The number of regions
    @param exclude_deletions: A list of region indices for which deletion states should not be included. Typically the CDRs.
    
    @return: A list of lists where each region has been numbered according to the scheme. Some regions will need renumbering. This should be taken care of after the function called.
    
    """

    
    _regions = [ [] for _ in xrange(n_regions) ]
    
    # Initialise the insertion index (-1 is a blank space) and the previous state.
    insertion = -1
    previous_state_id = 1 
    start_index, end_index  = None, None
    
    region = None
    # Iterate over the aligned state vector
    for (state_id, state_type ), si in state_vector:
       
        # Retrieve the region index
        if state_type != "i" or region is None: # BUG_FIX - JD 9/4/15 - do not allow a new region to start as an insertion.
            region = region_index_dict[region_string[state_id-1]] 
        
        # Check the state_types
        if state_type == "m": # It is a match
            
            # Check whether this position is in the scheme as an independent state
            if state_string[state_id-1]=="I": # No, it should be treated as an insertion 
                insertion +=1 # Increment the insertion annotation index
                rels[region] -= 1 # Update the relative numbering from the imgt states
            else: # Yes 
                insertion = -1 # Reset the insertions 
            
            # Add the numbering annotation to the appropriate region list            
            _regions[region].append( ( (state_id + rels[region], alphabet[insertion] ), sequence[si]  ) )
            previous_state_id = state_id # Record the previous state ID
            if start_index is None:
                start_index = si
            end_index = si
            
        elif state_type == "i": # It is an insertion
            insertion +=1 # Increment the insertion annotation index
            
            # Add the numbering annotation to the appropriate region list
            _regions[region].append( ( (previous_state_id + rels[region], alphabet[insertion]), sequence[si]  ) )
            if start_index is None:
                start_index = si
            end_index = si

        else: # It is a deletion
            
            # Check whether this position is in the scheme as an independent state
            if state_string[state_id-1]=="I": # No, therefore irrelevant to the scheme.
                rels[region] -= 1 # Update the relative numbering from the imgt states
                continue 
            
            # Check whether we should add the delete state to the numbering (don't if it needs ammending)  
            if region not in exclude_deletions: # Don't add deletions
                _regions[region].append( ( (state_id + rels[region], " "), "-"  ) )
            
            insertion = -1 # Reset the insertions
            previous_state_id = state_id # Record the previous state ID, should not be needed (no delete to insert state transition)
        
        assert insertion < 25, "Too many insertions for numbering scheme to handle" # We ran out of letters. Cows will not be happy :( 
            
    return _regions, start_index, end_index


# Functions to perform the numbering and the corrections for each of the implemented schemes.
# These have been written fairly verbosely so that the template of how to generate a function for a new scheme is more clear.
# They have two stages: Perform the mapping between imgt and the scheme; Renumber those regions that do not map nicely onto imgt (e.g. CDR insertions)
    
#######
# Aho #
#######

def number_aho(state_vector, sequence, chain_type):
    """    
    Apply the Aho numbering scheme
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in IMGT scheme, I is an insertion. (All X's for IMGT)

    XXXXXXX XXX XXXXXXXXXXXXXX XXXXXXXXXXXXXXXX XXXXXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXX XXXXX XXXX XXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXX
    AAAAAAA BBB CCCCCCCCCCCCCC DDDDDDDDDDDDDDDD EEEEEEEEEEEEEEE FFFFFFFFFFFFFFFFFFFF GGGGG HHHH IIIIIIIIIIIIIIIIIIII JJJJJJJJJJJJJ KKKKKKKKKKK


    Regions - (N.B These do not match up with any particular definition of CDR)

    A. 1-7 inclusive 
    B. 8-10 inclusive. Deletion occurs at 8
    C. 11-24 inclusive.
    D. 25-42 inclusive (deletion surround 28) 32-42 inclusive (deletions surround 36)
    E. 43-57 inclusive 
    F. 58-77 inclusive (deletions surround 63). Alpha chains have deletions at 74,75
    G. 78-82 inclusive 
    H. 83-86 inclusive  gaps on 86 then 85
    I. 87-106 inclusive
    J. 107-137 inclusive gaps on 123 symetrically.
    K. 
    
    """
    
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string =  'AAAAAAABBBCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFGGGGGHHHHIIIIIIIIIIIIIIIIIIIIJJJJJJJJJJJJJKKKKKKKKKKK'
#                     0      1  2             3               4              5                   6    7   8                   9            10
    region_index_dict = dict( zip( "ABCDEFGHIJK", range(11) ) )
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:0,
                         2:0,
                         3:0,
                         4:2,
                         5:2,
                         6:2,
                         7:2,
                         8:2,
                         9:2,
                         10:21}
    
    n_regions = 11
    
    exclude_deletions = [1,3,4,5,7,9]    
    
    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    ###############
    # Renumbering #
    ###############

    _numbering = [ _regions[0], _regions[1], _regions[2],[], _regions[4], [], _regions[6], [], _regions[8],_regions[9],_regions[10] ]

    #########################################
    # Move the deletion around 8-10 onto 8  #
    #########################################
    length=len( _regions[1] )
    if length < 3:
        annotations = [(8," "), (9," "), (10, " ")][3-length:]
        _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in xrange(length) ]
    else:
        _numbering[1] = _regions[1]


    #########
    # CDR 1 # - divided in two parts in the Aho scheme.
    ######### - gaps at 28 depending on the chain type (although not absolutely defined e.g. 8fab which this will not get correct).
    # The default number of gaps expected at 28 for each chain type is:
    pregaps = {"H":1,"K":2,"L":1,"A":1,"B":2}[chain_type]
    # This is a heuristic taken from the Aho paper.

    length = len( _regions[3] )
    # Only get rid of the default pregaps if there are insertions in the main loop.
    insertions = max( length-16, 0)
    for ins in xrange(max( length-16, 0)):
        if pregaps > 0:
            pregaps -= 1
        else:
            break
    # Pre-31. 
    positions  = [(25," "), (26," "), (27," "), (28," "), (29," "), (30," "), (31," ")]
    if length < 7-pregaps:
        _n = length
    else:
        _n = 7-pregaps
    annotations = [ positions[i] for i in sorted([0,1,6,5,4,2,3][:_n]) ]

    # Post 31. Delete in order 36,35,37,34,38,33,39,32,40,41,42. Insertions on 36 if present
    _n = max( length - _n, 0 )
    positions  = [(32," "),(33," "),(34," "),(35," "),(36," "),(37," "),(38," "),(39," "),(40," "),(41," "),(42," ")]
    if insertions:
        annotations += positions[:5] + [(36, alphabet[i]) for i in xrange(insertions)] + positions[5:]
    else:
        annotations += [ positions[i] for i in sorted([10,9,8,0,7,1,6,2,5,3,4][:_n]) ]
    _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in xrange(length) ]

    #########
    # CDR 2 #
    #########
    # Ignore all gaps in this region and number to the aho specification.
    # Gaps centred on 63 and also 74,75
    # Gaps only occur on 74,75 in the alpha chain according to the aho analysis.
    # This region describe 58 to 77 inclusive
        
    length = len(_regions[5])
    insertions = max( length - 20, 0 )
    positions = [ (_, " ") for _ in xrange(58,78) ]
    if insertions:
        annotations = positions[:6] + [(63, alphabet[i]) for i in xrange(insertions)] + positions[6:]          
    else:
        deletions = 20 -length
        if chain_type == "A":
            removal_order = [17,16,5,4,6,3,7,2,8,1,9,0]+range(10,16)+[18,19]
        else:
            removal_order = [5,4,6,3,7,2,8,1,9,0]+range(10,20)
        annotations = [ positions[i] for i in sorted( removal_order[deletions:] ) ]            
    _numbering[5] = [ (annotations[i], _regions[5][i][1]) for i in xrange(length) ]


    ###############################################
    # Move deletions around 83-86 onto 85 then 86 #
    ###############################################
    length=len( _regions[7] )
    if length < 4:
        annotations = [(83," "), (84," "), (85, " "), (86, " ")][:length]
        _numbering[7] = [ (annotations[i], _regions[7][i][1]) for i in xrange(length) ]
    else:
        _numbering[7] = _regions[7]

    #########
    # CDR 3 #
    #########
    # Deletions on 123. 
    # Point of the Aho scheme is that they have accounted for all possible positions.
    # Assumption is that no more insertions will occur.... 
    # We'll put insertions on 123 linearly.(i.e.ABCDEF...) if they ever do.
    # This does not happen as the HMM alignment does not recognise CDR3s that are that long.
    length = len(_regions[9])
    positions = [ (_, " ") for _ in xrange(107,139) ]
    length = len(_regions[9])
    insertions = max( length-32, 0)
    if insertions:
        annotations = positions[:17] + [(123, alphabet[i]) for i in xrange(insertions)] + positions[17:]
    else:
        deletions = 32 - length
        removal_order = [16,17,15,18,14,19,13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,0]
        annotations = [ positions[i] for i in sorted( removal_order[deletions:] ) ]    
    _numbering[9] = [ (annotations[i], _regions[9][i][1]) for i in xrange(length) ]

    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] + \
           _numbering[3] + \
           _numbering[4] + \
           _numbering[5] + \
           _numbering[6] + \
           _numbering[7] + \
           _numbering[8] + \
           _numbering[9] + \
           _numbering[10] , startindex, endindex
    
########
# IMGT #
########

def number_imgt(state_vector, sequence):
    """    
    Apply the IMGT numbering scheme for heavy or light chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in IMGT scheme, I is an insertion. (All X's for IMGT)
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXX
    11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 2222222222222 33333333333
    
    Regions - (N.B These do not match up with any particular definition of CDR)
    1. All positions before CDR3
    2. CDR positions 105 (inc) to 118 (exc)
    3. All positions after CDR3
    
    Region 2 is renumbered     
    
    """
    
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111222222222222233333333333'

    region_index_dict = {"1":0,"2":1,"3":2}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:0,
                         2:0}
    
    n_regions = 3
    
    exclude_deletions = [1]    
    
    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    ###############
    # Renumbering #
    ###############

    _numbering = [ _regions[0], [], _regions[2]  ]
    
    # CDR3 requires renumbering
    # CDR3 is numbered in a specific way for each CDRlength - i.e. we gap from the sides. Insertions are at 111 and 112 symetrically.    
    length = len( _regions[1] )
    si = 0
    
    previous_state_id = 104
    for annotation in get_cdr3_annotations(length, scheme="imgt"):
        if annotation is None:
            _numbering[1].append( ((previous_state_id+1, " "), "-"   ) )
            previous_state_id+=1
        else:
            _numbering[1].append( (annotation, _regions[1][si][1] ) )
            previous_state_id = annotation[0]
            si+=1
  
        # Return the full vector and the start and end indices of the numbered region of the sequence
    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] , startindex, endindex
  
###########
# Chothia #
###########

# Heavy chains
def number_chothia_heavy(state_vector, sequence):
    """
    Apply the Chothia numbering scheme for heavy chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Chothia scheme, I is an insertion.
    XXXXXXXXXI XXXXXXXXXXXXXXXXXXXX IIIIXX XXXXXXXXXXXXXXXXXXXX XIXII XXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXX XXXXXXIII XXXXXXXXXXXXX
    1111111111 22222222222222222222 333333 44444444444444444444 55555 666666666666666666666666666666666666666666666 777777777 8888888888888
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Put the insertions at Chothia position 6
     2  -  Simple mapping (treat "I" states as inserts and not own match states)
     3  -  CDRH1 - 30 (inc) to 32 (exc) put insertions on 31
     4  -  Simple mapping (treat "I" states as inserts and not own match states)
     5  -  CDRH2 - 52 (inc) 58 (exc) put insertions on 52 
     6  -  Simple mapping (treat "I" states as inserts and not own match states)
     7  -  CDRH3 93 (inc) to 103 (exc) put insertion on 100
     8  -  Simple mapping (treat "I" states as inserts and not own match states)


    Regions 1,3,5 and 7 are renumbered
    
    """
 
    # Set up the numbering 
 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111112222222222222333333333333344444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'
    
    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [0,2,4,6] # Don't put deletions in these regions

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    
    ###############
    # Renumbering #
    ###############

    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]
    
    # Chothia H region 1 (index 0)
    # Insertions are placed at Chothia position 6.
    # Count how many we recognised as insertion by the hmm
    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 
    # We will place all insertion in this region at Chothia position 6.
    if insertions:
        start = _regions[0][0][0][0] # The starting Chothia number as found by the HMM (could easily start from 2 for example)
        # I have a feeling this may be a source of a bug in very unusual cases. Can't break for now. Will catch mistakes in a validate function. 
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in xrange(start, 7) ] + [ (6, alphabet[_]) for _ in xrange(insertions) ] + [(7," "),(8," "),(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in xrange(length) ]
    else:
        _numbering[0] = _regions[0]

    
    # CDR1 
    # Chothia H region 3 (index 2)
    # put insertions onto 31
    length = len( _regions[2] )
    insertions = max(length - 9, 0) # Pulled back to the cysteine as heavily engineered cdr1's are not playing nicely
    annotations = [(_, " ") for _ in xrange(23,32)][:length-insertions] + [(31, alphabet[i]) for i in xrange(insertions) ]
    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in xrange(length) ]
 
    # CDR2
    # Chothia H region 5 (index 4) 
    # put insertions onto 52
    length = len( _regions[4] )
    # 50 to 57 inclusive
    insertions = max(length - 8, 0) # Eight positions can be accounted for, the remainder are insertions
    # Delete in the order, 52, 51, 50,53, 54 ,55, 56, 57
    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in xrange(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in xrange(length) ]
     
    # CDR3
    # Chothia H region 7 (index 6) 
    # put insertions onto 100
    length = len( _regions[6] )    
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="heavy")
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in xrange(length)  ]

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] + \
           _numbering[3] + \
           _numbering[4] + \
           _numbering[5] + \
           _numbering[6] + \
           _numbering[7]   , startindex, endindex                                        

# Light chains
def number_chothia_light(state_vector, sequence):
    """
    Apply the Chothia numbering scheme for light chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Chothia scheme, I is an insertion.
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX IIIIIIX XXXXXXXXXXXXXXXXXXXX XIIIIIIIXXX XXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXX XXXXXIIIIXX XXXXXXXXXXXXX
    11111111111111111111111111111 2222222 33333333333333333333 44444444444 5555555555555555555555555555555555555 66666666666 7777777777777
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRL1 - 24 (inc) to 35 (exc) put insertions on 30
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRL2 - 51 (inc) 55 (exc) put insertions on 54 
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  CDRL3 89 (inc) to 96 (exc) put insertion on 95
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
    
    Region 2, 3 and 5 are renumbered
    
    """
    
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666667777777777777'
    

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1: 0,
                         2:-6,
                         3:-6,
                         4:-13,
                         5:-16,
                         6:-20,
                         }    
    
    n_regions = 7
    
    exclude_deletions = [1,3,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    _numbering = [ _regions[0], [], _regions[2], [], _regions[4], [], _regions[6] ]
    
    
    ###############
    # Renumbering #
    ###############

    # CDR1 
    # Chothia L region 2 (index 1)
    # put insertions onto 30
    length = len( _regions[1] )
    insertions = max(length - 11, 0) # Eleven positions can be accounted for, the remainder are insertions
    # Delete forward from 31 
    annotations  =  [(24, " "),(25, " "), (26, " "), (27, " "), (28, " "),(29, " "),(30, " ")][:max(0,length)] 
    annotations += [(30, alphabet[i]) for i in xrange(insertions) ]
    annotations += [(31, " "),(32, " "),(33, " "),(34, " ")][ abs( min(0,length-11) ):] 
    _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in xrange(length) ]


    # CDR2
    # Chothia L region 4 (index 3) 
    # put insertions onto 54
    length = len( _regions[3] )
    insertions = length - 4
    annotations = [(51, " "),(52, " "),(53, " "),(54, " ")][:length-insertions] + [(54, alphabet[i]) for i in xrange(insertions) ]
    _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in xrange(length) ]
    
    # CDR3
    # Chothia L region 6 (index 5) 
    # put insertions onto 95
    length = len( _regions[5] )    
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="light")
    _numbering[5]  = [ (annotations[i], _regions[5][i][1]) for i in xrange(length)  ]
    

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] + \
           _numbering[3] + \
           _numbering[4] + \
           _numbering[5] + \
           _numbering[6]   , startindex, endindex      


#########
# Kabat #
#########

# Heavy chains
def number_kabat_heavy(state_vector, sequence):
    """
    Apply the Kabat numbering scheme for heavy chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Kabat scheme, I is an insertion.
    XXXXXXXXXI XXXXXXXXXXXXXXXXXXXX IIIIXXXXXX XXXXXXXXXXXXXXXX XIXII XXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXX XXXXXXIII XXXXXXXXXXXXX
    1111111111 22222222222222222222 3333333333 4444444444444444 55555 666666666666666666666666666666666666666666666 777777777 8888888888888
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Put the insertions at Chothia position 6
     2  -  Simple mapping (treat "I" states as inserts and not own match states)
     3  -  CDRH1 - 30 (inc) to 36 (exc) put insertions on 35
     4  -  Simple mapping (treat "I" states as inserts and not own match states)
     5  -  CDRH2 - 52 (inc) 58 (exc) put insertions on 52 
     6  -  Simple mapping (treat "I" states as inserts and not own match states)
     7  -  CDRH3 93 (inc) to 103 (exc) put insertion on 100
     8  -  Simple mapping (treat "I" states as inserts and not own match states)

    """
 
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [2,4,6]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    

    ###############
    # Renumbering #
    ###############
        
    # Renumbering required for 0, 2, 4, 6 regions in Chothia heavy
    
    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]


    # Kabat H region 1 (index 0)
    # Insertions are placed at Kabat position 6.
    # Count how many we recognised as insertion by the hmm
    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 
    # We will place all insertion in this region at Kabat position 6.
    if insertions:
        start = _regions[0][0][0][0] # The starting Kabat number as found by the HMM (could easily start from 2 for example)
        # I have a feeling this may be a source of a bug in very unusual cases. Can't break for now. Will catch mistakes in a validate function. 
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in xrange(start, 7) ] + [ (6, alphabet[_]) for _ in xrange(insertions) ] + [(7," "),(8," "),(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in xrange(length) ]
    else:
        _numbering[0] = _regions[0]
    
    
    # CDR1 
    # Kabat H region 3 (index 2)
    # Put insertions onto 35. Delete from 35 backwards
    length = len( _regions[2] )
    insertions = max(0,length - 13)
    annotations = [(_,' ') for _ in xrange(23, 36)][:length] # Pulled back to cysteine as some heavily engineered cdrs do not play nicely
    annotations += [(35, alphabet[i]) for i in xrange(insertions) ]
    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in xrange(length) ]
 
    # CDR2
    # Chothia H region 5 (index 4) 
    # put insertions onto 52
    length = len( _regions[4] )
    # 50 to 57 inclusive
    insertions = max(length - 8, 0) # Eight positions can be accounted for, the remainder are insertions
    # Delete in the order, 52, 51, 50,53, 54 ,55, 56, 57
    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in xrange(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in xrange(length) ]

     
    # CDR3
    # Chothia H region 7 (index 6) 
    # put insertions onto 100
    length = len( _regions[6] )    
    annotations = get_cdr3_annotations(length, scheme="kabat", chain_type="heavy") #  Chothia and Kabat the same here
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in xrange(length)  ]

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] + \
           _numbering[3] + \
           _numbering[4] + \
           _numbering[5] + \
           _numbering[6] + \
           _numbering[7]   , startindex, endindex                 
           
# Light chains    
def number_kabat_light(state_vector, sequence):
    """
    Apply the Kabat numbering scheme for light chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Kabat scheme, I is an insertion.
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX IIIIIIX XXXXXXXXXXXXXXXXXXXX XIIIIIIIXXX XXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXX XXXXXIIIIXX XXXXXXXXXXXXX
    11111111111111111111111111111 2222222 33333333333333333333 44444444444 5555555555555555555555555555555555555 66666666666 7777777777777
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRL1 - 24 (inc) to 35 (exc) put insertions on 27
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRL2 - 51 (inc) 55 (exc) put insertions on 54 
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  CDRL3 89 (inc) to 96 (exc) put insertion on 95
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
    
    """
    
    # Set up the numbering 
 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666667777777777777'
    
    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1: 0,
                         2:-6,
                         3:-6,
                         4:-13,
                         5:-16,
                         6:-20,
                         }    
    
    n_regions = 7
    
    exclude_deletions = [1,3,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    _numbering = [ _regions[0], [], _regions[2], [], _regions[4], [], _regions[6] ]
    
    
    ###############
    # Renumbering #
    ###############
    
    # CDR1 
    # Kabat L region 2 (index 1)
    # put insertions onto 27
    length = len( _regions[1] )
    insertions = max(length - 11, 0) # Eleven positions can be accounted for, the remainder are insertions
    # Delete forward from 28 
    annotations  =  [(24, " "),(25, " "), (26, " "), (27, " ")][:max(0,length)] 
    annotations += [(27, alphabet[i]) for i in xrange(insertions) ]
    annotations += [(28, " "),(29, " "),(30, " "),(31, " "),(32, " "),(33, " "),(34, " ")][ abs( min(0,length-11) ):] 
    _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in xrange(length) ]
  

    # CDR2
    # Kabat L region 4 (index 3) 
    # put insertions onto 54
    length = len( _regions[3] )
    insertions = max(length - 4, 0) # Four position can be accounted for, the remainder are insertions.
    annotations = [(51, " "),(52, " "),(53, " "),(54, " ")][:length-insertions] + [(54, alphabet[i]) for i in xrange(insertions) ]
    _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in xrange(length) ]
    
    # CDR3
    # Kabat L region 6 (index 5) 
    # put insertions onto 95
    length = len( _regions[5] )    
    insertions = max(length - 4, 0) # Four position can be accounted for, the remainder are insertions.
    annotations = get_cdr3_annotations(length, scheme="kabat", chain_type="light") 
    _numbering[5]  = [ (annotations[i], _regions[5][i][1]) for i in xrange(length)  ]
    

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] + \
           _numbering[3] + \
           _numbering[4] + \
           _numbering[5] + \
           _numbering[6]   , startindex, endindex      



#############################
# Martin (extended Chothia) #
#############################

# Heavy chains
def number_martin_heavy(state_vector, sequence):
    """
    Apply the Martin (extended Chothia) numbering scheme for heavy chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Martin scheme, I is an insertion.
    XXXXXXXXXI XXXXXXXXXXXXXXXXXXXX IIIIXX XXXXXXXXXXXXXXXXXXXX XIXII XXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXX XXXXXXIII XXXXXXXXXXXXX
    1111111111 22222222222222222222 333333 44444444444444444444 55555 666666666666666666666666666666666666666666666 777777777 8888888888888
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Put the insertions at Chothia position 8
     2  -  Simple mapping (treat "I" states as inserts and not own match states)
     3  -  CDRH1 - 30 (inc) to 32 (exc) put insertions on 31
     4  -  Simple mapping (treat "I" states as inserts and not own match states)
     5  -  CDRH2 - 52 (inc) 58 (exc) put insertions on 52 
     6  -  Simple mapping (treat "I" states as inserts and not own match states)
     7  -  CDRH3 93 (inc) to 103 (exc) put insertion on 100
     8  -  Simple mapping (treat "I" states as inserts and not own match states)


    Regions 1,3,5 and 7 are renumbered
    
    """
 
    # Set up the numbering 
 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111112222222222222333333333333344444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'
    
    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [2,4,6]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    

    ###############
    # Renumbering #
    ###############
        
    # Renumbering required for 0, 2, 4, 6 regions in Chothia heavy
    
    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]
    
    # Chothia H region 1 (index 0)
    # Insertions are placed at Chothia position 8.
    # Count how many we recognised as insertion by the hmm
    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 
    # We will place all insertion in this region at Chothia position 8.
    if insertions:
        start = _regions[0][0][0][0] # The starting Chothia number as found by the HMM (could easily start from 2 for example)
        # I have a feeling this may be a source of a bug in very unusual cases. Can't break for now. Will catch mistakes in a validate function. 
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in xrange(start, 9) ] + [ (8, alphabet[_]) for _ in xrange(insertions) ] + [(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in xrange(length) ]
    else:
        _numbering[0] = _regions[0]

    
    
    # CDR1 
    # Chothia H region 3 (index 2)
    # put insertions onto 31
    length = len( _regions[2] )
    insertions = max(length - 9, 0) # Pulled back to the cysteine as heavily engineered cdr1's are not playing nicely
    annotations = [(_, " ") for _ in xrange(23,32)][:length-insertions] + [(31, alphabet[i]) for i in xrange(insertions) ]
    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in xrange(length) ]
 
    # CDR2
    # Chothia H region 5 (index 4) 
    # put insertions onto 52
    length = len( _regions[4] )
    # 50 to 57 inclusive
    insertions = max(length - 8, 0) # Eight positions can be accounted for, the remainder are insertions
    # Delete in the order, 52, 51, 50,53, 54 ,55, 56, 57
    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in xrange(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in xrange(length) ]

     
    # CDR3
    # Chothia H region 7 (index 6) 
    # put insertions onto 100
    length = len( _regions[6] )    
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="heavy")
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in xrange(length)  ]

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] + \
           _numbering[3] + \
           _numbering[4] + \
           _numbering[5] + \
           _numbering[6] + \
           _numbering[7]   , startindex, endindex                                        

# Light chains
def number_martin_light(state_vector, sequence):
    """
    Apply the Martin numbering scheme for light chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Martin scheme, I is an insertion.
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX IIIIIIX XXXXXXXXXXXXXXXXXXXX XIIIIIIIXXX XXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXX XXXXXIIIIXX XXXXXXXXXXXXX
    11111111111111111111111111111 2222222 33333333333333333333 44444444444 5555555555555555555555555555555555555 66666666666 7777777777777
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRL1 - 30 (inc) to 31 (exc) put insertions on 30
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRL2 - 51 (inc) 55 (exc) put insertions on 54 
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  CDRL3 89 (inc) to 96 (exc) put insertion on 95
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
    
    Region 2, 3 and 5 are renumbered
    
    """
    
    # In terms of there seems very little difference (if any) between the martin scheme and chothia scheme.
    # For example the insertion at L52 in Martin already exists in the Chothia implementation (although the scientific description
    # is only discussed in the abhinandan & martin paper 2008. 
    
    # The only other difference we do not explicitely take care of is possible insertions at 68. However the alignment to the HMM
    # should take care of this. (insertions are allowed anywhere in ANARCI to allow flexibility) 
    return number_chothia_light(state_vector,sequence)


###########
# Wolfguy #
###########

def number_wolfguy_heavy(state_vector, sequence):
    """
    Apply the wolfguy numbering scheme for heavy chains 

    The scheme numbers the sequence using different segments so that the numbering tells you
    where in the antibody the sequence is describing. 

    XXXXXXXXXIXXXXXXXXXXXXXXXX XXXXXXXXXXXXXX XXXXXXXXXXXXXX XXXXXXXXXXXXXXXXXXIX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXXX XXXXXXXXXXX
    11111111111111111111111111 22222222222222 33333333333333 44444444444444444444 55555555555555555555555555555555 66666666666 77777777777
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRH1 - 155-199 (inc). Gap symmetrically about 175-176.
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRH2 - 251-299 (inc). Gap symmetrically about 271-272, then gap back from 294.
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  CDRH3 351-399 (inc). Gap according to the  
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
    
     Start gaps on rhs each time.
    """
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111111111111111111111222222222222223333333333333344444444444444444444555555555555555555555555555555556666666666677777777777'
    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:100, 
                         1:124,
                         2:160,
                         3:196,
                         4:226,
                         5:244,
                         6:283}    
    
    n_regions = 7
    
    exclude_deletions = [1,3,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    ###############
    # Renumbering #
    ###############
        
    # Renumbering required for 1, 3, 5 regions in wolfguy heavy    
    _numbering = [ _regions[0], [] , _regions[2], [], _regions[4] , [], _regions[6] ]


    # CDRH1
    # Delete symmetrically about 177. Delete right first.
    # May have to change this to reflect where the point of symmetry is
    ordered_deletions = [151]
    for p1,p2 in zip( range(152,176), range(199, 175,-1)): ordered_deletions += [ p1,p2 ]
    length = len( _regions[1] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[1]  = [ ((annotations[i]," "), _regions[1][i][1]) for i in xrange(length)  ]
    
    # CDRH2
    # Delete symmetrically about 271. Delete right first.
    # Then delete right from 288
    ordered_deletions = [251]
    for p1,p2 in zip( range(252,271), range(290, 271,-1)): ordered_deletions += [ p1,p2 ]
    ordered_deletions.append( 271 )
    ordered_deletions = range( 299, 290, -1 ) + ordered_deletions
    length = len( _regions[3] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[3]  = [ ((annotations[i]," "), _regions[3][i][1]) for i in xrange(length)  ]

    # CDRH3	
    # Delete symmetrically about 374. Delete right first. 
    # Scheme changes at length 8
    # Scheme changes at length 12
    ordered_deletions = []
    for p1,p2 in zip( range(356,374), range(391, 373,-1)): ordered_deletions += [ p1,p2 ]
    ordered_deletions = [ 354, 394, 355, 393, 392 ] + ordered_deletions
    ordered_deletions = [ 399, 398, 351, 352, 397, 353, 396, 395 ] + ordered_deletions
    length = len( _regions[5] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[5]  = [ ((annotations[i]," "), _regions[5][i][1]) for i in xrange(length)  ]    
  
    # Return the full vector and the start and end indices of the numbered region of the sequence
    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] + \
           _numbering[3] + \
           _numbering[4] + \
           _numbering[5] + \
           _numbering[6] , startindex, endindex   
            

def number_wolfguy_light(state_vector, sequence):
    """
    Apply the wolfguy numbering scheme for light chains 

    The scheme numbers the sequence using different segments so that the numbering tells you
    where in the antibody the sequence is describing. 

    XXXXXXX XXX XXXXXXXXXXXXX XXXXXXXXXXXXXXXXX XXXXXXXXXXXXXXX XXXXXXXXXXXXXX XXXIXXXXXXX XXXX XXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXX
    1111111 AAA BBBBBBBBBBBBB 22222222222222222 333333333333333 44444444444444 55555555555 6666 77777777777777777777 8888888888888 99999999999

    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     A  -  Move indels onto 508
     B  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRL1 - 551-599 (inc). Assign via the matching consensus sequence and length.
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRL2 - 651-699 (inc). Gap about 673 then right from 694
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  Move indels onto 713 and 714
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
     8  -  CDRL3 751-799 (inc). Gap symmetrically about 374-375 
     9  -  Simple mapping (treat "I" states as inserts and not own match states)
    
    """
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '1111111AAABBBBBBBBBBBBB222222222222222223333333333333334444444444444455555555555666677777777777777777777888888888888899999999999'
    region_index_dict = {"1":0,"A":1,"B":2,"2":3,"3":4,"4":5,"5":6,"6":7,"7":8,"8":9,"9":10}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:500,
                         1:500,
                         2:500,    
                         3:527,
                         4:560,
                         5:595,
                         6:631,
                         7:630,
                         8:630,                                                  
                         9:646,
                         10:683}    
    
    n_regions = 11
    
    exclude_deletions = [1,3,5,7,9]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    ###############
    # Renumbering #
    ###############
        
    # Renumbering required for 1, 3, 5 regions in wolfguy heavy    
    _numbering = [ _regions[0], [], _regions[2], [] , _regions[4], [], _regions[6], [], _regions[8], [], _regions[10] ]


    # Gaps in the first section go 508 instead of the imgt 510 equivalent
    length = len(_regions[1] )
    annotations = sorted([ (510,' '), (509, ' '), (508, ' ')][ :length ] + [(508,a) for a in alphabet[:max(0, length-3)]])  
    _numbering[1]  = [ (annotations[i], _regions[1][i][1]) for i in xrange(length)  ]
    
    # CDRL1
    # Number by predicting the canonical 
    length = len(_regions[3] )
    annotations = _get_wolfguy_L1( _regions[3], length)
    _numbering[3]  = [ ((annotations[i]," "), _regions[3][i][1]) for i in xrange(length)  ]
    
    # CDRL2
    # Delete about 673. Finally delete right from 694
    ordered_deletions = []
    for p1,p2 in zip( range(651,673), range(694, 672,-1)): ordered_deletions += [ p2,p1 ]
    ordered_deletions = range( 699, 694, -1 ) + ordered_deletions
    length = len( _regions[5] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[5]  = [ ((annotations[i]," "), _regions[5][i][1]) for i in xrange(length)  ]

    # The placement of the indel in wolfguy is different to that in imgt
    length = len( _regions[7] )
    insertions = max( 0, length - 4 )
    annotations = [(711, ' '), (712, ' '), (713, ' '), (714, ' ')][:length] + [ (714, a) for a in alphabet[:insertions] ]    
    _numbering[7]  = [ (annotations[i], _regions[7][i][1]) for i in xrange(length)  ]  
    
    # CDRL3
    # Delete symmetrically about 775. Delete right first. Finally delete 798 and 799
    ordered_deletions = []
    for p1,p2 in zip( range(751,775), range(799, 775,-1)): ordered_deletions += [ p1,p2 ]
    ordered_deletions.append( 775 )
  
    length = len( _regions[9] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[9]  = [ ((annotations[i]," "), _regions[9][i][1]) for i in xrange(length)  ]  
  
    # Return the full vector and the start and end indices of the numbered region of the sequence
    return _numbering[0] + \
           _numbering[1] + \
           _numbering[2] + \
           _numbering[3] + \
           _numbering[4] + \
           _numbering[5] + \
           _numbering[6] + \
           _numbering[7] + \
           _numbering[8] + \
           _numbering[9] + \
           _numbering[10] , startindex, endindex  


def _get_wolfguy_L1(seq, length):
    """
    Wolfguy's L1 annotation is based on recognising the length and the sequence pattern defined
    by a set of rules. If the length has not been characterised, we number symmetrically about the
    middle of the loop.
    """
    
    # These are the annotations for different lengths of L1 according to the wolfguy definitions.
    L1_sequences = {
    9: [['9',     'XXXXXXXXX', [551, 552, 554, 556, 563, 572, 597, 598, 599]]], 
    10: [['10',   'XXXXXXXXXX', [551, 552, 553, 556, 561, 562, 571, 597, 598, 599]]], 
    11: [['11a',  'RASQDISSYLA', [551, 552, 553, 556, 561, 562, 571, 596, 597, 598, 599]], 
         ['11b',  'GGNNIGSKSVH', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]], 
         ['11b.2','SGDQLPKKYAY', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]]], 
    12: [['12a',  'TLSSQHSTYTIE', [551, 552, 553, 554, 555, 556, 561, 563, 572, 597, 598, 599]], 
         ['12b',  'TASSSVSSSYLH', [551, 552, 553, 556, 561, 562, 571, 595, 596, 597, 598, 599]], 
         ['12c',  'RASQSVxNNYLA', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]], 
         ['12d',  'rSShSIrSrrVh', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]]], 
    13: [['13a',  'SGSSSNIGNNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 597, 598, 599]], 
         ['13b',  'TRSSGSLANYYVQ', [551, 552, 553, 554, 556, 561, 562, 563, 571, 572, 597, 598, 599]]], 
    14: [['14a',  'RSSTGAVTTSNYAN', [551, 552, 553, 554, 555, 561, 562, 563, 564, 571, 572, 597, 598, 599]], 
         ['14b',  'TGTSSDVGGYNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 596, 597, 598, 599]]], 
    15: [['15',   'XXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 594, 595, 596, 597, 598, 599]]], 
    16: [['16',   'XXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 594, 595, 596, 597, 598, 599]]], 
    17: [['17',   'XXXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 584, 594, 595, 596, 597, 598, 599]]]
    }    

    if length in L1_sequences: # Use the pre-defined motif 
        # Find the maximum scoring canonical form for this length. 
        curr_max = None, -10000
        for canonical in L1_sequences[length]:
            sub_score = 0
            for i in xrange( length ):
                try:
                    sub_score += blosum62[ (seq[i][1], canonical[1][i].upper() ) ]
                except KeyError:
                    sub_score += blosum62[ (canonical[1][i].upper(), seq[i][1] ) ]
            if sub_score > curr_max[1]:
                curr_max = canonical, sub_score

        # return the annotations
        return curr_max[0][2]
    else: # Use a symmetric numbering about the anchors.
        ordered_deletions = []
        for p1,p2 in zip( range(551,575), range(599, 575,-1)): ordered_deletions += [ p2,p1 ]
        ordered_deletions.append(575)
        return sorted( ordered_deletions[:length] )

######################
# Annotation of CDR3 #
######################
    
def get_cdr3_annotations(length, scheme="imgt", chain_type=""):
    """
    Given a length of a cdr3 give back a list of the annotations that should be applied to the sequence.
    
    This function should be depreciated
    """
    az = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" 
    za = "ZYXWVUTSRQPONMLKJIHGFEDCBA"
    
    if scheme=="imgt":
        start, end = 105, 118 # start (inclusive) end (exclusive)
        annotations = [None for _ in xrange(max(length,13))]
        front = 0
        back  = -1
        assert (length-13) < 50, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        for i in xrange(min(length,13)):
            if i%2:
                annotations[back] = (end+back, " ")
                back -= 1
            else:
                annotations[front] = (start+front, " ")
                front += 1
        for i in xrange(max(0,length-13)): # add insertions onto 111 and 112 in turn
            if i%2:
                annotations[back] = (112, za[back+6])
                back-=1
            else:
                annotations[front] = (111, az[front-7])
                front +=1        
        return annotations
    elif scheme in [ "chothia", "kabat"] and chain_type=="heavy": # For chothia and kabat
        # Number forwards from 93
        insertions = max(length - 10, 0)
        assert insertions < 25, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        ordered_deletions = [ (100, ' '), (99,' '), (98,' '), (97,' '), (96,' '), (95,' '), (101,' '),(102,' '),(94,' '), (93,' ') ]
        annotations = sorted( ordered_deletions[ max(0, 10-length): ] + [ (100,a) for a in az[:insertions ] ] )
        return annotations
    elif scheme in [ "chothia", "kabat"] and chain_type=="light":
        # Number forwards from 89
        insertions = length - 7
        assert insertions < 25, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        annotations = []
        start  = 89
        for _ in xrange(length - insertions):
            annotations.append( (start+_, " ") )
        for _ in xrange(insertions):
            annotations.append( (95, az[_]) )
        return annotations
    else:
        raise AssertionError("Unimplemented scheme")

