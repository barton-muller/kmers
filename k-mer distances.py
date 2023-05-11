#import modules

#for reading and manipulating sequence data
from Bio.Seq import Seq
from Bio import SeqIO

#for hash function
import mmh3

#for saving sketches
import json

#for finding intersections ignoring 14-mrs occur more than once
from collections import Counter

#for creating distance tree
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

#saves dictionary to file
def saveD(dict, name):

    with open(name, "w") as f:
        f.write(json.dumps(dict))

#reads dictionary saved to file and returns the dictionary object
def readD(name):
    # Read the dictionary from the text file
    with open(name, "r") as f:
        dict = json.loads(f.read())
    return dict

#returns the jaccard distanace between 2 input arrays of kmers
def Jaccard(sequence1, sequence2):

    #counter returns dictionary of unique items within input array with key of item and value no. of occurences
    # avoids double counting kmers which would make Jdistance between same kmer non 0
    counter1 = Counter(sequence1)
    counter2 = Counter(sequence2)
    
    #find intersection and union
    intersection = sum((counter1 & counter2).values())
    union = sum((counter1 | counter2).values())

    # return jaccard distance
    return 1 - intersection / union

#old method that didnt account for repeat kmers
"""def Jaccard(kmers1, kmers2):
    totaln = len(kmers1) + len(kmers2)
    
    overlap = len(set(kmers1).intersection(kmers2))
    union = totaln - overlap
    return 1-(overlap/union)"""

#returns sketch of array of kmers
def sketch(kmers):

    #iterate through kmers and hash
    hashed = []
    for i in kmers:
        hash_value = mmh3.hash(str(i))
        hashed.append(hash_value)

    hashed.sort(reverse = True) # sort low to high
    return hashed[0:1000] #return first 1000 hashes to create sketch

###########################   file loading and sketching    ###########################
"""
 this section loads all 4 files and creates 3 dictionaries representing the data
 file sequences containg all the sequences
 file kmers containg all the kmer arrays for each sequence
 sketches containg arrays of all the kmer array sketches
 all these dictionaries are in the same format with overarching dictionary having keys of the each file names and value of internal dictionary
 the internal dictionaries store the sequeces within the respective file
 the internal dictionaries have keys of sequence ID and values of the respective seuence/kmerarray/ sketch
 then saves the sketch dictioanry for retreval. 
 this proscess, particurly calculating sketches takes quite along time
 this proscess could all be done in one go in one dictionary 
"""
if input("load files (y?)") == "y": #can skip this long skep and just use saved sketch dictionary if already run
    
    #file names containg data were intrested in
    files = ["R6.fa", "TIGR4.fa", "14412_3#84.contigs_velvet.fa", "14412_3#82.contigs_velvet.fa"]

    #build dictionary with raw sequence data
    filesequences = {} #format: key=filename value= (dict of key=seqID value = sequence)

    #iterate though files
    for filename in files:

        seqdict = {} #create internal dictionary

        #read file and automaticaly iterate thorugh sequences within file as defined by filestructure
        for seq_record in SeqIO.parse(filename, "fasta"): #fasta is file type
            seqdict[seq_record.id] = seq_record.seq #add sequence to internal dict

        #add internal file dict to overarching dict
        filesequences[filename] = seqdict

    #iterating though the filesequce dictioanry create a simmilarly structed dict but convert sequences into kmer array
    filekmers = {}   #format [filename]:( [seqid]: array(kmers))
    #iterate though files in dict
    for filename in filesequences:

        kmerdict = {} #create new internal dict

        #iterate though sequences in each file
        for sequenceid in filesequences[filename]:

            sequence = filesequences[filename][sequenceid] #get sequece assosciated with id
            kmer_array = [] #create array to place kmers

            for k in range(len(sequence)-14): #run through array, -14 so dont go out of range
                kmer_array.append(sequence[0+k:14+k]) #splice array at 14 length intervals and add to array

            kmerdict[sequenceid] = kmer_array #add array to internal dict
        filekmers[filename] = kmerdict #ad internal kmer containg dict to overarching dict


    #same proscess but sketches kmerdicts
    filesketch = {}   #format [filename]:( [seqid]: sketch(array(kmers)))

    for filename in filekmers:
        sketchdict = {}
        for sequenceid in filekmers[filename]:
            kmers = filekmers[filename][sequenceid] 
            sketchdict[sequenceid] = sketch(kmers) #get sketch of kmer array using sketch function which
            
        filesketch[filename] = sketchdict

    #save dict of sketches to file named sketch.txt for retrival
    saveD(filesketch, "sketch.txt")

###########################    find jaccard distances and neighbour joining tree ###########################
"""
this section calculates Jaccard distances between every combination of sequence sketches
also takes quite a while
to create neighbor joining tree showing relations between sequences require Jaccard distances in lower triangle form
ie.
S1 [0]
S2 [0,1]
S3 [0,1,2]
   S1 S2 S3
 where each element is the jaccard distance between the sketches assoscited with the row and collum its in
 diagonals by defination are 0
"""

#read saved dictionary containg all sketch data from file
sketches = readD("sketch.txt") 

# jaques = [] #can create list format Jdistance, (1st file, 1stid), (2nd file,2ndid) 
#would be useful if you wanted Jdistance between 2 sequecnes

#initialise some lists to add to
#Jdistances are for each pair regardles of order so cut down on calcualations by only doing each pair once
combos = [] #keeps track of all pair calculations done

#will eventualy contain top triangle of jaccard distance
matrix = []
#will contain names of sequence corresponding to the axis of the matrix
names = []
#whilst it iterates though the first loop it will go though all file names in the right order 
firstloop = True

#iterate through files
for filename in sketches:

    #iterate though sequences in file
    for sequenceid in sketches[filename]:
        
        row = [] #will be row of matrix containg J distances

        #comapare selected sequence to all other sequences and find the Jdistances

        for otherfilename in sketches: #iterate thoug other files
                       
            for othersequenceid in sketches[otherfilename]: #iterate though sequences in each file
                
                # check if havent done comparison with this other sequence already in either order
                if [sequenceid, othersequenceid] not in combos and [othersequenceid, sequenceid] not in combos: 
                    
                    #find Jdist between the sketches of the 2 sequenes
                    row.append(Jaccard(sketches[filename][sequenceid],sketches[otherfilename][othersequenceid]))

                    #used if wanting entire jauques atrix
                    #jaques.append([Jaccard(sketches[filename][sequenceid],sketches[otherfilename][othersequenceid]), (filename, sequenceid), (otherfilename, othersequenceid)])

                    #add combination to combos so wont be calculated again reduces calcs from N^N to N!
                    combos.append([othersequenceid , sequenceid])

                    #if in first loop where evry sequence is iterated through add sequence Id to list for row collum headers
                    if firstloop == True:
                        names.append(othersequenceid)

        matrix.append(row) #add whole row of Jdistances to matrix, for every iteration this will be 1 smaller crewated upper triangle shape

    firstloop = False #set false for afterfirst loop 

# this section converts upper triangle shape to lower triangle for use in tree creation   
size = len(matrix[0])
distance_matrix = [[0] * i for i in range(1, size+1)] #create lower triangle of 0s correct size

#iterate though upper trig "matrix" to transpose into lower trig "distance matrix"
for i in range(len(matrix)):
    for j in range(len(matrix[i])):

        distance_matrix[j+i][i]= matrix[i][j]   #j+i becuase elements indexed by array position not matrix position             

# Create a DistanceMatrix object for use in tree construction
#names are the sequenceID coresponing to rows and collums of matrix from topleft
#distance matrix is a lower triangle matrix of distances
dm = DistanceMatrix(names=names, matrix=distance_matrix) 

# Construct the Neighbor-Joining tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm) #runs construction of tree on the distancematrix using neighbour joining method

# Print the tree in ascii
Phylo.draw_ascii(tree) 