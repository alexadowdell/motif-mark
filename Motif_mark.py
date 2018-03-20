#!/usr/bin/python3

#from IPython.core.display import SVG, display
import argparse
import cairo
import random

#####################################################################################################################
############################################ Set up Argparse ########################################################
#####################################################################################################################

parser = argparse.ArgumentParser(description="Motif Marker: Identifies and counts motif sequences that appear in gene sequences from a FASTA file and constructs an image of the motif locations in relation to intron/exons.")
parser.add_argument("-f", "--file", help="Required parameter. Please pass a FASTA file with gene names in the sequence ID and exons CAPITALIZED", required=True, type=str)
parser.add_argument("-m", "--motif", help="Require parameter. Please pass of text file of motifs (one per line) to be used to look for. Motifs must use IUPAC nucleotide code format", required=True, type=str)

args = parser.parse_args() #sets arguments as call-able 

f = args.file
m = args.motif

#f = "./sequences.fasta"
#m = "./motifs.txt"


#####################################################################################################################
############################################ Define main functions ##################################################
#####################################################################################################################

#Initiate and construct a dictionary with all IUPAC nucleotide codes to bases. Keys = nucleotide code, Value = Bases
iupac_dict = {"A": "[Aa]", "C":"[Cc]", "G":"[Gg]", "T":"[TtUu]","U":"[UuTt]", "R": "[AaGg]", "Y":"[TtCcUu]", 
              "S":"[CcGg]", "W":"[AaTtUu]","K":"[GgTtUu]", "M":"[AaCc]", "B":"[CcGgTtUu]", "D":"[AaGgTtUu]",
              "H":"[AaCcTtUu]", "V":"[AaCcGg]", "N":"[AaTtCcGgUu]", "a": "[Aa]", "c":"[Cc]", "g":"[Gg]", 
              "t":"[TtUu]","u":"[UuTt]", "r": "[AaGg]", "y":"[TtCcUu]", "s":"[CcGg]", "w":"[AaTtUu]","k":"[GgTtUu]", 
              "m":"[AaCc]", "b":"[CcGgTtUu]", "d":"[AaGgTtUu]", "h":"[AaCcTtUu]", "v":"[AaCcGg]", "n":"[AaTtCcGgUu]"}

def split_sequence(sequence):
    '''Parses fasta sequence and splits up a sequence in the form of 'preintronEXONpostintron', into seperated
    items in a list [intron_1, exon, intron_2] '''
    intron_1 = ""
    intron_2 = ""
    exon = ""
    exon_identified = False

    for char in sequence:
        if (char.islower() == False): # If lowercase letters stop being found, the exon was located. 
            exon +=  char #Exon identified
            exon_identified = True
        elif (exon_identified == False): #If the exon was not identified yet, the
            intron_1 = intron_1 + char
        else:
            intron_2 = intron_2 + char
    return [intron_1, exon, intron_2]


def fasta_conversion(fastaFile):
    ''' Takes in a fasta file. Construct a dictionary parsing the fasta with the key being the sequence header
    and the value being the full sequence in one string. Returns a fasta dictionary with entries in the form
    {Header:Sequence}'''
    with open(fastaFile) as fh:
        fastaDict = {}
        header = ""
        sequence = ""
        firstLine = True

        for line in fh:
            line = line.strip()
            if line[0] == ">": #if the line starts with a '>' its a header
                if (firstLine != True): #
                    fastaDict[header] = split_sequence(sequence) #for all additional lines, 
                header = line #if first line is true, set the header line
                sequence = "" #still first line, initialize the sequence 
            else:
                sequence += line
                firstLine = False

    fastaDict[header] = split_sequence(sequence) #special case
    return fastaDict

    
def parse_motifs(motifFile):
    '''Takes in a text file of motifs and creates a list of the motifs. For each motif, each nucleotide code is 
    evaluated and the corresponding IUPAC bases are associated with it. Dictionary is returned key=motif, value=list 
    of sets of bases corresponding to each character in the motif.'''
    motifList = []
    motifDict = {}
    with open(motifFile) as fh:
        for line in fh:
            motifList.append(line.strip()) #Create a list of all motifs in file

    for motif in motifList: #for each motif in list
        regex = []
        for char in motif: #look at each nucleotide individually
            regex.append(iupac_dict[char]) #finding corresponding base codes and return them as the value in the dictionary
        motifDict[motif] = regex #Dictionary: key=motif, value=list associated base combinations

    return motifDict

    
def extract_motif_pos(header, sequence, motifDict):
    '''Takes in a motif header, associated sequence and motif dictionary. Walks through a sequence looking for instances 
    of the associated motif by referencing the possible combinations of bases given the IUPAC nucleotide codes. Returns
    a list of positions the motif was found at in the sequence'''
    positionList = {} #Initialize a dictionary of motif positions
    for motif in motifDict.keys(): #For each motif in the dictionary
        motif_codes = motifDict[motif] #set the motif's associated regex nucleotide codes to a variable.
        frame = len(motif)
        for i in range(0,len(sequence) - frame):
            current_frame = sequence[i:i + frame] #all possible frames across the sequence the length of the motif, as bases
            counter = 0
            motif_match = True
            for base in current_frame: 
                if base not in motif_codes[counter]:#if base does not match one in the possible motif codes
                    motif_match = False
                counter += 1
            if (motif_match == True): #if motif match is found
                if header + "_" + motif in positionList:
                    positionList[header + "_" + motif].append(i) #extract the position and add it to the value of that header_motif in list form
                else:
                    positionList[header + "_" + motif] = [i] #if the motif was not found yet, add it as a new entry to the dictionary, 
                    #value is a new initialized list with the position

    return positionList


def identify_all_motifs(fastaFile, motifFile):
    '''Takes in a fasta file and a motif text file. Returns a positions list containing motifs and the occurences of those motifs
    for each fasta sequence, as well as the locations of the pre-intron, exon, and post-intron.'''
    fastaDict = fasta_conversion(fastaFile) #Create fasta dictionary {header:sequence}
    motifDict = parse_motifs(motifFile) #Create dictionary of motifs {motif:nucleotide code combinations in IUPAC format for the motif}
    positionList = []
    
    for entry in fastaDict.items(): 
        sequence = "".join(fastaDict[entry[0]])
        header = entry[0]
        content = extract_motif_pos(header, sequence, motifDict) #extract all positions of motifs
        content["intron_1"] = len(fastaDict[entry[0]][0]) #extract all positions of pre-introns
        content["exon"] = len(fastaDict[entry[0]][1]) #extract all positions of exons
        content["intron_2"] = len(fastaDict[entry[0]][2]) #extract all positions of post-introns
        positionList.append(content) #Create one global position list containing all motif positions, pre-introns, exons,
        # and post-introns for each sequence passed in.

    return positionList
    
#####################################################################################################################
######################################## Define plotting functions ##################################################
#####################################################################################################################

def add_gene_name(start_pos, gene):
    '''Add gene name to the left of gene drawn in black'''
    context.set_source_rgb(0.4,0.4,0.4)
    context.move_to(start_pos[0] - 100, start_pos[1]) #defines location of gene name placement
    context.show_text(gene) #prints gene name
    
    
def add_labels(r, g, b, alpha, start_pos, adjustment, motif):
    '''Add motif names in the random colors the markers were generated in'''
    context.set_source_rgba(r, g, b, alpha)
    context.move_to(start_pos[0] - 57, start_pos[1] + adjustment) #adjust motif label placement
    context.show_text(motif) #prints motif name

    
def draw_gene(intron_1 , intron_2, exon_len, context, start_pos):
    '''Draw the gene in grey/black intron, exon, intron, with the exon being a wider width rectangle'''
    context.set_line_width(1)
    context.set_source_rgb(0.4,0.4,0.4) #grey/black
    context.move_to(start_pos[0], start_pos[1])
    context.line_to(start_pos[0] + intron_1 + exon_len + intron_2, start_pos[1])
    context.stroke()
    context.set_line_width(10)
    context.move_to(start_pos[0] + intron_1, start_pos[1])
    context.line_to(start_pos[0] + intron_1 + exon_len, start_pos[1])
    context.stroke()

    
def draw_motifs(positionList, context, start_pos, motif):
    '''Draw markers at the position each motif was found in the sequence. Randomly generated 
    color corresponding to gene labels'''
    r = random.random()
    g = random.random()
    b = random.random()
    alpha = random.uniform(0.7, 1.0)

    context.set_line_width(10)
    context.set_source_rgba(r, g, b, alpha)
    for position in positionList:
        context.move_to(start_pos[0] + position, start_pos[1])
        context.line_to(start_pos[0] + position + len(motif), start_pos[1])
        context.stroke()
        context.move_to(start_pos[0] + position, start_pos[1] - 20)
        context.show_text(str(position))

    return [r, g, b, alpha]
    
    
#####################################################################################################################
################################################ Main code ##########################################################
#####################################################################################################################

#Main command, finds all motif positions, intron, and exon positions
global_positions = identify_all_motifs(f, m)

#Construct a drawing surface and identify a start position
surface = cairo.SVGSurface("motif_plot.svg", 1500, 800)
context = cairo.Context(surface)
context.set_line_width(10)
start_pos = [200,150]

for entry in global_positions:
    draw_gene(entry["intron_1"],entry["intron_2"],entry["exon"], context, start_pos)
    adjustment = -10 
    for item in entry.items():
        if item[0][0] == ">":
            motif = item[0].split("_")
            motif = motif[1]
            gene = item[0].split(" ")
            gene = gene[0][1:]
            color = draw_motifs(item[1], context, start_pos, motif)
            add_gene_name(start_pos, gene)
            add_labels(color[0], color[1], color[2], color[3], start_pos, adjustment, motif)
            adjustment = adjustment + 10
    start_pos = [start_pos[0], start_pos[1] + 100]

#display(SVG('example.svg'))

