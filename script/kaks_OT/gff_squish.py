#!/bin/python3

'''
GFF SQUISHER!
'''
import argparse, re

# my_gff_in_path = "./test.gff"
# my_gff_out_path = "./test_squished.gff"
# gff_header = re.compile("^##gff-version 3\n$")
# fasta_start = re.compile("^##FASTA$")
# in_header = True
# in_annot = False


def analyze_contig(gff_path):
    '''
    '''
    print("Analyzing contigs...")
    contig_dict = dict()
    add_to_seq = 0
    in_seq = False
    with open(gff_path, "r") as gff_in:
        for linus in gff_in:
            if re.compile("^>.*$").match(linus):
                in_seq = True
                curr_contig = linus.strip("\n")
                contig_dict[curr_contig] = add_to_seq #add cumsum of the contig size
            if not re.compile("^>").match(linus) and in_seq == True:
                add_to_seq+=len(linus.strip("\n"))
        return(contig_dict)

# print(analyze_contig(my_gff_path))

def squish_contigs(gff_in_path, gff_out_path):
    '''
    '''
    print("Squishing time!")
    contig_dict = analyze_contig(gff_in_path)
    contig_new_name = "squish"
    # print(contig_dict)
    header_passed = False
    annot_passed = False
    with open(gff_out_path, "w") as gff_out:
        with open(gff_in_path, "r") as gff_in:
            for linus in gff_in:
                # print(linus)
                if re.compile("^##").match(linus) and not header_passed:
                    continue
                if not re.compile("^##").match(linus) and not annot_passed: 
                    header_passed = True
                    annot = linus.split("\t")
                    curr_contig = ">" + annot[0]
                    # Shift the annotation positions
                    annot[3] = str(int(annot[3]) + contig_dict[curr_contig])
                    annot[4] = str(int(annot[4]) + contig_dict[curr_contig])
                    # Homogenize the contig names
                    annot[0] = contig_new_name
                    annot = "\t".join(annot)
                    # print(annot)
                    gff_out.write(annot)
                    # print(annot_passed)
                if re.compile("^##FASTA$").match(linus) and header_passed:
                    annot_passed = True
                    # print(linus)
                    gff_out.write(linus)
                    # print(">"+ contig_new_name + "\n")
                    gff_out.write(">" + contig_new_name + "\n")
                elif re.compile("^>").match(linus) and header_passed and annot_passed:
                    continue
                elif header_passed and annot_passed:
                    # print(linus)
                    gff_out.write(linus)
    return(0)

# squish_contigs(my_gff_in_path, my_gff_out_path)

def main():
    # Get command-line arguments
    parser = argparse.ArgumentParser(description = "Squishing a GFF like there's no tomorrow")
    parser.add_argument('--gff_in', type=str, help='goes in')
    parser.add_argument('--gff_out', default=10000, type=str, help='comes out')
    args = parser.parse_args()
    # parser.print_help()
    if args.gff_in == None or args.gff_out == None:
        parser.print_help()
    else: 
        squish_contigs(args.gff_in, args.gff_out)
    
if __name__ == '__main__':
    main()