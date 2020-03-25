import re

def readFASTA(file):
    """
    It reads FASTA files and returns them as a dictionary
    Works both for one-line and multi-line FASTA files
    """
    regex_name = re.compile("^>(\S+)")
    name = ""
    Protein_dic = {}
    
    file_fd = open(file,'r')
    for line in file_fd:
        line = line.strip()
        n = regex_name.match(line)
        if n:
            name = n.group(1)
            Protein_dic[name]=""
        else:
            Protein_dic[name] += line
            
    return Protein_dic

def printFASTA(dictionary, multiple=True, file=""):
    """
    Prints a fasta dictionary into a file
    if multiple is True it will print a single file (with a name specified in file)
    if multiple is False it will print a different file for each sequence,
    being the file name the sequence ID with .fa extension
    """
    if multiple:
        file_fd = open(file,'w')
        for prot in dictionary:
            file_fd.write(">" + prot + "\n" + dictionary[prot] + "\n")
        file_fd.close()
    else:
        for prot in dictionary:
            file_fd = open(prot+".fa",'w')
            file_fd.write(">" + prot + "\n" + dictionary[prot] + "\n")
            file_fd.close()
    