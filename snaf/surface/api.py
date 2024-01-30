import requests
import sys,os
import subprocess
import re

def run_ensembl(ens):
    # return the string for sequence, either full_length (using ENST) or peptide (using ENSP)
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/{}?".format(ens)
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    return r.text



def run_emboss(asequence,bsequence,python_executable):
    # return the string for alignment
    current_pwd = os.getcwd()
    os.chdir(os.path.dirname(__file__))
    command = '{} emboss.py --email whatever@gmail.com --stype protein --asequence {} --bsequence {}'.format(python_executable,asequence,bsequence)
    capture_pat = re.compile(r'Creating result file: (.*?.out.txt)')
    stdout = subprocess.run(command,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout
    outfile = re.search(capture_pat,stdout).group(1)
    with open('{}'.format(outfile),'r') as f:
        string_stream = ''
        f_list = f.readlines()
        for i,line in enumerate(f_list):
            if line.startswith('EMBOSS_001'):
                string_stream += f_list[i]
                string_stream += f_list[i+1]
    string_stream = '\n'.join(string_stream.split('\n')[:-1])
    clear_commands = ['rm *.txt','rm *.params']
    for clear_command in clear_commands:
        subprocess.run(clear_command,shell=True)
    os.chdir(current_pwd)
    return string_stream
    





