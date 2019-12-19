# test for PlasClass
import numpy as np
import os
import subprocess
import sys

def test_fasta_classifier():
    fpath = os.path.dirname(os.path.abspath(__file__))
    cmd = 'python ' + os.path.join(os.path.dirname(fpath), 'classify_fasta.py') + ' -f ' + os.path.join(fpath,'test.fa') + ' -o ' + os.path.join(fpath,'test.out')
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    test_dict = {}
    with open(os.path.join(fpath,'test.out')) as f:
        for line in f:
            splt = line.strip().split()
            test_dict[splt[0]] = np.around(float(splt[1]),8)
    gs_dict = {}
    with open(os.path.join(fpath,'test.gs')) as f:
        for line in f:
            splt = line.strip().split()
            gs_dict[splt[0]] = np.around(float(splt[1]),8)
    diff_pass = True
    for k,v in test_dict.items():
        if k not in gs_dict or v != gs_dict[k]:
            diff_pass = False
    for k,v in gs_dict.items():
        if k not in test_dict or v != test_dict[k]:
            diff_pass = False
    os.remove(os.path.join(fpath,'test.out'))
    if not diff_pass:
        print('TEST FAILED!')
        exit(-1)

def test_module():
    # suppress output
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    dn = open(os.devnull,'w')
    sys.stdout = dn
    sys.stderr = dn

    from plasclass import plasclass
    c = plasclass.plasclass(16)
    seq = 'GCTTGAACAAGTTTTTACCAGATGCACTCAAAGGTCTGAAACTCTGAAACTCTAGGACTAAAATAATAGGACACAGAAGATCTTGATCTGATTAGTTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGAATAATCTATTTTCTATATTTTATTGAAAGTTCTTTACTAAAAAGCTAGTAAACAAAAATAAAACCAAGATATGAACGAAAACAAAATGTGAGAAAATCCCGCTCATACTATGGTTTTCAAATCTGATACCGCATTGATTACATTAGCTAAAATCAATTCCAGGTGGCGATTGGGCAACGTTGTGTATGGTGGTTTGTATCGCAAGAAATACGATGCTGACTTAGAAAAAGCTGTTCAATACTATAGTTCG'
    s = c.classify(seq)

    # restore output
    dn.close()
    sys.stdout = old_stdout
    sys.stderr = old_stderr

    if np.around(s,6) != np.around(0.72923195,6):
        print('TEST FAILED!')
        exit(-1)

def all_tests():
    print("Testing...")
    test_fasta_classifier()
    test_module()
    print("Passed all tests")

if __name__=='__main__':
    all_tests()
