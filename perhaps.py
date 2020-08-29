import os
import getopt
import sys


description = "python perhaps.py -i <sample IID> -d <sample file location> -s <the chr and positions of SNPs>"


def arg(args):
    try:
        opts, temp = getopt.getopt(args, "hi:d:s:", ['help', 'iid=', 'dir=', 'snps='])
    except getopt.GetoptError:
        print(description)
        sys.exit(2)
    else:
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                print(description)
                sys.exit()
            elif opt in ('-i', '--iid'):
                iid = arg
            elif opt in ('-d', '--dir'):
                loc = arg
            elif opt in ('-s', '--snps'):
                snps = arg
        try:
            return [iid, loc, snps]
        except NameError:
            raise Exception('Not complete arguments!')


def cmd_res(cmd: str):
    with os.popen(cmd) as p:
        return p.read()


def find_file(IID: str, loc: str):
    loc = loc.rstrip('\\')
    file = '\\'.join([loc, IID])
    if os.path.isfile(file + '.bam'):
        os.system(rf".\windows_tools\samtools view {file + '.bam'} > results\{IID}.sam")
    elif os.path.isfile(file + '.sam'):
        os.system(rf"copy {file + '.sam'} results\{IID}.sam")
    else:
        raise Exception('File not find!')


def process(IID: str, loc: str, SNPs: str):
    [Chr, pos] = SNPs.split(':')

    if not os.path.isdir('results'):
        os.mkdir('results')

    find_file(IID, loc)
    readlen = cmd_res(r'.\windows_tools\awk "NR==1 {printf length($10)}" results\%s.sam' % IID).strip()

    os.system(r'.\windows_tools\cut -f 1-10 .\results\%s.sam | '
              r'.\windows_tools\awk "$6 !~/S/ {if ($1 in reads) print reads[$1]\" \"$0; reads[$1]=$0}" '
              r'> results\%s.sam.paired' % (IID, IID))

    os.system(r'.\windows_tools\awk "{if ($9<0) print -$9; else print $9}" results\%s.sam.paired'
              r' | .\windows_tools\uniq | .\windows_tools\sort -n | .\windows_tools\uniq  > results\hap.len' % IID)

    os.system(r'.\windows_tools\awk -v readlen=%s -v c=%s -v pos=%s -f .\windows_tools\main_process.awk'
              r' results\%s.sam.paired  > results\%s.hap' % (readlen, Chr, pos, IID, IID))

    os.system(r'.\windows_tools\awk "($0 ~/SNP1/ && $0~/SNP2/) {$1=$2=\"\"; print $0}" results\%s.hap '
              r'| .\windows_tools\sort | .\windows_tools\uniq -c' % IID)
    for i in [r'results\hap.len', rf'results\{IID}.hap', rf'results\{IID}.sam.paired']:
        os.remove(i)


if __name__ == '__main__':
    # IID = 'NA20525'  ## sample ID
    # rawfile = rf'.\test-data'  ## the location of the BAM or CRAM file, indexed
    # SNPs = '1:159205564-159205704-159205737'  ## the chr and positions of SNPs for directy haplotype detection.
    # process(IID, rawfile, SNPs)
    process(*arg(sys.argv[1:]))
