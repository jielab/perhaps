import os
import getopt
import sys
import subprocess

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


def find_file(IID: str, loc: str):
    loc = loc.rstrip('\\')
    file = '\\'.join([loc, IID])
    if os.path.isfile(file + '.bam'):
        run_cmd(rf".\windows_tools\samtools view {file + '.bam'} > results\{IID}.sam")
    elif os.path.isfile(file + '.sam'):
        run_cmd(rf"copy {file + '.sam'} results\{IID}.sam")
    else:
        raise Exception('File not find!')


def get_res_cmd(cmd: str):
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE) as p:
        result_f = p.stdout
        error_f = p.stderr
        errors = error_f.read()
        if errors:
            pass
        result_str = result_f.read()
        if result_f:
            result_f.close()
        if error_f:
            error_f.close()
        return result_str.decode('GBK')


def show_need(hap_res: str, snp_num: int):
    res = ''
    proper = True
    while True:
        for line in hap_res.splitlines():
            if int(line.split(' ')[-1]) == snp_num:
                res += line + '\n'
        if res != '':
            if proper:
                return 'All paired haplotype results:\n' + res, snp_num
            else:
                return f"Can't find the all paired haplotype, show maximum paired number ({snp_num}) " \
                       f"haplotype result:\n" + res, snp_num
        elif snp_num > 1:
            snp_num -= 1
            proper = False
        else:
            return 'None paired sequences find!', snp_num


def run_cmd(cmd: str):
    a = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    a.communicate(timeout=15)


def get_paired_sam(IID: str, paired_file='.\\results\\paired.temp'):
    with open(rf'.\\results\\{IID}.sam', 'r') as raw:
        with open(rf'.\\results\\{IID}.subset.sam', 'w') as new:
            with open(paired_file, 'r') as need:
                need_id = need.readlines()
                need_id = [ID.strip() for ID in need_id]
            for line in raw:
                if line.split('\t')[0] in need_id:
                    new.write(line)


def process(IID: str, loc: str, SNPs: str):
    [Chr, pos] = SNPs.split(':')
    pos_n = pos.count('-') + 1

    if not os.path.isdir('results'):
        os.mkdir('results')
    if not os.path.isdir('windows_tools'):
        raise Exception("Can't not find the windows_tools folder!")

    find_file(IID, loc)
    readlen = get_res_cmd(r'.\windows_tools\awk "NR==1 {printf length($10)}" results\%s.sam' % IID).strip()

    run_cmd(r'.\windows_tools\cut -f 1-10 .\results\%s.sam | '
            r'.\windows_tools\awk "$6 !~/S/ {if ($1 in reads) print reads[$1]\" \"$0; reads[$1]=$0}" '
            r'> results\%s.sam.paired' % (IID, IID))

    run_cmd(r'.\windows_tools\awk "{if ($9<0) print -$9; else print $9}" results\%s.sam.paired'
            r' | .\windows_tools\uniq | .\windows_tools\sort -n | .\windows_tools\uniq  > results\hap.len' % IID)

    run_cmd(r'.\windows_tools\awk -v readlen=%s -v c=%s -v pos=%s -f .\windows_tools\main_process.awk'
            r' results\%s.sam.paired | .\windows_tools\sort -k 4,4nr -k 1,1n '
            r'> results\%s.hap' % (readlen, Chr, pos, IID, IID))
    a = get_res_cmd(r'.\windows_tools\awk "($4 != \"\") {$1=$2=\"\"; print $0}" results\%s.hap '
                    r'| .\windows_tools\sort | .\windows_tools\uniq -c' % IID)
    res, num = show_need(a, pos_n)
    if num:
        run_cmd(r'.\windows_tools\awk -v num=%s "($4 == num) {print $2}" results\%s.hap '
                r'> results\paired.temp' % (num, IID))
        get_paired_sam(IID, paired_file=r'results\paired.temp')
    else:
        raise Exception('None paired haplotype find!')
    for i in [r'results\hap.len', rf'results\{IID}.hap', rf'results\{IID}.sam.paired', rf'results\{IID}.sam',
              r'results\paired.temp']:
        os.remove(i)

    return res


if __name__ == '__main__':
    process(*arg(sys.argv[1:]))
