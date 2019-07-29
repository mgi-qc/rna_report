import csv
import os
import sys
import glob
import argparse
import subprocess
from string import Template

parser = argparse.ArgumentParser()
parser.add_argument('-hu', help='Human metrics', action='store_true')
parser.add_argument('-m', help='Mouse metrics', action='store_true')
parser.add_argument('-e', help='ercc values', action='store_true')
args = parser.parse_args()

#test!
args.hu = True

if not (args.hu or args.m):
    sys.exit('Please specify -hu (human) or -m (mouse) flag.')

if args.hu:
    ar = 0.70

if args.m:
    ar = 0.60

if not glob.glob('*.picard.analysis.metrics.tsv'):
    sys.exit('picard.analysis.metrics.tsv file not found')
else:
    picard_file = glob.glob('*.picard.analysis.metrics.tsv')[0]

if args.e:
    ercc_file = picard_file.split('.')[0] + '.ercc.tsv'
    if not os.path.isfile(ercc_file):
        sys.exit('ercc file not found')

outfile = picard_file.split('.')[0] + '.picard.analysis.metrics.results.tsv'
report_outfile = picard_file.split('.')[0] + '.picard.analysis.metrics.report.txt'

results = {}
template_file_dict = {}

metrics_tracked = ['PASS_SAMPLES', 'ALN_FAIL', 'PCT_RIB_FAIL', 'PCT_SUM_FAIL', 'ERCC_FAIL']
for metric in metrics_tracked:
    template_file_dict[metric] = 0

if not os.path.isfile('/gscmnt/gc2783/qc/GMSworkorders/reports/RNA_report_template.txt'):
    sys.exit('Template file not found.')

with open('/gscmnt/gc2783/qc/GMSworkorders/reports/RNA_report_template.txt', 'r', encoding='utf-8') as fh:
    template = fh.read()
    template_file = Template(template)


with open(picard_file, 'r') as fh:

    picard_dict = csv.DictReader(fh, delimiter='\t')
    header = picard_dict.fieldnames

    last_succeeded_build_id = []
    tot_aln_rate = tot_ercc = tot_pct_cod = tot_int_bases = tot_rib_bases = 0
    tot_pct_sum = tot_pf_aln_bases = tot_cor_strnd_reads = 0

    for line in picard_dict:

        preferred_failures = []
        template_file_dict['WOID'] = line['WOID']

        line['ERCC_Value'] = 0
        if args.e:
            line['ERCC_Value'] = subprocess.check_output(['grep', '-m', '1', line['subject.name'],
                                                          ercc_file]).decode('utf-8').strip().split('\t')[1].split(':')[1]
            if float(line['ERCC_Value']) < 0.9:
                preferred_failures.append('ERCC_Value')
                template_file_dict['ERCC_FAIL'] += 1

        if float(line['ALIGNMENT_RATE']) < ar:
            line['QC_Status'] = 'Fail'
            template_file_dict['ALN_FAIL'] += 1
        else:
            line['QC_Status'] = 'Pass'
            template_file_dict['PASS_SAMPLES'] += 1

        if float(line['PCT_RIBOSOMAL_BASES']) > 0.05:
            preferred_failures.append('PCT_RIBOSOMAL_BASES')
            template_file_dict['PCT_RIB_FAIL'] += 1

        if float(line['PCT_SUM']) < 0.75:
            preferred_failures.append('PCT_SUM')
            template_file_dict['PCT_SUM_FAIL'] += 1

        line['Preferred_Metric_Failures'] = 'None'
        if len(preferred_failures) != 0:
            line['Preferred_Metric_Failures'] = ','.join(preferred_failures)

        tot_aln_rate += float(line['ALIGNMENT_RATE'])
        tot_ercc += float(line['ERCC_Value'])
        tot_pct_cod += float(line['PCT_CODING_BASES'])
        tot_int_bases += float(line['PCT_INTERGENIC_BASES'])
        tot_rib_bases += float(line['PCT_RIBOSOMAL_BASES'])
        tot_pct_sum += float(line['PCT_SUM'])
        tot_pf_aln_bases += float(line['PF_ALIGNED_BASES'])
        tot_cor_strnd_reads += float(line['PCT_CORRECT_STRAND_READS'])

        last_succeeded_build_id.append(line['last_succeeded_build.id'])

        results[line['subject.name']] = line

header.extend(['ERCC_Value', 'QC_Status', 'Preferred_Metric_Failures'])

seq_notes = []

print('Confluence Link:\nhttps://confluence.ris.wustl.edu/pages/viewpage.action?spaceKey=AD&title=WorkOrder+{}'.format(template_file_dict['WOID']))
seq_in = input('Would you like to add a sequencing note? (y/n)')

while True:
    if seq_in == 'y':
        #get seq note
        print('Enter "return q return" to exit when finished:')
        while True:
            note_line = input()
            if note_line != 'q':
                seq_notes.append(note_line)
            else:
                print('Skipping sequencing note')
                break
        break

    elif seq_in == 'n':
        break

    else:
        print('Please enter "y" or "n"')

seq_note = ''
for l in seq_notes:
    seq_note += l + '\n'

with open(report_outfile, 'w', encoding='utf-8') as fh:
    fh.write(template_file.substitute(WOID=template_file_dict['WOID'],
                                      SAMPLE_NUMBER=len(results),
                                      PASS_SAMPLES=template_file_dict['PASS_SAMPLES'],
                                      ALN_FAIL=template_file_dict['ALN_FAIL'],
                                      PCT_RIB_FAIL=template_file_dict['PCT_RIB_FAIL'],
                                      PCT_SUM_FAIL=template_file_dict['PCT_SUM_FAIL'],
                                      ERCC_FAIL=template_file_dict['ERCC_FAIL'],
                                      AVG_ALN_RATE=tot_aln_rate / len(results),
                                      AVG_ERCC=tot_ercc / len(results),
                                      AVG_CODING_BASES=tot_pct_cod / len(results),
                                      AVG_PCT_INT_BASES=tot_int_bases / len(results),
                                      AVG_RIB_BASES=tot_rib_bases / len(results),
                                      AVG_PCT_SUM=tot_pct_sum / len(results),
                                      AVG_PF_ALIGNED_BASES=tot_pf_aln_bases / len(results),
                                      AVG_PCT_CORRECT_STRAND_READS=tot_cor_strnd_reads / len(results),
                                      BUILDS=','.join(last_succeeded_build_id),
                                      RESULTS_SPREADSHEET=outfile,
                                      SEQ_NOTE=seq_note))

with open(outfile, 'w') as of:
    ofd = csv.DictWriter(of, fieldnames=header, delimiter='\t')
    ofd.writeheader()
    for v in results:
        ofd.writerow(results[v])



