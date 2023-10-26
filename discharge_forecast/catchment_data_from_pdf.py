import PyPDF2 as pypdf
from difflib import SequenceMatcher
import re
from glob import glob
import pandas as pd
import os

# import configs for the project:
from discharge_forecast import config

# location of the pdfs:
dpath = config.proj_base + '/data/catchment_characteristics'

#------------------DERIVE TABLE HEADERS BY COMPARING .pdfs------------------#
# load two example .pdfs to compare and read out common fields as parameter/column names

# obtain data:
with open(os.path.join(dpath,'Nedbørfeltparam-Bjørgum.pdf'), 'rb') as pdfFileObj:
    pdfReader = pypdf.PdfReader(pdfFileObj)
    pageObj = pdfReader.pages[0]
    pdftxt1 = pageObj.extract_text()

with open(os.path.join(dpath,'Nedbørfeltparam-Bjørgum.pdf'), 'rb') as pdfFileObj:
    pdfReader = pypdf.PdfReader(pdfFileObj)
    pageObj = pdfReader.pages[0]
    pdftxt2 = pageObj.extract_text()


# compare line by line:
param_header = []
for l1,l2 in zip (pdftxt1.split('\n'),pdftxt2.split('\n')):
    flag = 0

    l1_splt = l1.split(':')[0]
    l2_splt = l2.split(':')[0]
    
    match = SequenceMatcher(None, l1_splt, l2_splt).find_longest_match() 
    head = l1[match.a:match.a+match.size]
    
    try:
        end_parenth = re.search(r'\)',head).end()
        head = head[:end_parenth]
        flag = 1
    except:
        pass
    
    # find headers where (part of) a number/parameter value ended up in the header:
    if flag == 0:
        try:
            begin_numer = re.search(r' (\d+)',head).start()
            head = head[:begin_numer]
        except:
            pass

    try:
        begin_numer = re.search(r' -',head).start()
        head = head[:begin_numer]
    except:
        pass

    if head != l2:
        param_header.append(head)



#------------------READ DATA FROM ALL CATCHMENT AREA FILES------------------#
# go through each file to extract values matching the previously identified headers
pdf_ls = sorted(glob(os.path.join(dpath,'*.pdf')))

df = {}
unit = {}

for pdf_file in pdf_ls:
    catchment = pdf_file.split('.')[0].split('-')[1]

    # fix catchment name:
    catchment = re.sub(r'([a-z])(\d{1,2})', r'\1 \2', catchment.lower())

    df[catchment] = {}

    print(catchment,pdf_file)

    with open(pdf_file, 'rb') as pdfFileObj:
        pdfReader = pypdf.PdfReader(pdfFileObj)
        pageObj = pdfReader.pages[0]
        pdftxt = pageObj.extract_text()

    for line in pdftxt.split('\n'):
        # print(line)
        pheader = [ph for ph in param_header if ph in line]
        if len(pheader) == 1:
            phead = pheader[0]
            coll = []
            for st in phead:
                if st in ['(',')']:
                    coll.append('\\'+st)
                else:
                    coll.append(st)
            search_head = ''.join(coll)
            # print(phead)
            # read out part of line that does not equal header:
            endl = re.search(search_head,line).end()
            # print('{:}: {:}'.format(phead,line[endl:].split(':')[-1]))
            val_string = line[endl:].split(':')[-1].strip()
            
            # try extracting numerical value:
            phead = phead.replace('ø','o').replace('å','a').replace('æ','ae').strip()

            # exclude some headers that don't have number values but might contain numerals:
            if phead not in ['Kartdatum','Projeksjon','Beregn.punkt','Vassdragsnr.','Rapportdato','Kommune.','Fylke.','Vassdrag.','Norges vassdrags- og energidirektoratKartbakgrunn']:
                try:
                    pval_str = re.sub(r'[^0-9.\-]','',val_string)
                    unit[phead] = re.sub(r'[0-9.\-]','',val_string)
                    if '-' in pval_str.split('.')[-1]:
                        pval_str = pval_str[:re.search(r'-',pval_str).start()]
                        unit[phead] = val_string[re.search(pval_str,val_string).end():]
                    if 'Hypsografisk kurve' in unit[phead]:
                        unit[phead] = unit[phead][:1]
                    param_val = float(pval_str)
                except:
                    param_val = val_string
            else:
                param_val = val_string
            df[catchment][phead] = param_val



# Convert dictionary to pandas dataframe
DF = pd.DataFrame(df).drop(index='Norges vassdrags- og energidirektoratKartbakgrunn')

# fix some of the DF.stat_ids:


#DF.insert(0, 'unit', unit)
# put in placeholder for unit column where no unit could be derived:
#DF.unit= DF.unit.fillna('-')
# replace 
#DF = DF.replace('-m',pd.NA)

# Save to .csv file, choosing `encoding='utf-8-sig'` seems to be a good choice
trgt_path = config.proj_base + '/results/catchment_properties/smaakraft'
DF.to_csv(os.path.join(trgt_path,'feltparam_smaakraft.csv'),encoding='utf-8-sig')
# DF.to_csv(os.path.join(dpath,'nedborfeltparam-smaakraft.csv'))