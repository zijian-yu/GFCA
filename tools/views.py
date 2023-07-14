from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from GFCA.tools import *
from django.urls import reverse
from Bio.Seq import Seq


def index(request):
    return render(request, 'index.html')


def tools(request):
    return render(request, 'tools.html')


def gfca(request):
    return render(request, 'gfca.html')


def database(request):
    return render(request, 'database.html')


def about(request):
    return render(request, 'about.html')


def contact(request):
    return render(request, 'contact.html')


def help(request):
    return render(request, 'help.html')


def successful(request):
    return render(request, 'successful.html')


def blast_match(request):
    return render(request, 'blast_match.html')


def successful_2(request):
    return render(request, 'successful_2.html')

# 分割线----------------------------------------tools模块
def dotplot(request):
    return render(request, 'dotplot.html')


def file_merge(request):
    return render(request, 'file_merge.html')


def formatTofasta(request):
    return render(request, 'formatTofasta.html')


def find_replace(request):
    return render(request, 'find_replace.html')


def extract_row(request):
    return render(request, 'extract_row.html')


def duplicate_removal(request):
    return render(request, 'quchong.html')


# def cds_convert_pep(request):
#     return render(request, 'cds_convert_pep.html')


# index_email
def index_email(request):
    name = request.POST.get("name", None)
    email = request.POST.get("email", None)
    message = request.POST.get("message", None)

    if email is None or name is None or message is None:
        return render(request, '../templates/contact.html')
    else:
        c = name_set()
        keep_email_index(place=c, name=name, email=email, message=message)
        return render(request, '../templates/successful.html')


# dotplot
def dotplot_run(request):
    # 五个文件
    blast = request.FILES.get('blast', None)
    heng_len = request.FILES.get('heng_len', None)
    shu_len = request.FILES.get('shu_len', None)
    heng_gff = request.FILES.get('heng_gff', None)
    shu_gff = request.FILES.get('shu_gff', None)

    # 四种参数-值
    evalue = str(request.POST.get('evalue', None))
    score = request.POST.get('score', None)
    hitnum = request.POST.get('hitnum', None)
    repnum = request.POST.get('repnum', None)

    # 需要保存的染色体数
    heng_chromosomes = str(request.POST.get('heng_chromosomes', None))
    shu_chromosomes = str(request.POST.get('shu_chromosomes', None))

    # 需要保存的物种名称
    heng_species = str(request.POST.get('heng_species', None))
    shu_species = str(request.POST.get('shu_species', None))

    email = request.POST.get('email', None)

    if blast is None or heng_len is None or shu_len is None or heng_gff is None or shu_gff is None:
        return render(request, '../templates/dotplot.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        blast_name = blast.name
        shu_gff_name = shu_gff.name
        heng_gff_name = heng_gff.name
        shu_len_name = shu_len.name
        heng_len_name = heng_len.name

        keep_file(na=blast, c=c)
        keep_file(na=heng_len, c=c)
        keep_file(na=shu_len, c=c)
        keep_file(na=heng_gff, c=c)
        keep_file(na=shu_gff, c=c)

        dotplot_run_tools(place=c, email=email, evalue=evalue, score=score, hitnum=hitnum, repnum=repnum,
                          heng_chromosomes=heng_chromosomes, shu_chromosomes=shu_chromosomes,
                          heng_species=heng_species, shu_species=shu_species,
                          blast_name=blast_name,
                          shu_gff_name=shu_gff_name, heng_gff_name=heng_gff_name,
                          shu_len_name=shu_len_name, heng_len_name=heng_len_name
                          )
        return render(request, '../templates/successful.html')


# file_merge
def file_merge_run(request):

    file = request.FILES.getlist("file", None)
    email = request.POST.get("email", None)
    # file_name = request.POST.get("filename", None)

    if file is None:
        return render(request, '../templates/file_merge.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        keep_files(na=file, c=c)

        # file_merge_run_tools(place=c, email=email, file_name=file_name)
        file_merge_run_tools(place=c, email=email)
        return render(request, '../templates/successful.html')


# find_replace
def find_replace_run(request):
    file_up = request.FILES.get("file", None)
    old_name = request.POST.get("old_name", None)
    new_name = request.POST.get("new_name", None)
    email = request.POST.get("email", None)

    if file_up is None:
        return render(request, '../templates/find_replace.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        file_up.name = "upload.txt"
        keep_file(na=file_up, c=c)
        find_replace_run_tools(place=c, email=email, old_name=old_name, new_name=new_name)
        return render(request, '../templates/successful.html')


# Extract_column
def extract_column_run(request):
    row = request.FILES.get("row", None)
    row_num = request.POST.get("row_num", None)
    email = request.POST.get("email", None)

    if row is None:
        return render(request, '../templates/extract_row.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        row_name = row.name.split(".")[0]
        row_last = row.name.split(".")[1]
        keep_file(na=row, c=c)

        extr_row_run(place=c, email=email, row_num=row_num, row_name=row_name, row_last=row_last)
        return render(request, '../templates/successful.html')


# duplicate_removal_run
def duplicate_removal_run(request):
    file = request.FILES.get("file", None)
    email = request.POST.get('email', None)

    if file is None:
        return render(request, '../templates/quchong.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        file.name = '1.fasta'
        keep_file(na=file, c=c)

        quchong_run(place=c, email=email)
        return render(request, '../templates/successful.html')


# formatTofasta_run
def formatTofasta_run(request):
    format = request.FILES.get("file1", None)
    patten = request.POST.get("patten", None)
    patten_end = request.POST.get("patten_end", None)
    email = request.POST.get("email", None)

    if format is None:
        return render(request, '../templates/formatTofasta.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        file_name = format.name.split(".")[0]
        file_last = format.name.split(".")[1]
        keep_file(na=format, c=c)

        format_fasta_run(place=c, email=email, file_name=file_name, patten=patten, patten_end=patten_end, file_last=file_last)
        return render(request, '../templates/successful.html')


# blast_match_run
def blast_match_run(request):
    patten = request.POST.get("patten", None)
    species_data = request.POST.get("species_data", None)
    myFile1 = request.FILES.get("file1", None)

    e_value = str(request.POST.get('e_value', None))
    sco = request.POST.get('score_1', None)
    iden = request.POST.get('identity_1', None)
    email = request.POST.get("email", None)

    nucl_or_prot = "prot" if patten == 'blastp' else "nucl"

    if myFile1 is None:
        return render(request, 'blast_match.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        myFile1.name = '1.fasta'
        keep_file(na=myFile1, c=c)

        yishan_run_blast(place=c, patten=patten, species_data=species_data, e=e_value, email=email, score=sco, identity=iden, nucl_or_prot=nucl_or_prot)
        return render(request, '../templates/successful_2.html')


# -----------------------------------------开始-----------------------------------------
# cds 转化为 pep 序列，并将结果返回到页面：输入的字符串
def cds_convert_pep_sequences(request, cds_sequences):
    gene_sequences = {}
    current_gene_id = ""
    for line in cds_sequences.split("\n"):
        if line.startswith(">"):
            if current_gene_id:
                gene_sequences[current_gene_id] = current_sequence
            current_gene_id = line[1:]
            current_sequence = ""
        else:
            try:
                current_sequence += line.strip()
            except Exception as e:
                print("Error:", e)
                current_sequence = {}
                current_sequence[">"] = "Error: Please enter CDS sequences with gene ID."
                return render(request, 'cds_convert_pep.html', {'pep_sequences': current_sequence})

    if current_gene_id:
        gene_sequences[current_gene_id] = current_sequence

    # 将cds序列，翻译成pep序列，并返回结果
    pep_sequences = {}
    for gene_id, sequence in gene_sequences.items():
        try:
            pep_sequence = Seq(sequence).translate() + "\n"
            pep_sequences[">" + gene_id] = str(pep_sequence)
        except Exception as e:
            pep_sequences[">"] = "Error: Cannot translate sequence, Please check your cds sequence!"
    return render(request, 'cds_convert_pep.html', {'pep_sequences': pep_sequences})


# cds 转化为 pep 序列，并将结果返回到页面：文件
def cds_convert_pep_file(request, place, cds_file):
    pep_sequences = {}
    try:
        for line in SeqIO.parse(f"{path_get}/file_keep/{place}/{cds_file}", "fasta"):
            pep_sequences[">" + line.id + "\n"] = str(Seq(str(line.seq)).translate()) + "\n"
    except Exception as e:
        pep_sequences[">"] = "Error: Cannot translate sequence, Please check your cds file!"
    return render(request, 'cds_convert_pep.html', {'pep_sequences': pep_sequences})


# cds_convert_pep - 将 cds 序列转化为 pep 序列，并返回到页面，未使用 tools.py!
def cds_convert_pep(request):
    cds_sequences = request.POST.get("cds_sequences", None)
    cds_file = request.FILES.get("cds_file", None)

    # cds_sequences 分析
    if cds_file is None:
        if cds_sequences is None:
            return render(request, 'cds_convert_pep.html')
        else:
            return cds_convert_pep_sequences(request, cds_sequences)

    # 将 cds_file 文件，存在创建的文件夹中
    elif cds_sequences is '':
        if cds_file is None:
            return render(request, 'cds_convert_pep.html')
        else:
            c = name_set()
            keep_file(na=cds_file, c=c)
            return cds_convert_pep_file(request, c, cds_file.name)

    elif cds_file and cds_sequences:
        if cds_sequences is None:
            return render(request, 'cds_convert_pep.html')
        else:
            return cds_convert_pep_sequences(request, cds_sequences)
    
    elif cds_file is None and cds_sequences is None:
        return render(request, 'cds_convert_pep.html')
# -----------------------------------------结束-----------------------------------------



