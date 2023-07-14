from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from GFCA.tools import *
from django.urls import reverse


def blast(request):
    return render(request, 'blast.html')


def cd_hit(request):
    return render(request, 'cd_hit.html')


def circos(request):
    return render(request, 'circos.html')


def codonw(request):
    return render(request, 'codonw.html')


def colinearscan(request):
    return render(request, 'colinearscan.html')


def cpgfinder(request):
    return render(request, 'cpgfinder.html')


def diamond(request):
    return render(request, 'diamond.html')


def dupgen_finder(request):
    return render(request, 'dupgen_finder.html')


def fasttree(request):
    return render(request, 'fasttree.html')


def gsds(request):
    return render(request, 'gsds.html')


def heatmap(request):
    return render(request, 'heatmap.html')


def hmmer(request):
    return render(request, 'hmmer.html')


def iq_tree(request):
    return render(request, 'iq_tree.html')


def kaks_calculator(request):
    return render(request, 'kaks_calculator.html')


def mcscanx(request):
    return render(request, 'mcscanx.html')


def meme_suit(request):
    return render(request, 'meme_suit.html')


def orthofinder(request):
    return render(request, 'orthofinder.html')


def orthomcl(request):
    return render(request, 'orthomcl.html')


def paml(request):
    return render(request, 'paml.html')


def paraat(request):
    return render(request, 'paraat.html')


def pfam(request):
    return render(request, 'pfam.html')


def phyml(request):
    return render(request, 'phyml.html')


def sequence_alignment(request):
    return render(request, 'sequence_alignment.html')


def sequence_screening(request):
    return render(request, 'sequence_screening.html')


def sequenceserver(request):
    return render(request, 'sequenceserver.html')


def shinycircos(request):
    return render(request, 'shinycircos.html')


def string(request):
    return render(request, 'string.html')


# blast
def blast_run(request):
    patten = request.POST.get('patten', None)
    myFile1 = request.FILES.get('file1', None)  # abcd.fa
    myFile2 = request.FILES.get('file2', None)
    email = request.POST.get('email', None)
    e_value = str(request.POST.get('e_value', None))
    sco = request.POST.get('score_1', None)
    iden = request.POST.get('identity_1', None)
    nucl_or_prot = "prot" if patten == 'blastp' else "nucl"

    # 判断文件是否为空 如果为空则刷新页面
    if myFile1 is None or myFile2 is None:
        return render(request, 'blast.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        myFile1_name = myFile1.name.split(".")[0]  # abcd
        myFile2_name = myFile2.name.split(".")[0]

        myFile1.name = '1.fasta'
        myFile2.name = '2.fasta'
        keep_file(na=myFile1, c=c)
        keep_file(na=myFile2, c=c)

        blast_run_tools(place=c, patten=patten, e=e_value, email=email, score=sco, identity=iden, myFile1=myFile1_name,
                        myFile2=myFile2_name, nucl_or_prot=nucl_or_prot)

        return render(request, '../templates/successful.html')


# diamond
def diamond_run(request):
    patten = request.POST.get('patten', None)
    myFile1 = request.FILES.get('file1', None)
    myFile2 = request.FILES.get('file2', None)
    email = request.POST.get('email', None)
    e_value = str(request.POST.get('e_value', None))
    sco = request.POST.get('score_1', None)
    iden = request.POST.get('identity_1', None)
    nucl_or_prot = "prot" if patten == 'blastp' else "nucl"

    # 判断文件是否为空 如果为空则刷新页面
    if myFile1 is None or myFile2 is None:
        return render(request, 'blast.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        myFile1_name = myFile1.name.split(".")[0]
        myFile2_name = myFile2.name.split(".")[0]

        myFile1.name = '1.fasta'
        myFile2.name = '2.fasta'
        keep_file(na=myFile1, c=c)
        keep_file(na=myFile2, c=c)

        diamond_run_tools(place=c, patten=patten, e=e_value, email=email, score=sco, identity=iden,
                          myFile1=myFile1_name, myFile2=myFile2_name)
        return render(request, '../templates/successful.html')


# sequence_screening
def sequence_screening_run(request):
    email = request.POST.get('email', None)
    protein = request.FILES.get("protein", None)
    gene = request.FILES.get('gene', None)
    gff = request.FILES.get('gff', None)
    id = request.FILES.get('id', None)

    # 判断文件是否为空 如果为空则刷新页面
    if protein is None or gene is None or gff is None or id is None:
        return render(request, 'sequence_screening.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        protein.name = 'protein.fasta'
        gene.name = 'gene.fasta'
        gff.name = 'gff.fasta'
        id.name = 'id.fasta'

        keep_file(na=protein, c=c)
        keep_file(na=gene, c=c)
        keep_file(na=gff, c=c)
        keep_file(na=id, c=c)

        sequence_run(place=c, email=email)
        return render(request, '../templates/successful.html')


# hmmer
def hmmer_run(request):
    email = request.POST.get('email', None)
    hmm_model = request.FILES.get("hmm_model", None)
    hmm_fa = request.FILES.get('hmm_fa', None)
    e_value = request.POST.get('e_value', None)

    if hmm_fa is None or hmm_model is None:
        return render(request, '../templates/hmmer.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        hmm_model.name = "hmm_model.fasta"
        hmm_fa.name = "hmm_fa.fasta"
        keep_file(na=hmm_model, c=c)
        keep_file(na=hmm_fa, c=c)

        hmmer_run_tools(place=c, email=email, e_value=e_value)
        return render(request, '../templates/successful.html')


# pfam
def pfam_run(request):
    email = request.POST.get('email', None)
    myFile1 = request.FILES.get("file1", None)
    structure = request.POST.get('structure', None)

    if myFile1 is None:
        return render(request, '../templates/pfam.html')
    else:
        c = name_set()

        keep_email(e=c, address=email)

        myFile1.name = "1.fasta"
        keep_file(na=myFile1, c=c)

        pfam_run_tools(place=c, structure=structure, email=email)
        # bd_2(place=c, email=email, structure=structure)
        return render(request, '../templates/successful.html')


# MEME
def meme_suit_run(request):
    myFile1 = request.FILES.get("file1", None)
    patten = request.POST.get('patten', None)
    motif = request.POST.get('motif', None)
    motif_num = request.POST.get('motif_num', None)
    minw = request.POST.get('minw', None)
    maxw = request.POST.get('maxw', None)
    email = request.POST.get('email', None)

    # 判断文件是否为空 如果为空则刷新页面
    if myFile1 is None:
        return render(request, '../templates/meme_suit.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        myFile1.name = '1.fasta'
        keep_file(na=myFile1, c=c)
        meme_run_tools(place=c, patten=patten, email=email, motif=motif, motif_num=motif_num, minw=minw, maxw=maxw)

        return render(request, '../templates/successful.html')


# codonw
def codonw_run(request):
    email = request.POST.get('email', None)
    file = request.FILES.get("file", None)

    # 判断文件是否为空 如果为空则刷新页面
    if file is None:
        return render(request, '../templates/codonw.html')
    else:
        from .tools import codonw_run
        c = name_set()
        keep_email(e=c, address=email)

        file_name = file.name.split(".")[0]

        file.name = '1.fasta'
        keep_file(na=file, c=c)

        codonw_run_tools(place=c, email=email, file_name=file_name)
        return render(request, '../templates/successful.html')


# cpgfinder
def cpgfinder_run(request):
    myFile1 = request.FILES.get("file1", None)

    island = request.POST.get('island', None)
    cg_percent = request.POST.get('cg_percent', None)
    gc_ratio = request.POST.get('gc_ratio', None)
    email = request.POST.get('email', None)

    if myFile1 is None:
        return render(request, '../templates/cpgfinder.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        myFile_name = myFile1.name
        keep_file(na=myFile1, c=c)

        cpgfinder_run_tools(place=c, island=island, cg_percent=cg_percent, gc_ratio=gc_ratio, email=email,
                            myFile_name=myFile_name)
        return render(request, '../templates/successful.html')


# dupgen_finder
def dupgen_finder_run(request):
    patten = request.POST.get('patten', None)
    file1 = request.FILES.get("file1", None)
    file2 = request.FILES.get('file2', None)
    file3 = request.FILES.get('file3', None)
    file4 = request.FILES.get('file4', None)
    target = request.POST.get('target', None)
    outgroup = request.POST.get('outgroup', None)
    email = request.POST.get('email', None)

    if file1 is None or file2 is None or file3 is None or file4 is None:
        return render(request, '../templates/dupgen_finder.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        file1.name = f'{target}.gff'
        file2.name = f'{target}.blast'
        file3.name = f'{target}_{outgroup}.gff'
        file4.name = f'{target}_{outgroup}.blast'

        keep_file(na=file1, c=c)
        keep_file(na=file2, c=c)
        keep_file(na=file3, c=c)
        keep_file(na=file4, c=c)

        dupgen_finder_run_tools(email=email, place=c, patten=patten, target=target, outgroup=outgroup)
        return render(request, '../templates/successful.html')


# cd_hit
def cd_hit_run(request):
    file = request.FILES.get("file", None)
    similar = request.POST.get('similar', None)
    patten = request.POST.get('patten', None)
    email = request.POST.get('email', None)

    if file is None:
        return render(request, '../templates/cd_hit.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        file.name = '1.fasta'
        keep_file(na=file, c=c)

        cd_hit_run_tools(place=c, email=email, patten=patten, similar=similar)
        return render(request, '../templates/successful.html')


# sequence_alignment
def sequence_alignment_run(request):
    myFile1 = request.FILES.get("file1", None)
    patten = request.POST.get('patten', None)
    email = request.POST.get('email', None)

    if myFile1 is None:
        return render(request, '../templates/sequence_alignment.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        myFile1.name = "1.fasta"
        keep_file(na=myFile1, c=c)

        sequence_aligment_run_tools(place=c, email=email, patten=patten)
        return render(request, '../templates/successful.html')


# MCScanX
def mcscanx_run(request):
    myFile1 = request.FILES.get("file1", None)
    myFile2 = request.FILES.get("file2", None)
    email = request.POST.get('email', None)

    if myFile1 is None or myFile2 is None:
        return render(request, '../templates/mcscanx.html')
    else:
        c = name_set()
        file_name = myFile1.name.split(".")[0]

        keep_email(e=c, address=email)
        keep_file(na=myFile1, c=c)
        keep_file(na=myFile2, c=c)

        mcscanx_run_tools(place=c, email=email, file_name=file_name)
        return render(request, '../templates/successful.html')


# colinearscan
def colinearscan_run(request):
    blast = request.FILES.get("blast", None)
    gff1 = request.FILES.get("gff1", None)
    gff2 = request.FILES.get("gff2", None)
    e_value = str(request.POST.get('e_value', None))
    score = request.POST.get('score', None)
    hitnum = request.POST.get('hitnum', None)
    pos_order = request.POST.get('pos_order', None)
    cdsvchr = request.POST.get('cdsvchr', None)
    email = request.POST.get('email', None)

    if blast is None or gff1 is None or gff2 is None:
        return render(request, '../templates/colinearscan.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        blast_name = blast.name
        file_name = blast.name.split(".")[0]

        gff1_name = gff1.name
        gff2_name = gff2.name
        keep_file(na=blast, c=c)
        keep_file(na=gff1, c=c)
        keep_file(na=gff2, c=c)
        colinearscan_run_tools(place=c, score=score, e_value=e_value, hitnum=hitnum, email=email, blast_name=blast_name,
                               gff1_name=gff1_name, gff2_name=gff2_name, file_name=file_name, pos_order=pos_order,
                               cdsvchr=cdsvchr)
        return render(request, '../templates/successful.html')


# orthofinder
def orthofinder_run(request):
    myFile1 = request.FILES.getlist("file1", None)
    number = str(request.POST.get('number', None))
    gene_tree = request.POST.get('gene_tree', None)
    blast = request.POST.get('blast', None)
    alignment = request.POST.get('alignment', None)
    tree_way = request.POST.get('tree_way', None)
    email = request.POST.get('email', None)

    if myFile1 is None:
        return render(request, '../templates/orthofinder.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        keep_files(na=myFile1, c=c)

        orthfinder_run_tools(place=c, gene_tree=gene_tree, blast=blast, email=email, alignment=alignment, tree_way=tree_way, number=number)
        return render(request, '../templates/successful.html')


# paraat
def paraat_run(request):
    patten = request.POST.get('patten', None)
    homologs = request.FILES.get("homologs", None)
    cds = request.FILES.get('cds', None)
    pep = request.FILES.get('pep', None)
    email = request.POST.get('email', None)

    if homologs is None or cds is None or pep is None:
        return render(request, '../templates/paraat.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        homologs.name = 'hom.homologs'
        cds.name = 'cds.cds'
        pep.name = 'pep.pep'
        keep_file(na=homologs, c=c)
        keep_file(na=cds, c=c)
        keep_file(na=pep, c=c)

        paraat_run_tools(place=c, email=email, patten=patten)
        return render(request, '../templates/successful.html')


# kaks_calculator
def kaks_calculator_run(request):
    patten = request.POST.get('patten', None)
    kaks = request.FILES.get("kaks", None)
    email = request.POST.get('email', None)

    if kaks is None:
        return render(request, '../templates/kaks_calculator.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        kaks.name = 'kaks.axt'
        keep_file(na=kaks, c=c)

        kaks_calculator_run_tools(place=c, email=email, patten=patten)
        return render(request, '../templates/successful.html')


# fasttree
def fasttree_run(request):
    patten = request.POST.get('patten', None)
    myFile1 = request.FILES.get("file1", None)
    email = request.POST.get('email', None)

    if myFile1 is None:
        return render(request, '../templates/fasttree.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)

        myFile1.name = "1.phy"
        keep_file(na=myFile1, c=c)

        fasttree_run_tools(place=c, email=email, patten=patten)
        return render(request, '../templates/successful.html')


# iqtree
def iq_tree_run(request):
    email = request.POST.get('email', None)
    file = request.FILES.get("file", None)
    patten = request.POST.get('patten', None)
    bs_num = request.POST.get('bs_num', None)

    if file is None:
        return render(request, '../templates/iq_tree.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        file.name = 'result.phy'
        keep_file(na=file, c=c)

        iqtree_run_tools(place=c, email=email, patten=patten, bs_num=bs_num)
        return render(request, '../templates/successful.html')


# phyml
def phyml_run(request):
    myFile1 = request.FILES.get("file", None)
    patten = request.POST.get('patten', None)
    small = request.POST.get('small', None)
    model = request.POST.get('model', None)
    frequency = request.POST.get('frequency', None)
    nice = request.POST.get('nice', None)
    num = request.POST.get('num', None)
    email = request.POST.get('email', None)

    if myFile1 is None:
        return render(request, '../templates/phyml.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        myFile1.name = "1.phy"
        keep_file(na=myFile1, c=c)

        phyml_run_tools(place=c, email=email, patten=patten, model=model, frequency=frequency, nice=nice, num=num, small=small)
        return render(request, '../templates/successful.html')


# paml
def paml_run(request):
    aln = request.FILES.get("aln", None)
    nwk = request.FILES.get("nwk", None)
    text = request.POST.get('patten', None)
    text2 = request.POST.get('small', None)
    email = request.POST.get('email', None)

    if aln is None or nwk is None:
        return render(request, '../templates/paml.html')
    else:
        c = name_set()
        keep_email(e=c, address=email)
        aln.name = "1.pml"
        nwk.name = "1.tree"
        keep_file(na=aln, c=c)
        keep_file(na=nwk, c=c)

        paml_run_one(patten=text, small=text2, path=c, place=c, email=email)
        paml_run_two(place=c, email=email, patten=text, small=text2)
        return render(request, '../templates/successful.html')






