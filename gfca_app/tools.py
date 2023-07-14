#!/usr/bin/env python

import time
import random
import os
import re
import pymysql
from Bio import SeqIO
from GFCA.email import youjian

# 搭载到服务器时，改一下各个软件调用的线程数

path_get = os.getcwd()


# 在文件夹内 创建文件夹
def keep_in_file(na, c, fg=None):
    global path_get
    cmd = '''
        cd %s/file_keep/%s
        mkdir all_file
    ''' % (path_get, c)
    os.system(cmd)

    with open('%s' % os.path.join(f'{path_get}/file_keep/{c}/all_file/', na.name), 'wb') as f:
        for chunk1 in na.chunks():
            f.write(chunk1)
    f.close()


# 保存单个文件
def keep_file(na, c, fg=None):
    if fg is None:
        pass
    else:
        nam = na.name.split(".")
        na.name = nam[0]

    # 被我保存到本地啦，保存到哪里随意啦
    with open(f"{path_get}/file_keep/{c}/{na.name}", 'wb+') as f:
        for chunk1 in na.chunks():
            f.write(chunk1)
    f.close()


# 保存多个文件
def keep_files(na, c, fg=None, ret_files=None):
    global path_get
    cmd = '''
            cd %s/file_keep/%s
            mkdir orthomcl
            ''' % (path_get, c)
    os.system(cmd)

    for file in na:
        with open('%s' % os.path.join(f'{path_get}/file_keep/{c}/orthomcl/', file.name), 'wb') as f:
            for i in file.chunks():
                f.write(i)


# 解压文件
def zip_dp(e, id):
    global path_get

    cmd = '''
        cd %s/file_keep/%s
        gzip -d %s.fasta.gz
        ''' % (path_get, e, id)

    os.system(cmd)


# 此模块用于检查 文件是否存在 （b，c，d 参数未补全）
#     e ：文件夹名字
def name_check(e):
    """
    此模块用于检查 文件是否存在 （b，c，d 参数未补全）
    e ：文件夹名字
    """
    global path_get

    a = os.path.isfile('/home/a123123/PycharmProjects/untitled2/file_keep/%s/1.FASTA' % e)
    b = os.path.isfile('/home/a123123/PycharmProjects/untitled2/file_keep/%s/email.txt' % e)
    c = os.path.isfile('/home/a123123/PycharmProjects/untitled2/file_keep/%s/email.txt' % e)
    d = os.path.isfile('/home/a123123/PycharmProjects/untitled2/file_keep/%s/email.txt' % e)
    n = ''
    if b:
        n = '1.FASTA'
    elif a:
        n = 'd'
    elif a:
        n = 'd'
    elif a:
        n = 'd'
    return n


# 设定文件夹
def name_set():
    """
    此方法用于 文件保存文件夹的建立
    在 views 模块中被调用
    """
    global path_get

    a = time.asctime(time.localtime(time.time()))
    b = a[-15:] + str(random.randint(0, 1000))
    c = b.replace(" ", "")
    c = c.replace(":", "0")
    cmd = '''
        cd %s/file_keep
        mkdir %s
        ''' % (path_get, c)
    os.system(cmd)
    return c


# 保存邮件地址
def keep_email(e, address):
    global path_get

    with open(f'{path_get}/file_keep/{e}/email.txt', 'w') as f:
        f.write(address)


# index 保存邮件地址
def keep_email_index(e, name, email, message):
    global path_get

    with open(f'{path_get}/file_keep/{e}/advise.txt', 'w+') as f:
        f.write(name + '\n')
        f.write(email + '\n')
        f.write(message + '\n')

    g = f'{path_get}/file_keep/{e}/advise.txt'
    youjian(a="1457634627@qq.com", b=g, place=e)


# 创建conf配置文件
def text_create(name, msg, place):
    desktop_path = f"{path_get}/file_keep/{place}/"  # 新创建的txt文件的存放路径
    full_path = desktop_path + name + '.conf'  # 也可以创建一个.doc的word文档
    file = open(full_path, 'w')
    file.write(msg)  # msg也就是下面的Hello world!
    file.close()


def read_email(place):
    global path_get
    with open(f'{path_get}/file_keep/{place}/email.txt', 'r+') as f:
        return f.read()


def rawdata_run(place, email, name):
    global path_get

    cmd = f"""
        cd {path_get}/file_keep/{place}
        python ../programs/github_gff_deal.py -p {name} dna_top.fasta gff3.fasta
        sed -i -e 's/gene://' -e 's/transcript://' {name}.gff
    """
    os.system(cmd)

# set +o posix  在运行前加上这个。。。非常重要。。半天时间才试出来。。。哭了

    cmd2 = f"""
        cd {path_get}/file_keep/{place}
        set +o posix
        seqkit grep -f <(cut -f 7 {name}.gff) cds.fasta |
        seqkit seq --id-regexp "^(.*?)\\.\\d" -i > {name}.cds.fa
        seqkit grep -f <(cut -f 7 {name}.gff) pep.fasta |
        seqkit seq --id-regexp "^(.*?)\\.\\d" -i > {name}.pep.fa
        zip -r result.zip {name}*
    """
    os.system(cmd2)

    g = f'{path_get}/file_keep/{place}/result.zip'
    youjian(a=email, b=g, place=place)


# tools ---blast_dotplot---
def dotplot_run(place, email, evalue, score, heng_species, shu_species, patten):
    global path_get

    write_txt = (str("[dotplot]") + "\n" + str(f"blast = diamond.fasta") + "\n" +
                 str("gff1 =  11.gff") + "\n" + str("gff2 =  22.gff") + "\n" + str("lens1 = 11.len") + "\n" +
                 str("lens2 = 22.len") + "\n" + str(f"genome1_name =  {shu_species}") + "\n" + str(f"genome2_name =  {heng_species}") + "\n" +
                 str("multiple  = 1") + "\n" + str(f"score = {score}") + "\n" + str(f"evalue = {evalue}") + "\n" +
                 str("repeat_number = 10") + "\n" + str("position = order") + "\n" + str(
                "blast_reverse = false") + "\n" +
                 str("ancestor_left = none") + "\n" + str("ancestor_top = none") + "\n" + str(
                "markersize = 0.5") + "\n" +
                 str("figsize = 10,10") + "\n" + str(f"savefig = result.{patten}") + "\n" + "\n"
                 )
    text_create('ath', str(write_txt), place)

    cmd = f"""
        cd {path_get}/file_keep/{place}
        wgdi -d ath.conf
"""
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.{patten}'
    youjian(a=email, b=g, place=place)


# tools ---wgdi_共线性分析---
def wgdimcs_run(place, email, evalue, score):
    global path_get

    write_txt = (str("[collinearity]") + "\n" + str("gff1 =  11.gff") + "\n" + str("gff2 =  22.gff") + "\n" +
                 str("lens1 = 11.len") + "\n" + str("lens2 = 22.len") + "\n" + str(f"blast = diamond.fasta") + "\n" +
                 str("blast_reverse = false") + "\n" + str("multiple  = 1") + "\n" + str("process = 8") + "\n" +
                 str(f"evalue = {evalue}") + "\n" + str(f"score = {score}") + "\n" + str("grading = 50,40,25") + "\n" +
                 str("mg = 40,40") + "\n" + str("pvalue = 0.2") + "\n" + str("repeat_number = 10") + "\n" +
                 str("positon = order") + "\n" + str("savefile = result.collinearity.txt") + "\n" + "\n"
                 )
    text_create('ath', str(write_txt), place)

    cmd = f"""
        cd {path_get}/file_keep/{place}
        wgdi -icl ath.conf
"""
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.collinearity.txt'
    youjian(a=email, b=g, place=place)


# tools ---wgdi_ks分析---
def wgdikaks_run(place, email):
    global path_get

    write_txt = (str("[ks]") + "\n" + str("cds_file = cds.fasta") + "\n" + str("pep_file = 	pep.fasta") + "\n" +
                 str("align_software = muscle") + "\n" + str("pairs_file = colinearity.txt") + "\n" +
                 str("ks_file = ks_result.txt") + "\n"
                 )
    text_create('ath', str(write_txt), place)

    cmd = f"""
        cd {path_get}/file_keep/{place}
        wgdi -ks ath.conf
"""
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/ks_result.txt'
    youjian(a=email, b=g, place=place)


# tools ---绘制kaks点图---
def wgdikaks_dot_run(place, email, heng_species, shu_species):
    global path_get

    write_txt = (str("[blockks]") + "\n" + str(f"lens1 = 11.len") + "\n" +
                 str("lens2 = 22.len") + "\n" + str(f"genome1_name =  {heng_species}") + "\n" +
                 str(f"genome2_name =  {shu_species}") + "\n" + str("blockinfo =  blockinfo.csv") + "\n" +
                 str("pvalue = 0.05") + "\n" + str("tandem = true") + "\n" +
                 str("tandem_length = 200") + "\n" + str("markersize = 1") + "\n" + str("area = 0,3") + "\n" +
                 str("block_length =  5") + "\n" + str("figsize = 8,8") + "\n" +
                 str("savefig = ks.dotplot.pdf") + "\n" + "\n")
    text_create('ath', str(write_txt), place)

    cmd = f"""
        cd {path_get}/file_keep/{place}
        wgdi -bk ath.conf
"""
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/ks.dotplot.pdf'
    youjian(a=email, b=g, place=place)


# file_merge_run
def file_merge_run(place, email, file_name):
    global path_get

    cmd = f"""
            cd '{path_get}/file_keep/{place}'
            cat orthomcl/*.{file_name} >> result.{file_name}
        """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.{file_name}'
    youjian(a=email, b=g, place=place)


# find_replace_run
def find_replace_run(place, email, old_name, new_name):
    global path_get

    cmd = f"""
            cd '{path_get}/file_keep/{place}'
            sed -e 's/{old_name}/{new_name}/g' upload.txt > result.txt
        """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.txt'
    youjian(a=email, b=g, place=place)


# blast
def blast_run_tools(place, patten, e, email, score, identity, myFile1, myFile2):
    global path_get

    cmd = f"""
    source /etc/profile
    cd '{path_get}/file_keep/{place}'
    makeblastdb -in 1.fasta -dbtype prot -parse_seqids -out dbname
    {patten} -query 2.fasta -db dbname -outfmt 6 -evalue {e} -out result_2.fasta
"""
    os.system(cmd)

    # 初筛
    blast_file = open(f'{path_get}/file_keep/{place}/{myFile1}_{myFile2}.fasta', 'w')
    sequence_file = open(f'{path_get}/file_keep/{place}/sequence.txt', 'w')
    pbl_fasta = open(f'{path_get}/file_keep/{place}/result_2.fasta', 'r')
    list_id = []
    for line in pbl_fasta.readlines():
        ind = line.split()[2]
        sco = line.split()[11]
        # sequence = line.split()[1]
        # list_id.append(sequence)
        if float(ind) >= float(identity):
            if float(sco) >= float(score):
                blast_file.write(line)
                sequence = line.split()[1]
                list_id.append(sequence)

    # 去重
    lst2 = list(set(list_id))
    for line in lst2:
        sequence_file.write('>' + line + '\n')

    pbl_fasta.close()
    blast_file.close()
    sequence_file.close()

    cmd3 = f"""
        cd {path_get}/file_keep/{place}
        zip -r result.zip {myFile1}_{myFile2}.fasta sequence.txt
    """
    os.system(cmd3)

    g = f'{path_get}/file_keep/{place}/result.zip'
    youjian(a=email, b=g, place=place)


# blast---diamond
def diamond(place, patten, e, email, score, identity, myFile1, myFile2):
    global path_get

    cmd = f"""
    source /etc/profile
    cd {path_get}/file_keep/{place}
    diamond makedb --in 1.fasta -d nr
    diamond {patten} -d nr -q 2.fasta -o result_2.fasta -e {e} -b 2.0 --index-chunks 1
"""
    os.system(cmd)

    # 初筛
    blast_file = open(f'{path_get}/file_keep/{place}/{myFile1}_{myFile2}.fasta', 'w')
    sequence_file = open(f'{path_get}/file_keep/{place}/sequence.txt', 'w')
    pbl_fasta = open(f'{path_get}/file_keep/{place}/result_2.fasta', 'r')
    list_id = []
    for line in pbl_fasta.readlines():
        ind = line.split()[2]
        sco = line.split()[11]
        # sequence = line.split()[1]
        # list_id.append(sequence)
        if float(ind) >= float(identity):
            if float(sco) >= float(score):
                blast_file.write(line)
                sequence = line.split()[1]
                list_id.append(sequence)

    # 去重
    lst2 = list(set(list_id))
    for line in lst2:
        sequence_file.write('>' + line + '\n')

    pbl_fasta.close()
    blast_file.close()
    sequence_file.close()

    cmd3 = f"""
    cd {path_get}/file_keep/{place}
    zip -r result.zip {myFile1}_{myFile2}.fasta sequence.txt
"""
    os.system(cmd3)

    g = f'{path_get}/file_keep/{place}/result.zip'
    youjian(a=email, b=g, place=place)


# 序列筛选---sequence
def sequence_run(place):
    global path_get

    # 蛋白文件
    seq_file_1 = {}
    protein = open(f'{path_get}/file_keep/{place}/protein.fasta', 'r')
    for line_1 in protein.readlines():
        line_pro = line_1.strip()
        if line_pro[0] == '>':
            seq_id = line_pro.split()[0]
            seq_file_1[seq_id] = ''
        else:
            seq_file_1[seq_id] += line_pro
    # cds核酸文件
    seq_file_2 = {}
    gene = open(f'{path_get}/file_keep/{place}/gene.fasta', 'r')
    for line_2 in gene.readlines():
        line_gene = line_2.split()[0]
        if line_gene[0] == '>':
            seq_id = line_gene.split()[0]
            seq_file_2[seq_id] = ''
        else:
            seq_file_2[seq_id] += line_gene

    sequence_id = open(f'{path_get}/file_keep/{place}/id.fasta', 'r')
    gff_list = []
    pro_list = []
    cds_list = []
    for lst_id in sequence_id.readlines():
        id = lst_id.strip()
        # gff
        input_fasta = open(f'{path_get}/file_keep/{place}/gff.fasta', 'r')
        for line2 in input_fasta.readlines():
            line = '>' + line2.split()[1]
            if id == line:
                gff_list.append(line2)
        # 蛋白
        for pro_line in seq_file_1.keys():
            if id == pro_line:
                pro_list.append(pro_line + '\n' + seq_file_1[pro_line] + '\n')
        # 核酸
        for cds_line in seq_file_2.keys():
            if id == cds_line:
                cds_list.append(cds_line + '\n' + seq_file_2[cds_line] + '\n')

    str = ''
    new_gff = open(f'{path_get}/file_keep/{place}/newgff.txt', 'w')
    new_gff.write(str.join(gff_list))

    new_pro = open(f'{path_get}/file_keep/{place}/newpro.txt', 'w')
    new_pro.write(str.join(pro_list))

    new_cds = open(f'{path_get}/file_keep/{place}/newcds.txt', 'w')
    new_cds.write(str.join(cds_list))

    # 下面的程序会先执行
    # cmd3 = """
    #         cd '%s/file_keep/%s'
    #         zip -r result.zip newpro.txt newgff.txt newcds.txt
    #     """ % (path_get, place)
    # os.system(cmd3)
    #
    # g = f'{path_get}/file_keep/{place}/result.zip'
    # youjian(a=email, b=g, place=place)


# sequence-22
def sequence_2(place, email):
    global path_get
    cmd3 = f"""
    cd {path_get}/file_keep/{place}
    zip -r result.zip newpro.txt newgff.txt newcds.txt
"""
    os.system(cmd3)

    g = f'{path_get}/file_keep/{place}/result.zip'
    youjian(a=email, b=g, place=place)


# pfam
def bd(place, structure):
    # 张岚老师data文件夹使用命令
    # pfam_scan.pl -fasta 1.fasta -dir ~/../../../../www/tools/PfamScan/pfamdata -outfile pfam.txt
    global path_get
    cmd = f"""
    source /etc/profile
    cd '{path_get}/file_keep/{place}'
    pfam_scan.pl -fasta 1.fasta -dir ~/../../../../www/tools/PfamScan/pfamdata -outfile pfam.txt
"""
    os.system(cmd)
    # 根据用户设置的结构域进行筛选---有问题
    # hh = structure
    # pfam_file = open(f'{path_get}/file_keep/{place}/pfam.txt', 'r')
    # pfam_new = open(f'{path_get}/file_keep/{place}/pfam_new.txt', 'w')
    # for pfam_line in pfam_file.readlines():
    #     if pfam_line[0] != '#' and pfam_line[0] != '\n':
    #         struture = pfam_line.split()[6]
    #         result = hh.find(struture)
    #         if result != -1:
    #             pfam_new.write(pfam_line)
    # pfam_file.close()
    # pfam_new.close()
    #
    # g = '%s/file_keep/%s/pfam_new.txt' % (path_get, place)
    # youjian(a=email, b=g, place=place)


def bd_2(place, email, structure):
    # hh = structure
    # pfam_file = open(f'{path_get}/file_keep/{place}/pfam.txt', 'r')
    # pfam_new = open(f'{path_get}/file_keep/{place}/pfam_new.txt', 'w')
    # for pfam_line in pfam_file.readlines():
    #     if pfam_line[0] != '#' and pfam_line[0] != '\n':
    #         struture = pfam_line.split()[6]
    #         result = hh.find(struture)
    #         if result != -1:
    #             pfam_new.write(pfam_line)

    # cmd3 = """
    #         cd '%s/file_keep/%s'
    #         zip -r result.zip pfam_new.txt
    # """ % (path_get, place)
    # os.system(cmd3)

    g = '%s/file_keep/%s/pfam.txt' % (path_get, place)
    youjian(a=email, b=g, place=place)


# meme
def meme_run(place, email, patten, motif, motif_num, minw, maxw):
    global path_get

    cmd = f"""
    source /etc/profile
    cd {path_get}/file_keep/{place}
    meme 1.fasta -{patten} -oc . -mod {motif} -nmotifs {motif_num} -minw {minw} -maxw {maxw}
    zip -r result.zip logo*.png meme.*
"""
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.zip'
    youjian(a=email, b=g, place=place)


# paml_deal
def paml_sy(patten, small, path):
    global path_get
    file_data = ""
    codeml_file = open(f'{path_get}/file_keep/paml/{patten}/{small}/codeml.ctl', 'r')
    for line in codeml_file.readlines():
        hh = re.sub(r'\d{9,}', f'{path}', line)
        file_data += hh
    co_new = open(f'{path_get}/file_keep/paml/{patten}/{small}/codeml.ctl', 'w')
    co_new.write(file_data)


# paml
def sy(place, email, patten, small):
    cmd = f"""
    cd {path_get}/file_keep/{place}
    codeml ../paml/{patten}/{small}/codeml.ctl
"""
    os.system(cmd)

    g = '%s/file_keep/%s/resoult.txt' % (path_get, place)
    youjian(a=email, b=g, place=place)


# mcscanx
def gxx(place, email, file_name):
    cmd = f"""
    cd {path_get}/file_keep/{place}
    MCScanX {file_name}
"""
    os.system(cmd)

    cmd2 = f"""
    cd {path_get}/file_keep/{place}
    duplicate_gene_classifier {file_name} > gene_classifier.txt
    zip -r result.zip {file_name}.html {file_name}.collinearity {file_name}.tandem gene_classifier.txt
"""
    os.system(cmd2)

    g = f"{path_get}/file_keep/{place}/result.zip"
    youjian(a=email, b=g, place=place)


# 多序列比对
def jh(place, email, patten):
    if patten == "clustalw2":
        cmd2 = f'{path_get}/file_keep/{place}'
        os.chdir(cmd2)
        cmd = "clustalw2 1.fasta"
        os.system(cmd)
        g = '%s/file_keep/%s/1.aln' % (path_get, place)
        youjian(a=email, b=g, place=place)
    elif patten == "muscle":
        cmd2 = f'{path_get}/file_keep/{place}'
        os.chdir(cmd2)
        cmd = "muscle3.8.31_i86linux64 -in 1.fasta -clwout 1.aln"
        os.system(cmd)
        g = '%s/file_keep/%s/1.aln' % (path_get, place)
        youjian(a=email, b=g, place=place)
    else:
        cmd2 = f'{path_get}/file_keep/{place}'
        os.chdir(cmd2)
        cmd = "mafft --auto 1.fasta > mafft.fasta"
        os.system(cmd)
        g = '%s/file_keep/%s/mafft.fasta' % (path_get, place)
        youjian(a=email, b=g, place=place)


# hmmer
def hmm(place, email, e_value):
    """
    hmm_model.fasta  多序列比对文件
    hmm_fa.fasta  未比对的蛋白文件
    db_model.hmm  建立模型的结果文件
    result.out  hmmer最终结果文件
    """
    cmd = f"""
        cd {path_get}/file_keep/{place}
        hmmbuild db_model.hmm hmm_model.fasta
        hmmsearch -o result.out --noali -E {e_value} db_model.hmm hmm_fa.fasta
    """
    os.system(cmd)

    g = '%s/file_keep/%s/result.out' % (path_get, place)
    youjian(a=email, b=g, place=place)


# fasttree
def xt(place, email, patten):
    if patten == "JTT+CAT":
        cmd2 = f'{path_get}/file_keep/{place}'
        os.chdir(cmd2)
        cmd = 'FastTree 1.phy > 1.nwk'
        os.system(cmd)
    else:
        cmd2 = f'{path_get}/file_keep/{place}'
        os.chdir(cmd2)
        cmd = "FastTree -gtr -nt  1.phy > 1.nwk"
        os.system(cmd)

    g = '%s/file_keep/%s/1.nwk' % (path_get, place)
    youjian(a=email, b=g, place=place)


# phyml
def phyml_run(email, place, patten, nice, model, frequency, num, small):
    global path_get

    cmd = f"""
        cd {path_get}/file_keep/{place}
        phyml -i 1.phy -d {patten} -b {num} -m {small} -f {frequency} -v e -a e -o {nice}
        zip -r result.zip 1.phy_*
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.zip'
    youjian(a=email, b=g, place=place)


# heat_maps
def map1(place, email, gbd, gds, dbds, jh, bzh):
    dir = {
        "line": "-r NA",
        "row": "-c NA",
        "NULL": ""
    }

    cmd = f"""
        cd {path_get}/file_keep/{place}
        Rscript ../programs/3.R -in heat.txt {dir[jh]} -s {bzh} -l {gbd} -m {gds} -g {dbds}
    """

    os.system(cmd)
    g = '%s/file_keep/%s/Rplots.pdf' % (path_get, place)
    youjian(a=email, b=g, place=place)


# paraAT_run
def para_run(place, email, patten):
    cmd = f"""
            cd '{path_get}/file_keep/{place}'
            ParaAT.pl -h hom.homologs -n cds.cds -a pep.pep -p ../programs/proc -o result -f {patten}
            cat result/*.{patten} > result.{patten}
        """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.{patten}'
    youjian(a=email, b=g, place=place)


# kaks_cal_run
def kaks_cal_run(place, email, patten):
    cmd = f"""
        cd {path_get}/file_keep/{place}
        KaKs_Calculator -i kaks.axt -o result_kaks.txt -m {patten}
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result_kaks.txt'
    youjian(a=email, b=g, place=place)


# iqtree
def iqtree_run(place, email, patten, bs_num):
    cmd = f"""
            cd {path_get}/file_keep/{place}
            iqtree2 -s result.phy -m {patten} -bb {bs_num} -bnni -redo
        """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.phy.treefile'
    youjian(a=email, b=g, place=place)


# extr_row_run
def extr_row_run(place, email, row_num, row_name, row_last):
    cmd = f"""
            cd {path_get}/file_keep/{place}
            cut -f {row_num} {row_name}.{row_last} > {row_name}.new.{row_last}
        """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/{row_name}.new.{row_last}'
    youjian(a=email, b=g, place=place)


# formatTofasta
def format_fasta_run(place, email, file_name, patten, patten_end, file_last):
    SeqIO.convert(f"{path_get}/file_keep/{place}/{file_name}.{file_last}", f"{patten}",
                  f"{path_get}/file_keep/{place}/{file_name}.{patten_end}", f"{patten_end}")

    g = f'{path_get}/file_keep/{place}/{file_name}.{patten_end}'
    youjian(a=email, b=g, place=place)


# codonw_run
def codonw_run(place, email, file_name):
    cmd = f"""
        cd {path_get}/file_keep/{place}
        codonw 1.fasta -all_indices -nomenu {file_name}.out {file_name}.blk
        zip -r codonw_result.zip {file_name}.out {file_name}.blk
    """
    os.system(cmd)

    g = '%s/file_keep/%s/codonw_result.zip' % (path_get, place)
    youjian(a=email, b=g, place=place)


# fastqc_run
def fastqc_run(place, email):
    cmd = f"""
        cd {path_get}/file_keep/{place}

        fastqc -o output dir -f fastq|bam|sam seqfile1 .. seqfileN

        -o：用来指定输出文件的所在目录，注意是不能自动新建目录的。
        输出的结果：是.zip文件，默认自动解压缩，命令里加上--noextract则不解压缩。
        -f：用来强制指定输入文件格式，默认会自动检测。
        -c：用来指定一个contaminant文件，fastqc会把overrepresented sequences往这个contaminant文件里搜索。
        contaminant文件的格式是"Name\tSequences"，#开头的行是注释。

        
        fastqc -t 1 -o {path_get}/file_keep/{place} {path_get}/file_keep/{place}/orthomcl/*.*
    """
    os.system(cmd)

    cmd3 = f'''
        cd {path_get}/file_keep/{place}
        zip -r result.zip *.fq_fastqc.html *.fq_fastqc.zip
    '''
    os.system(cmd3)

    g = '%s/file_keep/%s/result.zip' % (path_get, place)
    youjian(a=email, b=g, place=place)


# cpgfinder_run
def cpgfinder_run(place, email, island, cg_percent, gc_ratio, myFile_name):
    cmd = f"""
        cd {path_get}/file_keep/{place}
        cpgfinder -f {myFile_name} -l {island} -p {cg_percent} -r {gc_ratio} > cpgfinder_result.txt
    """
    os.system(cmd)

    g = '%s/file_keep/%s/cpgfinder_result.txt' % (path_get, place)
    youjian(a=email, b=g, place=place)


# 去重
def quchong_run(place, email):
    id_list = []
    new_file = open(f'{path_get}/file_keep/{place}/new.fasta', 'w')
    for line in SeqIO.parse(f'{path_get}/file_keep/{place}/1.fasta', 'fasta'):
        if line.id not in id_list:
            id_list.append(line.id)
            new_file.write(">" + str(line.id) + "\n" + str(line.seq) + "\n")
    new_file.close()

    g = '%s/file_keep/%s/new.fasta' % (path_get, place)
    youjian(a=email, b=g, place=place)


# orthofinder
def orthfinder_run(place, email, gene_tree, blast, alignment, tree_way, number):
    # -M    基因树推断方法（默认为dendroblast）可选：dendroblast ，msa
    # -S    序列比对使用的程序 （默认为Blast）可选：blast, mmseqs, blast_gz, diamond（推荐使用diamond，比对速度快）
    # # -A    多序列联配方式，该选项仅当 - M    msa    选项时才有效（默认为mafft）可选：muscle, mafft
    # # -T    建树方式，该选项仅当 - M    msa    选项时才有效 （默认为fasttree）可选：iqtree, raxml - ng, fasttree, raxml
    # -I    设定MCL的膨胀系数（默认为1.5）
    if gene_tree == str("msa"):
        cmd = f"""
            cd {path_get}/file_keep/{place}
            orthofinder -f orthomcl -M {gene_tree} -S {blast} -A {alignment} -T {tree_way} -I {number} 
        """
        os.system(cmd)

        cmd2 = f"""
            cd {path_get}/file_keep/{place}/orthomcl/
            zip -r orthofinder_result.zip OrthoFinder/
        """
        os.system(cmd2)

        g = '%s/file_keep/%s/orthomcl/orthofinder_result.zip' % (path_get, place)
        youjian(a=email, b=g, place=place)
    else:
        cmd = f"""
            cd {path_get}/file_keep/{place}
            orthofinder -f orthomcl/ -M {gene_tree} -S {blast} -I {number} 
        """
        os.system(cmd)

        cmd2 = f"""
            cd {path_get}/file_keep/{place}/orthomcl/
            zip -r orthofinder_result.zip OrthoFinder/
        """
        os.system(cmd2)

        g = f'{path_get}/file_keep/{place}/orthomcl/orthofinder_result.zip'
        youjian(a=email, b=g, place=place)


# cd_hit_run 蛋白质序列或核酸序列聚类的工具
def cd_hit_run(place, email, similar, patten):
    """
    cd-hit -i db -o db90 -c 0.9 -n 5 -M 16000 -d 0 -T 8
        -i 输入文件，fasta格式的序列
        -o 输出文件路径和名字
        -c 相似性（clustering threshold），0.9表示相似性大于等于90%的为一类
        -n 两两序列进行序列比对时选择的 word size
        -d 0表示使用 fasta 标题中第一个空格前的字段作为序列名字
        -M 16000，16GB RAM
        -T 使用的线程数
        Choose of word size:
        -n 5 for thresholds 0.7 ~ 1.0
        -n 4 for thresholds 0.6 ~ 0.7
        -n 3 for thresholds 0.5 ~ 0.6
        -n 2 for thresholds 0.4 ~ 0.5
        -aL：控制代表序列比对严格程度的参数，默认为0，若设为0.8则表示比对区间要占到代表（长）序列的80%。
        -AL：控制代表序列比对严格程度的参数，默认为99999999，若设为40则表示代表序列的非比对区间要短于40bp。
        -aS：控制短序列比对严格程度的参数，默认为0，若设为0.8则表示比对区间要占到短序列的80%。
        -AS：控制短序列比对严格程度的参数，默认为99999999，若设为40则表示短序列的非比对区间要短于40bp。
    """
    cmd = f"""
        cd {path_get}/file_keep/{place}
        cd-hit -i 1.fasta -o cd_hit_result.fasta -c {similar} -n {patten} -d 0 -T 1
        zip -r cd_hit_result.zip cd_hit_result.fasta.clstr cd_hit_result.fasta
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/cd_hit_result.zip'
    youjian(a=email, b=g, place=place)


# colinearscan_run
def colinearscan_run(place, email, score, e_value, hitnum, pos_order, cdsvchr, blast_name, gff1_name, gff2_name,
                     file_name):
    global path_get
    # perl ../programs/run_colinearscan.pl ZM-OS.filter.blast 1e-5 0 30 ZM.gff OS.gff order n ZM_OS

    cmd = f"""
        cd '{path_get}/file_keep/{place}'
        perl ../programs/run_colinearscan.pl {blast_name} {e_value} {score} {hitnum} {gff1_name} {gff2_name} {pos_order} {cdsvchr} {file_name}
    """
    os.system(cmd)

    cmd2 = f"""
        cd '{path_get}/file_keep/{place}'
        zip -r colinearscan_result.zip block pair
    """
    os.system(cmd2)

    g = f'{path_get}/file_keep/{place}/colinearscan_result.zip'
    youjian(a=email, b=g, place=place)


# circos
def circos_run1(place, email, colour1, colour2, colour3, colour4, colour5, name, ran_long_name,
                gff_file_name, kaks_name):
    global path_get
    # python3.6 ../programs/circos_tab_none.py OS
    cmd = f"""
cd {path_get}/file_keep/{place}
python3.6 ../programs/circos_tab_none.py {name} {ran_long_name} {gff_file_name} {kaks_name} {colour1} {colour2} {colour3} {colour4} {colour5}
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/OS.genefam.4-10.proximal.211111.pdf'
    youjian(a=email, b=g, place=place)


# circos
def circos_run2(place, email, in_colour, tab_name, colour1, colour2, colour3, colour4, colour5,
                name, ran_long_name, gff_file_name, kaks_name):
    global path_get
    # python3.6 ../programs/circos_tab_have.py PV
    # python3.6 ../programs/circos_tab_have.py PV pv.lens pv.new.gff pv.ks.txt Pheseolus_vulgaris.beta.gene_pair.tab.txt red green yellow blue red pink

    cmd = f"""
cd {path_get}/file_keep/{place}
python3.6 ../programs/circos_tab_have.py {name} {ran_long_name} {gff_file_name} {kaks_name} {tab_name} {in_colour} {colour1} {colour2} {colour3} {colour4} {colour5}
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/{name}_{ran_long_name}_{gff_file_name}_{kaks_name}_{tab_name}_{in_colour}_{colour1}_{colour2}_{colour3}_{colour4}_{colour5}.genefam.beta.gene.pdf'
    youjian(a=email, b=g, place=place)


# 用户名
user_data = 'root'
# 密码
pd_data = '131477'


class Orthomcl(object):
    def __init__(self):
        global user_data
        global pd_data

        self.E = 1e-10
        self.name = re.sub(r'[ :]*', '', time.asctime(time.localtime(time.time())))[-15:]
        self.percentMatchCutoff = 50
        self.Inflation = 1.5

    def databasrs_set(self):
        global user_data
        global pd_data
        conn = pymysql.connect(host='localhost', user=user_data, password=pd_data, charset='utf8')
        cursor = conn.cursor()
        sql = f'create database {self.name}'
        cursor.execute(sql)
        cursor.close()
        conn.close()
        # print(f'{self.name} 数据库建立完成')

    def del_database(self):
        global user_data
        global pd_data
        conn = pymysql.connect(host='localhost', user=user_data, password=pd_data, charset='utf8')
        cursor = conn.cursor()
        sql = f'drop database {self.name}'
        cursor.execute(sql)
        cursor.close()
        conn.close()

    def handle_fasta(self):
        dir = os.getcwd()
        fl_list = os.listdir(dir)
        config = f'''# this config assumes a mysql database named 'orthomcl'.  adjust according
# to your situation.
dbVendor=mysql 
dbConnectString=dbi:mysql:{self.name}
dbLogin={user_data}
dbPassword={pd_data}
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff={self.percentMatchCutoff}
evalueExponentCutoff={self.E}
oracleIndexTblSpc=NONE'''

        with open('orthomcl.config.template', 'w+') as f:
            f.write(config)
        os.mkdir('compliantFasta')

        for fl_name in fl_list:
            if fl_name[-2:] == 'py':
                continue
            if '.' in fl_name[:2]:
                set_name = fl_name[:1] + 'oo'
            elif '.' in fl_name[:3]:
                set_name = fl_name[:2] + 'o'
            else:
                set_name = fl_name[:3]
            st = f'orthomclAdjustFasta {set_name} {fl_name} 1'
            os.system(st)
            os.system(f'mv ./{set_name}.fasta  compliantFasta/')
            # print('格式化完成')

    def handle_data(self):
        cmd = [
            'orthomclInstallSchema orthomcl.config.template orthomcl.config.log',
            'orthomclFilterFasta compliantFasta/ 10 20',
            'makeblastdb -in goodProteins.fasta -dbtype prot -out orthomc',
            'blastp -db orthomc -query goodProteins.fasta -out orthomcl.blastout -evalue %s -outfmt 7 -num_threads 24 -num_alignments 5' % self.E,
            'grep -P "^[^#]" orthomcl.blastout > blastresult',
            'orthomclBlastParser blastresult compliantFasta/> similarSequences.txt',
            'orthomclLoadBlast orthomcl.config.template similarSequences.txt',
            'orthomclPairs orthomcl.config.template pairs.log cleanup=yes',
            'orthomclDumpPairsFiles orthomcl.config.template',
            f'mcl mclInput --abc -I {self.Inflation} -o mclOutput',
            'orthomclMclToGroups G_ 1 < mclOutput > groups.txt'
        ]

        for i in cmd:
            # time.sleep(1)
            os.system(i)
            # print('%s -> ok' % i)


# orthomcl
def orth(place, email, percentMatchCutoff, inflation, evalue):
    from file_keep.orthomcl import run
    # run(percentMatchCutoff=percentMatchCutoff, Inflation=inflation, E=evalue)

    # 改变工作目录
    os.chdir(f"{path_get}/file_keep/{place}")

    run_r = Orthomcl()

    run_r.Inflation = inflation
    # run_r.name = name
    run_r.name = re.sub(r'[ :]*', '', time.asctime(time.localtime(time.time())))[-15:]
    run_r.percentMatchCutoff = percentMatchCutoff
    run_r.E = evalue

    run_r.databasrs_set()
    run_r.handle_fasta()
    run_r.handle_data()

    run_r.del_database()
    os.mkdir('result')
    # 将文件处理结果放入同一个文件夹 可以注释掉
    os.system('mv ./pairs  result/')
    os.system('mv ./groups.txt  ./result/')
    print('工作结束')

    cmd = f"""
        cd {path_get}/file_keep/{place}
        run(percentMatchCutoff=percentMatchCutoff, Inflation=inflation, E=evalue)
        python3 ../orthomcl.py
    """
    os.system(cmd)

    cmd = '''
        cd %s/file_keep/%s
        zip -r result.zip result
    ''' % (path_get, place)
    os.system(cmd)

    g = '%s/file_keep/%s/result.zip' % (path_get, place)
    youjian(a=email, b=g, place=place)
