import time
import random
import os
import re
import pymysql
from .email import youjian
import string
from Bio import SeqIO
import shutil
import schedule


path_get = os.getcwd()


# 在文件夹内 创建文件夹
def keep_in_file(na, c, fg=None):
    cmd = f'''
        cd {path_get}/file_keep/{c}
        mkdir all_file
    '''
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
    cmd = f'''
        cd {path_get}/file_keep/{c}
        mkdir orthomcl
    '''
    os.system(cmd)

    for file in na:
        with open('%s' % os.path.join(f'{path_get}/file_keep/{c}/orthomcl/', file.name), 'wb') as f:
            for i in file.chunks():
                f.write(i)


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


def delete_folder(c):
    """删除创建的文件夹"""
    folder_path = f"{path_get}/file_keep{c}"
    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)


def name_set():
    """ 在file_keep中创建文件夹
    :return: c 随机文件夹名字
    """
    c = ''.join(random.sample(string.digits, 9))
    os.mkdir(f"file_keep/{c}")
    schedule.every(3).days.do(delete_folder, c)
    return c


# 保存邮件地址
def keep_email(e, address):
    with open(f'{path_get}/file_keep/{e}/email.txt', 'w') as f:
        f.write(address)


def read_email(place):
    global path_get
    with open(f'{path_get}/file_keep/{place}/email.txt', 'r+') as f:
        return f.read()


def zip_file(place, email, files):
    cmd = f"""
        cd {path_get}/file_keep/{place}
        zip -r result.zip {files}
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.zip'
    youjian(a=email, b=g, place=place)


# ================================= gfca start  =============================================


# index 保存邮件地址
def keep_email_index(place, name, email, message):

    with open(f'{path_get}/file_keep/{place}/advise.txt', 'w+') as f:
        f.write(name + '\n')
        f.write(email + '\n')
        f.write(message + '\n')

    g = f'{path_get}/file_keep/{place}/advise.txt'
    youjian(a="1457634627@qq.com", b=g, place=place)


# blast_初筛
def blast_choose(place, myFile1, myFile2, identity, score):
    blast_file = open(f'{path_get}/file_keep/{place}/{myFile1}_{myFile2}.fasta', 'w')
    sequence_file = open(f'{path_get}/file_keep/{place}/sequence.txt', 'w')
    pbl_fasta = open(f'{path_get}/file_keep/{place}/result_2.fasta', 'r')
    gather = set()
    for line in pbl_fasta.readlines():
        ind, sco = line.split()[2], line.split()[11]
        if float(ind) >= float(identity):
            if float(sco) >= float(score):
                blast_file.write(line)
                # sequence = line.split()[1]
                gather.add(line.split()[1])
    for line in gather:
        sequence_file.write(line + '\n')


# blast
def blast_run_tools(place, patten, e, email, score, identity, myFile1, myFile2, nucl_or_prot):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        makeblastdb -in 1.fasta -dbtype {nucl_or_prot} -parse_seqids -out dbname
        {patten} -query 2.fasta -db dbname -outfmt 6 -evalue {e} -out result_2.fasta -num_threads 8
    """
    os.system(cmd)

    blast_choose(place=place, myFile1=myFile1, myFile2=myFile2, identity=identity, score=score)

    files = f"{myFile1}_{myFile2}.fasta sequence.txt"
    zip_file(place=place, email=email, files=files)


# yishan_blast
def yishan_run_blast(place, patten, species_data, e, email, score, identity, nucl_or_prot):

    cmd = f"""
source /etc/profile
cd {path_get}/file_keep/{place}
{patten} -query 1.fasta -db ../../../templates/Nyssaceae/blast_database/{nucl_or_prot}/{species_data} -outfmt 6 -evalue {e} -out result.txt -num_threads 8 -max_target_seqs 1
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.txt'
    youjian(a=email, b=g, place=place)


# diamond
def diamond_run_tools(place, patten, e, email, score, identity, myFile1, myFile2):
    cmd = f"""
    source /etc/profile
    cd {path_get}/file_keep/{place}
    diamond makedb --in 1.fasta -d nr
    diamond {patten} -d nr -q 2.fasta -o result_2.fasta -e {e} -b 2.0 --index-chunks 1
"""
    os.system(cmd)

    blast_choose(place=place, myFile1=myFile1, myFile2=myFile2, identity=identity, score=score)
    files = f"{myFile1}_{myFile2}.fasta sequence.txt"
    zip_file(place=place, email=email, files=files)


# sequence_screening_clas
class sequence_run2(object):
    def __init__(self, place):
        self.id_list = []
        self.place = place

    def new_prot(self):
        new_prot = open(f'{path_get}/file_keep/{self.place}/new_pro.txt', 'w')
        for line in SeqIO.parse(f'{path_get}/file_keep/{self.place}/protein.fasta', 'fasta'):
            if line.id in self.id_list:
                new_prot.write('>' + str(line.id) + '\n' + str(line.seq) + '\n')

    def new_cds(self):
        new_cds = open(f'{path_get}/file_keep/{self.place}/new_cds.txt', 'w')
        for line in SeqIO.parse(f'{path_get}/file_keep/{self.place}/gene.fasta', 'fasta'):
            if line.id in self.id_list:
                new_cds.write('>' + str(line.id) + '\n' + str(line.seq) + '\n')

    def new_gff(self):
        new_gff = open(f'{path_get}/file_keep/{self.place}/new_gff.txt', 'w')
        gff_file = open(f'{path_get}/file_keep/{self.place}/gff.fasta', 'r')
        for line in gff_file:
            gff_id = line.split()[1]
            if gff_id in self.id_list:
                new_gff.write(line)

    def main(self):
        id_file = open(f'{path_get}/file_keep/{self.place}/id.fasta', 'r')
        for line in id_file.readlines():
            self.id_list.append(line.split()[0])

        self.new_prot()
        self.new_cds()
        self.new_gff()


# sequence_screening
def sequence_run(place, email):
    run = sequence_run2(place)
    run.main()

    files = 'new_pro.txt new_cds.txt new_gff.txt'
    zip_file(place=place, email=email, files=files)


# hmmer
def hmmer_run_tools(place, email, e_value):
    cmd = f"""
        cd {path_get}/file_keep/{place}
        hmmbuild db_model.hmm hmm_model.fasta
        hmmsearch -o result.out --noali -E {e_value} db_model.hmm hmm_fa.fasta
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.out'
    youjian(a=email, b=g, place=place)


# pfam
def pfam_run_tools(place, structure, email):
    # pfam_scan.pl -fasta 1.fasta -dir ~/../../../../www/tools/PfamScan/pfamdata -outfile pfam.txt
    cmd = f"""
    source /etc/profile
    cd {path_get}/file_keep/{place}
    pfam_scan.pl -fasta 1.fasta -dir ~/../../../../www/tools/PfamScan/pfamdata -outfile pfam.txt
"""
    os.system(cmd)
    structure_list = list(structure.split(';'))

    # 根据用户设置的结构域进行筛选
    file = open(f'{path_get}/file_keep/{place}/pfam.txt', 'r')
    pfam_new = open(f'{path_get}/file_keep/{place}/pfam_new.txt', 'w')
    deal = ('#', ' ', '\n')
    for line in file:
        if line.startswith(deal):
            continue
        structure = line.split()[6]
        if structure in structure_list:
            pfam_new.write(line)

    files = 'pfam.txt pfam_new.txt'
    zip_file(email=email, place=place, files=files)


# meme
def meme_run_tools(place, email, patten, motif, motif_num, minw, maxw):
    cmd = f"""
    cd {path_get}/file_keep/{place}
    meme 1.fasta -{patten} -oc . -mod {motif} -nmotifs {motif_num} -minw {minw} -maxw {maxw}
"""
    os.system(cmd)

    files = 'logo*.png meme.*'
    zip_file(email=email, place=place, files=files)


# codonw_run
def codonw_run_tools(place, email, file_name):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        codonw 1.fasta -all_indices -nomenu {file_name}.out {file_name}.blk
    """
    os.system(cmd)

    files = f'{file_name}.out {file_name}.blk'
    zip_file(email=email, files=files, place=place)


# cpgfinder_run
def cpgfinder_run_tools(place, email, island, cg_percent, gc_ratio, myFile_name):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        cpgfinder -f {myFile_name} -l {island} -p {cg_percent} -r {gc_ratio} > cpgfinder_result.txt
    """
    os.system(cmd)

    files = 'cpgfinder_result.txt'
    zip_file(email=email, place=place, files=files)


# dupgen_finder_run
def dupgen_finder_run_tools(email, place, patten, target, outgroup):
    if patten == 'general':
        cmd = f"""
source /etc/profile
cd {path_get}/file_keep/{place}
mv {target}* ../programs/dupgen_finder/data/
cd ../programs/dupgen_finder
perl DupGen_finder.pl -i data -t {target} -c {outgroup} -o results
zip -r result.zip results/
mv result.zip {path_get}/file_keep/{place}
rm -rf data/*
rm -rf result*
        """
        os.system(cmd)

        g = f'{path_get}/file_keep/{place}/result.zip'
        youjian(a=email, b=g, place=place)

    else:
        cmd = f"""
source /etc/profile
cd {path_get}/file_keep/{place}
mv {target}* ../programs/dupgen_finder/data/
cd ../programs/dupgen_finder
perl DupGen_finder-unique.pl -i data -t {target} -c {outgroup} -o results
zip -r result.zip results/
mv result.zip {path_get}/file_keep/{place}
rm -rf data/*
rm -rf result*
        """
        os.system(cmd)

        g = f'{path_get}/file_keep/{place}/result.zip'
        youjian(a=email, b=g, place=place)


# cd_hit_run
def cd_hit_run_tools(place, email, similar, patten):
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
        source /etc/profile
        cd {path_get}/file_keep/{place}
        cd-hit -i 1.fasta -o cd_hit_result.fasta -c {similar} -n {patten} -d 0 -T 1
    """
    os.system(cmd)

    files = 'cd_hit_result.fasta.clstr cd_hit_result.fasta'
    zip_file(email=email, place=place, files=files)


# sequence_aligment
def sequence_aligment_run_tools(place, email, patten):
    if patten == "clustalw2":
        cmd = f"""
            source /etc/profile
            cd {path_get}/file_keep/{place}
            clustalw2 1.fasta
        """
        os.system(cmd)
        g = f'{path_get}/file_keep/{place}/1.aln'
        youjian(a=email, b=g, place=place)

    elif patten == "muscle":
        cmd = f"""
            source /etc/profile
            cd {path_get}/file_keep/{place}
            muscle -in 1.fasta -clwout 1.aln
        """
        os.system(cmd)
        g = f'{path_get}/file_keep/{place}/1.aln'
        youjian(a=email, b=g, place=place)
    elif patten == "mafft":
        cmd = f"""
            source /etc/profile
            cd {path_get}/file_keep/{place}
            mafft --auto 1.fasta > mafft.fasta
        """
        os.system(cmd)
        g = f'{path_get}/file_keep/{place}/mafft.fasta'
        youjian(a=email, b=g, place=place)


# mcscanx
def mcscanx_run_tools(place, email, file_name):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        MCScanX {file_name}
        duplicate_gene_classifier {file_name} > gene_classifier.txt
    """
    os.system(cmd)

    files = f'{file_name}.html {file_name}.collinearity {file_name}.tandem gene_classifier.txt'
    zip_file(email=email, place=place, files=files)


# colinearscan
def colinearscan_run_tools(place, email, score, e_value, hitnum, pos_order, cdsvchr, blast_name, gff1_name, gff2_name,
                           file_name):
    # perl ../programs/run_colinearscan.pl ZM-OS.filter.blast 1e-5 0 30 ZM.gff OS.gff order n ZM_OS

    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        perl ../programs/run_colinearscan.pl {blast_name} {e_value} {score} {hitnum} {gff1_name} {gff2_name} {pos_order} {cdsvchr} {file_name}
    """
    os.system(cmd)

    files = 'block pair'
    zip_file(email=email, place=place, files=files)


# orthofinder
def orthfinder_run_tools(place, email, gene_tree, blast, alignment, tree_way, number):
    # -M    基因树推断方法（默认为dendroblast）可选：dendroblast ，msa
    # -S    序列比对使用的程序 （默认为Blast）可选：blast, mmseqs, blast_gz, diamond（推荐使用diamond，比对速度快）
    # # -A    多序列联配方式，该选项仅当 - M    msa    选项时才有效（默认为mafft）可选：muscle, mafft
    # # -T    建树方式，该选项仅当 - M    msa    选项时才有效 （默认为fasttree）可选：iqtree, raxml - ng, fasttree, raxml
    # -I    设定MCL的膨胀系数（默认为1.5）
    if gene_tree == str("msa"):
        cmd = f"""
            source /etc/profile
            cd {path_get}/file_keep/{place}
            orthofinder -f orthomcl -M {gene_tree} -S {blast} -A {alignment} -T {tree_way} -I {number} 
        """
        os.system(cmd)
    else:
        cmd = f"""
            source /etc/profile
            cd {path_get}/file_keep/{place}
            orthofinder -f orthomcl/ -M {gene_tree} -S {blast} -I {number} 
        """
        os.system(cmd)

    cmd2 = f"""
        cd {path_get}/file_keep/{place}/orthomcl/OrthoFinder
        zip -r result.zip *
        mv result.zip ../../
    """
    os.system(cmd2)

    g = f'{path_get}/file_keep/{place}/result.zip'
    youjian(a=email, b=g, place=place)


# paraAT
def paraat_run_tools(place, email, patten):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        ParaAT.pl -h hom.homologs -n cds.cds -a pep.pep -p ../programs/proc -o result -f {patten}
        cat result/*.{patten} > result.{patten}
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.{patten}'
    youjian(a=email, b=g, place=place)


# kaks_calculator
def kaks_calculator_run_tools(place, email, patten):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        KaKs_Calculator -i kaks.axt -o result_kaks.txt -m {patten}
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result_kaks.txt'
    youjian(a=email, b=g, place=place)


# fasttree
def fasttree_run_tools(place, email, patten):
    if patten == "JTT+CAT":
        cmd = f"""
            source /etc/profile
            cd {path_get}/file_keep/{place}
            FastTree 1.phy > 1.nwk
        """
        os.system(cmd)
    else:
        cmd = f"""
            source /etc/profile
             cd {path_get}/file_keep/{place}
            FastTree -gtr -nt  1.phy > 1.nwk
        """
        os.system(cmd)

    g = f'{path_get}/file_keep/{place}/1.nwk'
    youjian(a=email, b=g, place=place)


# iqtree
def iqtree_run_tools(place, email, patten, bs_num):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        iqtree2 -s result.phy -m {patten} -bb {bs_num} -bnni -redo
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.phy.treefile'
    youjian(a=email, b=g, place=place)


# phyml
def phyml_run_tools(email, place, patten, nice, model, frequency, num, small):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        phyml -i 1.phy -d {patten} -b {num} -m {small} -f {frequency} -v e -a e -o {nice}
    """
    os.system(cmd)

    files = "1.phy_*"
    zip_file(email=email, place=place, files=files)


# paml_deal
def paml_run_one(patten, small, path, email, place):
    file_data = ""
    codeml_file = open(f'{path_get}/file_keep/paml/{patten}/{small}/codeml.ctl', 'r')
    for line in codeml_file.readlines():
        hh = re.sub(r'\d{9,}', f'{path}', line)
        file_data += hh
    co_new = open(f'{path_get}/file_keep/paml/{patten}/{small}/codeml.ctl', 'w')
    co_new.write(file_data)


# paml
def paml_run_two(place, email, patten, small):
    cmd = f"""
        source /etc/profile
        cd {path_get}/file_keep/{place}
        codeml ../paml/{patten}/{small}/codeml.ctl
    """
    os.system(cmd)

    g = '%s/file_keep/%s/resoult.txt' % (path_get, place)
    youjian(a=email, b=g, place=place)


# ================================= gfca end  ==============================================


# ================================= tools start =============================================


# dotplot
def dotplot_run_tools(place, email, evalue, score, hitnum, repnum, heng_chromosomes, shu_chromosomes, heng_species,
                      shu_species, blast_name, shu_gff_name, heng_gff_name, shu_len_name, heng_len_name):
    global path_get
    # perl ../programs/dotplot.pl 1e-5 100 5 20 1_2_3 1_2_3 si.lens zm.lens si.gff zm.gff zm_si.blast zm si
    cmd = f"""
        cd {path_get}/file_keep/{place}
        perl ../programs/dotplot.pl {evalue} {score} {hitnum} {repnum} {heng_chromosomes} {shu_chromosomes} {shu_len_name} {heng_len_name} {shu_gff_name} {heng_gff_name} {blast_name} {heng_species} {shu_species} 
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.png'
    youjian(a=email, b=g, place=place)


# file_merge_run
# def file_merge_run_tools(place, email, file_name):
def file_merge_run_tools(place, email):
    cmd = f"""
        cd '{path_get}/file_keep/{place}'
        cat orthomcl/* >> result.txt
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.txt'
    youjian(a=email, b=g, place=place)


# find_replace_run
def find_replace_run_tools(place, email, old_name, new_name):
    cmd = f"""
        cd '{path_get}/file_keep/{place}'
        sed -e 's/{old_name}/{new_name}/g' upload.txt > result.txt
    """
    os.system(cmd)

    g = f'{path_get}/file_keep/{place}/result.txt'
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


# quchong_run
def quchong_run(place, email):
    id_list = []
    new_file = open(f'{path_get}/file_keep/{place}/new.fasta', 'w')
    for line in SeqIO.parse(f'{path_get}/file_keep/{place}/1.fasta', 'fasta'):
        if line.id not in id_list:
            id_list.append(line.id)
            new_file.write(">" + str(line.id) + "\n" + str(line.seq) + "\n")
    new_file.close()

    g = f'{path_get}/file_keep/{place}/new.fasta'
    youjian(a=email, b=g, place=place)


# formatTofasta
def format_fasta_run(place, email, file_name, patten, patten_end, file_last):
    SeqIO.convert(f"{path_get}/file_keep/{place}/{file_name}.{file_last}", f"{patten}",
                  f"{path_get}/file_keep/{place}/{file_name}.{patten_end}", f"{patten_end}")

    g = f'{path_get}/file_keep/{place}/{file_name}.{patten_end}'
    youjian(a=email, b=g, place=place)


# ================================= tools end ================================================


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




