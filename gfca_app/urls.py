#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yzj
# @FileName :urls.py
# @Time :2021/10/27 20:40


from django.contrib import admin
from django.urls import path, include
from tools import urls
from . import views


urlpatterns = [
    path('blast', views.blast, name='blast_url'),
    path('cd_hit', views.cd_hit, name='cd_hit_url'),
    path('circos', views.circos, name='circos_url'),
    path('codonw', views.codonw, name='codonw_url'),
    path('colinearscan', views.colinearscan, name='colinearscan_url'),
    path('cpgfinder', views.cpgfinder, name='cpgfinder_url'),
    path('diamond', views.diamond, name='diamond_url'),
    path('dupgen_finder', views.dupgen_finder, name='dupgen_finder_url'),
    path('fasttree', views.fasttree, name='fasttree_url'),
    path('gsds', views.gsds, name='gsds_url'),
    path('heatmap', views.heatmap, name='heatmap_url'),
    path('hmmer', views.hmmer, name='hmmer_url'),
    path('iq_tree', views.iq_tree, name='iq_tree_url'),
    path('kaks_calculator', views.kaks_calculator, name='kaks_calculator_url'),
    path('mcscanx', views.mcscanx, name='mcscanx_url'),
    path('meme_suit', views.meme_suit, name='meme_suit_url'),
    path('orthofinder', views.orthofinder, name='orthofinder_url'),
    path('orthomcl', views.orthomcl, name='orthomcl_url'),
    path('paml', views.paml, name='paml_url'),
    path('paraat', views.paraat, name='paraat_url'),
    path('pfam', views.pfam, name='pfam_url'),
    path('phyml', views.phyml, name='phyml_url'),
    path('sequence_alignment', views.sequence_alignment, name='sequence_alignment_url'),
    path('sequence_screening', views.sequence_screening, name='sequence_screening_url'),
    path('sequenceserver', views.sequenceserver, name='sequenceserver_url'),
    path('shinycircos', views.shinycircos, name='shinycircos_url'),
    path('string', views.string, name='string_url'),


    path('blast_run', views.blast_run, name='blast_run_url'),
    path('diamond_run', views.diamond_run, name='diamond_run_url'),
    path('sequence_screening_run', views.sequence_screening_run, name='sequence_screening_run_url'),
    path('hmmer_run', views.hmmer_run, name='hmmer_run_url'),
    path('pfam_run', views.pfam_run, name='pfam_run_url'),
    path('meme_suit_run', views.meme_suit_run, name='meme_suit_run_url'),
    path('codonw_run', views.codonw_run, name='codonw_run_url'),
    path('cpgfinder_run', views.cpgfinder_run, name='cpgfinder_run_url'),
    # gsds、sequenceserver 不需要run - 已经配置
    path('dupgen_finder_run', views.dupgen_finder_run, name='dupgen_finder_run_url'),
    path('cd_hit_run', views.cd_hit_run, name='cd_hit_run_url'),
    # string、heatmap 没有配置
    path('sequence_alignment_run', views.sequence_alignment_run, name='sequence_alignment_run_url'),
    path('mcscanx_run', views.mcscanx_run, name='mcscanx_run_url'),
    path('colinearscan_run', views.colinearscan_run, name='colinearscan_run_url'),
    # ortgomcl 没有搭建
    path('orthofinder_run', views.orthofinder_run, name='orthofinder_run_url'),
    path('paraat_run', views.paraat_run, name='paraat_run_url'),
    path('kaks_calculator_run', views.kaks_calculator_run, name='kaks_calculator_run_url'),
    path('fasttree_run', views.fasttree_run, name='fasttree_run_url'),
    path('iq_tree_run', views.iq_tree_run, name='iq_tree_run_url'),
    path('phyml_run', views.phyml_run, name='phyml_run_url'),
    path('paml_run', views.paml_run, name='paml_run_url'),
]

