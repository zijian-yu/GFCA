#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yzj
# @FileName :urls.py
# @Time :2021/10/24 16:04

from django.contrib import admin
from django.urls import path, include
from tools import urls
from . import views


urlpatterns = [
    path('dotplot', views.dotplot, name='dotplot_url'),
    path('file_merge', views.file_merge, name='file_merge_url'),
    path('find_replace', views.find_replace, name='find_replace_url'),
    path('extract_row', views.extract_row, name='extract_row_url'),
    path('duplicate_removal', views.duplicate_removal, name='quchong_url'),
    path('formatTofasta', views.formatTofasta, name='formatTofasta_url'),
    path('cds_convert_pep', views.cds_convert_pep, name='cds_convert_pep_url'),

    path('dotplot_run', views.dotplot_run, name='dotplot_run_url'),
    path('file_merge_run', views.file_merge_run, name='file_merge_run_url'),
    path('find_replace_run', views.find_replace_run, name='find_replace_run_url'),
    path('extract_column_run', views.extract_column_run, name='extract_column_run_url'),
    path('duplicate_removal_run', views.duplicate_removal_run, name='duplicate_removal_run_url'),
    path('formatTofasta_run', views.formatTofasta_run, name='formatTofasta_run_url'),
    # path('cds_convert_pep_run', views.cds_convert_pep_run, name='cds_convert_pep_run_url'),
    path('index_email', views.index_email, name='index_email_url'),
]









