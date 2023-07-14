"""GFCA URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from tools import urls, views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.index),

    path('index', views.index, name='index_url'),
    path('tools', views.tools, name='tools_url'),
    path('gfca', views.gfca, name='gfca_url'),
    path('database', views.database, name='database_url'),
    path('help', views.help, name='help_url'),
    path('contact', views.contact, name='contact_url'),

    path('successful', views.successful, name='successful_url'),
    path('successful_2', views.successful_2, name='successful_url'),
    path('help', views.help, name='help_url'),
    path('index_email', views.index_email, name='index_email_url'),
    path('blast_match', views.blast_match, name=''),
    path('blast_match_run', views.blast_match_run, name=''),

    path('tools/', include('tools.urls')),
    path('gfca/', include('gfca_app.urls')),

]


