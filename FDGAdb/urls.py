"""FDGAdb URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin
from search import views as view

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^$', view.home),
    url(r'^FDRdb/', view.home),
    url(r'^home/', view.home_sub),
    url(r'^home_search/(?P<tissue>[\w-]+)/', view.home_tissue),
    url(r'^home_species/(?P<species>[\w-]+)/', view.home_species),
    url(r'^browse/', view.browse),
    url(r'^help/', view.help),
    url(r'^browse_data/', view.browse_data),
    url(r'^browse_miRNA/(?P<miRNA>[\w-]+)/', view.browse_miRNA),
    url(r'^browse_lncRNA/(?P<lncRNA>[\w-]+)/', view.browse_lncRNA),
    url(r'^browse_mRNA/(?P<mRNA>[\w-]+)/', view.browse_mRNA),
    url(r'^browse_circRNA/(?P<circRNA>[\w-]+)/', view.browse_circRNA),
    url(r'^browse_siRNA/(?P<siRNA>[\w-]+)/', view.browse_siRNA),
    url(r'^browse_snRNA/(?P<snRNA>[\w-]+)/', view.browse_snRNA),
    url(r'^browse_tissue/(?P<tissue>[\w-]+)/', view.browse_tissue),
    url(r'^browse_species/(?P<species>[\w-]+)/', view.browse_species),
    url(r'^browse_disease/(?P<disease>[\w-]+)/', view.browse_disease),
    url(r'^browse_detail/(?P<RNA>[\w\s-]+)/(?P<disease>[\w\s-]+)/(?P<species>[\w\s-]+)/(?P<pmid>[\w\s-]+)/', view.browse_detail),
    url(r'^browse_detail_data/(?P<data>[\w\s-]+)/(?P<species>[\w\s-]+)/(?P<disease>[\w\s-]+)/', view.browse_detail_data),
    url(r'^browse_data_type/(?P<type>[\w-]+)/', view.browse_data_type),
    url(r'^browse_data_species/(?P<species>[\w-]+)/', view.browse_data_species),
    url(r'^browse_data_disease/(?P<disease>[\w-]+)/', view.browse_data_disease),
    url(r'^browse_data_tissue/(?P<tissue>[\w-]+)/', view.browse_data_tissue),
    url(r'^download_csv/(?P<typename>[\w-]+)/', view.download_csv),
    url(r'^download_txt/(?P<typename>[\w-]+)/', view.download_txt),
    url(r'^search/', view.search),
    url(r'^download/', view.download),
    url(r'^link/', view.link),
    url(r'^help_overview/', view.help_overview),
    url(r'^help_home/', view.help_home),
    url(r'^results/', view.results),
    url(r'^analysis/', view.analysis),
    url(r'^Apply/', view.upload),
    url(r'^diff_download/(?P<dataname>[\w\s-]+)/', view.diff_download),
    url(r'^expression_example/', view.expression_example),
    url(r'^class_example/', view.class_example),
]
