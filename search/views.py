

import django.conf
from django.shortcuts import render
from django.http import HttpResponse
from django.template import loader
from .models import fibrosis_gene, fibrosis_data
from .forms import SearchForm, SearchspeciesForm
from wsgiref.util import FileWrapper
from django.db.models import Q
import re, os, json, pandas as pd, numpy as np
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import math
# Create your views here.
import zipfile



mol_info = pd.read_table(os.path.join(django.conf.settings.DATA_DIR,"molecular_info.txt"), header=0, index_col=0)
symbol_id = pd.read_table(os.path.join(django.conf.settings.DATA_DIR,"id_symbol_anno.txt"), header=0, index_col=0)
go_info = pd.read_table(os.path.join(django.conf.settings.DATA_DIR,"go_anno.txt"), header=None, index_col=0)
kegg_info = pd.read_table(os.path.join(django.conf.settings.DATA_DIR,"kegg_anno.txt"), header=None, index_col=0)
go_all = pd.read_table(os.path.join(django.conf.settings.DATA_DIR,"go_all.txt"), header=None, index_col=None)
kegg_all = pd.read_table(os.path.join(django.conf.settings.DATA_DIR,"kegg_all.txt"), header=None, index_col=None)

def home(request):
    return HttpResponse(loader.get_template('FDRdb.html').render({}, request))

def browse(request):
    template = loader.get_template("Browse.html")
    return HttpResponse(template.render(request))

def browse_data(request):
    template = loader.get_template("Browse_data.html")
    return HttpResponse(template.render(request))
def help(request):
    template = loader.get_template("Help.html")
    return HttpResponse(template.render(request))

def home_sub(request):
    template = loader.get_template("Home.html")
    return HttpResponse(template.render(request))

def help_overview(request):
    template = loader.get_template("help_overview.html")
    return HttpResponse(template.render(request))

def help_home(request):
    template = loader.get_template("help_home.html")
    return HttpResponse(template.render(request))


def browse_miRNA(request,miRNA):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(gene_symbol=miRNA)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")

    return HttpResponse(template.render({"result_info":result_info,"file_name":miRNA},request))

def browse_lncRNA(request,lncRNA):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(gene_symbol=lncRNA)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":lncRNA},request))

def browse_mRNA(request,mRNA):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(gene_symbol=mRNA)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":mRNA},request))

def browse_siRNA(request,siRNA):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(gene_symbol=siRNA)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":siRNA},request))

def browse_snRNA(request,snRNA):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(gene_symbol=snRNA)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":snRNA},request))

def browse_circRNA(request,circRNA):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(gene_symbol=circRNA)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":circRNA},request))

def browse_tissue(request,tissue):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(tissue=tissue)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_namel":tissue},request))

def browse_disease(request,disease):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(disease=disease)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_namel":disease},request))


def browse_species(request,species):
    template = loader.get_template("Browse_table.html")
    result_info = fibrosis_gene.objects.filter(species=species)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":species},request))

def browse_detail(request,RNA,disease,species,pmid):
    template = loader.get_template("Detail.html")
    species = species.replace(' ', '_')
    disease = disease.replace(' ', '_')
    result_info = fibrosis_gene.objects.filter(Q(gene_symbol=RNA) & Q(species=species) & Q(disease=disease) & Q(pmid=pmid))
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":RNA},request))

def browse_detail_data(request,data,species,disease):
    template = loader.get_template("Detail_data.html")
    data = data.replace(' ', '_')
    species = species.replace(' ', '_')
    disease = disease.replace(' ', '_')
    result_info = fibrosis_data.objects.filter(Q(data_id=data) & Q(species=species) & Q(disease=disease))
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":data},request))

def browse_data_type(request,type):
    template = loader.get_template("Browse_data_table.html")
    result_info = fibrosis_data.objects.filter(type=type)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":type},request))

def browse_data_species(request,species):
    template = loader.get_template("Browse_data_table.html")
    result_info = fibrosis_data.objects.filter(species=species)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":species},request))

def browse_data_disease(request,disease):
    template = loader.get_template("Browse_data_table.html")
    result_info = fibrosis_data.objects.filter(disease=disease)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":disease},request))

def browse_data_tissue(request,tissue):
    template = loader.get_template("Browse_data_table.html")
    result_info = fibrosis_data.objects.filter(tissue=tissue)
    for s in result_info:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.strip().split(";")
    return HttpResponse(template.render({"result_info":result_info,"file_name":tissue},request))


@csrf_exempt
def results(request):
    template = loader.get_template("Results.html")
    if request.method == "POST":
        form = SearchForm(request.POST)
        if form.is_valid():
            search_species = form.cleaned_data["search_species"]
            search_tissues = form.cleaned_data["search_tissues"]
            gf = form.cleaned_data["filter"]
            gf = gf.replace(' ', '_')
            gf = [i for i in re.split("[,;\s]\s*", gf) if not i == ""]
            if search_species == 'ALL':
                if search_tissues == 'ALL':
                    result_fibrosis = [fibrosis_gene.objects.filter(
                        Q(gene_symbol=i) | Q(gene_id=i) | Q(disease=i)) for i in gf]
                else:
                    result_fibrosis = [fibrosis_gene.objects.filter(
                        Q(gene_symbol=i) | Q(gene_id=i) | Q(disease=i)).filter(Q(tissue=search_tissues)) for i in gf]
            else:
                if search_tissues == 'ALL':
                    result_fibrosis = [fibrosis_gene.objects.filter(
                        Q(gene_symbol=i) | Q(gene_id=i) | Q(disease=i)).filter(Q(species=search_species)) for i in gf]
                else:
                    result_fibrosis = [fibrosis_gene.objects.filter(
                        Q(gene_symbol=i) | Q(gene_id=i) | Q(disease=i)).filter(Q(tissue=search_tissues)).filter(Q(species=search_species)) for i in gf]
            result_fibrosis = result_fibrosis[0]
            result_fibrosis_final = {result_fibrosis[i] for i in range(0, len(result_fibrosis))}
            for s in result_fibrosis_final:
                s.species = s.species.replace('_', ' ')
                s.tissue = s.tissue.replace('_', ' ')
                s.expression = s.expression.replace('-', ' ')
                s.disease = s.disease.replace('_', ' ')
                s.pmid = s.pmid.replace(';', ' ')
                s.pmid = s.pmid.strip().split(' ')
            number = len(result_fibrosis_final)
            return HttpResponse(template.render({"result_fibrosis": result_fibrosis_final,"number":number}, request))
    else:
        form = SearchForm()

    # return HttpResponse(loader.get_template('GI/search.html').render({"form":form},request))
    return render(request, "Results.html", {'form': form})


def home_tissue(request,tissue):
    template = loader.get_template("home_search.html")
    result_fibrosis = fibrosis_gene.objects.filter(tissue=tissue)
    for s in result_fibrosis:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.replace(';', ' ')
        s.pmid = s.pmid.strip().split(' ')

    return HttpResponse(template.render({"result_fibrosis": result_fibrosis,"tissue_name": tissue}, request))

def home_species(request,species):
    template = loader.get_template("home_search.html")
    result_fibrosis = fibrosis_gene.objects.filter(species=species)
    for s in result_fibrosis:
        s.species = s.species.replace('_', ' ')
        s.tissue = s.tissue.replace('_', ' ')
        s.expression = s.expression.replace('-', ' ')
        s.disease = s.disease.replace('_', ' ')
        s.pmid = s.pmid.replace(';', ' ')
        s.pmid = s.pmid.strip().split(' ')

    return HttpResponse(template.render({"result_fibrosis": result_fibrosis}, request))



def search(request):
    template = loader.get_template("Search.html")
    return HttpResponse(template.render(request))

def download(request):
    template = loader.get_template("Download.html")
    return HttpResponse(template.render(request))

def link(request):
    template = loader.get_template("Links.html")
    return HttpResponse(template.render(request))


def download_csv(request,typename):

    filename = "/alidata/database/FDGAdb/data/download/"+typename+"_FDs_info.csv" # Select your file here.
    wrapper = FileWrapper(open(filename, 'rb'))
    response = HttpResponse(wrapper, content_type='application/octet-stream')
    response['Content-Length'] = os.path.getsize(filename)
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format(typename+"_FDs_info.csv")

    return response

def download_txt(request,typename):

    filename = "/alidata/database/FDGAdb/data/download/"+typename+"_FDs_info.txt" # Select your file here.
    wrapper = FileWrapper(open(filename, 'rb'))
    response = HttpResponse(wrapper, content_type='application/octet-stream')
    response['Content-Length'] = os.path.getsize(filename)
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format(typename+"_FDs_info.txt")

    return response

def analysis(request):
    template = loader.get_template("Analysis.html")
    return HttpResponse(template.render(request))

@csrf_exempt
def upload(request):
            if request.method == 'POST':
                data = request.FILES['analysis']
                classtype = request.FILES['class']
                dataname = '%s/%s' % (settings.MEDIA_ROOT, data.name)
                classname = '%s/%s' % (settings.MEDIA_ROOT, classtype.name)
                ff = open(dataname, 'w+')
                ff.close()
                with open(dataname, 'wb') as f:
                    for fdata in data.chunks():
                        f.write(fdata)
                tt = open(classname, 'w+')
                tt.close()
                with open(classname, 'wb') as t:
                    for tclass in classtype.chunks():
                        t.write(tclass)
                data1 = pd.read_table(dataname, header=0, index_col=0)
                class1 = pd.read_table(classname, header=0, index_col=0)
                data2 = data1
                which = lambda lst: list(np.where(lst)[0])
                lst1 = map(lambda x: x == 1, class1.iloc[0, :])
                lst2 = map(lambda x: x == 2, class1.iloc[0, :])
                samples1 = which(lst1)
                samples2 = which(lst2)
                samples_all = len(samples1)+len(samples2)
                data2 = data2.iloc[:, samples1+samples2]
                type1 = data2.iloc[:, samples1].mean(axis=1)
                type2 = data2.iloc[:, samples2].mean(axis=1)
                # fold_change
                fold = type2 - type1
                #plt.hist(fold)
                #plt.title("Histogram of fold-change")
                #fdimg_name = '%s/%s' % (settings.MEDIA_ROOT, data.name) + '.histogram_fd.png'
                #plt.savefig(fdimg_name)
                #plt.close()
                # p
                pvalue = []
                filtered_ids = []
                filtered_ids_heat = []
                for i in range(0, data2.shape[0]):
                    ttest = stats.ttest_ind(data2.iloc[i, samples1], data2.iloc[i, samples2])
                    pvalue.append(ttest[1])
                    if (abs(list(fold)[i]) >= 1) and (pvalue[i] <= 0.05):
                        filtered_ids.append(i)
                    if pvalue[i] <= 0.05:
                        filtered_ids_heat.append(i)
                filtered = data2.iloc[filtered_ids, :]
                filtered_heat = data2.iloc[filtered_ids_heat, :]

                # Histogram of the p-values
                #plt.hist(-np.log(pvalue))
                #plt.title("Histogram of -log(p-value)")
                #pimg_name = '%s/%s' % (settings.MEDIA_ROOT, data.name) + '.histogram_p.png'
                #plt.savefig(pimg_name)
                #plt.close()
                # volcanoimg_name
                genearray = np.asarray(pvalue)
                result = pd.DataFrame({'pvalue': genearray, 'FoldChange': fold})
                result['log(pvalue)'] = -np.log10(result['pvalue'])
                result['GROUP'] = 'none'
                result['fib_relate'] = 'NA'
                result['size'] = np.abs(result['FoldChange']) / 10
                result.loc[(result.FoldChange > 1) & (result.pvalue < 0.05), 'GROUP'] = 'up'
                result.loc[(result.FoldChange < -1) & (result.pvalue < 0.05), 'GROUP'] = 'down'
                form = SearchspeciesForm(request.POST)
                up_fibrosis = 0
                down_fibrosis = 0
                kegg_index = 0
                kegg_address = '/alidata/database/EMTRegulome/emtMotif/search/static/FDRdb/static/analysis_data/no_result.kegg_result.csv'
                keggimg_address = '/static/FDRdb/static/analysis_data/no_result.kegg.png'
                go_index = 0
                go_address = '/alidata/database/EMTRegulome/emtMotif/search/static/FDRdb/static/analysis_data/no_result.go_result.csv'
                goimg_address = '/static/FDRdb/static/analysis_data/no_result.go.png'
                if form.is_valid():
                    species_select = form.cleaned_data["species_select"]
                    if species_select == 'human':
                        A = data2._stat_axis.values
                        B = mol_info.loc[mol_info.species == 'homo_sapiens', 'symbol'].tolist()
                        C = mol_info.loc[mol_info.species == 'homo_sapiens', 'id'].tolist()
                        result.loc[[i in B for i in result._stat_axis.values], 'fib_relate'] = 'fibrosis'
                        result.loc[[str(i) in C for i in result._stat_axis.values], 'fib_relate'] = 'fibrosis'
                        up_fibrosis = len(which(map(lambda x,y: (x == 'fibrosis') and (y == 'up'), result.loc[:,'fib_relate'],result.loc[:,'GROUP'])))
                        down_fibrosis = len(which(map(lambda x,y: (x == 'fibrosis') and (y == 'down'), result.loc[:,'fib_relate'],result.loc[:,'GROUP'])))
                        kegg_bg = kegg_all.shape[0]
                        go_bg = go_all.shape[0]
                        sig_list = result.loc[(result.pvalue < 0.05),]._stat_axis.values
                        name_all = set(sig_list)
                        sig_list_matrix = pd.DataFrame({'name': list(name_all)})
                        id_list_1 = pd.merge(sig_list_matrix, symbol_id, left_on='name', right_on='symbol').loc[:, 'id']
                        id_list_2 = pd.merge(sig_list_matrix, symbol_id, left_on='name', right_on='id').loc[:, 'id']
                        id_list = id_list_1.tolist() + id_list_2.tolist()
                        inbg = len(list(set(id_list) & set(kegg_all.iloc[:, 0])))
                        kegg_result = pd.DataFrame({'pathway_kegg': kegg_info.iloc[:,0],'p_value': 'NA','-ln p': 'NA', 'set': 'NA','background': 'NA', 'in background': 'NA','in set': 'NA'})
                        kegg_result.index = kegg_result.loc[:,'pathway_kegg']
                        kegg_result = kegg_result.drop(columns='pathway_kegg')
                        goinbg = len(list(set(id_list) & set(go_all.iloc[:, 0])))
                        go_result = pd.DataFrame({'pathway_go': go_info.iloc[:,0],'p_value': 'NA','-ln p': 'NA', 'set': 'NA','background': 'NA', 'in background': 'NA','in set': 'NA'})
                        go_result.index = go_result.loc[:,'pathway_go']
                        go_result = go_result.drop(columns='pathway_go')
                        if len(id_list) > 0:
                            kegg_index = 1
                            for i in range(0,kegg_info.shape[0]):
                                set1 = len(kegg_info.iloc[i,:].dropna())
                                inset_1 = len(list(set(id_list) & set(kegg_info.iloc[i,:].dropna())))-1
                                kegg_result.loc[kegg_info.iloc[i,0],'p_value'] = stats.hypergeom.sf(inset_1, kegg_bg, inbg, set1)
                                kegg_result.loc[kegg_info.iloc[i,0],'-ln p'] = -math.log(kegg_result.loc[kegg_info.iloc[i,0],'p_value'],10)
                                kegg_result.loc[kegg_info.iloc[i,0],'set'] = set1
                                kegg_result.loc[kegg_info.iloc[i,0],'in set'] = inset_1+1
                                kegg_result.loc[kegg_info.iloc[i,0],'background'] = kegg_bg
                                kegg_result.loc[kegg_info.iloc[i,0],'in background'] = inbg
                            kegg_result_final_orgin = kegg_result.loc[(kegg_result.p_value < 0.05),:]
                            kegg_result_final = kegg_result_final_orgin.copy()
                            if kegg_result_final.shape[0] > 0:
                                kegg_result_final.sort_values(by="p_value", axis=0, ascending=True, inplace=True)
                                kegg_address = '/alidata/database/EMTRegulome/emtMotif/search/static/FDRdb/static/analysis_data/' + '%s' % (
                                    data.name) + '.kegg_result.csv'
                                kegg_result_final.to_csv(kegg_address)
                                if kegg_result_final.shape[0] > 20:
                                    kegg_result_final_final = kegg_result_final.iloc[0:20,:]
                                    KEGG_PATHWAY = list(kegg_result_final_final._stat_axis.values)
                                    kegg_data = list(kegg_result_final_final.loc[:,'-ln p'])
                                    plt.barh(range(len(kegg_data)), kegg_data, tick_label=KEGG_PATHWAY)
                                    plt.subplots_adjust(left=0.5)
                                    plt.yticks(fontsize=7)
                                    plt.title("Pathway enrichment analysis (KEGG)")
                                    plt.xlabel('-log10(p value)', fontsize=12)
                                    keggimg_name = '%s/%s' % (settings.MEDIA_ROOT, data.name) + '.kegg.png'
                                    plt.savefig(keggimg_name)
                                    plt.close()
                                else:
                                    KEGG_PATHWAY = list(kegg_result_final._stat_axis.values)
                                    kegg_data = list(kegg_result_final.loc[:,'-ln p'])
                                    plt.barh(range(len(kegg_data)), kegg_data, tick_label=KEGG_PATHWAY)
                                    plt.subplots_adjust(left=0.5)
                                    plt.yticks(fontsize=7)
                                    plt.title("Pathway enrichment analysis (KEGG)")
                                    plt.xlabel('-log10(p value)', fontsize=12)
                                    keggimg_name = '%s/%s' % (settings.MEDIA_ROOT, data.name) + '.kegg.png'
                                    plt.savefig(keggimg_name)
                                    plt.close()
                                keggimg_address = '/static/FDRdb/static/analysis_data/' + '%s' % (data.name) + '.kegg.png'
                            go_index = 1
                            for i in range(0,go_info.shape[0]):
                                set1 = len(go_info.iloc[i,:].dropna())
                                inset_1 = len(list(set(id_list) & set(go_info.iloc[i,:].dropna())))-1
                                go_result.loc[go_info.iloc[i,0],'p_value'] = stats.hypergeom.sf(inset_1, go_bg, goinbg, set1)
                                go_result.loc[go_info.iloc[i,0],'-ln p'] = -math.log(go_result.loc[go_info.iloc[i,0],'p_value'],10)
                                go_result.loc[go_info.iloc[i,0],'set'] = set1
                                go_result.loc[go_info.iloc[i,0],'in set'] = inset_1+1
                                go_result.loc[go_info.iloc[i,0],'background'] = go_bg
                                go_result.loc[go_info.iloc[i,0],'in background'] = goinbg
                            go_result_final_orgin = go_result.loc[(go_result.p_value < 0.05),:]
                            go_result_final = go_result_final_orgin.copy()
                            if go_result_final.shape[0] > 0:
                                go_result_final.sort_values(by="p_value", axis=0, ascending=True, inplace=True)
                                go_address = '/alidata/database/EMTRegulome/emtMotif/search/static/FDRdb/static/analysis_data/' + '%s' % (
                                    data.name) + '.go_result.csv'
                                go_result_final.to_csv(go_address)
                                if go_result_final.shape[0] > 20:
                                    go_result_final_final = go_result_final.iloc[0:20,:]
                                    GO_PATHWAY = list(go_result_final_final._stat_axis.values)
                                    go_data = list(go_result_final_final.loc[:,'-ln p'])
                                    plt.barh(range(len(go_data)), go_data, tick_label=GO_PATHWAY)
                                    plt.subplots_adjust(left=0.5)
                                    plt.yticks(fontsize=7)
                                    plt.title("BP enrichment analysis (GO)")
                                    plt.xlabel('-log10(p value)', fontsize=12)
                                    goimg_name = '%s/%s' % (settings.MEDIA_ROOT, data.name) + '.go.png'
                                    plt.savefig(goimg_name)
                                    plt.close()
                                else:
                                    GO_PATHWAY = list(go_result_final._stat_axis.values)
                                    go_data = list(go_result_final.loc[:,'-ln p'])
                                    plt.barh(range(len(go_data)), go_data, tick_label=GO_PATHWAY)
                                    plt.subplots_adjust(left=0.5)
                                    plt.yticks(fontsize=7)
                                    plt.title("BP enrichment analysis (GO)")
                                    plt.xlabel('-log10(p value)', fontsize=12)
                                    goimg_name = '%s/%s' % (settings.MEDIA_ROOT, data.name) + '.go.png'
                                    plt.savefig(goimg_name)
                                    plt.close()
                                goimg_address = '/static/FDRdb/static/analysis_data/' + '%s' % (data.name) + '.go.png'
                    elif species_select == 'mice':
                        A = data2._stat_axis.values
                        B = mol_info.loc[mol_info.species == 'mus_musculus', 'symbol']
                        C = mol_info.loc[mol_info.species == 'mus_musculus', 'id']
                        result.loc[[i in B for i in result._stat_axis.values], 'fib_relate'] = 'fibrosis'
                        result.loc[[i in C for i in result._stat_axis.values], 'fib_relate'] = 'fibrosis'
                        up_fibrosis = len(which(map(lambda x,y: (x == 'fibrosis') and (y == 'up'), result.loc[:,'fib_relate'],result.loc[:,'GROUP'])))
                        down_fibrosis = len(which(map(lambda x,y: (x == 'fibrosis') and (y == 'down'), result.loc[:,'fib_relate'],result.loc[:,'GROUP'])))
                    elif species_select == 'rat':
                        A = data2._stat_axis.values
                        B = mol_info.loc[mol_info.species == 'rattus_norvegicus', 'symbol']
                        C = mol_info.loc[mol_info.species == 'rattus_norvegicus', 'id']
                        result.loc[[i in B for i in result._stat_axis.values], 'fib_relate'] = 'fibrosis'
                        result.loc[[i in C for i in result._stat_axis.values], 'fib_relate'] = 'fibrosis'
                        up_fibrosis = len(which(map(lambda x,y: (x == 'fibrosis') and (y == 'up'), result.loc[:,'fib_relate'],result.loc[:,'GROUP'])))
                        down_fibrosis = len(which(map(lambda x,y: (x == 'fibrosis') and (y == 'down'), result.loc[:,'fib_relate'],result.loc[:,'GROUP'])))
                up_gene = len(which(result.GROUP == 'up'))
                down_gene = len(which(result.GROUP == 'down'))
                if len(filtered) != 0:
                    ax = sns.scatterplot(x="FoldChange", y="log(pvalue)", hue='GROUP',hue_order=('down', 'none', 'up'),palette=("#377EB8", "grey", "#E41A1C"),data=result)
                    ax.set_ylabel('-log10(p value)', fontweight='bold')
                    ax.set_xlabel('FoldChange', fontweight='bold')
                    volcano_fig = ax.get_figure()
                    volcanoimg_name = '%s/%s' % (settings.MEDIA_ROOT, data.name) + '.histogram_volcano.png'
                    volcano_fig.savefig(volcanoimg_name)
                    plt.close()
                    result_diff_bool = 0
                else:
                    result_diff_bool = 1
                if len(filtered_heat) != 0:
                    # heatmap
                    sns.clustermap(filtered_heat, cmap='RdYlGn_r', standard_scale=0)
                    heatimg_name = '%s/%s' % (settings.MEDIA_ROOT, data.name) + '.histogram_heat.png'
                    plt.savefig(heatimg_name)
                    plt.close()
                    result_diff_heat = 0
                else:
                    result_diff_heat = 1

                #fdimg_address = '/static/FDRdb/static/analysis_data/' + '%s' % (data.name) + '.histogram_fd.png'
                #pimg_address = '/static/FDRdb/static/analysis_data/' + '%s' % (data.name) + '.histogram_p.png'
                volcanoimg_address = '/static/FDRdb/static/analysis_data/' + '%s' % (data.name) + '.histogram_volcano.png'
                heatimg_address = '/static/FDRdb/static/analysis_data/' + '%s' % (data.name) + '.histogram_heat.png'

                result_address = '/alidata/database/EMTRegulome/emtMotif/search/static/FDRdb/static/analysis_data/' + '%s' % (data.name) + '.diff_result.csv'
                result.to_csv(result_address)
                dataname = data.name.split(".")[0]
                all_result_address = '/alidata/database/EMTRegulome/emtMotif/search/static/FDRdb/static/analysis_data/' + '%s' % (dataname) + '.all_result.zip'
                f = zipfile.ZipFile(all_result_address, 'w', zipfile.ZIP_DEFLATED)
                if len(filtered) != 0:
                    for i in ['/alidata/database/EMTRegulome/emtMotif/search'+volcanoimg_address,'/alidata/database/EMTRegulome/emtMotif/search'+heatimg_address,'/alidata/database/EMTRegulome/emtMotif/search'+keggimg_address,kegg_address,'/alidata/database/EMTRegulome/emtMotif/search'+goimg_address,go_address,result_address]:
                        file = i.split('/')[-1]
                        f.write(i, file)
                    f.close()
                else:
                    for i in [result_address]:
                        file = i.split('/')[-1]
                        f.write(i, file)
                    f.close()
                return HttpResponse(loader.get_template("result_diff.html").render(
                    {"result_diff_bool":result_diff_bool,"result_diff_heat":result_diff_heat,"dataname": dataname, "up_fibrosis":up_fibrosis,"down_fibrosis":down_fibrosis,"All_gene_num": data2.shape[0],"Diff_genes_up": up_gene,
                     "Diff_genes_down": down_gene,"number_sample": samples_all,"class1_sample": len(samples1),"class2_sample": len(samples2),
                     "volcanoimg_name": volcanoimg_address, "heatimg_name": heatimg_address, "keggimg_name": keggimg_address, "kegg_index": kegg_index, "goimg_name": goimg_address, "go_index": go_index},request))


def diff_download(request,dataname):
    filename_sub = dataname + '.all_result.zip'
    filename = '/alidata/database/EMTRegulome/emtMotif/search/static/FDRdb/static/analysis_data/' + filename_sub
    wrapper = FileWrapper(open(filename, 'rb'))
    response = HttpResponse(wrapper, content_type='application/octet-stream')
    response['Content-Length'] = os.path.getsize(filename)
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format(filename_sub)

    return response

def expression_example(request):

    filename = "/alidata/database/FDGAdb/data/download/ESCA_mRNA.txt" # Select your file here.
    wrapper = FileWrapper(open(filename, 'rb'))
    response = HttpResponse(wrapper, content_type='application/octet-stream')
    response['Content-Length'] = os.path.getsize(filename)
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format("ESCA_mRNA.txt")

    return response

def class_example(request):

    filename = "/alidata/database/FDGAdb/data/download/ESCA_class.txt" # Select your file here.
    wrapper = FileWrapper(open(filename, 'rb'))
    response = HttpResponse(wrapper, content_type='application/octet-stream')
    response['Content-Length'] = os.path.getsize(filename)
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format("ESCA_class.txt")

    return response