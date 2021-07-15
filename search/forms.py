from django import forms


class SearchForm(forms.Form):
    search_species = forms.CharField(label="search_species", max_length=100,required=True)
    search_tissues = forms.CharField(label="search_tissues", max_length=100,required=True)
    filter = forms.CharField(label="filter", max_length=100, required=True)


class SearchspeciesForm(forms.Form):
    species_select = forms.CharField(label="species_select",max_length=100,required=False)
