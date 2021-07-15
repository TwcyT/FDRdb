from __future__ import unicode_literals

from django.db import models

# Create your models here.

class fibrosis_gene(models.Model):
	tissue = models.CharField(max_length=100)
	disease = models.CharField(max_length=100)
	species = models.CharField(max_length=100)
	gene_symbol = models.CharField(max_length=100)
	gene_id = models.CharField(max_length=100)
	source = models.CharField(max_length=100,null=True)
	type = models.CharField(max_length=100)
	expression = models.CharField(max_length=100)
	target = models.CharField(max_length=1000,null=True)
	pathway = models.CharField(max_length=100,null=True)
	method = models.CharField(max_length=1000,null=True)
	drug = models.CharField(max_length=100,null=True)
	function = models.CharField(max_length=10000,null=True)
	pmid = models.CharField(max_length=100,null=True)
	title = models.CharField(max_length=500,null=True)
	link = models.CharField(max_length=200,null=True)
	year = models.CharField(max_length=100,null=True)

class fibrosis_data(models.Model):
	database = models.CharField(max_length=100)
	data_id = models.CharField(max_length=100)
	type = models.CharField(max_length=100)
	species = models.CharField(max_length=100)
	disease = models.CharField(max_length=100)
	platform = models.CharField(max_length=100)
	title = models.CharField(max_length=500)
	link = models.CharField(max_length=200,null=True)
	pmid = models.CharField(max_length=100,null=True)
	pmid_link = models.CharField(max_length=1000,null=True)
	tissue = models.CharField(max_length=100,null=True)
	sample_num = models.CharField(max_length=100,null=True)
	describ = models.CharField(max_length=10000,null=True)
	keyword = models.CharField(max_length=100,null=True)

