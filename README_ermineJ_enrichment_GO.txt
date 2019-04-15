# Enriquecimento GO
## considerações a levar-se em conta
- pegar somente os genes considerados expressos (tabela do edgeR)
- e os que possuem GO atribuído
- aí sim, utilizar o ermineJ
- Score.txt -> atribuir 0 aos genes DE, e definir "quanto menor, melhor"

##############################################################################
## preparo dos arquivos para utilizar no ermineJ ##

- o script (/enrichment_GO/script.py)transforma o arquivo contendo a anotação GO gerada pelo Blast2GO (Ss_full_GO_annot.csv) em outro contendo apenas os genes que apresentam anotação GO (Annotation_file.csv).

# após na pasta /"enrichment prepare ermineJ"/, estão contidos scripts para o preparo da anotação para o enriquecimento GO:

- o script (annotation_expresses.py) gera o arquivo (new_annotation_file.csv) a partir de outro arquivo contendo os genes considerados expressos pela análise no edgeR (new_expresses_genes.csv), sendo este o arquivo da anotação dos genes expressos (new_annotation_file.csv) a ser utilizado no ermineJ;
-o script (score_file.py) gera o arquivo (score_file.txt) que contém os scores para serem utilizados no ermineJ (ORA analysis), onde é atribuído 0 para os genes DE e 1 para os genes expressos, porém não DE.

##################################################################################
## utilizando o ermineJ ##
####
primeiramente é preciso adicionar exceção de segurança para o ermineJ:
-abrir o Control Panel do Java 8:

cd /usr/java/jre1.8.0_201/bin/
./ControlPanel

--aba Security
- marcar 'Very High'
- Edit Site List...
https://home.pavlab.msl.ubc.ca
ou
https://home.pavlab.msl.ubc.ca:443
-Add
-OK
####
executar o arquivo ermineJ.jnlp para abrir a interfácie gráfica do programa:
- Gene Ontology file (downloadable here: archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz) -> baixar e inserir no programa
- Gene annotation file -> inserir o new_annotation_file.csv gerado pelo script python
- Start
#####

-aba Analysis
-- Run Analysis
--- ORA - Over-representation analysis
---- inserir o score_file.txt -> Column: 2
---- Selecionar as 3 categorias GO
---- Definir o tamanho máximo e mínimo do gene para ser considerado ????????????
---- Minimum gene size -> 5; maximum gene size -> 1000
---- desmarcar as duas caixas e definir um Gene score threshold para qualquer valor entre 0 e 1 (0.5 p.e.)
---- Finish

-aba Analysis
-- Save Analysis

Finalizado!
