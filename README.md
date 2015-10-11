# seq_explorer

Web application to allow creation of projects to explore and compare DNA and/or protein data.

Comprised of two parts:
* polymer based web application
* python rest api with postgres back end

REST API:
```
/entrez/nucleotide/summary?term=<search_term>&retmax=<retmax> - search for DNA info 
/entrez/nucleotide/summary/<id> - view specific DNA info
/entrez/protein/summary?term=<search_term>&retmax=<retmax> - search for protein info
/entrez/protein/summary/<id> - view specific protein info
/project - list/create project
/project/<project_id> - view/create specific project
/project/<project_id>/dna/<entrez_id> - view genbank dna data part of this project
/project/<project_id>/protein/<entrez_id> - view genbank protein data part of this project
/fetch_request - list/create requests to add dna or protein information to a project
/fetch_request/<request_id> - view request status
```
