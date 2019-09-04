First_using=FALSE
if(First_using==TRUE) {
        install.packages("BiocManager")
        library(BiocManager)
        #BiocManager::install("STRINGdb")
        install.packages("igraph")
        install.packages(c("boot", "cluster", "MASS", "Matrix", "mgcv", "nlme", "rpart", "survival"))
        BiocManager::install("biomaRt")
}

#########  Loading required database ##############
library(BiocManager)
library(igraph)

#########  Prepare id converting dictionary ##############
library(gdata)
ID_dic<-read.xls("mart_export.xlsx",sheet = 1)
ID_dic2<-read.xls("mart_export2.xlsx",sheet = 1)


#########  Loading potential Targets  #########

# Load predicted targets and targets realted to macrophage polarization
M_targets <- read.xls("Network_new.xlsx",sheet = 2)
C_targets <- read.xls("Network_new.xlsx",sheet = 1)

M_target<-as.matrix(M_targets$stringId)
M_target<-M_target[-which(duplicated(M_target))]
C_target<-as.matrix(C_targets$stringId)
C_target<-C_target[-which(duplicated(C_target))]

M_targets2<-sub("9606.","",M_target)
C_targets2<-sub("9606.","",C_target)
map<-match(c(C_targets2,M_targets2),ID_dic[,4])
MC_target<-cbind(ID_dic)[map,2]
MC_target<-MC_target[-which(is.na(MC_target))]
MC_target_no_dup<-unique(MC_target)

#########  Prepare database for protein-protein interactions  #########

#Prepare stringDB and remove duplications
stringdb<-read.table(gzfile("./9606.protein.links.v11.0.txt.gz"))
stringdb<-as.matrix(stringdb)
dup4<-which(as.numeric(stringdb[,3]) < 950)
head(dup4) # do not remove anyitem because no duplication was found
stringdb3<-stringdb[-dup4,c(1,2)]

stringdb2_node1<-sub("9606.","",stringdb3[,1])
stringdb2_node2<-sub("9606.","",stringdb3[,2])

map_node1<-match(stringdb2_node1,ID_dic[,4])
map_node2<-match(stringdb2_node2,ID_dic[,4])
stringdb2<-cbind(ID_dic[map_node1,2],ID_dic[map_node2,2])


#Prepare HPRD 
hprddb<-read.table("HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt",sep="\t")
hprddb2<-as.matrix(hprddb[,c(2,5)])
dup3<-which(duplicated(hprddb2))
hprddb2<-hprddb2[-dup3,]

#Prepare HIPPIE database
hippiedb<-read.table("hippie_current.txt",sep="\t")
hippiedb2<-as.matrix(hippiedb[,c(2,4)])
dup3<-which(duplicated(hippiedb2))
hippiedb2<-hippiedb2[-dup3,]


#prepare reactomeDB and remove duplications
reactomedb<-read.table("reactome.homo_sapiens.interactions.tab-delimited.txt",sep = "\t")
dup<-which(sapply(c(1:dim(reactomedb)[1]), FUN = function(x){match(reactomedb[x,3],reactomedb[x,6])})==1)
reactomedb<-reactomedb[-dup,]
dup2<-which(sapply(c(1:dim(reactomedb)[1]), FUN = function(x){!all(is.na(match(as.matrix(reactomedb[x,]),"-")))}))
reactomedb<-reactomedb[-dup2,]
reactomedb2<-sub("entrezgene/locuslink:","",as.matrix(reactomedb[,c(3,6)]))
a<-grep("entrezgene/locuslink:",c(reactomedb2[,1]))
b<-grep("entrezgene/locuslink:",c(reactomedb2[,2]))
reactomedb2<-reactomedb2[-unique(c(a,b)),]
dup3<-which(duplicated(reactomedb2))
reactomedb2<-reactomedb2[-dup3,]

Interaction_db<-rbind(reactomedb2,stringdb2,hippiedb2,hprddb2)

#########  Construct protein-protein interaction network by database  #########


# Consruct the network using other databases
network_constructor<-function(database,target) {
        map<-which(!is.na(match(database[,1],target)))
        col1<-rep(target,length(map))
        col2<-database[map,2]
        return(cbind(col1,col2))
}

MC_network2<-sapply(MC_target_no_dup,network_constructor,database=Interaction_db)
MC_network2<-do.call(rbind, MC_network2)
dup5<-which(duplicated(MC_network2))
MC_network_nodup<-MC_network2[-dup5,]
remove_na<-which(is.na(MC_network_nodup[,2]))
MC_network_nodup<-MC_network_nodup[-remove_na,]
rownames(MC_network_nodup)<-NULL
nodes<-c(MC_network_nodup[,2],MC_network_nodup[,1])
nodes<-nodes[-which(duplicated(nodes))]
dim(MC_network_nodup)[1]

#########  Export the network  #########
rownames(MC_network_nodup)<-NULL
export_node1<-MC_network_nodup[,1]
export_node2<-MC_network_nodup[,2]


map_node1<-as.character(ID_dic2[match(export_node1,ID_dic2[,2]),5])
map_node2<-as.character(ID_dic2[match(export_node2,ID_dic2[,2]),5])
export_network<-cbind(map_node1,map_node2)
rownames(export_network)
export_network<-export_network[-which(is.na(export_network[,2])),]

write.table(export_network,file = "initial_network.csv",quote = FALSE,row.names = FALSE)


#########  Topology analysis (performed on  Cytoscape) #########

# The results are listed below:
# 1. Initial Network: 5920 nodes and 15068 Interactions
# 2. Hub node network: 597 nodes and 5516 Interactions
# 3. Main candidiate network: 189 nodes and 186 Interactions
# The main candidate network is exported as main_nodes.csv
# Enrichment analysis using String DB

#########  Connecting compounds and identified main candidates  #########
main_nodes<-read.csv("main_candidate_node.csv",sep = ",")
compounds_list=as.character(unique(C_targets$Compounds))

CR_list<-as.matrix(C_targets[,c(1,3)])
map<-match(sub(pattern = "9606.","",CR_list[,2]),ID_dic2[,4])
CR_list[,2]<-cbind(as.character(ID_dic2[map,5]))
CR_list<-CR_list[-c(48,84),]


make_node_list<-function(node_list=main_nodes,target_list=CR_list, list_caption=compounds_list[1]) {
        id<-match("name",colnames(main_nodes))
        node_list<-as.matrix(node_list)
        colnames(node_list)<-NULL
        R_list<-target_list[which(!is.na(match(target_list[,1], list_caption))),2]
        matched_node<-node_list[match(R_list,node_list[,id],nomatch = 0),id]
        compound_node<-rep(list_caption,length(matched_node))
        return(cbind(compound_node,matched_node))
}
node_list<-lapply(compounds_list, make_node_list, node_list=main_nodes, target_list=CR_list)
CR_network<-do.call(rbind, node_list)
colnames(CR_network)<-NULL
CR_network<-cbind(CR_network,rep(x = "Compounds",dim(CR_network)[1]))

#########  Construct Compound-main candidate-enrichment network #########
main_edge<-read.csv("main_candidate_edge.csv",sep = ",",stringsAsFactors = FALSE)
main_network<-lapply(main_edge$name, function(x){return(t(strsplit(x," ")[[1]][c(1,4)]))})
main_network<-do.call(rbind,main_network)
main_network<-cbind(main_network,rep(x = "Macrophages",dim(main_network)[1]))

enrich_nodes<-read.csv("enrich_node.csv",sep = ",", stringsAsFactors = FALSE)
enrich_netwrok_construct<-function(x=enrich_nodes,y=1){
        node_asso2<-gsub("\\[|\\]", "", x$Associated.Genes.Found)[y]
        node_asso<-strsplit(node_asso2,", ")[[1]]
        return_vec<-cbind(node_asso,rep(x$GOTerm[y],length(node_asso)))
        colnames(return_vec)<-NULL
        return(return_vec)
}
enrich_network<-lapply(c(1:dim(enrich_nodes)[1]), enrich_netwrok_construct,x=enrich_nodes)
enrich_network<-do.call(rbind,enrich_network)
enrich_network<-cbind(enrich_network,rep(x = "Enrichment",dim(enrich_network)[1]))



# Output whole network

whole_network<-rbind(CR_network,main_network,enrich_network)
write.csv(whole_network,file = "CRE_network.csv",quote = FALSE,row.names = FALSE,col.names = FALSE)

#########  Analysis and Visualizatiion (performed on  Cytoscape) #########
# Using Cluego for enrichment analysis
# Visualizing the network using different style in Cytoscape --> Final.cys