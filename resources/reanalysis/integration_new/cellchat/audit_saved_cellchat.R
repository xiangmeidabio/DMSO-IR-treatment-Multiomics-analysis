#!/usr/bin/env Rscript
# Read saved original networks; do not install/load CellChat or rerun inference.
main<-function(){
 suppressPackageStartupMessages({library(Matrix);library(jsonlite)})
 args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.',mustWork=TRUE)
 src<-file.path(root,'new','找到一个基因');out<-file.path(root,'reanalysis','integration_new','cellchat')
 td<-file.path(out,'tables');rd<-file.path(out,'resources');dir.create(td,recursive=TRUE,showWarnings=FALSE);dir.create(rd,showWarnings=FALSE)
 groups<-c('Ctrl','DMSO','IR','DMSO-IR');filenames<-setNames(c('Ctrl','DMSO','IR','DMSO_IR'),groups)
 final<-read.csv(file.path(root,'reanalysis','Result1','tables','metadata_harmonized.csv'),row.names=1,check.names=FALSE)
 rna<-read.csv(gzfile(file.path(root,'reanalysis','Result5','communication','tables','candidate_genes_per_cell.csv.gz')),check.names=FALSE)
 genes<-c('Ulbp1','Klrk1','Hcst','Tnf','Tnfrsf1a','Tnfrsf1b','Sell','Cd34','Icam1','Spn')
 ids<-c('ULBP1_KLRK1_HCST','ULBP1_KLRK1','TNF_TNFRSF1A','TNF_TNFRSF1B','SELL_CD34','ICAM1_SPN')
 audits<-maps<-flows<-weights<-positive<-paths<-candidates<-expression<-list();params<-list();dbs<-list();seen<-character()
 for(g in groups){
  file<-file.path(src,paste0(filenames[g],'.cellchat.rds'));message('Reading saved ',g)
  x<-readRDS(file);stopifnot(isS4(x),identical(class(x)[[1]],'CellChat'))
  # Reading S4 attributes directly avoids inventing or registering any CellChat class.
  meta<-attr(x,'meta');data<-attr(x,'data');net<-attr(x,'net');np<-attr(x,'netP');db<-attr(x,'DB');lr<-attr(x,'LR')$LRsig
  barcodes<-rownames(meta);stopifnot(!anyDuplicated(barcodes),identical(colnames(data),barcodes))
  fm<-final[match(barcodes,rownames(final)),];stopifnot(!anyNA(fm$celltype),all(fm$orig.ident==g))
  old<-as.character(meta$labels);new<-as.character(fm$celltype)
  cross<-as.data.frame(table(original_label=old,harmonized_label=new));cross<-cross[cross$Freq>0,];cross$condition<-g
  maps[[g]]<-cross
  mapping<-lapply(split(new,old),unique);stopifnot(all(lengths(mapping)==1L));mapping<-unlist(mapping)
  stopifnot(!anyDuplicated(unname(mapping)))
  dn<-dimnames(net$prob);stopifnot(identical(dim(net$pval),dim(net$prob)),identical(dn[[1]],dn[[2]]))
  cnt<-table(old)
  p_lt<-net$prob;p_lt[net$pval>=.05]<-0
  p_le<-net$prob;p_le[net$pval>.05]<-0
  w_lt<-apply(p_lt,c(1,2),sum);w_le<-apply(p_le,c(1,2),sum)
  e_lt<-max(abs(net$weight-w_lt));e_le<-max(abs(net$weight-w_le))
  c_lt<-apply(p_lt>0,c(1,2),sum);c_le<-apply(p_le>0,c(1,2),sum)
  # Exact aggregate semantics are checked against the stored weight/count.
  rule<-if(e_lt<1e-10 && max(abs(net$count-c_lt))==0) 'pval < 0.05' else if(e_le<1e-10 && max(abs(net$count-c_le))==0) 'pval <= 0.05' else 'unresolved'
  used<-if(rule=='pval < 0.05') p_lt else p_le
  if(rule=='unresolved') stop('Saved aggregate does not match candidate threshold rules: ',g)
  audit<-data.frame(condition=g,source_file=file,n_cells=nrow(meta),n_genes=nrow(data),n_labels=length(cnt),
    exact_final_condition_barcodes=setequal(barcodes,rownames(final)[final$orig.ident==g]),
    original_label_matches_final=sum(old==new),n_label_disagreements=sum(old!=new),
    n_LR_candidates=dim(net$prob)[3],n_saved_pathways=dim(np$prob)[3],
    total_saved_weight=sum(net$weight),total_saved_count=sum(net$count),
    aggregate_threshold=rule,weight_error_p_lt_005=e_lt,weight_error_p_le_005=e_le,
    count_error_p_lt_005=max(abs(net$count-c_lt)),count_error_p_le_005=max(abs(net$count-c_le)))
  audits[[g]]<-audit;params[[g]]<-attr(x,'options');dbs[[g]]<-db;seen<-c(seen,barcodes)
  a<-expand.grid(source_original=dn[[1]],target_original=dn[[2]],stringsAsFactors=FALSE)
  a$condition<-g;a$source<-unname(mapping[a$source_original]);a$target<-unname(mapping[a$target_original]);a$source_n<-as.integer(cnt[a$source_original]);a$target_n<-as.integer(cnt[a$target_original]);a$saved_weight<-as.vector(net$weight);a$saved_count<-as.vector(net$count)
  weights[[g]]<-a
  ix<-which(net$prob>0,arr.ind=TRUE)
  a<-data.frame(condition=g,source_original=dn[[1]][ix[,1]],target_original=dn[[2]][ix[,2]],interaction=dn[[3]][ix[,3]],saved_probability=net$prob[ix],saved_permutation_p=net$pval[ix],included_in_saved_aggregate=used[ix]>0)
  a$source<-unname(mapping[a$source_original]);a$target<-unname(mapping[a$target_original]);m<-match(a$interaction,lr$interaction_name);a$pathway<-lr$pathway_name[m];a$ligand<-lr$ligand[m];a$receptor<-lr$receptor[m];positive[[g]]<-a
  pnames<-dimnames(np$prob)[[3]]
  a<-expand.grid(source_original=dimnames(np$prob)[[1]],target_original=dimnames(np$prob)[[2]],pathway=pnames,stringsAsFactors=FALSE)
  a$condition<-g;a$source<-unname(mapping[a$source_original]);a$target<-unname(mapping[a$target_original]);a$saved_pathway_probability<-as.vector(np$prob);paths[[g]]<-a
  flows[[g]]<-data.frame(condition=g,pathway=pnames,saved_information_flow=as.numeric(apply(np$prob,3,sum)))
  # Full candidate arrays include zeros and candidates not selected in this group's model.
  can<-list()
  for(id in ids)for(s in dn[[1]])for(t in dn[[2]]){
   present<-id %in% dn[[3]];pr<-if(present)net$prob[s,t,id] else NA_real_;pv<-if(present)net$pval[s,t,id] else NA_real_
   d<-db$interaction[match(id,db$interaction$interaction_name),]
   can[[length(can)+1L]]<-data.frame(condition=g,source_original=s,target_original=t,source=unname(mapping[s]),target=unname(mapping[t]),
    source_n=as.integer(cnt[s]),target_n=as.integer(cnt[t]),interaction=id,pathway=d$pathway_name,
    ligand=d$ligand,receptor=d$receptor,selected_in_saved_LR_candidates=present,saved_probability=pr,saved_permutation_p=pv,
    included_in_saved_aggregate=if(present)used[s,t,id]>0 else FALSE,
    source_or_target_below10=cnt[s]<10 | cnt[t]<10)
  }
  candidates[[g]]<-do.call(rbind,can)
  rr<-rna[match(barcodes,rna$barcode),];stopifnot(identical(rr$barcode,barcodes))
  er<-list()
  for(ge in genes){
   z<-as.numeric(data[ge,]);current<-log1p(rr[[paste0(ge,'_CP10k')]])
   for(ty in unique(new)){
    ii<-which(new==ty);q<-unname(quantile(z[ii],c(.25,.5,.75),type=7));
    er[[length(er)+1L]]<-data.frame(condition=g,celltype=ty,gene=ge,n_cells=length(ii),saved_data_mean=mean(z[ii]),
     saved_data_percent_nonzero=100*mean(z[ii]>0),saved_data_Q1=q[1],saved_data_median=q[2],saved_data_Q3=q[3],
     saved_data_triMean=.25*q[1]+.5*q[2]+.25*q[3],max_abs_difference_from_current_log1p_CP10k=max(abs(z[ii]-current[ii])))
   }
  }
  expression[[g]]<-do.call(rbind,er)
  write.csv(lr,file.path(td,paste0(filenames[g],'_saved_LR_definitions.csv')),row.names=FALSE)
  rm(x,data);invisible(gc())
  message(g,': weight=',format(audit$total_saved_weight,digits=7),' count=',audit$total_saved_count,' label mismatches=',audit$n_label_disagreements)
 }
 stopifnot(length(seen)==23006,!anyDuplicated(seen),setequal(seen,rownames(final)))
 write_tab<-function(x,name)write.csv(x,file.path(td,name),row.names=FALSE,na='')
 combine<-function(z)do.call(rbind,z)
 au<-combine(audits);fl<-combine(flows);we<-combine(weights);ca<-combine(candidates);ex<-combine(expression)
 write_tab(au,'saved_object_audit.csv');write_tab(combine(maps),'original_to_harmonized_label_mapping.csv');write_tab(we,'saved_celltype_network_weight_count.csv');write_tab(fl,'saved_pathway_information_flow_long.csv');write_tab(ca,'saved_candidate_LR_all_sources_targets.csv');write_tab(ca[ca$source=='HSC'|ca$target=='HSC',],'saved_HSC_candidate_LR.csv');write_tab(ex,'saved_input_expression_vs_current.csv')
 for(pair in list(list(combine(positive),'saved_LR_positive_probabilities.csv.gz'),list(combine(paths),'saved_pathway_probabilities_all_edges.csv.gz'))){con<-gzfile(file.path(td,pair[[2]]),'wt');write.csv(pair[[1]],con,row.names=FALSE);close(con)}
 allpaths<-sort(unique(fl$pathway));wide<-data.frame(pathway=allpaths)
 for(g in groups){v<-fl$saved_information_flow[match(paste(allpaths,g),paste(fl$pathway,fl$condition))];v[is.na(v)]<-0;wide[[g]]<-v}
 wide$IR_minus_Ctrl<-wide$IR-wide$Ctrl;wide$DMSO_IR_minus_IR<-wide[['DMSO-IR']]-wide$IR;wide$DMSO_IR_minus_Ctrl<-wide[['DMSO-IR']]-wide$Ctrl;wide$DMSO_minus_Ctrl<-wide$DMSO-wide$Ctrl
 wide$original_direction_rule<-wide$IR_minus_Ctrl>0 & wide$DMSO_IR_minus_IR<0
 wide$DMSO_IR_closer_to_Ctrl<-abs(wide$DMSO_IR_minus_Ctrl)<abs(wide$IR_minus_Ctrl)
 wide$proximity_gain<-abs(wide$IR_minus_Ctrl)-abs(wide$DMSO_IR_minus_Ctrl)
 write_tab(wide,'saved_pathway_information_flow_four_conditions.csv');write_tab(wide[wide$original_direction_rule,],'pathways_matching_original_direction_rule.csv')
 wb<-read.csv(file.path(rd,'original_Overlap_workbook_rows.csv'),check.names=FALSE);overlap<-wb$E[!wb$E %in% c('','#N/A','Overlap')]
 ov<-data.frame(pathway=sort(unique(c(overlap,wide$pathway[wide$original_direction_rule]))));ov$in_original_workbook_overlap<-ov$pathway %in% overlap;ov$in_saved_network_direction_rule<-ov$pathway %in% wide$pathway[wide$original_direction_rule];write_tab(ov,'original_workbook_overlap_vs_saved_network.csv')
 write_tab(dbs[[1]]$interaction[dbs[[1]]$interaction$interaction_name %in% ids,],'saved_database_candidate_definitions.csv')
 complex<-dbs[[1]]$complex;write_tab(data.frame(complex=rownames(complex),complex)[rownames(complex) %in% c('KLRK1_HCST','NKG2D_HCST'),],'saved_NKG2D_complex_definition.csv')
 saveRDS(dbs[[1]],file.path(rd,'original_saved_CellChatDB.mouse.rds'))
 write_json(params,file.path(rd,'saved_object_options.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
 current<-new.env();load(file.path(root,'reanalysis','Result5','communication','resources','CellChatDB.mouse.rda'),envir=current)
 audit<-list(analysis='Extraction and descriptive comparison of original saved CellChat outputs; inference was not rerun.',
  all_cells_match_final_23006=TRUE,conditions=groups,DB_identical_all_four=vapply(dbs,function(d)identical(d,dbs[[1]]),logical(1)),
  DB_identical_to_previously_downloaded_snapshot=identical(dbs[[1]],current$CellChatDB.mouse),
  saved_DB_interactions=nrow(dbs[[1]]$interaction),saved_DB_complexes=nrow(dbs[[1]]$complex),
  original_direction_rule_n=sum(wide$original_direction_rule),original_workbook_overlap_n=length(overlap),
  original_rule_and_closer_n=sum(wide$original_direction_rule & wide$DMSO_IR_closer_to_Ctrl),
  expression_max_abs_difference_from_current=max(ex$max_abs_difference_from_current_log1p_CP10k),
  CellChat_version='Not encoded in saved object; manuscript reports v2.2.0 and original script installs GitHub without a fixed version.',
  no_new_CellChat_inference=TRUE,no_animal_level_tests=TRUE,completed_at=format(Sys.time(),tz='UTC',usetz=TRUE))
 write_json(audit,file.path(out,'audit_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
 writeLines(capture.output(sessionInfo()),file.path(out,'extraction_sessionInfo.txt'))
 print(au);print(wide[wide$pathway %in% c('TNF','NKG2D','SELL','ICAM','VCAM','NOTCH','COMPLEMENT'),]);print(audit)
}
main()
