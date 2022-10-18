# Temporary directory.
tmp.dir <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
tmp.file <- normalizePath(tempfile(), winslash='/', mustWork=FALSE)

# Bitmap extensions accepted in uploaded images.
raster.ext <- c('.jpg', '.JPG', '.png', '.PNG')

# Confirm button labels.
lab.sgl <- 'Search by single gene ID (e.g. ENSMUSG00000000031) or symbols'
lab.mul <- 'Search by single or multiple gene IDs (e.g. ENSMUSG00000000031 ENSMUSG00000000093)'

# Search portion of URLs on the landing page.
brain.hum.url <- '?_inputs_&deg-datDEG-P=0&deg-edg.lim.nor=%22CNF-TMM%22&dat-CV1=-10000&shmAll-net-dpwColNet=0&shmAll-shmMhNet=%22shm1%22&shmAll-togDrop=1&upl-svgInpath1=null&shmAll-col.n=2&deg-ssg.update=0&deg-datDEG-CV1=-10000&shmAll-net-dpwNetTar=0&shmAll-match=0&upl-dimName=%22None%22&shmAll-net-measure=%22correlation%22&shmAll-net-ds=%223%22&shmAll-genCon=%22gene%22&deg-rok.dis.nor=%22CNF-TMM%22&shmAll-val.lgd.text=10&shmAll-color=%22yellow%2Corange%2Cred%22&upl-target=null&shmAll-net-dpwNetType=0&shmAll-t=2&shmAll-lgd.incld=%22Yes%22&shmAll-disDrop=0&sear-ids.but=1&deg-ssg.fdr=0.05&shmAll-net-mhmNav=%22mhmPlot%22&shmAll-scale.shm=1.4&dat-dtSel_cell_clicked=%7B%7D&deg-deg-ids.but=0&shmAll-net-dpwNetAdj=0&dat-scale=%22Row%22&shmAll-ggly.but=0&shmAll-vdo.itvl=1&sidebarItemExpanded=null&dat-CV2=10000&shmAll-val.lgd.feat=%22No%22&upl-svgInpath2=null&shmAll-scaleDrop=0&shmAll-net-mhm.but=0&shmAll-net-dpwModSize=0&shmAll-net-cpt.nw=0&shmAll-scroDrop=0&shmAll-net-color.net=%22yellow%2Corange%2Cred%22&dat-log=%22No%22&deg-ssg.fc=1&dat-dtSel_state=null&shmAll-interNav=%22interPlot%22&shmAll-title.size=12&deg-datDEG-CV2=10000&shmAll-net-thr=%22p%22&shmAll-net-col.but.net=0&shmAll-togDrop_state=false&shmAll-lgd.size=0.5&dat-dtSel_search=%22%22&shmAll-scale.ly=1&shmAll-val.lgd.row=1&shmAll-vdo.val.lgd=%22No%22&shmAll-net-cor.abs=%22No%22&upl-cusHelp=0&shmAll-vdo.dim=%22640x480%22&upl-fileIn=%22brain_Prudencio%22&shmAll-net-mhm.v=0.2&shmAll-val.lgd.key=0.03&shmAll-vdo.key.row=2&shmAll-lgd.row=2&shmAll-vdo.key.size=0.04&deg-datDEG-fil.but=0&shmAll-dropdown=0&dat-fil.but=0&sear-ids.in=%22ENSG00000000971%20CFH%3A%20complement%20factor%20H%22&shmAll-cs.v=%22Selected%20rows%22&shmAll-net-mat.scale=%22Row%22&shmAll-colDrop=0&shmAll-ext=%22NA%22&shmAll-fs=0&shmAll-vdo.but=0&shmAll-line.size=0.1&shmAll-net-netNav=%22netPlot%22&shmAll-net-min.size=15&shmAll-net-gen.sel=%22ENSG00000000971%22&upl-geneInpath=null&shmAll-line.color=%22grey70%22&deg-sam.con=%22feature%22&dat-A=0&dat-dtSel_columns_selected=null&upl-tar=null&dat-P=0&shmAll-titDrop=0&deg-ssg.meth=%5B%22edgeR%22%2C%22limma%22%5D&shm.sup=%22shmPanelAll%22&upl-met=null&shmAll-tis=null&shmAll-net-net.type=%22signed%22&shmAll-net-adj.in=%220.906%22&shmAll-vdo.label=%22No%22&deg-deg-ids.in=null&deg-datDEG-A=0&dat-dtSel_cells_selected=%5B%5D&shmAll-lgd.key.size=0.04&shmAll-lgd.lab.size=2.5&shmAll-togSld=0.7&shmAll-vdo.res=400&shmAll-net-max.edg=50&shmAll-col.but=0&shmAll-lgd.label=%22No%22&shmAll-val.lgd=0&dat-dtSel_rows_selected=null&sidebarCollapsed=true&shmAll-vdo.lab.size=2&shmAll-scrollH=450&upl-config=null&shmAll-relaSize=null&shmAll-res=300&right.bar=true'

mouse.url <- '?_inputs_&upl-fileIn=%22mouse_Merkin%22&shmAll-net-mhm.v=0.2&deg-datDEG-P=0&shmAll-val.lgd.key=0.03&shmAll-vdo.key.row=2&deg-edg.lim.nor=%22CNF-TMM%22&shmAll-lgd.row=3&dat-CV1=-10000&shmAll-vdo.key.size=0.04&shmAll-togDrop=1&deg-datDEG-fil.but=0&upl-svgInpath1=null&shmAll-col.n=3&shmAll-dropdown=0&dat-fil.but=0&deg-ssg.update=0&shmAll-cs.v=%22Selected%20rows%22&shmAll-net-mat.scale=%22Row%22&shmAll-colDrop=4&shmAll-ext=%22NA%22&sear-ids.in=%22ENSMUSG00000000031%20H19%3A%20H19%2C%20imprinted%20maternally%20expressed%20transcript%22&shmAll-match=0&shmAll-fs=0&upl-dimName=%22None%22&shmAll-vdo.but=0&shmAll-genCon=%22gene%22&shmAll-net-measure=%22correlation%22&shmAll-net-ds=%223%22&deg-rok.dis.nor=%22CNF-TMM%22&upl-target=null&shmAll-val.lgd.text=10&shmAll-line.size=0.1&shmAll-t=2&shmAll-net-min.size=15&shmAll-lgd.incld=%22Yes%22&deg-datDEG-CV1=-10000&shmAll-color=%22yellow%2Corange%2Cred%22&shmAll-disDrop=0&sear-ids.but=0&shmAll-net-gen.sel=%22ENSMUSG00000000031%22&upl-geneInpath=null&shmAll-line.color=%22grey70%22&deg-sam.con=%22feature%22&dat-A=0&shmAll-scale.shm=1&upl-tar=null&deg-ssg.fdr=0.05&dat-dtSel_columns_selected=null&dat-dtSel_cell_clicked=%7B%7D&dat-scale=%22Row%22&shmAll-ggly.but=0&shmAll-vdo.itvl=1&sidebarItemExpanded=null&dat-P=0&shmAll-titDrop=0&deg-ssg.meth=%5B%22edgeR%22%2C%22limma%22%5D&dat-CV2=10000&shm.sup=%22shmPanelAll%22&upl-met=null&upl-svgInpath2=null&shmAll-scaleDrop=5&shmAll-net-mhm.but=0&shmAll-net-cpt.nw=0&shmAll-scroDrop=0&shmAll-tis=null&shmAll-val.lgd.feat=%22No%22&dat-log=%22No%22&shmAll-vdo.label=%22No%22&shmAll-net-net.type=%22signed%22&shmAll-net-adj.in=%220.896%22&shmAll-togSld=0.75&shmAll-vdo.res=400&deg-datDEG-A=0&shmAll-net-max.edg=50&shmAll-col.but=0&deg-ssg.fc=1&shmAll-title.size=12&shmAll-net-color.net=%22yellow%2Corange%2Cred%22&deg-datDEG-CV2=10000&shmAll-lgd.key.size=0.04&shmAll-lgd.lab.size=2.5&shmAll-lgd.label=%22No%22&dat-dtSel_cells_selected=%5B%5D&dat-dtSel_state=null&shmAll-scaleDrop_state=false&shmMhNet=%22shm1%22&shmAll-val.lgd=0&shmAll-net-col.but.net=0&shmAll-net-thr=%22p%22&dat-dtSel_rows_selected=null&shmAll-lgd.size=0.5&dat-dtSel_search=%22%22&shmAll-colDrop_state=false&shmAll-togDrop_state=false&sidebarCollapsed=true&shmAll-scale.ly=1&shmAll-val.lgd.row=1&shmAll-vdo.val.lgd=%22No%22&shmAll-scrollH=450&upl-config=null&shmAll-relaSize=null&shmAll-net-cor.abs=%22No%22&shmAll-res=300&upl-cusHelp=0&shmAll-vdo.dim=%22640x480%22&right.bar=true&shmAll-vdo.lab.size=2'

chicken.url <- '?_inputs_&upl-fileIn=%22chicken_Cardoso.Moreira%22&shmAll-net-mhm.v=0.2&deg-datDEG-P=0&shmAll-val.lgd.key=0.03&shmAll-vdo.key.row=2&deg-edg.lim.nor=%22CNF-TMM%22&shmAll-lgd.row=3&dat-CV1=-10000&shmAll-vdo.key.size=0.04&shmAll-togDrop=1&deg-datDEG-fil.but=0&upl-svgInpath1=null&shmAll-col.n=5&shmAll-dropdown=0&dat-fil.but=0&shmAll-disDrop_state=false&deg-ssg.update=0&shmAll-cs.v=%22Selected%20rows%22&shmAll-net-mat.scale=%22Row%22&shmAll-colDrop=2&shmAll-ext=%22NA%22&sear-ids.in=%22ENSGALG00000000059%20%22&shmAll-match=0&shmAll-fs=0&upl-dimName=%22None%22&shmAll-vdo.but=0&shmAll-genCon=%22gene%22&shmAll-net-measure=%22correlation%22&shmAll-net-ds=%223%22&deg-rok.dis.nor=%22CNF-TMM%22&upl-target=null&shmAll-val.lgd.text=10&shmAll-line.size=0.1&shmAll-t=2&shmAll-net-min.size=15&shmAll-lgd.incld=%22Yes%22&deg-datDEG-CV1=-10000&shmAll-color=%22yellow%2Corange%2Cred%22&shmAll-disDrop=1&sear-ids.but=0&shmAll-net-gen.sel=%22ENSGALG00000000059%22&upl-geneInpath=null&shmAll-line.color=%22grey70%22&deg-sam.con=%22feature%22&dat-A=0&shmAll-scale.shm=0.6&upl-tar=null&deg-ssg.fdr=0.05&dat-dtSel_columns_selected=null&dat-dtSel_cell_clicked=%7B%7D&deg-deg-ids.but=0&dat-scale=%22Row%22&shmAll-ggly.but=0&shmAll-vdo.itvl=1&sidebarItemExpanded=null&dat-P=0&shmAll-titDrop=0&deg-ssg.meth=%5B%22edgeR%22%2C%22limma%22%5D&dat-CV2=10000&shm.sup=%22shmPanelAll%22&upl-met=null&upl-svgInpath2=null&shmAll-scaleDrop=1&shmAll-net-mhm.but=0&shmAll-net-cpt.nw=0&shmAll-scroDrop=0&shmAll-tis=null&shmAll-val.lgd.feat=%22No%22&dat-log=%22No%22&shmAll-vdo.label=%22No%22&shmAll-net-net.type=%22signed%22&shmAll-net-adj.in=%220.896%22&deg-deg-ids.in=null&shmAll-togSld=0.75&shmAll-vdo.res=400&shmAll-net-max.edg=50&shmAll-col.but=0&deg-datDEG-A=0&shmAll-title.size=12&deg-ssg.fc=1&deg-datDEG-CV2=10000&shmAll-net-color.net=%22yellow%2Corange%2Cred%22&shmAll-lgd.key.size=0.04&shmAll-lgd.label=%22No%22&shmAll-lgd.lab.size=2.5&dat-dtSel_cells_selected=%5B%5D&dat-dtSel_state=null&shmMhNet=%22shm1%22&shmAll-val.lgd=0&shmAll-net-col.but.net=0&shmAll-net-thr=%22p%22&dat-dtSel_rows_selected=null&shmAll-lgd.size=0.5&dat-dtSel_search=%22%22&shmAll-colDrop_state=false&shmAll-togDrop_state=false&sidebarCollapsed=true&shmAll-scale.ly=1&shmAll-val.lgd.row=1&shmAll-scaleDrop_state=false&shmAll-vdo.val.lgd=%22No%22&shmAll-scrollH=450&upl-config=null&shmAll-relaSize=null&shmAll-net-cor.abs=%22No%22&shmAll-res=300&upl-cusHelp=0&shmAll-vdo.dim=%22640x480%22&right.bar=true&shmAll-vdo.lab.size=2'

organ.arab.url <- '?_inputs_&deg-datDEG-P=0&deg-edg.lim.nor=%22CNF-TMM%22&dat-CV1=-10000&shmAll-net-dpwColNet=0&shmAll-shmMhNet=%22shm1%22&shmAll-togDrop=1&upl-svgInpath1=null&shmAll-col.n=2&deg-ssg.update=0&deg-datDEG-CV1=-10000&shmAll-net-dpwNetTar=0&shmAll-match=0&upl-dimName=%22None%22&shmAll-net-measure=%22correlation%22&shmAll-net-ds=%223%22&shmAll-genCon=%22gene%22&deg-rok.dis.nor=%22CNF-TMM%22&shmAll-val.lgd.text=10&shmAll-color=%22yellow%2Corange%2Cred%22&upl-target=null&shmAll-net-dpwNetType=0&shmAll-t=2&shmAll-lgd.incld=%22Yes%22&shmAll-disDrop=0&sear-ids.but=1&deg-ssg.fdr=0.05&shmAll-net-mhmNav=%22mhmPlot%22&shmAll-scale.shm=1.4&dat-dtSel_cell_clicked=%7B%7D&deg-deg-ids.but=0&shmAll-net-dpwNetAdj=0&dat-scale=%22Row%22&shmAll-ggly.but=0&shmAll-vdo.itvl=1&sidebarItemExpanded=null&dat-CV2=10000&shmAll-val.lgd.feat=%22No%22&upl-svgInpath2=null&shmAll-scaleDrop=0&shmAll-net-mhm.but=0&shmAll-net-dpwModSize=0&shmAll-net-cpt.nw=0&shmAll-scroDrop=0&shmAll-net-color.net=%22yellow%2Corange%2Cred%22&dat-log=%22No%22&deg-ssg.fc=1&dat-dtSel_state=null&shmAll-interNav=%22interPlot%22&shmAll-title.size=12&deg-datDEG-CV2=10000&shmAll-net-thr=%22p%22&shmAll-net-col.but.net=0&shmAll-togDrop_state=true&shmAll-lgd.size=0.5&dat-dtSel_search=%22%22&shmAll-scale.ly=1&shmAll-val.lgd.row=1&shmAll-vdo.val.lgd=%22No%22&shmAll-net-cor.abs=%22No%22&upl-cusHelp=0&shmAll-vdo.dim=%22640x480%22&upl-fileIn=%22organ_Mustroph%22&shmAll-net-mhm.v=0.2&shmAll-val.lgd.key=0.03&shmAll-vdo.key.row=2&shmAll-lgd.row=2&shmAll-vdo.key.size=0.04&deg-datDEG-fil.but=0&shmAll-dropdown=0&dat-fil.but=0&sear-ids.in=%22NDHE%20NADH%20dehydrogenase%20ND4L%22&shmAll-cs.v=%22Selected%20rows%22&shmAll-net-mat.scale=%22Row%22&shmAll-colDrop=0&shmAll-ext=%22NA%22&shmAll-fs=0&shmAll-vdo.but=0&shmAll-line.size=0.1&shmAll-net-netNav=%22netPlot%22&shmAll-net-min.size=15&shmAll-net-gen.sel=%22NDHE%22&upl-geneInpath=null&shmAll-line.color=%22grey70%22&deg-sam.con=%22feature%22&dat-A=0&dat-dtSel_columns_selected=null&upl-tar=null&dat-P=0&shmAll-titDrop=0&deg-ssg.meth=%5B%22edgeR%22%2C%22limma%22%5D&shm.sup=%22shmPanelAll%22&upl-met=null&shmAll-tis=null&shmAll-net-net.type=%22signed%22&shmAll-net-adj.in=%220.906%22&shmAll-vdo.label=%22No%22&deg-deg-ids.in=null&deg-datDEG-A=0&shmAll-lgd.key.size=0.04&shmAll-lgd.lab.size=2.5&dat-dtSel_cells_selected=%5B%5D&shmAll-togSld=0.7&shmAll-vdo.res=400&shmAll-net-max.edg=50&shmAll-col.but=0&shmAll-lgd.label=%22No%22&shmAll-val.lgd=0&dat-dtSel_rows_selected=null&sidebarCollapsed=true&shmAll-vdo.lab.size=2&shmAll-scrollH=450&upl-config=null&shmAll-relaSize=null&shmAll-res=300&right.bar=true'

shoot.arab.url <- '?_inputs_&deg-datDEG-P=0&deg-edg.lim.nor=%22CNF-TMM%22&dat-CV1=-10000&shmAll-net-dpwColNet=0&shmAll-shmMhNet=%22shm1%22&shmAll-togDrop=1&upl-svgInpath1=null&shmAll-col.n=2&deg-ssg.update=0&deg-datDEG-CV1=-10000&shmAll-net-dpwNetTar=0&shmAll-match=0&upl-dimName=%22None%22&shmAll-net-measure=%22correlation%22&shmAll-net-ds=%223%22&shmAll-genCon=%22gene%22&deg-rok.dis.nor=%22CNF-TMM%22&shmAll-val.lgd.text=10&shmAll-color=%22yellow%2Corange%2Cred%22&upl-target=null&shmAll-net-dpwNetType=0&shmAll-t=2&shmAll-lgd.incld=%22Yes%22&shmAll-disDrop=0&sear-ids.but=0&deg-ssg.fdr=0.05&shmAll-net-mhmNav=%22mhmPlot%22&shmAll-scale.shm=1.4&dat-dtSel_cell_clicked=%7B%7D&deg-deg-ids.but=0&shmAll-net-dpwNetAdj=0&dat-scale=%22Row%22&shmAll-ggly.but=0&shmAll-vdo.itvl=1&sidebarItemExpanded=null&dat-CV2=10000&shmAll-val.lgd.feat=%22No%22&upl-svgInpath2=null&shmAll-scaleDrop=0&shmAll-net-mhm.but=0&shmAll-net-dpwModSize=0&shmAll-net-cpt.nw=0&shmAll-scroDrop=0&shmAll-net-color.net=%22yellow%2Corange%2Cred%22&dat-log=%22No%22&deg-ssg.fc=1&dat-dtSel_state=null&shmAll-interNav=%22interPlot%22&shmAll-title.size=12&deg-datDEG-CV2=10000&shmAll-net-thr=%22p%22&shmAll-net-col.but.net=0&shmAll-togDrop_state=false&shmAll-lgd.size=0.5&dat-dtSel_search=%22%22&shmAll-scale.ly=1&shmAll-val.lgd.row=1&shmAll-vdo.val.lgd=%22No%22&shmAll-net-cor.abs=%22No%22&upl-cusHelp=0&shmAll-vdo.dim=%22640x480%22&upl-fileIn=%22shoot_Mustroph%22&shmAll-net-mhm.v=0.2&shmAll-val.lgd.key=0.03&shmAll-vdo.key.row=2&shmAll-lgd.row=2&shmAll-vdo.key.size=0.04&deg-datDEG-fil.but=0&shmAll-dropdown=0&dat-fil.but=0&sear-ids.in=%22NDHE%20NADH%20dehydrogenase%20ND4L%22&shmAll-cs.v=%22Selected%20rows%22&shmAll-net-mat.scale=%22Row%22&shmAll-colDrop=0&shmAll-ext=%22NA%22&shmAll-fs=0&shmAll-vdo.but=0&shmAll-line.size=0.1&shmAll-net-netNav=%22netPlot%22&shmAll-net-min.size=15&shmAll-net-gen.sel=%22NDHE%22&upl-geneInpath=null&shmAll-line.color=%22grey70%22&deg-sam.con=%22feature%22&dat-A=0&dat-dtSel_columns_selected=null&upl-tar=null&dat-P=0&shmAll-titDrop=0&deg-ssg.meth=%5B%22edgeR%22%2C%22limma%22%5D&shm.sup=%22shmPanelAll%22&upl-met=null&shmAll-tis=null&shmAll-net-net.type=%22signed%22&shmAll-net-adj.in=%220.906%22&shmAll-vdo.label=%22No%22&deg-deg-ids.in=null&deg-datDEG-A=0&shmAll-lgd.key.size=0.04&shmAll-lgd.lab.size=2.5&dat-dtSel_cells_selected=%5B%5D&shmAll-togSld=0.7&shmAll-vdo.res=400&shmAll-net-max.edg=50&shmAll-col.but=0&shmAll-lgd.label=%22No%22&shmAll-val.lgd=0&dat-dtSel_rows_selected=null&sidebarCollapsed=true&shmAll-vdo.lab.size=2&shmAll-scrollH=450&upl-config=null&shmAll-relaSize=null&shmAll-res=300&right.bar=true'

root.arab.url <- '?_inputs_&deg-datDEG-P=0&deg-edg.lim.nor=%22CNF-TMM%22&dat-CV1=-10000&shmAll-net-dpwColNet=0&shmAll-shmMhNet=%22shm1%22&shmAll-togDrop=1&upl-svgInpath1=null&shmAll-col.n=2&deg-ssg.update=0&deg-datDEG-CV1=-10000&shmAll-net-dpwNetTar=0&shmAll-match=0&upl-dimName=%22None%22&shmAll-net-measure=%22correlation%22&shmAll-net-ds=%223%22&shmAll-genCon=%22gene%22&deg-rok.dis.nor=%22CNF-TMM%22&shmAll-val.lgd.text=10&shmAll-color=%22yellow%2Corange%2Cred%22&upl-target=null&shmAll-net-dpwNetType=0&shmAll-t=2&shmAll-lgd.incld=%22Yes%22&shmAll-disDrop=0&sear-ids.but=0&deg-ssg.fdr=0.05&shmAll-net-mhmNav=%22mhmPlot%22&shmAll-scale.shm=1.4&dat-dtSel_cell_clicked=%7B%7D&deg-deg-ids.but=0&shmAll-net-dpwNetAdj=0&dat-scale=%22Row%22&shmAll-ggly.but=0&shmAll-vdo.itvl=1&sidebarItemExpanded=null&dat-CV2=10000&shmAll-val.lgd.feat=%22No%22&upl-svgInpath2=null&shmAll-scaleDrop=0&shmAll-net-mhm.but=0&shmAll-net-dpwModSize=0&shmAll-net-cpt.nw=0&shmAll-scroDrop=0&shmAll-net-color.net=%22yellow%2Corange%2Cred%22&dat-log=%22No%22&deg-ssg.fc=1&dat-dtSel_state=null&shmAll-interNav=%22interPlot%22&shmAll-title.size=12&deg-datDEG-CV2=10000&shmAll-net-thr=%22p%22&shmAll-net-col.but.net=0&shmAll-togDrop_state=false&shmAll-lgd.size=0.5&dat-dtSel_search=%22%22&shmAll-scale.ly=1&shmAll-val.lgd.row=1&shmAll-vdo.val.lgd=%22No%22&shmAll-net-cor.abs=%22No%22&upl-cusHelp=0&shmAll-vdo.dim=%22640x480%22&upl-fileIn=%22root_Mustroph%22&shmAll-net-mhm.v=0.2&shmAll-val.lgd.key=0.03&shmAll-vdo.key.row=2&shmAll-lgd.row=2&shmAll-vdo.key.size=0.04&deg-datDEG-fil.but=0&shmAll-dropdown=0&dat-fil.but=0&sear-ids.in=%22NDHE%20NADH%20dehydrogenase%20ND4L%22&shmAll-cs.v=%22Selected%20rows%22&shmAll-net-mat.scale=%22Row%22&shmAll-colDrop=0&shmAll-ext=%22NA%22&shmAll-fs=0&shmAll-vdo.but=0&shmAll-line.size=0.1&shmAll-net-netNav=%22netPlot%22&shmAll-net-min.size=15&shmAll-net-gen.sel=%22NDHE%22&upl-geneInpath=null&shmAll-line.color=%22grey70%22&deg-sam.con=%22feature%22&dat-A=0&dat-dtSel_columns_selected=null&upl-tar=null&dat-P=0&shmAll-titDrop=0&deg-ssg.meth=%5B%22edgeR%22%2C%22limma%22%5D&shm.sup=%22shmPanelAll%22&upl-met=null&shmAll-tis=null&shmAll-net-net.type=%22signed%22&shmAll-net-adj.in=%220.906%22&shmAll-vdo.label=%22No%22&deg-deg-ids.in=null&deg-datDEG-A=0&shmAll-lgd.key.size=0.04&shmAll-lgd.lab.size=2.5&dat-dtSel_cells_selected=%5B%5D&shmAll-togSld=0.7&shmAll-vdo.res=400&shmAll-net-max.edg=50&shmAll-col.but=0&shmAll-lgd.label=%22No%22&shmAll-val.lgd=0&dat-dtSel_rows_selected=null&sidebarCollapsed=true&shmAll-vdo.lab.size=2&shmAll-scrollH=450&upl-config=null&shmAll-relaSize=null&shmAll-res=300&right.bar=true'

stage.arab.url <- '?_inputs_&shmAll-vdoNav=%22video%22&deg-datDEG-P=0&deg-edg.lim.nor=%22CNF-TMM%22&dat-CV1=-10000&shmAll-net-dpwColNet=0&shmAll-shmMhNet=%22shm1%22&shmAll-togDrop=0&upl-svgInpath1=null&shmAll-col.n=2&deg-ssg.update=0&deg-datDEG-CV1=-10000&shmAll-net-dpwNetTar=0&shmAll-match=0&upl-dimName=%22None%22&shmAll-net-measure=%22Correlation%22&shmAll-net-ds=%223%22&shmAll-genCon=%22gene%22&deg-rok.dis.nor=%22CNF-TMM%22&shmAll-val.lgd.text=10&shmAll-color=%22yellow%2Corange%2Cred%22&upl-target=null&shmAll-net-dpwNetType=0&shmAll-t=2&shmAll-lgd.incld=%22Yes%22&shmAll-disDrop=0&deg-ssg.fdr=0.05&shmAll-scale.shm=0.9&dat-dtSel_cell_clicked=%7B%7D&shmAll-net-dpwNetAdj=0&shmAll-ggly.but=0&shmAll-vdo.itvl=1&sidebarItemExpanded=null&dat-CV2=10000&shmAll-val.lgd.feat=%22No%22&upl-svgInpath2=null&shmAll-scaleDrop=0&shmAll-net-mhm.but=0&shmAll-net-dpwModSize=0&shmAll-net-cpt.nw=0&shmAll-scroDrop=0&shmAll-net-color.net=%22yellow%2Corange%2Cred%22&dat-log=%22No%22&deg-ssg.fc=1&dat-dtSel_state=null&shmAll-interNav=%22interPlot%22&shmAll-title.size=12&deg-datDEG-CV2=10000&shmAll-net-thr=%22p%22&shmAll-net-col.but.net=0&shmAll-lgd.size=0.5&dat-dtSel_search=%22%22&shmAll-scale.ly=1&shmAll-val.lgd.row=1&shmAll-vdo.val.lgd=%22No%22&shmAll-net-cor.abs=%22No%22&upl-cusHelp=0&shmAll-vdo.dim=%22640x480%22&upl-fileIn=%22growthStage_Mustroph%22&shmAll-net-mhm.v=0.2&shmAll-val.lgd.key=0.03&shmAll-vdo.key.row=2&shmAll-lgd.row=2&shmAll-vdo.key.size=0.04&deg-datDEG-fil.but=0&shmAll-dropdown=0&dat-fil.but=0&shmAll-cs.v=%22Selected%20rows%22&shmAll-net-mat.scale=%22Row%22&shmAll-colDrop=0&shmAll-ext=%22jpg%22&shmAll-fs=0&shmAll-vdo.but=0&shmAll-line.size=0.1&shmAll-net-min.size=15&shmAll-dld.but=0&shmAll-net-gen.sel=%22gene1%22&upl-geneInpath=null&shmAll-line.color=%22grey70%22&deg-sam.con=%22feature%22&dat-A=0&dat-dtSel_columns_selected=null&upl-tar=null&shmAll-net-clusNav=%22mhmPlot%22&dat-P=0&shmAll-titDrop=0&deg-deg-sch.mode=%22Multiple%22&deg-ssg.meth=%5B%22edgeR%22%2C%22limma%22%5D&shm.sup=%22shmPanelAll%22&upl-met=null&shmAll-tis=null&shmAll-net-net.type=%22signed%22&shmAll-net-adj.in=%221%22&shmAll-vdo.label=%22No%22&deg-datDEG-A=0&sear-sch.mul=%22gene1%22&shmAll-lgd.key.size=0.04&shmAll-lgd.lab.size=2.5&sear-sch.mode=%22Multiple%22&dat-dtSel_cells_selected=%5B%5D&shmAll-togSld=0.67&shmAll-vdo.res=400&shmAll-net-max.edg=50&shmAll-col.but=0&shmAll-shms.in=%22arabidopsis.thaliana_organ_shm1.svg%22&dat-dat.all.but=0&shmAll-lgd.label=%22No%22&shmAll-net-dpbThr=0&shmAll-val.lgd=0&sear-sch.mul.but=0&dat-dtSel_rows_selected=null&sidebarCollapsed=true&shmAll-vdo.lab.size=2&shmAll-scrollH=550&upl-config=null&shmAll-net-dpbMea=0&dat-scaleDat=%22Row%22&shmAll-res=300&shmAll-relaSize=1&right.bar=true'

mus.multi.dim.url <- '?_inputs_&shmAll-vdoNav=%22video%22&deg-datDEG-P=0&deg-edg.lim.nor=%22CNF-TMM%22&scell-covisMan-ncomT=2&dat-CV1=-10000&shmAll-net-dpwColNet=0&shmAll-shmMhNet=%22shm1%22&shmAll-togDrop=0&shmAll-col.n=2&dat-dtAll_cells_selected=%5B%5D&scell-covisAuto-normCoclus=%22fct%22&scell-covisMan-minSize=100&deg-ssg.update=0&shmAll-shmPar=%22basic%22&deg-datDEG-CV1=-10000&shmAll-net-dpwNetTar=0&shmAll-raster=%22Yes%22&shmAll-net-ds=%223%22&deg-rok.dis.nor=%22CNF-TMM%22&shmAll-net-measure=%22Correlation%22&scell-covisAuto-tailor-dimCell=%22UMAP%22&shmAll-genCon=%22gene%22&shmAll-val.lgd.text=10&shmAll-color=%22yellow%2Corange%2Cred%22&shmAll-net-dpwNetType=0&shmAll-t=2&scell-covisAuto-filBlkCV1=0.1&shmAll-lgd.incld=%22Yes%22&scell-covisAuto-tailor-tailorHelp=0&shmAll-disDrop=0&dat-normDat=%22none%22&deg-ssg.fdr=0.05&shmAll-scale.shm=1.1&shmAll-net-dpwNetAdj=0&shmAll-ggly.but=0&scell-covisMan-nn.graph=%22buildSNNGraph%22&shmAll-dims=%22UMAP%22&scell-covisMan-ncomU=2&shmAll-vdo.itvl=1&sidebarItemExpanded=null&dat-CV2=10000&scell-covisAuto-clusMeth=%22wt%22&shmAll-val.lgd.feat=%22No%22&shmAll-scaleDrop=0&shmAll-net-mhm.but=0&shmAll-net-dpwModSize=0&shmAll-dimLgdRows=2&dat-dtab.shm=%22dTabAll%22&shmAll-scroDrop=0&shmAll-net-cpt.nw=0&shmAll-profile=%22Yes%22&dat-log=%22No%22&dat-page=300&deg-ssg.fc=1&dat-tran.scale.but.prof=0&dat-sig.but=0&scell-covisMan-maxSize=3000&scell-covisAuto-maxRank=50&shmAll-net-color.net=%22yellow%2Corange%2Cred%22&shmAll-interNav=%22interPlot%22&shmAll-title.size=12&deg-datDEG-CV2=10000&shmAll-net-thr=%22p%22&shmAll-net-col.but.net=0&shmAll-lgd.size=0.8&bulk=%7B%22collapsible%22%3Atrue%2C%22collapsed%22%3Afalse%2C%22closable%22%3Afalse%2C%22visible%22%3Atrue%2C%22status%22%3A%22primary%22%2C%22solidHeader%22%3Atrue%2C%22width%22%3A12%7D&shmAll-val.lgd.row=1&shmAll-scale.ly=1&shmAll-vdo.val.lgd=%22No%22&scell-covisAuto-filBlkCV2=200&shmAll-net-cor.abs=%22No%22&scell-covisMan-maxRank=50&scell-covisMan-hvgN=3000&shmAll-vdo.dim=%22640x480%22&scell-covisAuto-filPGen=0.01&upl-fileIn=%22multiDimensions_Attilio%22&shmAll-net-mhm.v=0.2&shmAll-val.lgd.key=0.03&scell-methCovis=%22auto%22&shmAll-vdo.key.row=2&shmAll-lgd.row=3&shmAll-vdo.key.size=0.04&shmAll-vdo.bar.width=0.1&deg-datDEG-fil.but=0&scell-covisMan-parManBut=0&shmAll-dropdown=0&dat-fil.but=0&dat-tran.scale.but.sel=0&shmAll-cs.v=%22Selected%20rows%22&shmAll-colDrop=0&shmAll-transBut=0&shmAll-rematch-matHelp=0&shmAll-alpOverBut=0&shmAll-fs=0&scell-covisAuto-tailor-selBlkBut=0&shmAll-vdo.but=0&shmAll-ext=%22jpg%22&shmAll-coal=%22No%22&shmAll-net-mat.scale=%22Row%22&shmAll-line.size=0.1&shmAll-net-min.size=15&scell-covisHelp=0&scell-direc=%22toBulk%22&scell-covisAuto-parAutoBut=0&shmAll-dld.but=0&shmAll-net-gen.sel=%22X1190002F15Rik%22&dat-sig.max=%22%22&shmAll-line.color=%22grey70%22&deg-sam.con=%22feature%22&shmAll-alpOver=1&scell-covisAuto-dimSel=%22PCA%22&upl-tar=null&scell-covisAuto-tailor-coclusPlotBut=0&dat-A=0&scell-covisAuto-filBlkP=0.1&dat-dtAll_rows_selected=null&scell-covisMan-tabSetCell=%22datCell%22&scell-covisAuto-graphMeth=%22knn%22&scell-covisMan-cntThr=0&scell-covisMan-scell.cluster=%22cluster_walktrap%22&scell-covisMan-hvgP=0.1&dat-dtAll_columns_selected=null&scell-covisMan-ntopT=500&shmAll-net-clusNav=%22mhmPlot%22&dat-P=0&shmAll-titDrop=0&deg-ssg.meth=%5B%22edgeR%22%2C%22limma%22%5D&deg-deg-sch.mode=%22Multiple%22&shm.sup=%22shmPanelAll%22&shmAll-net-adj.in=%221%22&shmAll-tis=null&shmAll-net-net.type=%22signed%22&deg-datDEG-A=0&shmAll-vdo.label=%22No%22&scell-covisMan-ntopU=500&scell-covisAuto-filBlkA=1&scell-covisAuto-minRank=5&dat-sig.min=%22%22&sear-sch.mode=%22Multiple%22&sear-sch.mul=%22X1190002F15Rik%22&shmAll-togSld=0.67&shmAll-vdo.res=400&shmAll-lgd.key.size=0.04&shmAll-lgd.lab.size=2.5&shmAll-net-max.edg=50&shmAll-col.but=0&scell-covisMan-nmads=3&scell-covisMan-dimredNav=%22Plot%22&dat-dat.all.but=0&shmAll-lgd.label=%22No%22&shmAll-net-dpbThr=0&cell=%7B%22collapsible%22%3Atrue%2C%22collapsed%22%3Afalse%2C%22closable%22%3Afalse%2C%22visible%22%3Atrue%2C%22status%22%3A%22primary%22%2C%22solidHeader%22%3Atrue%2C%22width%22%3A12%7D&scell-covisAuto-tabSetCellAuto=%22datCell%22&shmAll-val.lgd=0&scell-covisMan-normBlk=%22VST%22&scell-covisMan-minRank=5&sear-sch.mul.but=0&dat-dtAll_search=%22%22&dat-dtAll_cell_clicked=%7B%7D&sidebarCollapsed=true&scell-covisAuto-tailor-selBlkCancel=0&dat-dtAll_state=%7B%22time%22%3A1665217894198%2C%22start%22%3A0%2C%22length%22%3A81%2C%22order%22%3A%5B%5D%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Afalse%2C%22regex%22%3Atrue%2C%22caseInsensitive%22%3Atrue%7D%2C%22columns%22%3A%5B%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%2C%7B%22visible%22%3Atrue%2C%22search%22%3A%7B%22search%22%3A%22%22%2C%22smart%22%3Atrue%2C%22regex%22%3Afalse%2C%22caseInsensitive%22%3Atrue%7D%7D%5D%2C%22scroller%22%3A%7B%22topRow%22%3A0%2C%22baseScrollTop%22%3A0%2C%22baseRowTop%22%3A0%2C%22scrollTop%22%3A0%7D%7D&shmAll-dpwAlpOver=0&shmAll-scrollH=450&upl-config=null&shmAll-net-dpbMea=0&dat-scaleDat=%22Row%22&shmAll-res=300&shmAll-relaSize=1&shmAll-vdo.lab.size=2&right.bar=true&scell-covisMan-rematchCell-matHelp=0&scell-covisMan-pcs=50&scell-covisAuto-filPCell=0.1&scell-covisAuto-asgThr=0'

# Extract parameter values from url.
url_val <- function(na, lis.url) {
  # if (!exists('lis.url')) return('null')
  if (!na %in% names(lis.url$par)) return('null')
  if (length(lis.url$par)==0) val <- 'null' else val <- lis.url$par[[na]]
  # In "ifelse", the length of returned value is same with the first argument.
  # val <- ifelse(length(lis.url$par)==0, 'null', lis.url$par[[na]])
  gsub('\\"', '', val)
}

# Import internal functions.

com_roc <- get('com_roc', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
qc_cell <- get('qc_cell', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
check_obj <- get('check_obj', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
img_pa_na <- get('img_pa_na', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
covis_trans <- get('covis_trans', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
svg_separ <- get('svg_separ', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
sc_qc_plot <- get('sc_qc_plot', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
detect_cluster <- get('detect_cluster', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
dim_color_coclus <- get('dim_color_coclus', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
dim_color2cell <- get('dim_color2cell', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
dim_color <- get('dim_color', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
nn_graph <- get('nn_graph', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
svg_raster <- get('svg_raster', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

raster_path <- get('raster_path', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# dat_fun <- get('dat_fun', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

df_is_as <- get('df_is_as', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

scale_all <- get('scale_all', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

venn_inter <- get('venn_inter', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

deg_lis <- get('deg_lis', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

edgeR <- get('edgeR', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

limma <- get('limma', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

deseq2 <- get('deseq2', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

distt <- get('distt', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

up_dn <- get('up_dn', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

deg_ovl_mat  <- get('deg_ovl_mat', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

deter_core <- get('deter_core', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

rela_size <- get('rela_size', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

norm_data <- get('norm_data', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

cord_parent <- get('cord_parent', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

use <- get('use', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

cord <- get('cord', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

xy0 <- get('xy0', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

xy <- get('xy', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

tit_id <- get('tit_id', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

out_ply <- get('out_ply', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

sort_gen_con <- get('sort_gen_con', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

test_ffm <- get('test_ffm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

read_hdf5 <- get('read_hdf5', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

matrix_hm <- get('matrix_hm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Function to extract nearest genes.
sub_na <- get('sub_na', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

adj_mod <- get('adj_mod', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

filter_data <- get('filter_data', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

col_bar <- get('col_bar', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

lay_shm <- get('lay_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

nod_lin <- get('nod_lin', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Break combined path to a group (g=TRUE) or siblings (g=FALSE).
path_br <- get('path_br', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# The outline or tissue nodes are checked for combines paths. If combined paths are detected, those outside a group are broken to a group while those inside a group are broken as siblings.  
path_br_all <- get('path_br_all', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# 'a' nodes are not removed.
svg_attr <- get('svg_attr', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

svg_df <- get('svg_df', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Separate SHMs of grob and ggplot. Different SHMs (under different SVGs) of same 'gene_condition' are indexed with suffixed of '_1', '_2', ...
gg_shm <- get('gg_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
grob_shm <- get('grob_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Subset data matrix by correlation or distance measure.
submatrix <- get('submatrix', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Adjust legend key size and rows in ggplot.
gg_lgd <- get('gg_lgd', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Add value keys SHMs.
gg_2lgd <- get('gg_2lgd', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Prepare interactive SHMs in html.
html_ly <- get('html_ly', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Make videos.
video <- get('video', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Shown popup window. 
modal <- get('modal', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Data table in manual covis.
dat_covis_man <- function(sce, nr=1000, nc=100) {
  dat <- round(assay(sce), 3)
  if (nrow(dat) >= nr) r.idx <- nr else r.idx <- nrow(dat)
  if (ncol(dat) >= nc) c.idx <- nc else c.idx <- ncol(dat)
  datatable(as.matrix(dat[seq_len(r.idx), seq_len(c.idx)]), selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
  options=list(pageLength=20, lengthMenu=c(10, 20, 50, 100), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE,  scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, fixedColumns=list(leftColumns=1)) 
  ) 
}


# Extract a 1-column data frame of URLs. If no column of URL is present, the default google-search URLs are composed.
link_dat <- function(df.met) {
  cna <- colnames(df.met); rna <- rownames(df.met)
  link.idx <- grep('link|links', cna, ignore.case=TRUE)[1]
  if (is.na(link.idx)) {
    # Iterative operation on data frame: vectorization is faster than for/lapply loop.
    link <- paste0('<a href=\"https://www.google.com/search?q=', rna, '" target="_blank">link</a>')
    # link <- lapply(rownames(df.met), function(x) a("link", href=paste0('https://www.google.com/search?q=', x), target="_blank"))
    # link <- unlist(lapply(link, as.character))
  } else { link <- df.met[, link.idx]; link <- gsub('\\\\"', '\'', link) }
  return(data.frame(link=link, row.names=rownames(df.met)))
}


# Import input matrix, able to deal with/separate numeric matrix, character matrix, and mixture of both.
fread_df <- function(input, isRowGene=TRUE, header=TRUE, sep='auto', fill=TRUE, rep.aggr='mean', check.names=FALSE) {
  
  if (is(input, 'dgCMatrix')|is(input, 'matrix')) input <- as.data.frame(as.matrix(input))
  if (!is(input, 'data.frame')) {
  df0 <- tryCatch({
    fread(input=input, header=header, sep=sep, fill=fill, check.names=check.names)
  }, error = function(error_condition) {
    # Deals with only one column with row names.
    fread(input = input, header = FALSE, sep = sep, fill = FALSE, check.names = check.names)
  }) 
    cna <- make.names(colnames(df0))
    if (cna[1]=='V1') cna <- cna[-1] else cna <- cna[-ncol(df0)] 
    df1 <- as.data.frame(df0); rownames(df1) <- make.names(df1[, 1])
    df1 <- df1[, -1, drop = FALSE]; colnames(df1) <- cna
    if(isRowGene==FALSE) df1 <- t(df1)
    cna <- colnames(df1); rna <- rownames(df1) 
  } else { df1 <- input; rna <- rownames(df1); cna <- colnames(df1) }
  # Covert factors to character. Only data.frame works not matrix.
  fct.idx <- vapply(df1, is.factor, logical(1))
  df1[fct.idx] <- lapply(df1[fct.idx], as.character) 
  # Subsetting identical column names in a matrix will not trigger appending numbers.
  df1 <- as.matrix(df1)
  # Isolate data and row metadata.
  na <- vapply(seq_len(ncol(df1)), function(i) { tryCatch({ as.numeric(df1[, i]) }, warning=function(w) { return(rep(NA, nrow(df1))) }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(df1)) )
  if (nrow(df1)==1) na <- matrix(na, byrow=TRUE, ncol=ncol(df1))
  na <- as.data.frame(na); rownames(na) <- rna; colnames(na) <- cna
  vap <- df_is_as(na, is.na); idx <- colSums(vap)!=0
  df.num <- na[!idx]; colnames(df.num) <- cna <- cna[!idx]
  df.met.all <- as.data.frame(df1)[idx]
  cat('Preparing URLs .. \n')
  df.link <- link_dat(df.met.all) # Works if ncol(df.met.all) is 0.
  if (ncol(df.met.all) > 0) {
  cat('Preparing metadata .. \n')
    met.idx <- grep('^metadata$', colnames(df.met.all), ignore.case = TRUE)[1]
    if (!is.na(met.idx)) { 
      df.met <- df.met.all[, met.idx, drop = FALSE] 
      colnames(df.met) <- 'metadata'
      df.met <- cbind(df.met, df.link)
    } else df.met <- df.link
  } else df.met <- df.link

  # Only row metadata.
  if (ncol(df.num) == 0) {
    return(list(df.aggr = NULL, df.met=as.data.frame(df.met), df.rep = NULL, con.na = FALSE))
  }
  form <- grepl("__", cna); if (sum(form)==0) { cna <- colnames(df.num) <- paste0(cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE
  if(sum(is.na(as.numeric(as.matrix(df.num))))>=1) return('Make sure all values in data matrix are numeric.')
  
  df.rep <- df.num; rna <- rownames(df.rep)
  df.rep <- df_is_as(df.rep, as.numeric)
  # Aggregate replicates.
  if (any(duplicated(cna)) & !is.null(rep.aggr)) {

    # To keep colnames, "X" should be a character, not a factor.
    if (rep.aggr=='mean') df.num <- sapply(X=unique(cna), function(x) rowMeans(df.num[, cna==x, drop=FALSE]))
    if (rep.aggr=='median') {
      df.num <- sapply(X=unique(cna), function(x) Biobase::rowMedians(df.num[, cna==x, drop=FALSE]))
      rownames(df.num) <- rna
    }

  }; df.aggr <- df_is_as(df.num, as.numeric)
  return(list(df.aggr=as.data.frame(df.aggr), df.met=as.data.frame(df.met), df.rep=as.data.frame(df.rep), con.na=con.na))

}

# Separate colour ingredients.
col_sep <- function(color) {

  color <- gsub(' |\\.|-|;|,|/', '_', color)
  color <- strsplit(color, '_')[[1]]
  color <- color[color!='']; return(color)

}

# Check/process "sample__condition" in the se extracted from tar.
se_from_db <- function(se) {
  dat <- assay(se); cold <- colData(se)
  form <- grepl("__", colnames(dat))
  if (sum(form)==0) {
    if (all(c('sample', 'condition') %in% colnames(cold))) { 
      cna <- colnames(dat) <- paste0(cold$sample, '__', cold$condition)
      if (any(duplicated(cna))) return('Duplicated "sample__condition" replicates are detected in the selected dataset!')
    } else if ('sample' %in% colnames(cold)) {
      if (any(duplicated(cold$sample))) return('The "sample" should not be duplicated in the absence of "condition"!')
      colnames(dat) <- cold$sample
    }
  }; return(dat)
}


# Extract target svgs in tar into tmp folder, and return the paths. 
extr_svg <- function(file, name) {
  dir <- paste0(tempdir(check=TRUE), '/svg_shm')
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  untar(file, exdir=dir, tar='tar')
  pa <- paste0(dir, '/', name) 
  if (file.exists(pa)) return(pa) else return()
}

# Extract svg path/na from uploaded or internal tar files if not found in 'example' folder.
svg_pa_na <- function(svg.path, pa.svg.upl, raster.ext) {
  svg.na <- NULL; for (i in seq_along(svg.path)) {
    # Extract svg names. 
    str <- strsplit(svg.path[[i]], '/')[[1]]
    na0 <- str[length(str)]
    if (!grepl(paste0('\\', c('.svg', raster.ext), '$', collapse='|'), na0)) return('No aSVG/template file is detected! Solution: 1) select another aSVG and rematch it to data; 2) add an aSVG/template file for the selected data in the backend aSVG tar file or uploaded aSVG tar file.')
    svg.na <- c(svg.na, na0)
    # Complete uploaded svg paths.
    if (!grepl('example/', svg.path[[i]])) {
      # The data/svg precedence: uploaded tar > internal tar > default examples. The duplicated data/svgs are removed according to this precedence when processing data/svg upstream.
      pa0 <- NULL; if (!is.null(pa.svg.upl)) pa0 <- extr_svg(file=pa.svg.upl, name=na0)
      if (is.null(pa.svg.upl)|is.null(pa0)) {
        tar.all <- list.files('example', pattern='\\.tar$', full.names=TRUE)
        tar.svg <- tar.all[!grepl('data_shm.tar$', tar.all)][1]
        pa0 <- extr_svg(file=tar.svg, name=na0)
      }; if (is.null(pa0)) return(paste0("This aSVG/template file is not detected: ", na0, "!")) else svg.path[[i]] <- pa0
    }
  }; return(list(svg.path=svg.path, svg.na=svg.na))
}

# Check suffixes if multiple svgs.
svg_suffix <- function(svg.path, svg.na, raster.ext) {
  if (length(svg.na)>1) {
    ext <- paste0('_shm\\d+\\', c('.svg', raster.ext), '$', collapse='|')
    if (!all(grepl(ext, svg.na, perl=TRUE))) return("Suffixes of aSVGs and templates should be indexed as '_shm1.svg', '_shm1.png', '_shm2.svg', '_shm2.png', '_shm3.svg', '_shm3.png', ...")
    ord <- order(gsub('.*_(shm.*)$', '\\1', svg.na))
    svg.path <- svg.path[ord]; svg.na <- svg.na[ord]  
  }; return(list(svg.path=svg.path, svg.na=svg.na))
}

# Convert aSVG features to draggable items.
# ns() is the namespace in shiny modules.
ft2tag <- function(ft){
  lapply(ft, function(i) { tag("span", list(class = class(i), tags$span(class = "glyphicon glyphicon-move"), i)) }
  )
}


## Rematch features.
# Create a panel for each data feature, where aSVG features can be dropped.
# ns() is the namespace in shiny modules.
ft_dat <- function(x, ns) {
 span(class = "panel panel-default",
   div(class = "panel-heading", x), 
   div(class = "panel-body", id = ns(x))
  )
}

ft_lis_dat <- function(x, ns) {
 span(class = "panel panel-default",
   div(class = "panel-heading", names(x)), 
   div(class = "panel-body", id = ns(names(x)), ft2tag(x[[1]]))
  )
}

# Allow features are draggable across panels.
# ns() is the namespace in shiny modules.
ft_js <- function(x, ns) {
  sortable_js(css_id = ns(x),
    options = sortable_options(
      multiDrag = NULL, sort = FALSE, animation = 1000, direction = NULL, 
      group = list(name = "sortGroup1", put = TRUE),
      onSort = sortable_js_capture_input(ns(x))
    )
  )
}

# Complete matching interface.
match_interface <- function(to.ft, to.div.id='ftSVG', to.div.tit='Features in aSVG', from.ft, from.div.tit='Features in data',   ns) {
  to.ft <- sort(to.ft); from.ft <- sort(from.ft) 
  frow <- fluidRow(
    span(class = "panel panel-default", style = 'margin-left:0px',
      div(class = "panel-heading", strong(to.div.tit)), 
      div(class = "panel-body", id = ns(to.div.id), ft2tag(to.ft)) 
      ),
    div(class = "panel panel-default", 
      div(class = "panel-heading", strong(from.div.tit)),  
      lapply(from.ft, ft_dat, ns = ns)  
    ), lapply(c(to.div.id, from.ft), ft_js, ns = ns) # Items are interchangeable across ftSVG and sam.all.
  ); return(frow) 
} 


# Clean trash files in animation and video.
ggly_rm <- function() {
  if (dir.exists('www/ggly/')) {
    cat("Removing animation files in 'www/ggly/' ... \n")
    unlink('www/ggly/lib', recursive=TRUE)
    file.remove(list.files('www/ggly/', '*.html$', full.names=TRUE))
  } else dir.create('www/ggly', recursive=TRUE)
  if (dir.exists('R/www')) {
    cat("Removing animation files in 'R/www/ggly/' ... \n")
    unlink('R/www', recursive=TRUE)
  }
}
vdo_rm <- function() {
  if (dir.exists('www/video/')) {
    cat("Removing video file in 'www/video/' ... \n")
    file.remove(list.files('www/video/', '*.mp4$', full.names=TRUE))
  } else dir.create('www/video/', recursive=TRUE)
 if (dir.exists('R/www')) {
    cat("Removing video file in 'R/www/video' ... \n")
    unlink('R/www', recursive=TRUE)
 }
}






