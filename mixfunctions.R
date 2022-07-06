##########
## This is a "model based" approach in which we build a graph of the history of populations as splits, mixtures, and drift.
## This code defines a graph, a simulation procedure, and an inference procedure based on "optim".

defaultcontrol=list(maxit=1000,trace=1,factr=1e10,pgtol=0)

cnode<-function(id,pl=NA,pr=NA,cl=NA,cr=NA,d=NA,t=NA,p=NA,w=NA,type=NA){
    ## This is the generic function to create or update a node
    ## Pass an id to create a new node
    ## Pass a cnode to update it
    ## Returns a cnode
    if(is(id,"cnode")){
        if(!is.na(pl)) id$pl=pl
        if(!is.na(pr)) id$pr=pr
        if(!is.na(cl)) id$cl=cl
        if(!is.na(cr)) id$cr=cr
        if(!is.na(d)) id$d=d
        if(!is.na(t)) id$t=t
        if(!is.na(p)) id$p=p
        if(!is.na(w)) id$w=w
        if(!is.na(type)) id$type=type
        return(id)
    }else{
        if(is.na(type)) type="split"
        r=list(id=id,
               pl=pl,
               pr=pr,
               cl=cl,
               cr=cr,
               d=d,
               t=t,
               p=p,
               w=w,
               type=type)
        class(r)="cnode"
        return(r)
    }
}

simCoal<-function(n,times="coal",labels=paste0("t",1:n),outgroup=numeric()){
    ## Simulate a coalescent tree in the "clarity graph" framework
    ## Returns a "clarity graph" object
    nodes=list()
    for(i in 1:n) nodes[[i]]=cnode(i,t=0,d=0,w=1)
    alive=1:n
    outgroupnums=alive[labels%in%outgroup]
    alive=alive[!labels%in%outgroup]
    time=0
    for(i in n:2){
        ii=length(nodes)+1
        if(times=="coal"){
            t=rexp(1,i)
        }else if(times=="test"){
            t=(n-i+1) 
        }else{
            stop("Invalid times argument!")
        }
        time=time+t
        if((length(alive)<2)&&(length(outgroup)>0)){ ## Add in the outgroups last
            alive=c(alive,outgroupnums)
        }
        pair=sample(alive,2)
##        print(paste("i=",i,"pair = ",pair[1],pair[2]))
        nodes[[ii]]=cnode(ii,cl=pair[1],cr=pair[2],t=time,d=0)
        nodes[[pair[1]]]=cnode(nodes[[pair[1]]],pl=ii,d=time - nodes[[pair[1]]]$t,w=1)
        nodes[[pair[2]]]=cnode(nodes[[pair[2]]],pl=ii,d=time - nodes[[pair[2]]]$t,w=1)
        nodes[[pair[1]]]$pl=ii
        nodes[[pair[2]]]$pl=ii
        alive=c(alive[!alive%in%pair],ii)
    }
    g=list(nl=nodes,
           tips=1:n,
           root=length(nodes),
           tip.label=labels,
           internal=(n+1):length(nodes),
           mix=numeric(),
           spare=numeric(),
           mixparmap=numeric(),
           outgroup=outgroup,
           outmix=outgroup,
           n=n)
    class(g)="cg"
    g
}

mixedge<-function(g,source,target,alpha,w){
    ## Add a mixture edge to a clarity graph object

    if(length(g$spare)>0){
        ti=tail(g$spare,1)
        g$spare=g$spare[-(g$spare%in%ti)]
    }else{
        ti=length(g$nl)+1
    }
    ##    if(target%in%tipsunder(source)) stop("Invalid mixture edge?")
    parsource=g$nl[[source]]$pl
##    if(is.na(parsource)) stop("ERROR: tried to make the root a mixture!")
    ## Add the new node
    g$nl[[ti]]=cnode(ti,
                     pl=parsource,
                     cl=source,
                     cr=target,
                     d=g$nl[[source]]$d*(1-alpha),
                     t=(g$nl[[parsource]]$t - g$nl[[source]]$t)*alpha +  g$nl[[source]]$t,
                     w=1,
                     type="mixture")
################
## TODO: If the target already has a right parent, we need to create a new node
    if(!is.na(g$nl[[target]]$pr)){
        stop("Unimplemented exception: target would have three parents!")
    }
    ## Update the target
    g$nl[[target]]=cnode(g$nl[[target]],pr=ti,w=g$nl[[target]]$w * w)
    ## Update the source node (called source)
    g$nl[[source]]=cnode(g$nl[[source]],pl=ti,d=g$nl[[source]]$d*alpha)
    ## Update the original source's parent node (called parsource)
    if(!is.na(parsource)) { ## It was not the root.
        if(g$nl[[parsource]]$cl==source){ ## Check whether we rewire the left or right child
            g$nl[[parsource]]=cnode(g$nl[[parsource]],cl=ti)
        }else{
            g$nl[[parsource]]=cnode(g$nl[[parsource]],cr=ti)
        }
    }else{ ## It was not the root
        g$root=ti
        g$nl[[source]]=cnode(g$nl[[source]],d=0.5,w=1)
    }
    ## Update the list of nodes
    g$internal=c(g$internal,ti)
    g$mix=c(g$mix,ti)
    g
}

removemixedge<-function(g,i,careful=TRUE){
    ## Remove the node at index i
    ## If we are careful, we expect a proper graph
    ## Otherwise we accept missing right children for "dangling" mixture nodes with no right child
    if(g$nl[[i]]$type!="mixture") stop(paste("ERROR: Invalid request to remove node",i,"which is not a mixture edge"))
    cl = g$nl[[i]]$cl
    cr = g$nl[[i]]$cr
    pl = g$nl[[i]]$pl
    pr = g$nl[[i]]$pr
    if(is.na(cl)) stop(paste("ERROR: node",i,"has no left child!"))
    if(is.na(cr)&&careful) stop(paste("ERROR: node",i,"has no right child!"))
    if(is.na(pl) && (i != g$root)) stop(paste("ERROR: node",i,"has no left parent!"))
    ## Update left child
    g$nl[[cl]]$pl = g$nl[[i]]$pl
    g$nl[[cl]]$d  = g$nl[[i]]$d + g$nl[[cl]]$d
    ## Update right child
    if(!is.na(cr)) g$nl[[cr]]$pr= NA
    if(!is.na(cr)) g$nl[[cr]]$w = 1
    ## Update parent if present
    if(g$root!=i){
        isleftchild=(g$nl[[pl]]$cl==i)
        if(isleftchild) g$nl[[pl]]$cl = cl
        else g$nl[[pl]]$cr = cl
    }else{ # Its the root; have to make the left child the new root
        g$root=cl
    }
    ## Update right parent if present
    if((!is.na(pr))&&(!is.na(cr))){
        ## Must have been the right child
        ## Now point to the right child of the removed node
        ## Which is always safe
        g$nl[[pr]]$cr = cr
        g$nl[[cr]]$pr = pr
    }
    ## update graph information
    if(!is.na(cr)) g$nl[[i]] =cnode(id=cr,w=1,d=0,type="spare")
    g$mix=g$mix[-which(g$mix==i)]
    g$internal=g$internal[-which(g$internal==i)]
    g$spare=c(i,g$spare)
    g
}

myregraftmixedge<-function(g,i,source,target,alpha,w){
    if(is(g,"cglist")){
        ret=lapply(g,function(x){
            myregraftmixedge(x,i,source,target,alpha,w)
        })
        class(ret)="cglist"
        return(ret)
    }
    ## Does a complete removal and replacement of a mixture edge, with specified parameters
    if(!is.na(i)) g=removemixedge(g,i)
    g=mixedge(g,source,target,alpha,w)
    g
}


tipsunder<-function(g,i,w=FALSE){
    ## Return the set of tips underneath a specified node
    ## if w (weighted) return this as a matrix of the weight for tip t (rows) under node i
    ## if !w (unweighted) return this as a list of tip ids
    if(w){
        ## Weighted version
        r=rep(0,g$n)
        if(!is.na(g$nl[[i]]$cl)) {
            if(g$nl[[g$nl[[i]]$cl]]$pl==i){ # Is this the left parent?
                w = 1 - g$nl[[g$nl[[i]]$cl]]$w
            }else w = g$nl[[g$nl[[i]]$cl]]$w
            print(paste("node",i,"left child",g$nl[[i]]$cl,"weight",w))
            r = r + w * tipsunder(g,g$nl[[i]]$cl,w=w)
        }
        if(!is.na(g$nl[[i]]$cr)){
            if(g$nl[[g$nl[[i]]$cr]]$pl==i){ # Is this the left parent?
                w = 1 - g$nl[[g$nl[[i]]$cr]]$w
            }else w= g$nl[[g$nl[[i]]$cr]]$w
            print(paste("node",i,"right child",g$nl[[i]]$cr,"weight",w))
            r=r + w * tipsunder(g,g$nl[[i]]$cr,w=w)
        }
        if((is.na(g$nl[[i]]$cl))&(is.na(g$nl[[i]]$cr)) & (i%in%g$tips) ) {
            r[i]=r[i] + 1
        }
        r
    }else{
        ## Unweighted version
        r=c()
        if(!is.na(g$nl[[i]]$cl)) {
            r=c(r,tipsunder(g,g$nl[[i]]$cl,w=w))
        }
        if(!is.na(g$nl[[i]]$cr)){
            r=c(r,tipsunder(g,g$nl[[i]]$cr,w=w))
        }
        if((is.na(g$nl[[i]]$cl))&(is.na(g$nl[[i]]$cr)) & (i%in%g$tips)) r=c(r,i)
        r
    }
}

nodesunder<-function(g,i,visited=numeric()){
    r=c()
    if(i %in% visited) stop(paste("ERROR in nodesunder: node",i,"has already been visited!"))
    if(!is.na(g$nl[[i]]$cl)) {
        r=c(r,nodesunder(g,g$nl[[i]]$cl,visited=c(visited,i)))
    }
    if(!is.na(g$nl[[i]]$cr)){
        r=c(r,nodesunder(g,g$nl[[i]]$cr,visited=c(visited,i)))
    }
    r=c(r,i)
    unique(r)
}
c_get<-function(g,n,getleft,what="d"){
    ## Get edge properties from a parent node n to either left or right child
    if(getleft & is.na(n$cl)) return(NULL)
    if((!getleft) & is.na(n$cr)) return(NULL)

####### Special cases for mixture nodes
    if((!getleft) & (n$type=="mixture") & (what=="d")) return(0) # Special case: mixture edges have d=0 on the right
    ## Get the child node
    if(getleft) {tn=g$nl[[n$cl]]
    }else tn=g$nl[[n$cr]]    
    if((what=="w") & (!is.na(tn$pr))) {
        ## Special case: mixture edges have 1-w to the left parent and w to the right parent
        if(tn$pr==n$i) {return(tn[[what]])
        }else return(1-tn[[what]])
    }
    
    ## Otherwise normal
    if(getleft) return(g$nl[[n$cl]][[what]])
    return(g$nl[[n$cr]][[what]])
}

edges.cg<-function(g){
    ## Extract properties of all edges as a data frame
    edge.length=do.call("c",
                        lapply(rev(g$nl[g$internal]),function(n){
                            c(c_get(g,n,TRUE,"d"),
                              c_get(g,n,FALSE,"d"))
                        })
                        )
    edge.w=do.call("c",
                        lapply(rev(g$nl[g$internal]),function(n){
                            c(c_get(g,n,TRUE,"w"),
                              c_get(g,n,FALSE,"w"))
                        })
                        )
    edge=do.call("rbind",
                 lapply(rev(g$nl[g$internal]),function(n){
                     tr=c()
                     if(!is.na(n$cl)) tr= rbind(tr,c(n$id,n$cl))
                     if(!is.na(n$cr)) tr= rbind(tr,c(n$id,n$cr))
                     tr
                 })
                 )
    edge.ismix=do.call("c",
                        lapply(rev(g$nl[g$internal]),function(n){
                            pr=c(c_get(g,n,TRUE,"pr"),
                                 c_get(g,n,FALSE,"pr"))
                            pr[is.na(pr)]=FALSE
                            as.numeric(pr==n$id)
                        })
                        )
    colnames(edge)=c("parent","child")
    as.data.frame(cbind(edge,
                        weight=edge.w,
                        length=edge.length,
                        mix=edge.ismix))
}

locationasmatrix=function(g){
    ## Location of each node 
    tlayout=cbind(t=sapply(g$nl,function(x)x$t),
                  p=sapply(g$nl,function(x)x$p),
                  index=sapply(g$nl,function(x)x$id))
    rownames(tlayout)=tlayout[,"index"]
    tlayout
}

plot.cg=function(g,ref=g,arrows.length=0.1,edges=NULL,
                 arrows.col=c("grey","red"),
                 text.col=c("darkgrey","darkred"),
                 digits=1,mindepth=1e-2,maxdepth=Inf,rightalign=FALSE,
                 tips=NULL,cex.edge.text=0.5,
                 label.mixture=TRUE,
                 label.nonmixture=TRUE,
                 label.internal=TRUE,
                 labels.col="black",
                 showedges=NA,
                 keeplocation=FALSE,
                 lwd=1,textdelta=0,
                 format=c("triangular","rectangular"),rdelta=0.1,rdelta2=0.1,vadj.edge=0.2,
                 adj.node=0,adj.edge=0,show=TRUE,showaxis=TRUE,cex.labels=1,
                 ...){
    ## Plot without requiring igraph
    if(!keeplocation) g=assignlocation(g,mindepth,maxdepth,rightalign)
    if(all(is.null(edges))){
        edges=edges.cg(g)
        edges[,"weight"]=format(edges[,"weight"],digits=digits)
    }
    tlayout=locationasmatrix(g)
    labels=tlayout[,"index"]
    if(all(is.null(tips))) labels[1:length(g$tip.label)]=g$tip.label
    if(!all(is.null(tips))) labels[1:length(tips)]=tips
    if(!label.internal) labels[(length(g$tip.label)+1):length(labels)]=""
    if(show){
    plot(tlayout,type="n",xlab="",ylab="",axes=FALSE,...)
    if(is.na(showedges)&& any(edges[,"weight"]<0.99)) showedges=TRUE
    else if(is.na(showedges)) showedges=FALSE
    if(!label.mixture){
        edges[edges[,"mix"]==1,"weight"]=""
    }
    if(!label.nonmixture){
        edges[edges[,"mix"]==0,"weight"]=""
    }
    if(length(lwd)<dim(edges)[1])lwd=rep(lwd,dim(edges)[1])
    for(i in 1:dim(edges)[1]) {
        tt=tlayout[as.numeric(edges[i,1:2]),"t"]
        tp=tlayout[as.numeric(edges[i,1:2]),"p"]
        if(format[1]=="triangular"){
            arrows(tt[1],
                   tp[1],
                   tt[2],
                   tp[2],
                   col=arrows.col[edges[i,5]+1],lwd=lwd[i],
                   length=arrows.length)
            if(showedges) text(mean(tt),mean(tp),edges[i,"weight"],adj=adj.edge,
                               cex=cex.edge.text,col=text.col[edges[i,5]+1])
        }else{
            arrows(tt[1]-rdelta2*edges[i,5],
                   tp[1],
                   tt[1]-rdelta2*edges[i,5],
                   tp[2]+rdelta*edges[i,5],
                   col=arrows.col[edges[i,5]+1],lwd=lwd[i],
                   length=0)
            arrows(tt[1]-rdelta2*edges[i,5],
                   tp[2]+rdelta*edges[i,5],
                   tt[2],
                   tp[2]+rdelta*edges[i,5],
                   col=arrows.col[edges[i,5]+1],lwd=lwd[i],
                   length=arrows.length)
            if(showedges) text(tt[1],tp[2]+rdelta*edges[i,5]+vadj.edge,
                               edges[i,"weight"],adj=adj.edge,
                               cex=cex.edge.text,col=text.col[edges[i,5]+1])
        }
    }
    text(tlayout[,"t"]+textdelta,tlayout[,"p"],labels=labels,col=labels.col,adj=adj.node,cex=cex.labels)
    if(showaxis) axis(1)
    }
    tlayout=cbind(as.data.frame(tlayout),label=labels)
    order=tlayout[order(tlayout[1:g$n,"p"],decreasing=F),"label"]
    invisible(list(layout=tlayout,edges=edges,order=order))
}

plot.cglist=function(g,main="",...){
    if(length(main)!=length(g)) main=rep(main[1],length(g))
    r=lapply(1:length(g),function(i){
        plot(g[[i]],main=main[i], ...)
    })
    invisible(r)
}

plotold.cg<-function(g,digits=2,...){
    ## Plot using igraph
##    g=assignlocation(g)
    require(igraph)
##    require(visNetwork)
    edges=edges.cg(g)
    edges[,"weight"]=format(edges[,"weight"],digits=digits)
    edges[,"length"]=format(edges[,"length"],digits=digits)
##    tu=lapply(1:length(g$nl),function(x)tipsunder(g,x))
##    ntipsunder=sapply(tu,length)
    times=sapply(rev(g$nl),function(x)x$t)
##    graph.empty(n=length(g$nl), directed=TRUE)
    graph=graph_from_data_frame(edges, directed = TRUE)
    ## graph=set_vertex_attr(graph, "time", index = V(graph), times)
    ## topo_sort(graph)
##    torder=as.numeric(V(graph)$name)
##    map=cbind(i=1:length(torder),j=torder)
##    sapply(1:length(torder),function(x)
    lm=layout_with_sugiyama(graph)$layout
##    lm[,2]=times
    #lm[,1]=topo_sort(graph)
##    igv <- visIgraph(igraph = graph, layout = "layout_with_sugiyama")
    ##    visNodes(igv, shape = "circle", value = 30)
##    tlayout=locationasmatrix(g)
    ##    lm=lm[map[order(map[,2]),1],]
##    tlayout=tlayout[map[order(map[,1]),2],]
    plot(graph,edge.label = edges$length,layout = lm)#,...)
    return(invisible(graph))
}

checkgraph=function(g){
    ## Check that every node thinks that its parent is the same as the parent thinks about its child...
    if(is(g,"cglist")) return(checkgraph(g[[1]]))
    ok=TRUE
    ## Make sure there are no cycles
    tmp<-try(nodesunder(g,g$root),silent=TRUE)
    if(is(tmp,"try-error")){
        ok=FALSE
        print(paste("ERROR: circularity for node (?)"))
        print(attr(tmp,"condition"))
    }
    ## make sure that the root is not a mixture node
    ## if(g$nl[[g$root]]$type=="mixture"){
    ##     ok=FALSE
    ##     print(paste("ERROR: The root node",g$root,"is a mixture node!"))
    ## }
    ## Make sure that all the edges have their attributes
    tedges=edges.cg(g)
    if(any(is.na(tedges))){
        ok=FALSE
        print(paste("ERROR: some edges are missing attributes!"))
        print(tedges)
    }
    ## Checks for local structure
    usednodes=1:length(g$nl)
    usednodes=usednodes[!usednodes%in%g$spare]
    for(i in usednodes){
        ## Checks for targets of mixture nodes
        ## if((!is.na(g$nl[[i]]$pr)) && ( g$nl [[g$nl[[i]]$pl]]$type=="mixture" )){
        ##     ok=FALSE
        ##     print(paste("ERROR: node",i,"has two mixture edge parents!"))
        ## }
        ## Checks of parent/child matches
        if(!is.na(g$nl[[i]]$cl)){
            ii=g$nl[[i]]$cl
            if(!(i %in% na.omit(c(g$nl[[ii]]$pl,g$nl[[ii]]$pr))) ){
                ok=FALSE
                print(paste("ERROR: node",i,"has left child",ii,"but that has left parent",g$nl[[ii]]$pl))
            }
        }
        if(!is.na(g$nl[[i]]$cr)){
            ii=g$nl[[i]]$cr
            if(!(i %in% na.omit(c(g$nl[[ii]]$pl,g$nl[[ii]]$pr))) ){
                ok=FALSE
                print(paste("ERROR: node",i,"has right child",ii,"but that has left parent",g$nl[[ii]]$pl))
            }
        }
        if(!is.na(g$nl[[i]]$pl)){
            ii=g$nl[[i]]$pl
            if(!(i %in% c(g$nl[[ii]]$cl,g$nl[[ii]]$cr))){
                ok=FALSE
                print(paste("ERROR: node",i,"has left parent",ii,"but that has children ",c(g$nl[[ii]]$cl,g$nl[[ii]]$cr)))
            }
            iii=na.omit(c(g$nl[[i]]$cl,g$nl[[i]]$cr))
            if(ii %in% iii){
                ok=FALSE
                print(paste("ERROR: node",i,"has left parent",ii,"which is also its child!"))
            }
        }
        if(!is.na(g$nl[[i]]$pr)){
            ii=g$nl[[i]]$pr
            if(!(i %in% c(g$nl[[ii]]$cl,g$nl[[ii]]$cr))){
                ok=FALSE
                print(paste("ERROR: node",i,"has right parent",ii,"but that has children ",c(g$nl[[ii]]$cl,g$nl[[ii]]$cr)))
            }
        }
    }
    return(invisible(ok))
}


##############################################
## Everything to do with computing induced covarianecs
ccov_tree<-function(g){
    ## Compute the induced covariance for a tree-like g
    ## Weightings are ignored
    c=matrix(0,ncol=g$n,nrow=g$n)
    for(i in g$internal){
        if(!is.na(g$nl[[i]]$d)){
            tu=tipsunder(g,i)
            c[tu,tu]=c[tu,tu]+g$nl[[i]]$d
        }
    }
    for(i in g$tips){
        c[i,i]=c[i,i]+g$nl[[i]]$d
    }
    c
}
cenumerate_trees<-function(g){
    ## Enumerate all trees that can be induced by different switches of mixture edges
    allvals=lapply(g$mix,function(x){c(0,1)})
    names(allvals)=g$mix
    expand.grid(allvals)
}
cremovechild<-function(pnode,child){
    ## Provide a NODE pnode and an INDEX child
    ## Returns pnode with the child removed
    isleftchild = (pnode$cl==child)
    if(isleftchild){ # Remove and make the right child into the left child
        pnode$cl = pnode$cr
        pnode$cr=NA
    }else{
        pnode$cr=NA
    }
    pnode
}

csetgraphto_tree<-function(g,nfix){
    ## Force mixture edges to be fixed, inducing a tree topology and removing tree weightings
    tf=as.numeric(colnames(nfix)) # extract the target
    fval=nfix[1,] # Remove the matrix info
    for(i in 1:length(tf)){
        tfi = g$nl[[ tf[i]  ]] # The node we are examining
        if(is.na(tfi$cr)){ ## if a child has two mixture parents, the right may have been moved to the left by a previous edge
            tfic = tfi$cl # The target node into which the mixture goes
        }else{
            tfic =  tfi$cr # The target node into which the mixture goes
        }
        targeti=g$nl[[ tfic ]]
        
        ## tf[[i]] is the index
        ## tfi is a copy so can only be used for extracting indices and data
        ## targeti is a copy of the child with index tfic
        
        isleftparent=(targeti$pl==tf[i]) # Check if it's the left parent
        if(is.na(isleftparent)) stop("!!")
        if(fval[i]==0) { ## Disable this mixture edge, ie disable cr and keep cl
            if(isleftparent){
                g$nl[[tfic]]$pl=NA ## t$cr is not a thing!!! XXXXXXXXXXXXXX
            }else {
                g$nl[[tfic]]$pr=NA
            }
            g$nl[[ tfic ]]$w=1 # remove any reweighting
            g$nl[[ tf[i] ]]$cr=NA
        }else{ # Follow this mixture edge, i.e. keep cr as cl, and disable cl
            op=targeti$pl
            if(isleftparent){ # Get rid of the other parent
                ## Remove the node as a child of the other parent
                g$nl[[op]] = cremovechild(g$nl[[op]],tfic)
                g$nl[[tfic]]$pr=NA
            }else { # swap to the left
                ## Remove the node as a child of the other parent
                g$nl[[op]] = cremovechild(g$nl[[op]],tfic)
                g$nl[[tfic]]$pl=g$nl[[tfic]]$pr
                g$nl[[tfic]]$pr=NA
            }
            g$nl[[tf[i]]]$type="split" # No mixture edges anymore
            g$nl[[ tfic ]]$d=0 # mixture edges have distance zero
            g$nl[[ tfic ]]$w=1 # remove any reweighting
        }
    }
    g               
}
cweight<-function(g,mnode){
    ## Extract the weight associated with a mixture node mnode
    ## This is stored in its child
    tchild=g$nl[[mnode]]$cr
    g$nl[[tchild]]$w
}

clistoftrees<-function(g){
    ## Returns the list of all induced trees
    alltrees=cenumerate_trees(g) # tree ids
    clist=lapply(1:dim(alltrees)[1],function(i){
        tg<-try(csetgraphto_tree(g,alltrees[i,,drop=FALSE]))
        if(class(tg)=="try-error")browser()
        tg
    })
    clist
}

cmultmatrix<-function(g,mixsetting){
    ## Takes a mixture DAG g
    ## And a configuration of mixture on/off settings, mixsetting
    ## which is a row vector of length #mixtures with column names giving the node labels of the mixture edges
    ## with values 0 or 1 for whether we follow them or not
    ## We then return the mixture weighting matrix
    cmult=matrix(1,g$n,g$n)
    allw=sapply(g$mix,cweight,g=g) # weights
    for(j in 1:length(allw)){
        ttips=tipsunder(g,g$nl[[g$mix[j] ]]$cr) ## tips under the mixture target
        tw=allw[j]
        if(mixsetting[j]==0) tw = 1 - tw
        cmult[ttips,ttips]= cmult[ttips,ttips] * tw
        cmult= cmult * tw
    }
    cmult 
}

cmultweight<-function(g,mixsetting){
    ## Takes a mixture DAG g
    ## And a configuration of mixture on/off settings, mixsetting
    ## which is a row vector of length #mixtures with column names giving the node labels of the mixture edges
    ## with values 0 or 1 for whether we follow them or not
    ## We then return the weight
    cmult=1
    allw=sapply(g$mix,cweight,g=g) # weights
    for(j in 1:length(allw)){
        tw=allw[j]
        if(mixsetting[j]==0) tw = 1 - tw
        cmult= cmult * tw
    }
    cmult 
}

cenumerate_weights<-function(g,alltrees=NULL){
    ## Enumerate all tree weighting matrices and return them as a list
    if(all(is.null(alltrees))) alltrees = cenumerate_trees(g)
    sapply(1:dim(alltrees)[1],function(i)
        cmultweight(g,alltrees[i,,drop=FALSE])
        )
}

cenumerate_weightmatrices<-function(g,alltrees=NULL){
    ## Enumerate all tree weighting matrices and return them as a list
    if(all(is.null(alltrees))) alltrees = cenumerate_trees(g)
    lapply(1:dim(alltrees)[1],function(i)
        cmultmatrix(g,alltrees[i,,drop=FALSE])
        )
}

cenumerate_covmatrices<-function(g,alltrees=NULL){
    ## Enumerate all tree covariance matrices and return them as a list
    if(all(is.null(alltrees))) alltrees = cenumerate_trees(g)
    
    lapply(1:dim(alltrees)[1],function(i){
        tg<-try(csetgraphto_tree(g,alltrees[i,,drop=FALSE]))
        if(class(tg)=="try-error")browser()
        ccov_tree(tg)
    })
}

ccov_dag<-function(g){
    ## Enumerate all trees implied by a DAG and compute the covariance matrix for each.
    if(is(g,"cglist")){
        ret = lapply(g,ccov_dag)
        return(ret)
    }
    if(length(g$mix)==0){
        csum=ccov_tree(g)
        rownames(csum)=colnames(csum)=g$tip.label
        return(csum)
    }
    alltrees=cenumerate_trees(g) # tree ids
    allmatrices=cenumerate_weightmatrices(g,alltrees)
    allcovs=cenumerate_covmatrices(g,alltrees)
    ## allw=sapply(g$mix,cweight,g=g) # weights
    clist=lapply(1:dim(alltrees)[1],function(i){
        c=allcovs[[i]]
        cmult=allmatrices[[i]]
        c * cmult 
    })
    csum=Reduce('+',clist)
    rownames(csum)=colnames(csum)=g$tip.label
    csum
}



assignlocation<-function(g,mindepth=0,maxdepth=Inf,rightalign=FALSE){
    ret=matrix(NA,nrow=length(g$nl),ncol=2)
    nleft=0
    nnodes=length(g$nl)
    depth=0
    i=g$root
    node=g$root
    for(i in 1:length(g$nl)) g$nl[[i]]$tmp=0 # tracks completed
    mydepth=function(x){
        ifelse(is.na(x),0,min(max(mindepth,x),maxdepth))
    }
    maxtime=0
    while(TRUE){
        ## Situations:
        ## 1. We have a left and need to go there
        lefttodo = ((!is.na(g$nl[[node]]$cl)) && (g$nl[[g$nl[[node]]$cl]]$tmp==0))
        ## 2. We place ourself
        selftodo = (g$nl[[node]]$tmp==0)
        ## 3. Right to do
        righttodo = ((!is.na(g$nl[[node]]$cr)) &&
                     (g$nl[[g$nl[[node]]$cr]]$tmp==0) &&
                     (g$nl[[node]]$type!="mixture") )
        ## 4. Return to parent if we have one
        atroot = (g$root==node)
        ## 5. Stop if -back- at root
        if(lefttodo){
            node=g$nl[[node]]$cl
            depth=depth+mydepth(g$nl[[node]]$d)
            next;
        }else if(selftodo) {
#            print(paste("doing node",node,"to",nleft,",",depth))
            g$nl[[node]]$tmp=1
            g$nl[[node]]$t=depth
            g$nl[[node]]$p=nleft #/nnodes
            maxtime=max(maxtime,depth)
            if(is.na(g$nl[[node]]$pl)){ nleft=nleft + 1
            }else if(g$nl[[ g$nl[[node]]$pl ]]$type=="split") {  nleft=nleft + 1
##            }else if(g$nl[[ node ]]$type=="split") {  nleft=nleft + 1
            }
            next;
        }else if(righttodo){
            node=g$nl[[node]]$cr
            depth=depth+mydepth(g$nl[[node]]$d)
            next;            
        }else if(!atroot){
            depth=depth - mydepth(g$nl[[node]]$d)
            node=g$nl[[node]]$pl
            next;
        }else{
            break;
        }
    }
    if(rightalign){
        for(i in g$tips){
            g$nl[[i]]$t=maxtime
        }
    }
    g
}


###########
## Parameters:
## For every branch there is a drift parameter
## For every mixture node there is an alpha (just another drift parameter) and a weight w
## So for every node there is a drift parameter plus for every mixture edge a weight
mixscale<-function(x,inv=FALSE){
    if(!inv){ ## Take R -> [0,1]
##        r=1/(1+exp(-x)) ## Take R -> [0,1]
        return(1/(1+exp(-x)))
#        if(r<0.001) r=0.001
#        if(r>0.999) r=0.999
        return(r)
    }else{ ## Take [0,1] -> R
#        if(x<0.001) x=0.001
#        if(x>0.999) x=0.499
##        r = -log(1/(x)-1) ## Take [0,1] -> R
        return(-log(1/(x)-1))
    }
}
driftscale<-function(x,inv=FALSE){
    if(!inv){ ## Take R to R+
        return(exp(x))
    }else{ ## Take R+ to R
        return(log(x))
    }    
}
getp<-function(g){
    if(is(g,"cglist")) return(getp(g[[1]]))
    ## Extracts the indices of parameters of each type
    driftp=(1:length(g$nl))
    driftp=driftp[!driftp%in%g$root]
    if(length(g$mix)>0)
        for(i in g$mix) {
            driftp=driftp[!driftp%in%g$nl[[i]]$cl]
        }
    mixp=g$mix
    list(drift=driftp,mix=mixp)
}
npars<-function(g){
    ## Extract the number of each class of parameter, plus the total
    if(is(g,"cglist")){
        np=npars(g[[1]])
        return(np)
    }
    p<-getp(g)
    ret=c(nd=length(p$drift),
          nm=length(p$mix),
          tot=0)
    ret["tot"]=ret["nd"]+ret["nm"]
    return(ret)
}
transformpars<-function(g,pars,inv=F){
    ## Transform parameters from R to their required ranges
    if(is(g,"cglist")){
        np=npars(g)
        res=lapply(1:length(g),function(i){
            tpars=pars[(i-1)*np["tot"]+(1:np["tot"])]
            transformpars(g[[i]],tpars,inv)
        })
        return(do.call("c",res))
    }
    np=npars(g)
    if(np["tot"] != length(pars)) stop("Invalid parameterisation")
    for(i in 1:np["nd"]) pars[i] = driftscale(pars[i],inv=inv)
    if(np["nm"]>0) for(i in 1:np["nm"]) pars[np["nd"]+i] = mixscale(pars[np["nd"] + i],inv=inv)
    pars
}

parameterise<-function(g,pars,transform=TRUE,n=NA){
    ## Take parameters and put them in their correct place in g
    ## Optionally, treat these parameters as coming from R and transform them to their required ranges
    ## If we provide enough pars to parameterise multiple graphs, we return a list of them, or optionally via n, a specific one
    p=getp(g)
    np=npars(g)
    ng=length(pars)/np["tot"]
    if(floor(ng)!=ng) stop("Invalid parameterisation")
    if(!is(g,"cglist")){
        g=list(g)
        class(g)="cglist"
    }
    ret=lapply(1:ng,function(i){
        tg=g[[i]]
        tpars=pars[(i-1)*np["tot"] + (1:np["tot"]) ]
        if(transform) tpars=transformpars(tg,tpars)
        for(j in 1:np["nd"])  tg$nl[[ p$drift[j] ]]$d = tpars[j]
        if(np["nm"]>0) for(j in 1:np["nm"]) tg$nl[[ tg$nl[[p$mix[j]]]$cr ]]$w = tpars[j+np["nd"]]
        tg
    })
    for(i in 1:ng) if(length(g[[i]]$mixparmap)>0){
                       ii= g[[i]]$mixparmap
                       for(j in 1:np["nm"]) {
                           ret[[i]]$nl[[ ret[[i]]$nl[[p$mix[j]]]$cr ]]$w =
                               ret[[ii]]$nl[[ ret[[ii]]$nl[[p$mix[j]]]$cr ]]$w
                       }
                   }
    if(length(ret)==1) {
        if(is.na(n)) n=1
        return(ret[[n]])
    }
    class(ret)="cglist"
    return(ret)
}
gparvec<-function(g,invtrans=FALSE){
    ## Extract the parameter vector from g
    ## Optionally transform this into R
    if(is(g,"cglist")){
        r=lapply(g,gparvec,invtrans=invtrans)
        return(do.call("c",r))
    }
    p=getp(g)
    np=npars(g)
    pars=numeric(np["tot"])
    for(i in 1:length(p$drift)) pars[i] = g$nl[[ p$drift[i] ]]$d 
    if(np["nm"]>0) for(i in 1:np["nm"])  pars[i+np["nd"]] = g$nl[[ g$nl[[p$mix[i]]]$cr ]]$w
    if(invtrans) pars=transformpars(g,pars,T)
    pars
    
}
c_Center<- function (Y, col = TRUE, row = TRUE, method = "mean") 
{
    ret = Y
    if (col) {
        Cn1 = c_getC(dim(Y)[1])
        if (method == "matrix") {
            ret = Cn1 %*% Y
        }
        else if (method == "mean") {
            ret = ret - rowMeans(ret)
        }
        else stop("Invalid method: must be \"matrix\" or \"mean\"")
    }
    if (row) {
        Cn2 = c_getC(dim(Y)[2])
        if (method == "matrix") {
            ret = ret %*% Cn2
        }
        else if (method == "mean") {
            ret = t(t(ret) - colMeans(ret))
        }
        else stop("Invalid method: must be \"matrix\" or \"mean\"")
    }
    rownames(ret) = rownames(Y)
    colnames(ret) = colnames(Y)
    ret
}

ctree_loss2<-function(gref,dataref,center=FALSE){
    ## Evaluate the loss for a parameterised graph of class cg, NOT a cglist
    pred<-ccov_dag(gref)
    if(center){
        pred=c_Center(pred)
        dataref=c_Center(dataref)
    }
    dataref=dataref[gref$tip.label,gref$tip.label]
    loss<-sum(na.omit(as.numeric((pred-dataref)^2)))
}

ctree_loss_raw<-function(gref,dataref,center=FALSE){
    ## Report the loss for a PARAMETERISED graph
    if(is(dataref,"list")){
        losslist<-lapply(1:length(dataref),function(i){
            ctree_loss2(gref[[i]],dataref[[i]],center)
        })
        loss = Reduce("+",losslist)
    }else{
        loss = ctree_loss2(gref,dataref,center)
    }
    loss
}

ctree_loss<-function(pars,gref,dataref,transform=TRUE,center=FALSE){
    ## Evaluate the loss for a set of parameters pars
    ## Which are optionally treated as coming from R and hence transformed
    ## The loss is the mean square error to a reference dataset
    if(is(gref,"cglist")){
        checkgraph(gref[[1]])
    }else{
        checkgraph(gref)
    }
    gref<-parameterise(gref,pars,transform)
    if(is(gref,"cglist")){
        if((!is(dataref,"list")) || (length(gref) != length(dataref))){
            stop("Mismatch betweek number of graphs and parameters!")
        }
    }
    loss = ctree_loss_raw(gref,dataref,center)
    return(loss)
}


############################################################
############################################################
############################################################
## Everything about inferring topology
mypruneregraft=function(g,source,target,careful=TRUE){
    ## Run this ONLY on type=="split" nodes!
    gstart=g
    psourceroot=FALSE
    ptargetroot=FALSE
    ## Parent of source needs removing
    ptarget=g$nl[[target]]$pl ## 
    psource=g$nl[[source]]$pl ## This node is going to be moved
    ocsource=c(g$nl[[psource]]$cl,g$nl[[psource]]$cr) # other child of this parent
    ocsource=ocsource[ocsource!=source]
    ppsource=g$nl[[psource]]$pl # parent of the parent of the source, will become parent to other child of this parent
    if(ocsource==target) {
        print(paste("source",source,"and target",target,"have same parent",psource))
        return(g)
    }
    if(is.na(ppsource) ) {
        print(paste("NOTE: parent of source is the root ( node",psource,")"))
        psourceroot=TRUE
        ## psource will become internal
        g$nl[[psource]]$w = 1
        ## ocsource will become the root
        g$nl[[ocsource]]$d = NA
        g$nl[[ocsource]]$pl = NA
        g$root=ocsource
    }else{
        ## Update the other child of the parent of the source
        g$nl[[ocsource]]$d = g$nl[[ocsource]]$d + g$nl[[psource]]$d
        g$nl[[ocsource]]$pl = ppsource
        ## Update the parent of parent of the source
        psourceisleftchild = (g$nl[[ppsource]]$cl==psource)
        if(psourceisleftchild){
            g$nl[[ppsource]]$cl=ocsource
        }else{
            g$nl[[ppsource]]$cr=ocsource
        }
    }
    ## Update the parent of the target
    if(is.na(ptarget)) print(paste("parent of target is the root ( node",ptarget),")")
    targetisleftchild = (g$nl[[ptarget]]$cl==target)
    if(targetisleftchild){
        g$nl[[ptarget]]$cl=psource
    }else{
        g$nl[[ptarget]]$cr=psource
    }
    ## Update the split node, the old parent node
    g$nl[[psource]]$d=g$nl[[target]]$d/2
    g$nl[[psource]]$cl=target
    g$nl[[psource]]$cr=source
    g$nl[[psource]]$pl=ptarget
    ## Update the target and source
    g$nl[[target]]$d=g$nl[[target]]$d/2
    g$nl[[target]]$pl=psource    
    g$nl[[source]]$pl=psource
    ## Done!
    if(careful) {
        if(!checkgraph(g)){
            print(paste("Moving",source,"to",target))
            print(g)
            stop("Created invalid graph!")
        }
    }
    g
}

randomnodes<-function(g,n=2,root=FALSE,mix=FALSE,tips=TRUE,internal=TRUE){
    ## Return n random nodes of g satisfying the criterion
    valid=numeric()
    if(tips) valid=c(valid,g$tips)
    if(internal) valid=c(valid,g$internal)
    if(mix) valid=c(valid,g$mix)
    if(!root) valid=valid[!(valid%in%g$root)]
    if(length(valid)<n) stop(paste("Error: requested",n,"nodes but only",length(valid),"meet the criterion!"))
    if(length(valid)==n) return(valid) # sample behaves badly if length(valid)==1
    sample(valid,n)
}

isvalidregraftpair<-function(g,swap){
    ## Asks: is swap[1] -> swap[2] a valid prune pair?
    done=TRUE
    ## Reject if:
    ## nodes share a parent
    if((!is.na(g$nl[[swap[1]]]$pl))&&(!is.na(g$nl[[swap[2]]]$pl))) if(g$nl[[swap[1]]]$pl==g$nl[[swap[2]]]$pl) done=FALSE
    ## The target is the parent of the source
    if((!is.na(g$nl[[swap[1]]]$pl)) && (g$nl[[swap[1]]]$pl==swap[2])) done=FALSE
    ## The target is a descendent of the source
    if(swap[2] %in% nodesunder(g,swap[[1]])) done=FALSE
    ## The source is a descendent of the target
    if(swap[1] %in% nodesunder(g,swap[[2]])) done=FALSE
    ## Either is an outgroup
    if(length(g$outgroup)>0){
        tog=which(g$tip.label==g$outgroup)
        if(any(swap %in%tog)) done=FALSE
        ## The root does not have the outgroup
        if((g$nl[[g$root]]$cl !=tog) &&
           (g$nl[[g$root]]$cr != tog)) done=FALSE
    }
    return(done)
}

randomregraftpair<-function(g,maxtries=400,...){
    ## Return random nodes that are not siblings
    ## Careful not to go crazy if there are no valid options...
    done=FALSE
    ntries=0
    while(!done){
#        print(paste("... random prune pair try",ntries)) ## DEBUG
        ret=randomnodes(g,...)
#        print(paste(ret,collapse=",")) ## DEBUG
        for(i in 1:length(ret)) {
            while(TRUE){
                ii=ret[i]
                n=g$nl[[ii]]
                targetsparentismixture= (!is.na(n$pl)) && (g$nl[[n$pl]]$type=="mixture")
                if(targetsparentismixture) ret[i]=n$pl
                else break;
            }
        }
        done=isvalidregraftpair(g,ret)
        ntries=ntries+1
        
        if(ntries==maxtries) stop("Error! Reached maximum number of attempts to find valid tree move!")
    }
    return(ret)
}

isvalidregraftmixture<-function(g,rem,ret){
    ## Asks if removing a node rem (can be NA for no removal) and then adding a mixture edge from ret[1] to ret[2] is valid?
    done=TRUE
    ## Is the proposal sound? We need:
    parsource=g$nl[[ret[1]]]$pl
    partarget=g$nl[[ret[2]]]$pl
    ## clsource=g$nl[[ret[1]]]$cl
    ## crsource=g$nl[[ret[1]]]$cr
    ## cltarget=g$nl[[ret[2]]]$cl
    ## crtarget=g$nl[[ret[2]]]$cr

    add=all(is.na(rem))
    if(!add){
        if(parsource==rem) parsource=g$nl[[rem]]$pl
        if(partarget==rem) partarget=g$nl[[rem]]$pl
        ## Reject if (after removal of the mixture node):
        if(rem %in% ret) done=FALSE
        g=removemixedge(g,rem)
        rem=NA
    }
    ## nodes share a parent
    if((!is.na(parsource)) && (!is.na(partarget)) && (parsource==partarget) ) done=FALSE
    ## The target is the parent of the source
    if((!is.na(parsource)) && (parsource==ret[2]) ) done=FALSE
#######  CARE HERE:
    ## The target is the child of the source
    if((!is.na(partarget)) && (partarget==ret[2]) ) done=FALSE
#######  CARE HERE:
    ## The target is the child of the source (via a mixture edge)
    if((!is.na(g$nl[[ret[2]]]$pr)) &&
       (g$nl[[ret[2]]]$pr==ret[1])) done=FALSE
#######        
    ## The target is not already a mixture node target
    if((!is.na(g$nl[[ret[2]]]$pr))){
        if(is.na(rem) ||  (g$nl[[ret[2]]]$pr!=rem))
            done=FALSE
    }
    ## The targets parent is not a mixture node
##    if((!is.na(partarget)) && (g$nl[[partarget]]$type=="mixture")) done=FALSE
    ## The target is a descendent of the source
    ##        if(ret[2] %in% nodesunder(g,ret[[1]])) done=FALSE 
    ## The source is a descendent of the target
    if(ret[1] %in% nodesunder(g,ret[[2]])) done=FALSE

    ## The source or the target is a forbidden outgroup
    ## Either is an outgroup
    if(length(g$outmix)>0){
        tog=which(g$tip.label==g$outmix)
        if(any(ret%in%tog)) done=FALSE
    }    
    ## Done...
    return(done)
}

randomregraftmixture<-function(g,add=FALSE,maxtries=200,...){
    ntries=0
    done=FALSE
    while(!done){
        if(add){
            rem=NA
        }else{
        ## Which node to remove?
            if(length(g$mix)==1) rem=g$mix
            else rem=sample(g$mix,1)
        }
        ## Which mixture to propose?
        ret=randomnodes(g,...)
        print(c(rem,ret))
        done=isvalidregraftmixture(g,rem,ret)
        if(ntries==maxtries) stop("Error! Reached maximum number of attempts to find valid tree move!")
    }
    
    alpha=0 # runif(1,0.01,0.5)
    w=runif(1,0.01,0.5)
    ret=list(rem=rem,source=ret[1],target=ret[2],alpha=alpha,w=w)
    return(ret)
}

reversemixture<-function(g,mixrev){
    ####### INCOMPLETE
    if(is(g,"cglist")){
        ret=lapply(g,reversemixture,mixrev=mixrev)
        class(ret)="cglist"
        return(ret)
    }
    ## Swap the role of mixrev (a mixture node) and the other parent of the target
    target=g$nl[[mixrev]]$cr
    opar= g$nl[[target]]$pl
    if(opar==target) stop("reversemixture problem")
    ## Swap at the parents of the target
    isleftchild=(g$nl[[opar]]$cl==target)
    if(isleftchild){
        ochild=g$nl[[opar]]$cr
    }else{
        ochild=g$nl[[opar]]$cl
    }
    ## Find out the parent of the other child
    paropar=g$nl[[opar]]$pl
    ## Find out the parent of the mixrev
    parmix=g$nl[[mixrev]]$pl
    g$nl[[ochild]]$pl=paropar
    parisleftchild=(g$nl[[paropar]]$cl==opar)
    mixisleftchild=(g$nl[[parmix]]$cl==mixrev)
    ## Reverse the children of the parents
    if(parisleftchild){
        g$nl[[paropar]]$cl=mixrev
    }else{
        g$nl[[paropar]]$cr=mixrev
    }
    g$nl[[mixrev]]$pl=paropar
    if(mixisleftchild){
        g$nl[[parmix]]$cl=opar
    }else{
        g$nl[[parmix]]$cr=opar
    }
    g$nl[[opar]]$pl=paropar
    
    ## tpl=g$nl[[target]]$pl
    ## tpr=g$nl[[target]]$pr
    ## ## Swap at the target
    ## g$nl[[target]]$pl=tpr
    ## g$nl[[target]]$pr=tpl
    g
}

infergraphpar<-function(g,ctree_loss,
                       data,
                       lower=-5,
                       method="L-BFGS-B",
                       control=defaultcontrol,
                       ...){
    ## Infer a graph's best parameters
    p=gparvec(g,invtrans=TRUE)
    lower=rep(lower,length(p))
    opt=optim(p,ctree_loss,
                gref=g,
                method=method,
                dataref=data,lower=lower,
                control=control,...)
    g=parameterise(g,opt$par)
    list(g=g,par=opt$par,loss=opt$value)
}

mypruneregraftstep<-function(g){
    ## Do a complete prune/regraft step
    if(is(g,"cglist")){
        swap=randomregraftpair(g[[1]])
        gtest=lapply(g,function(x)mypruneregraft(x,swap[1],swap[2]))
        class(gtest)="cglist"
    }else{
        swap=randomregraftpair(g)
        gtest=mypruneregraft(g,swap[1],swap[2])
    }
    list(g=gtest,swap=swap,mixswap=NA)
}

myregraftmixedgestep<-function(g,add=FALSE){
    if(is(g,"cglist")){
        swap=randomregraftmixture(g[[1]],add=add)
        gtest=lapply(g,function(x){
            myregraftmixedge(x,swap$rem,swap$source,swap$target,swap$alpha,swap$w)})
        class(gtest)="cglist"
    }else{
        swap=randomregraftmixture(g,add=add)
        gtest=myregraftmixedge(g,swap$rem,swap$source,swap$target,swap$alpha,swap$w)
    }
    list(g=gtest,swap=NA,mixswap=swap)
}

dagstep<-function(g,data,control=defaultcontrol,freqs=c(1/3,1/3,1/3),verbose=FALSE,...){
    ## Do one iteration of the graph
    movetype=sample(1:3,1,prob=freqs)
    if((is(g,"cglist")&&(length(g[[1]]$mix)==0)) ||
       (is(g,"cg")&&(length(g$mix)==0)))movetype=1 # No mixture edges to worry about
    if(verbose) print(paste("Proposing move of type",movetype,"..."))
    if(movetype==1){
        proposal=mypruneregraftstep(g)
    }else if(movetype==2){        
        proposal=myregraftmixedgestep(g)
    }else if(movetype==3){
        proposal0=mypruneregraftstep(g)
        proposal=myregraftmixedgestep(proposal0$g)
        proposal$swap=proposal0$swap
    }else stop("Invalid move type in dagstep?!")
    if(verbose){
        print(paste("proposal: movetype",movetype,
              "swap:",paste(proposal$swap,collapse=","),
              "mixswap:",paste(proposal$mixswap,collapse=",")))
    }
    
    inf=infergraphpar(proposal$g,ctree_loss,data,control=control,...)
    proposal$inf=inf
    proposal
}

infer_dag<-function(g0,
                     data,
                     init="graph",
                     maxiter=100,
                     losstol=0.01,
                     control=defaultcontrol,
                     verbose=TRUE,
                     ...){
    ## Infer a dag by using prune/regraft moves
    if(init=="random") {
        p=runif(length(gparvec(g0)))
    }else if(init=="graph") {
        p=gparvec(g0,invtrans=TRUE)
    }else{
        stop("Invalid init provided")
    }
    if(verbose) print(paste("DAG Iteration 1 paramterising..."))
    last=infergraphpar(g0,ctree_loss,data,control=control,...)
    lastloss=last$loss
    losses=rep(NA,maxiter)
    losses[1]=lastloss
    if(maxiter>1) for(i in 2:maxiter){
        proposal=dagstep(last$g,data,control,verbose=verbose,...)
        if(!all(is.na(proposal$swap))){
            if(verbose) print(paste("DAG Iteration",i,"repruning",proposal$swap[1],
                                    "to",proposal$swap[2],"..."))
        }
        if(!all(is.na(proposal$mixswap))){
            if(verbose) print(paste("DAG Iteration",i,"replacing mixture node",proposal$mixswap$rem,"with",proposal$mixswap$source,"->",proposal$mixswap$target,
                                    "..."))
        }
        thisloss=proposal$inf$loss
        if(thisloss<lastloss){
            if(verbose) print(paste("DAG Iteration",i,"accepting swap, from loss",lastloss,"to loss",thisloss))
            lastloss=thisloss
            last=proposal$inf
        }else{
            if(verbose) print(paste("DAG Iteration",i,"rejected swap with loss",thisloss,", retaining loss",lastloss))
        }
        losses[i]=lastloss
        if(lastloss<losstol) break;
    }
    return(list(g=last$g,last=last,losses=losses))
}


randomgraphlist=function(ntips,nmix=0,ngraphs=1,...){
    ## Return a valid starting point for infer_dag
    ## With a specified number of tips, mixture edges, and replicates of the graph
    g=simCoal(ntips,...)
    if(nmix>0){
        for(i in 1:nmix){
            swap=randomregraftmixture(g,add=TRUE)
            g=myregraftmixedge(g,swap$rem,swap$source,swap$target,swap$alpha,swap$w)
        }
    }
    if(ngraphs==1) return(g)
    glist=list()
    for(i in 1:ngraphs) glist[[i]]=g
    class(glist)="cglist"
    return(glist)
}

######################
######################
######################
######################
######################
#############
readtm<-function(stem){
    ## Read a TreeMix graph
    v=read.table(paste0(stem,".vertices.gz"))[,1:10]
    for(i in c(1,6,7,8,9))v[,i]=as.character(v[,i])
    e=read.table(paste0(stem,".edges.gz"))
    for(i in c(1:2))e[,i]=as.character(e[,i])
    list(v=v,e=e)
 }

tm2cg<-function(x){
    ## Convert a TreeMix model to a ClarityGraph
    if(is(x,"list")){
        e=x[["e"]]
        v=x[["v"]]
    }
    if(is(x,"character")){
        x=readtm(stem)
        v=x[["v"]]
        e=x[["e"]]
    }
    if(any(is.null(e))) stop("Require e")
    if(any(is.null(v))) stop("Require v")

    ## Set up the connection between treemix node labels and ours
    rownames(v)=v[,1]
    ttip=which(v[,5]=="TIP")
    tint=which((v[,5]=="NOT_TIP")&(v[,4]=="NOT_MIG"))
    tmig=which(v[,4]=="MIG")
    troot=which(v[,3]=="ROOT")
    n=length(ttip)
    l=length(tint)
    m=length(tmig)
    link=as.data.frame(
        rbind(cbind(1:n,v[ttip,1]),
                    cbind(n+(1:l),v[tint,1]),
                    cbind(n+l+(1:m),v[tmig,1])
              ))
    link[,1]=as.numeric(link[,1])
    link[,2]=as.character(link[,2])
    rownames(link)=link[,2]
    
    #########
    ## Add all the nodes
    nodes=list()
    extraedges=c()
    typelink=c("MIG"="mixture","NOT_MIG"="split")
    for(i in 1:dim(link)[1]) nodes[[i]]=cnode(link[i,1],w=1,d=0,
                                              type=as.character(typelink[v[link[i,2],4]]))
    for(i in 1:dim(link)[1]) {
        ii=link[i,1]
        j=link[i,2]
        ## Get my parent
        te=which(e[,2]==v[j,1])
        tp=link[e[te,1],1]
        if(j %in% v[troot,1]){ # This is the root
            root=ii
        }else{
            ## Order not migration edge as the first one
            to=order(e[te,5],decreasing=TRUE)
            te=te[to]
            tp=tp[to]
            ## Deal with the left parent (always present)
            nodes[[ii]]$pl=tp[1] # Set this nodes parent
            nodes[[ii]]$d=e[te[1],3] # Set this nodes distance from parent
            isleftchild=is.na(nodes[[tp[1]]]$cl) ## (v[v[j,6],7]==j)
            if(isleftchild){
                nodes[[tp[1]]]$cl=ii
            }else{
                nodes[[tp[1]]]$cr=ii
            }
            ## Deal with the right parent (for migration edges)
            if(length(tp)>1){
                nodes[[ii]]$pr=tp[2] # Set this nodes parent
                nodes[[tp[2]]]$cr=ii # Set this nodes parent right child to this node
                nodes[[ii]]$w=e[te[2],4] # Set the weight
##                oc=nodes[[tp[2]]]$cl
##                nodes[[oc]]$w=1-nodes[[ii]]$w # Set the weight for the other child
            }
            ## They allow more than two parents
            if(length(tp)>2){
                extraedges=c(extraedges,te[-(1:2)])
            }
        }
    }
    
    #######
    tip.label=v[ttip,2]
    g=list(nl=nodes,
           tips=1:n,
           root=root,
           tip.label=tip.label,
           internal=(n+1):length(nodes),
           mix=numeric(),
           spare=numeric(),
           outgroup=numeric(),
           n=n)
    g$mix=which(sapply(g$nl,function(x)x$type)=="mixture")
    class(g)="cg"
    ## Fix any edges that couldn't be placed
    if(length(extraedges)>0){
        ## Remove mixture edges that are doing nothing
        te=e[extraedges,]
        te=cbind(te,oldparent=link[te[,1],1])
        te=cbind(te,newparent=sapply(te[,8],function(x)nodes[[x]]$cl))
        te=cbind(te,target=link[te[,2],1])
        te=cbind(te,alpha=sapply(1:length(extraedges),function(i){
            dp=g$nl[[te[i,"oldparent"]]]$d
            dc=g$nl[[te[i,"newparent"]]]$d
            ifelse(dc>0,dc/(dp+dc),0)
        }))
        trem=sort(te$oldparent,TRUE)
        for(rem in trem) {
            g=removemixedge(g,rem,FALSE)
        }
        warning(paste("Had to remove",length(trem),
                      "edges that create multi-furcating graphs!"))
        ## for(i in 1:length(extraedges)){
        ##     g=mixedge(g,te[i,"newparent"],te[i,"target"],te[i,"alpha"],te[i,4])
        ## }
    }

    g
}

############################################
############################################
