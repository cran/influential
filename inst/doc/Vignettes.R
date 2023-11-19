## ----include = FALSE--------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(tibble.print_min = 4L, tibble.print_max = 4L)
options(width=90)
options(rmarkdown.html_vignette.check_title = FALSE)
set.seed(60)

## ----echo=FALSE-------------------------------------------------------------------------
spaces <- function (n) {
  paste(rep("&nbsp;", n), collapse = "")
}

## ----fig.align="right", echo=FALSE, out.width="25%", out.extra='style="float:right; padding:10px"'----
knitr::include_graphics(path = "../man/figures/Symbol.png", error = FALSE)

## ----setup------------------------------------------------------------------------------
library(influential)

## ----exptl_data_fcor, eval=FALSE--------------------------------------------------------
#  
#  # Prepare a sample dataset
#  set.seed(60)
#  my_data <- matrix(data = runif(n = 10000, min = 2, max = 300),
#                         nrow = 50, ncol = 200,
#                         dimnames = list(c(paste("sample", c(1:50), sep = "_")),
#                                         c(paste("gene", c(1:200), sep = "_")))
#  )

## ----fcor_calc, eval=FALSE--------------------------------------------------------------
#  
#  # Calculate correlations between all pairs of genes
#  
#  correlation_tbl <- fcor(data = my_data,
#                          method = "spearman",
#                          mutualRank = TRUE,
#                          pvalue = "TRUE", adjust = "BH",
#                          flat = TRUE)

## ----echo=FALSE-------------------------------------------------------------------------
knitr::kable(head(coexpression.data))

## ----g_dataframe------------------------------------------------------------------------
# Preparing the data
MyData <- coexpression.data

# Reconstructing the graph
My_graph <- graph_from_data_frame(d=MyData)

## ---------------------------------------------------------------------------------------
class(My_graph)

## ----echo=FALSE-------------------------------------------------------------------------
knitr::kable(head(coexpression.adjacency, n=15)[10:15,10:15])

## ----g_adj, eval=FALSE------------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.adjacency
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_adjacency_matrix(MyData)

## ----echo=FALSE-------------------------------------------------------------------------
set.seed(60)
My_Data <- matrix(data = sample(c(0,1), replace = TRUE, size = 20), 
                  nrow = 4, ncol = 5,
                  dimnames = list(c(paste("cell", c(1:4), sep = "_")),
                                  c(paste("Gene", c(1:5), sep = "_"))))

knitr::kable(My_Data)

## ----g_inc, eval=FALSE------------------------------------------------------------------
#  # Reconstructing the graph
#  My_graph <- graph_from_adjacency_matrix(MyData)

## ----g_sif, eval=FALSE------------------------------------------------------------------
#  # Reconstructing the graph
#  My_graph <- sif2igraph(Path = "Sample_SIF.sif")
#  
#  class(My_graph)
#  #> [1] "igraph"

## ----Vertices, eval=FALSE---------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  My_graph_vertices <- V(My_graph)
#  
#  head(My_graph_vertices)
#  #> + 6/794 vertices, named, from 775cff6:
#  #> [1] ADAMTS9-AS2 C8orf34-AS1 CADM3-AS1   FAM83A-AS1  FENDRR      LANCL1-AS1

## ----DC, eval=FALSE---------------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculating degree centrality
#  My_graph_degree <- degree(My_graph, v = GraphVertices, normalized = FALSE)
#  
#  head(My_graph_degree)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>         172         121         168          26         189         176

## ----BC, eval=FALSE---------------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculating betweenness centrality
#  My_graph_betweenness <- betweenness(My_graph, v = GraphVertices,
#                                      directed = FALSE, normalized = FALSE)
#  
#  head(My_graph_betweenness)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>   21719.857   28185.199   26946.625    2940.467   33333.369   21830.511

## ----NC, eval=FALSE---------------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculating neighborhood connectivity
#  neighrhood.co <- neighborhood.connectivity(graph = My_graph,
#                                             vertices = GraphVertices,
#                                             mode = "all")
#  
#  head(neighrhood.co)
#  #>  ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>   11.290698    4.983471    7.970238    3.000000   15.153439   13.465909

## ----H_index, eval=FALSE----------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculating H-index
#  h.index <- h_index(graph = My_graph,
#                     vertices = GraphVertices,
#                     mode = "all")
#  
#  head(h.index)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>          11           9          11           2          12          12

## ----LH_index, eval=FALSE---------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculating Local H-index
#  lh.index <- lh_index(graph = My_graph,
#                     vertices = GraphVertices,
#                     mode = "all")
#  
#  head(lh.index)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>        1165         446         994          34        1289        1265

## ----CI, eval=FALSE---------------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculating Collective Influence
#  ci <- collective.influence(graph = My_graph,
#                            vertices = GraphVertices,
#                            mode = "all", d=3)
#  
#  head(ci)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>        9918       70560       39078         675       10716        7350

## ----CR, eval=FALSE---------------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculating ClusterRank
#  cr <- clusterRank(graph = My_graph,
#                    vids = GraphVertices,
#                    directed = FALSE, loops = TRUE)
#  
#  head(cr)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>   63.459812    5.185675   21.111776    1.280000  135.098278   81.255195

## ----cond.prob--------------------------------------------------------------------------
# Preparing the data
MyData <- centrality.measures        

# Assessing the conditional probability
My.conditional.prob <- cond.prob.analysis(data = MyData,       
                                          nodes.colname = rownames(MyData),
                                          Desired.colname = "BC",
                                          Condition.colname = "NC")

print(My.conditional.prob)

## ----double.cent.assess, eval=FALSE-----------------------------------------------------
#  # Preparing the data
#  MyData <- centrality.measures
#  
#  # Association assessment
#  My.metrics.assessment <- double.cent.assess(data = MyData,
#                                              nodes.colname = rownames(MyData),
#                                              dependent.colname = "BC",
#                                              independent.colname = "NC")
#  
#  print(My.metrics.assessment)
#  #> $Summary_statistics
#  #>         BC NC
#  #> Min.              0.000000000                   1.2000
#  #> 1st Qu.           0.000000000                  66.0000
#  #> Median            0.000000000                 156.0000
#  #> Mean              0.005813357                 132.3443
#  #> 3rd Qu.           0.000340000                 179.3214
#  #> Max.              0.529464720                 192.0000
#  #>
#  #> $Normality_results
#  #>                               p.value
#  #> BC    1.415450e-50
#  #> NC 9.411737e-30
#  #>
#  #> $Dependent_Normality
#  #> [1] "Non-normally distributed"
#  #>
#  #> $Independent_Normality
#  #> [1] "Non-normally distributed"
#  #>
#  #> $GAM_nonlinear.nonmonotonic.results
#  #>      edf  p-value
#  #> 8.992406 0.000000
#  #>
#  #> $Association_type
#  #> [1] "nonlinear-nonmonotonic"
#  #>
#  #> $HoeffdingD_Statistic
#  #>         D_statistic P_value
#  #> Results  0.01770279   1e-08
#  #>
#  #> $Dependence_Significance
#  #>                       Hoeffding
#  #> Results Significantly dependent
#  #>
#  #> $NNS_dep_results
#  #>         Correlation Dependence
#  #> Results  -0.7948106  0.8647164
#  #>
#  #> $ConditionalProbability
#  #> [1] 55.35386
#  #>
#  #> $ConditionalProbability_split.half.sample
#  #> [1] 55.90331

## ----double.cent.assess.noRegr., eval=FALSE---------------------------------------------
#  # Preparing the data
#  MyData <- centrality.measures
#  
#  # Association assessment
#  My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,
#                                                           nodes.colname = rownames(MyData),
#                                                           centrality1.colname = "BC",
#                                                           centrality2.colname = "NC")
#  
#  print(My.metrics.assessment)
#  #> $Summary_statistics
#  #>         BC NC
#  #> Min.              0.000000000                   1.2000
#  #> 1st Qu.           0.000000000                  66.0000
#  #> Median            0.000000000                 156.0000
#  #> Mean              0.005813357                 132.3443
#  #> 3rd Qu.           0.000340000                 179.3214
#  #> Max.              0.529464720                 192.0000
#  #>
#  #> $Normality_results
#  #>                               p.value
#  #> BC    1.415450e-50
#  #> NC 9.411737e-30
#  #>
#  #> $Centrality1_Normality
#  #> [1] "Non-normally distributed"
#  #>
#  #> $Centrality2_Normality
#  #> [1] "Non-normally distributed"
#  #>
#  #> $HoeffdingD_Statistic
#  #>         D_statistic P_value
#  #> Results  0.01770279   1e-08
#  #>
#  #> $Dependence_Significance
#  #>                       Hoeffding
#  #> Results Significantly dependent
#  #>
#  #> $NNS_dep_results
#  #>         Correlation Dependence
#  #> Results  -0.7948106  0.8647164
#  #>
#  #> $ConditionalProbability
#  #> [1] 55.35386
#  #>
#  #> $ConditionalProbability_split.half.sample
#  #> [1] 55.68163

## ----IVI.from.indices, eval=FALSE-------------------------------------------------------
#  # Preparing the data
#  MyData <- centrality.measures
#  
#  # Calculation of IVI
#  My.vertices.IVI <- ivi.from.indices(DC = MyData$DC,
#                                     CR = MyData$CR,
#                                     NC = MyData$NC,
#                                     LH_index = MyData$LH_index,
#                                     BC = MyData$BC,
#                                     CI = MyData$CI)
#  
#  head(My.vertices.IVI)
#  #> [1] 24.670056  8.344337 18.621049  1.017768 29.437028 33.512598

## ----IVI, eval=FALSE--------------------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculation of IVI
#  My.vertices.IVI <- ivi(graph = My_graph, vertices = GraphVertices,
#                         weights = NULL, directed = FALSE, mode = "all",
#                         loops = TRUE, d = 3, scale = "range")
#  
#  head(My.vertices.IVI)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>    39.53878    19.94999    38.20524     1.12371   100.00000    47.49356

## ----net.for.vis, eval=FALSE------------------------------------------------------------
#  # Reconstructing the graph
#  set.seed(70)
#  My_graph <-  igraph::sample_gnm(n = 50, m = 120, directed = TRUE)
#  
#  # Calculating the IVI values
#  My_graph_IVI <- ivi(My_graph, directed = TRUE)
#  
#  # Visualizing the graph based on IVI values
#  My_graph_IVI_Vis <- cent_network.vis(graph = My_graph,
#                                       cent.metric = My_graph_IVI,
#                                       directed = TRUE,
#                                       plot.title = "IVI-based Network",
#                                       legend.title = "IVI value")
#  
#  My_graph_IVI_Vis

## ----Spreading.score, eval=FALSE--------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculation of Spreading score
#  Spreading.score <- spreading.score(graph = My_graph,
#                                     vertices = GraphVertices,
#                                     weights = NULL, directed = FALSE, mode = "all",
#                                     loops = TRUE, d = 3, scale = "range")
#  
#  head(Spreading.score)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>   42.932497   38.094111   45.114648    1.587262  100.000000   49.193292

## ----Hubness.score, eval=FALSE----------------------------------------------------------
#  # Preparing the data
#  MyData <- coexpression.data
#  
#  # Reconstructing the graph
#  My_graph <- graph_from_data_frame(MyData)
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculation of Hubness score
#  Hubness.score <- hubness.score(graph = My_graph,
#                                     vertices = GraphVertices,
#                                     directed = FALSE, mode = "all",
#                                     loops = TRUE, scale = "range")
#  
#  head(Hubness.score)
#  #> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1
#  #>   84.299719   46.741660   77.441514    8.437142   92.870451   88.734131

## ----SIRIR, eval=FALSE------------------------------------------------------------------
#  # Reconstructing the graph
#  My_graph <-  sif2igraph(Path = "Sample_SIF.sif")
#  
#  # Extracting the vertices
#  GraphVertices <- V(My_graph)
#  
#  # Calculation of influence rank
#  Influence.Ranks <- sirir(graph = My_graph,
#                                     vertices = GraphVertices,
#                                     beta = 0.5, gamma = 1, no.sim = 10, seed = 1234)
#  

## ----exir.data, eval=FALSE--------------------------------------------------------------
#  # Prepare sample data
#  gene.names <- paste("gene", c(1:2000), sep = "_")
#  
#  set.seed(60)
#  tp2.vs.tp1.DEGs <- data.frame(logFC = rnorm(n = 700, mean = 2, sd = 4),
#                                FDR = runif(n = 700, min = 0.0001, max = 0.049))
#  
#  set.seed(60)
#  rownames(tp2.vs.tp1.DEGs) <- sample(gene.names, size = 700)
#  
#  set.seed(70)
#  tp3.vs.tp2.DEGs <- data.frame(logFC = rnorm(n = 1300, mean = -1, sd = 5),
#                                FDR = runif(n = 1300, min = 0.0011, max = 0.039))
#  
#  set.seed(70)
#  rownames(tp3.vs.tp2.DEGs) <- sample(gene.names, size = 1300)
#  
#  set.seed(80)
#  regression.data <- data.frame(R_squared = runif(n = 800, min = 0.1, max = 0.85))
#  
#  set.seed(80)
#  rownames(regression.data) <- sample(gene.names, size = 800)

## ----diff_data_assembl, eval=FALSE------------------------------------------------------
#  my_Diff_data <- diff_data.assembly(tp2.vs.tp1.DEGs,
#                                     tp3.vs.tp2.DEGs,
#                                     regression.data)
#  
#  my_Diff_data[c(1:10),]

## ----exptl_data, eval=FALSE-------------------------------------------------------------
#  set.seed(60)
#  MyExptl_data <- matrix(data = runif(n = 100000, min = 2, max = 300),
#                         nrow = 50, ncol = 2000,
#                         dimnames = list(c(paste("cancer_sample", c(1:25), sep = "_"),
#                                           paste("normal_sample", c(1:25), sep = "_")),
#                                         gene.names))
#  
#  # Log transform the data to bring them closer to normal distribution
#  MyExptl_data <- log2(MyExptl_data)
#  
#  MyExptl_data[c(1:5, 45:50),c(1:5)]

## ----condition.col, eval=FALSE----------------------------------------------------------
#  MyExptl_data <- as.data.frame(MyExptl_data)
#  MyExptl_data$condition <- c(rep("C", 25), rep("N", 25))

## ----ExIR, eval=FALSE-------------------------------------------------------------------
#  
#  #The table of differential/regression previously prepared
#  my_Diff_data
#  
#  #The column indices of differential values in the Diff_data table
#  Diff_value <- c(1,3)
#  
#  #The column indices of regression values in the Diff_data table
#  Regr_value <- 5
#  
#  #The column indices of significance (P-value/FDR) values in
#  # the Diff_data table
#  Sig_value <- c(2,4)
#  
#  #The matrix/data frame of normalized experimental
#  # data previously prepared
#  MyExptl_data
#  
#  #The name of the column delineating the conditions of
#  # samples in the Exptl_data matrix
#  Condition_colname <- "condition"
#  
#  #The desired list of features
#  set.seed(60)
#  MyDesired_list <- sample(gene.names, size = 500)  #Optional
#  
#  #Running the ExIR model
#  My.exir <- exir(Desired_list = MyDesired_list,
#                  cor_thresh_method = "mr", mr = 100,
#                  Diff_data = my_Diff_data, Diff_value = Diff_value,
#                  Regr_value = Regr_value, Sig_value = Sig_value,
#                  Exptl_data = MyExptl_data, Condition_colname = Condition_colname,
#                  seed = 60, verbose = FALSE)
#  
#  names(My.exir)
#  #> [1] "Driver table"         "DE-mediator table"     "Biomarker table"      "Graph"
#  
#  class(My.exir)
#  #> [1] "ExIR_Result"

## ----exir.vis, eval=FALSE---------------------------------------------------------------
#  My.exir.Vis <- exir.vis(exir.results = My.exir,
#                          n = 5,
#                          y.axis.title = "Gene")
#  
#  My.exir.Vis

## ----comp_manipulate, eval=FALSE--------------------------------------------------------
#  # Select which genes to knockout
#  set.seed(60)
#  ko_vertices <- sample(igraph::as_ids(V(My.exir$Graph)), size = 5)
#  
#  # Select which genes to up-regulate
#  set.seed(1234)
#  upregulate_vertices <- sample(igraph::as_ids(V(My.exir$Graph)), size = 5)
#  
#  Computational_manipulation <- comp_manipulate(exir_output = My.exir,
#                                                ko_vertices = ko_vertices,
#                                                upregulate_vertices = upregulate_vertices,
#                                                beta = 0.5, gamma = 1, no.sim = 100, seed = 1234)

