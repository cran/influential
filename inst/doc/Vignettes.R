## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(tibble.print_min = 4L, tibble.print_max = 4L)
set.seed(60)

## ----setup--------------------------------------------------------------------
library(influential)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(coexpression.data))

## ----g_dataframe--------------------------------------------------------------
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(d=MyData)        # Reconstructing the graph

## -----------------------------------------------------------------------------
class(My_graph)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(coexpression.adjacency)[,1:10])

## ----g_adj, eval=FALSE--------------------------------------------------------
#  MyData <- coexpression.adjacency        # Preparing the data
#  
#  My_graph <- graph_from_adjacency_matrix(MyData)        # Reconstructing the graph

## ----Vertices, eval=FALSE-----------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  My_graph_vertices <- V(My_graph)        # Extracting the vertices

## ----DC, eval=FALSE-----------------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  My_graph_degree <- degree(My_graph, v = GraphVertices, normalized = FALSE) # Calculating degree centrality

## ----BC, eval=FALSE-----------------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  My_graph_betweenness <- betweenness(My_graph, v = GraphVertices,    # Calculating betweenness centrality
#                                      directed = FALSE, normalized = FALSE)

## ----NC, eval=FALSE-----------------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  neighrhood.co <- neighborhood.connectivity(graph = My_graph,    # Calculating neighborhood connectivity
#                                             vertices = GraphVertices,
#                                             mode = "all")

## ----H_index, eval=FALSE------------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  h.index <- h_index(graph = My_graph,    # Calculating H-index
#                     vertices = GraphVertices,
#                     mode = "all")

## ----LH_index, eval=FALSE-----------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  lh.index <- lh_index(graph = My_graph,    # Calculating Local H-index
#                     vertices = GraphVertices,
#                     mode = "all")

## ----CI, eval=FALSE-----------------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  ci <- collective.influence(graph = My_graph,    # Calculating Collective Influence
#                            vertices = GraphVertices,
#                            mode = "all", d=3)

## ----CR, eval=FALSE-----------------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  cr <- clusterrank(graph = My_graph,    # Calculating ClusterRank
#                    vids = GraphVertices,
#                    directed = FALSE, loops = TRUE)

## ----cond.prob----------------------------------------------------------------
MyData <- centrality.measures        # Preparing the data

My.conditional.prob <- cond.prob.analysis(data = MyData,       # Assessing the conditional probability
                                          nodes.colname = rownames(MyData),
                                          Desired.colname = "BC",
                                          Condition.colname = "NC")

print(My.conditional.prob)

## ----double.cent.assess, eval=FALSE-------------------------------------------
#  MyData <- centrality.measures        # Preparing the data
#  
#  My.metrics.assessment <- double.cent.assess(data = MyData,       # Association assessment
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

## ----double.cent.assess.noRegr., eval=FALSE-----------------------------------
#  MyData <- centrality.measures        # Preparing the data
#  
#  My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,       # Association assessment
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

## ----IVI.from.indices---------------------------------------------------------
MyData <- centrality.measures        # Preparing the data

My.vertices.IVI <- ivi.from.indices(DC = centrality.measures$DC,       # Calculation of IVI
                                   CR = centrality.measures$CR,
                                   NC = centrality.measures$NC,
                                   LH_index = centrality.measures$LH_index,
                                   BC = centrality.measures$BC,
                                   CI = centrality.measures$CI)

## ----IVI, eval=FALSE----------------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  My.vertices.IVI <- ivi(graph = My_graph, vertices = GraphVertices, # Calculation of IVI
#                         weights = NULL, directed = FALSE, mode = "all",
#                         loops = TRUE, d = 3, scaled = TRUE)

## ----Spreading.score, eval=FALSE----------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  Spreading.score <- spreading.score(graph = My_graph,     # Calculation of Spreading score
#                                     vertices = GraphVertices,
#                                     weights = NULL, directed = FALSE, mode = "all",
#                                     loops = TRUE, d = 3, scaled = TRUE)

## ----Hubness.score, eval=FALSE------------------------------------------------
#  MyData <- coexpression.data        # Preparing the data
#  
#  My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph
#  
#  GraphVertices <- V(My_graph)        # Extracting the vertices
#  
#  Hubness.score <- hubness.score(graph = My_graph,     # Calculation of Hubness score
#                                     vertices = GraphVertices,
#                                     directed = FALSE, mode = "all",
#                                     loops = TRUE, scaled = TRUE)

## ----SIRIR--------------------------------------------------------------------
set.seed(1234)
My_graph <- igraph::sample_gnp(n=50, p=0.05)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

Influence.Ranks <- sirir(graph = My_graph,     # Calculation of influence rank
                                   vertices = GraphVertices, 
                                   beta = 0.5, gamma = 1, no.sim = 10, seed = 1234)

knitr::kable(Influence.Ranks[c(order(Influence.Ranks$rank)[1:10]),])

