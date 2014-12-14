(* ::Package:: *)

(* Mathematica Package *)

BeginPackage["Genealogy`CoalescentTree`",{"Combinatorica`"}]
(* Exported symbols added here with SymbolName::usage *)  

CoalescentTree::usage = "CoalescentTree[tlist,elist] is a data structure for a coalescent tree. Time is measured backward in unit of  2 Ne, where Ne is the effective population size. We use same lables for tree leaves and  sequences, and they are set as integers from 1 to n (the number of sequences).  The internal nodes are labled by integers from n+1 to 2 n-1, according the occurent order of coalescent events. tlist is a sequence of coalescent times  starting from zero at tree leaves, and elist is a sequence of branch pars between which coalescence occurs."

RandomCoalescentTree::usage = "RandomCoalescentTree[n] "

ToCombinatoricaGraph::usage = "ToCombinatoricaGraph[tree]  "

SameCoalescentTreeQ::usage = "SameCoalescentTreeQ [tree1,tree2] returns True if tree1 is same as tree2, False otherwise."

TreeTopology::usage = "TreeTopology[tree]  "

TreeToAdjacencyList::usage = "TreeToAdjacencyList  "

TreeNodeParent::usage = "TreeNodeParent  "

SubCoalescentTree::usage = "SubCoalescentTree  "

TotalBranchHeight::usage = "TotalBranchHeight[tree]  "

BranchHeights::usage = "BranchHeights [tree]"

TreeHeight::usage = "TreeHeight[tree]  "

LeafNumber::usage = "LeafNumber[tree]  "

TreeCorrelation::usage = "TreeCorrelation  "

TreeRFDistance::usage="TreeRFDistance "

Begin["`Private`"] (* Begin Private Context *) 

RandomCoalescentEvents[n_Integer] :=
    Module[ {set = Range[n], subset,i},
        Table[
         subset = RandomSample[set, 2];
         set = Append[Complement[set, subset], n + i];
         Thread[{subset, n + i}], {i, n - 1}]
    ]

RandomCoalescentTimes[n_Integer] :=
    Module[ {k},
        PadLeft[Accumulate[
          Table[RandomReal[ExponentialDistribution[(k (k - 1))/2]], {k, n, 2, -1}]], 2 n - 1]
    ];

RandomCoalescentTree[n_Integer] :=
    CoalescentTree[RandomCoalescentTimes[n], RandomCoalescentEvents[n]]

SubCoalescentTree[tree_CoalescentTree, leaves_List] :=
    Module[ {tlist, elist, bool, i = 1, nodes},
        {tlist, elist} = List @@ tree;
        bool = Table[False, {Length[tlist]}];
        bool[[leaves]] = True;
        While[i <= Length[elist],
         Which[
           ! bool[[elist[[i, 1, 1]]]], 
           elist = Delete[elist, i] /. elist[[i, 1, 2]] -> elist[[i, 2, 1]],
           ! bool[[elist[[i, 2, 1]]]], 
           elist = Delete[elist, i] /. elist[[i, 1, 2]] -> elist[[i, 1, 1]],
           True, 
           bool[[elist[[i, 1, 2]]]] = True;
           i++
         ];
    	];
        nodes = Union[Flatten[elist]];
        CoalescentTree[tlist[[nodes]], elist /. Thread[nodes -> Range[Length[nodes]]]]
    ]

SameCoalescentTreeQ[CoalescentTree[tlist1_List, elist1_List], 
  CoalescentTree[tlist2_List, elist2_List]] :=
    TrueQ[tlist1 == tlist2 && ((Sort[#] & /@ elist1) == (Sort[#] & /@ elist2))]
    
(*Root denoted by -1*)
TreeToAdjacencyList[CoalescentTree[_List, elist_List]] :=
    Append[SortBy[Flatten[elist, 1], First][[All, 2]], -1]
 
BranchHeights[CoalescentTree[tlist_List, elist_List]] :=
    Subtract @@ Reverse[Transpose[Partition[tlist[[Flatten[elist]]], 2]]]
 
TreeHeight[trees : {_CoalescentTree ..}] :=
    trees[[All, 1, -1]] - trees[[All, 1, 1]]
TreeHeight[CoalescentTree[tlist_List, _List]] :=
    Last[tlist] - First[tlist]
 
TotalBranchHeight[trees : {_CoalescentTree ..}] :=
    Module[ {tlists, elists, temp},
        tlists = trees[[All, 1]];
        elists = trees[[All, 2]];
        temp = MapThread[Partition[Part[#1, #2], 2] &, {tlists, Flatten[#] & /@ elists}];
        Total[temp[[All, All, 2]] - temp[[All, All, 1]], {2}]
    ]
  
TotalBranchHeight[CoalescentTree[tlist_List, elist_List]] :=
    Abs[Subtract @@ Total[Partition[tlist[[Flatten[elist]]], 2]]]

LeafNumber[CoalescentTree[_List, elist_List]] :=
    Length[elist] + 1;
    
TreeTopology[CoalescentTree[_List, elist_List]] :=
    elist[[-1, -1, -1]] //. Thread[elist[[All, 1, 2]] -> elist[[All, {1, 2}, 1]]]
     
ToCombinatoricaGraph[CoalescentTree[tlist_List, elist_List]] :=
    Module[ {rule, nestedlist, leaves, n, vls},
        rule = Thread[elist[[All, 1, 2]] -> elist[[All, {1, 2}, 1]]];
        nestedlist = elist[[-1, -1, -1]] //. rule;
        leaves = Level[nestedlist, {-1}];
        n = Length[leaves];
        vls = Thread[{0, tlist}];
        vls[[leaves, 1]] = Range[n];
        Scan[(vls[[#[[1]], 1]] = N[Mean[vls[[#[[2]], 1]]]]) &, rule];
        Graph[Partition[Flatten[elist, 1], 1], Partition[vls, 1], 
         EdgeDirection -> True]
    ]  
   
treebipartition[CoalescentTree[tlist_List, elist_List]] :=
    Module[ {rule, ls, temp, n},
        n = Length[elist] + 1;
        rule = Thread[elist[[All, 1, 2]] -> elist[[All, {1, 2}, 1]]];
        ls = Join[Partition[Range[n], 1], 
          Sort[Flatten[#]] & /@ (elist[[;; -2, All, 1]] //. rule)];
        temp = SortBy[Flatten[elist, 1], First];
        {tlist[[temp[[All, 2]]]] - tlist[[temp[[All, 1]]]], ls}
    ]
  
TreeCorrelation[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {ls1, ls2, d},
        ls1 = treebipartition[tree1];
        ls2 = treebipartition[tree2];
        d = Total[Flatten[Pick[KroneckerProduct[ls1[[1]], ls2[[1]]],Outer[SameQ, ls1[[2]], ls2[[2]], 1]]]];
        d/Sqrt[Total[ls1[[1]]^2] Total[ls2[[1]]^2]]
    ]   
    
treebipartition2[CoalescentTree[tlist_List, elist_List]] :=
 Module[{n, rule},
  n = Length[elist] + 1;
  rule = Thread[elist[[All, 1, 2]] -> elist[[All, {1, 2}, 1]]];
  Sort[Flatten[#]] & /@ (elist[[;; -3, All, 1]] //. rule)
  ] 
  
UnrootedRFDistance[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {ls1, ls2},
        ls1 = treebipartition2[tree1];
        ls2 = treebipartition2[tree2];
        (*Length[Intersection[ls1, ls2]]/Length[ls1]*)
        2 (Length[ls1] - Length[Intersection[ls1, ls2]])
    ] 
    
TreeClades[CoalescentTree[tlist_List, elist_List]] := 
 Module[{rule, temp},
  rule = Thread[elist[[All, 1, 2]] -> elist[[All, {1, 2}, 1]]];
  temp = Most[elist[[All, 1, 2]]] //. rule;
  Sort[Flatten[#]] & /@ temp
  ]

TreeRFDistance[tree1_CoalescentTree, tree2_CoalescentTree] := 
 2 Length[Complement[TreeClades[tree1], TreeClades[tree2]]]    
    
TreeNodeParent[tree_CoalescentTree, v_Integer?Positive] :=
    If[ v >= 2 LeafNumber[tree] - 1, 
        -1,
        Cases[Flatten[tree[[2]], 1], {v, _}][[1, 2]]
    ]        
   
End[] (* End Private Context *)

EndPackage[]
