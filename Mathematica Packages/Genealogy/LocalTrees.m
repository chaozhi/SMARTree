(* Mathematica Package *)

BeginPackage["Genealogy`LocalTrees`", { "Genealogy`CoalescentTree`"}]
(* Exported symbols added here with SymbolName::usage *) 

mSMCCoalescent::usage = "mSMCCoalescent "

mSMCRecombination::usage = "mSMCRecombination "

mSMCTransition::usage = "mSMCTransition "

NextLocalTree::usage = "NextLocalTree "

RandomNextLocalTree::usage = "RandomNextLocalTree[tree] "

RandomSequentialTrees::usage = "RandomSequentialTrees [tree,r,x0,x1] "

RandomSNPs::usage = "RandomSNPs [trees, \[Theta]]"

Begin["`Private`"] (* Begin Private Context *)

RandomWaitingTime[sq : {{_?NumericQ, _?NumericQ} ..}] :=
    Module[ {diff, wt, res,i},
        diff = Differences[Append[sq[[All, 1]], Infinity]];
        Do[wt = RandomReal[ExponentialDistribution[sq[[i, 2]]]];
           If[ wt <= diff[[i]],
               res = sq[[i]] + {wt, 0};
               Break[]
           ], {i, Length[diff]}];
        res
    ]
    

LineageNumber[CoalescentTree[tlist_List, elist_List], t0_?NonNegative] :=
    Module[ {sq, n},
        n = Length[elist] + 1;
        sq = Transpose[{Drop[tlist, n - 1], Range[n, 1, -1]}];
        If[ t0 < sq[[-1, 1]],
            sq = Select[sq, First[#] > t0 &];
            Prepend[sq, {t0, sq[[1, -1]] + 1}],
            {{t0, 1}}
        ]
    ]    
        
        
mSMCCoalescent[tree_CoalescentTree, t0_?NonNegative] :=
    Module[ {tlist, elist, t1, k, et, i1},
        {tlist, elist} = List @@ tree;
        {t1, k} = RandomWaitingTime[LineageNumber[tree, t0]];
        If[ k == 1,
            {t1, {2 LeafNumber[tree]-1, -1}},
            et = Partition[tlist[[Flatten[elist]]], 2];
            i1 = Flatten[Position[et, _?(#[[1]] <= t1 < #[[2]] &), {1}, Heads -> False]];
            i1 = RandomChoice[i1];
            i1 = {Quotient[i1, 2, 1] + 1, Mod[i1, 2, 1]};
            {t1, Extract[elist,i1]}
        ]
    ]    
    
mSMCRecombination[CoalescentTree[tlist_List, elist_List]] :=
    Module[ {et, w, i0, t0},
        et = Partition[tlist[[Flatten[elist]]], 2];
        w = et[[All, 2]] - et[[All, 1]];
        i0 = RandomChoice[w -> Range[Length[w]]];
        t0 = RandomReal[et[[i0]]];
        i0 = {Quotient[i0, 2, 1] + 1, Mod[i0, 2, 1]};
        {t0, Extract[elist,i0]}
    ]   

mSMCTransition[tree_CoalescentTree] :=
    Module[ {t0,e0,t1,e1},
        {t0, e0} = mSMCRecombination[tree];
        {t1,e1} = mSMCCoalescent[tree, t0];
        {e0,t1,e1}
    ]
    
mSMCTransition[tree_CoalescentTree, 
  erls : {{_Integer, _Integer} ..}] :=
    Module[ {t0, t1, tls, e1,pos},
        If[ Complement[erls, Flatten[tree[[2]], 1]] =!= {},
            Print["Branches do not belong to the tree!"];
            Return[$Failed]
        ];
        tls = Partition[tree[[1, Flatten[erls]]], 2];
        pos = RandomChoice[(tls[[All, 2]] - tls[[All, 1]]) ->Range[Length[tls]]];
        t0 = RandomReal[tls[[pos]]];
        {t1, e1} = mSMCCoalescent[tree, t0];
        {erls[[pos]], t1, e1}
    ]
      
mSMCTransition[tree_CoalescentTree,e0:{_Integer,_Integer}] :=
    Module[ {t0,t1,e1,tlist},
        tlist = tree[[1]];
        t0 = RandomReal[tlist[[e0]]];
        {t1,e1} = mSMCCoalescent[tree, t0];
        {e0,t1,e1}
    ]/;MemberQ[Flatten[tree[[2]],1],e0]
    
mSMCTransition[tree_CoalescentTree, 
  transform : {{_Integer, _Integer}, -1, {_Integer, _Integer}}] :=
    Module[ {e0, t1, e1, tlist, t1max, t1min, n, tsq, nlin, linmax, x},
        {e0, t1, e1} = transform;
        tlist = tree[[1]];
        t1max = If[ e1[[2]] == -1,
                    Infinity,
                    tlist[[e1[[2]]]]
                ];
        t1min = Max[tlist[[{e0[[1]], e1[[1]]}]]];
        n = LeafNumber[tree];
        tsq = Append[Drop[tlist, n - 1], Infinity];
        nlin[t_] :=
            If[ t == Infinity,
                1,
                n - Count[tsq, _?(# <= t &)] + 1
            ];
        linmax = nlin[t1max];
        x = First[RandomWaitingTime[LineageNumber[tree, RandomReal[tlist[[e0]]]]]];
        While[! (RandomReal[] < linmax/nlin[x] && t1min < x < t1max),
         x = First[RandomWaitingTime[LineageNumber[tree, RandomReal[tlist[[e0]]]]]]
         ];
        {e0, x, e1}
    ]    
          
RandomNextLocalTree[tree_CoalescentTree] :=
    Module[ {transform},
        transform = mSMCTransition[tree];
        If[ First[transform] ===Last[transform],
            Return[{None, tree}]
        ];
        {transform,NextLocalTree[tree,transform]}
    ]     

NextLocalTree[tree_CoalescentTree,transform : {{_Integer, _Integer}, _?NonNegative, {_Integer, _Integer}}] :=
    Module[ {n, tlist, elist, e0, t1, e1, p, k, q, w, condition, i0, i1, 
      index0, index1, rg},
        {e0, t1, e1} = transform;
        If[ e0 == e1,
            Return[tree]
        ];
        n = LeafNumber[tree];
        {tlist, elist} = List @@ tree;
        i0 = First[Position[elist, e0]];
        {p, k} = e0;
        q = First[Extract[elist, ReplacePart[i0, 2 -> (-i0[[2]] + 3)]]];
        w = First[e1];
        index0 = First[i0];
        (*the newly emerged sequences coaleses on a point just before the index1 coalescent events of tree*)
        index1 = Count[Drop[tlist, n], _?(# < t1 &)] + 1;
        condition = (Intersection[e0, e1] =!= {});
        If[ e1 == {2 n - 1, -1},
            If[ condition,
                tlist[[-1]] = t1,
                elist = 
                 Append[Delete[elist, index0] /. {k -> q}, {{w, 2 n}, {p, 2 n}}];
                elist = 
                 elist /. 
                  Dispatch[Thread[Range[k + 1, 2 n] -> Range[k, 2 n - 1]]];
                tlist = Append[Delete[tlist, n + index0], t1]
            ],
            (*replace w by 2n*)
            i1 = First[Position[elist, e1]];
            elist[[Sequence @@ i1, 1]] = 2 n;
            If[ index1 > index0,
                tlist = Delete[Insert[tlist, t1, n + index1], n + index0];
                elist = 
                 Delete[Insert[elist, {{w, 2 n}, {p, 2 n}}, index1], index0],
                tlist = Insert[Delete[tlist, n + index0], t1, n + index1];
                elist = 
                 Insert[Delete[elist, index0], {{w, 2 n}, {p, 2 n}}, index1]
            ];
            elist = If[ Last[e1] == k,
                        elist /. {k -> 2 n},
                        elist /. {k -> q}
                    ];
            rg = Range[Min[index0, index1], Max[index0, index1]];
            elist = elist /. Dispatch[Thread[elist[[rg, 1, 2]] -> n + rg]];
        ];
        CoalescentTree[tlist, elist]
    ]     
     
    
(*RandomNextTree[tree_CoalescentTree] :=
    NestWhile[RandomNextLocalTree[Last[#]] &, {0, tree}, 
     SameCoalescentTreeQ[Last[#], tree] &] 
*)     
          
RandomSequentialTrees[tree_CoalescentTree, r_?NonNegative, 
  x0_/;0<=x0<=1, x1_/;0<=x1<=1] :=
    Module[ {f},
        f = Function[{xx},
          Module[ {x, tran, tr},
              {x, tran, tr} = xx;
              x += 
               RandomReal[ExponentialDistribution[r TotalBranchHeight[tr]/2]];
              {tran, tr} = RandomNextLocalTree[tr];
              {x, tran, tr}
          ]
          ];
        NestWhileList[f, {x0, "Initial", tree}, First[#] <= x1 &]
    ]/;x1>x0     
     
RandomSequentialTrees[tree_CoalescentTree, r_?NonNegative, 
  x0_Integer, x1_Integer] :=
    Module[ {f},
        f = Function[{xx},
              Module[ {x, tran, tr,p},
                  {x, tran, tr} = xx;
                  p = 1 - Exp[-r TotalBranchHeight[tr]/2];
                  x += RandomInteger[GeometricDistribution[p]] + 1;
                  {tran, tr} = RandomNextLocalTree[tr];
                  {x, tran, tr}
              ]
          ];
        NestWhileList[f, {x0, "Initial", tree}, First[#] <= x1 &]
    ] /; x1 > x0         

RandomSNPs[trees : {{_?NumericQ, _CoalescentTree} ..}, theta_?NonNegative] :=
    Module[ {i,res, x, tr, xnext, lam, heights,elist, subres, deltx, index, set, n,y,ends},
        n=LeafNumber[trees[[1,2]]];
        res = Table[
          {x, tr} = trees[[i]];
          xnext = If[ i == Length[trees],
                      1,
                      First[trees[[i + 1]]]
                  ];
          lam = TotalBranchHeight[tr] theta/2;
          heights = BranchHeights[tr];
          elist = tr[[2]];
          subres = {};
          While[True,
           deltx = RandomReal[ExponentialDistribution[lam]];
           If[ x + deltx > xnext,
               Break[]
           ];
           index = RandomChoice[heights -> Range[Length[heights]]];
           set = Flatten[Most[Flatten[elist, 1][[index]]] //. 
              Thread[elist[[All, 1, 2]] -> elist[[All, All, 1]]]];
           y = Table[0, {n}];
           y[[set]] = 1;
           x += deltx;
           AppendTo[subres, {x, y}]
           ];
          subres, {i, Length[trees]}];
        res = Join @@ res;
        ends = Table[
            {x, tr} = trees[[i]];
            lam = TotalBranchHeight[tr] theta/2;
            heights = BranchHeights[tr];
            elist = tr[[2]];
            index = RandomChoice[heights -> Range[Length[heights]]];
            set = Flatten[Most[Flatten[elist, 1][[index]]] //.Thread[elist[[All, 1, 2]] -> elist[[All, All, 1]]]];
            y = Table[0, {n}];
            y[[set]] = 1;
            {x, y}, {i, {1, -1}}];
        res = Prepend[res,First[ends]];
        If[ res[[-1,1]]!=1,
            res = Append[res,Last[ends]]
        ];
        res
    ]
    
    
End[] (* End Private Context *)

EndPackage[]