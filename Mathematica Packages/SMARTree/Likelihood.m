(* Mathematica Package *)

BeginPackage["SMARTree`Likelihood`", { "Genealogy`CoalescentTree`"}]
(* Exported symbols added here with SymbolName::usage *)  

TreeTransitionProb::usage = "TreeTransitionProb [tree]"

TreeLogLikelihood::usage = "TreeLogLikelihood[yls,treels,allelefreq,\[Theta],\[Epsilon]]  "

TreeLogPriorProb::usage = "TreeLogPriorProb [tree]"

PolymorphicProb::usage="PolymorphicProb"

TreeTransitionTypeI::usage = "TreeTransitionTypeI  "

TreeTransitionTypeII::usage = "TreeTransitionTypeII  "

(*TreeTransitionTypeII2::usage = "TreeTransitionTypeII2  "*)

Begin["`Private`"] (* Begin Private Context *)
   
(*
TreeLogLikelihood::trees = 
  "The trees are not a list of coalescent trees!";
(*Assuming alleles at SNP sites starting from 0*)  
TreeLogLikelihood[data_?VectorQ, tree_CoalescentTree, 
  allelefreq_?VectorQ, theta_?NumericQ, epsilon_?NumericQ] :=
    First[TreeLogLikelihood[{data},{tree},{allelefreq},theta,epsilon]];
*)
TreeLogLikelihood[data_?VectorQ, tree_CoalescentTree, 
  allelefreq_?VectorQ, theta_?NumericQ, epsilon_?NumericQ,modelid_?StringQ] :=
    First[TreeLogLikelihood[{data}, {tree}, {allelefreq}, theta,epsilon,modelid]]; 
TreeLogLikelihood[snps_?MatrixQ, trees : {_CoalescentTree ..}, afreq_?MatrixQ, theta_?NumericQ, epsilon_?NumericQ,modelid_?StringQ] :=
    Module[ {nsite, nleaf, na, like, p, q, tt, pt, qt, kt, ep, eq, y,i,temp,forest = trees,res},
        {nsite, nleaf} = Dimensions[snps];
        na = Last[Dimensions[afreq]];
        like = ConstantArray[0, {nsite, 2 nleaf - 1}];
        (*Assume each snp site has two possible alleles: 0 and 1;
        like[[All, ;;nleaf]] to be modified for >2 alleles*)
        like[[All, ;;nleaf]] = Replace[snps, {1 -> {epsilon, 1 - epsilon},0 -> {1 - epsilon, epsilon}}, {2}];
        p = forest[[All, 2, All, 1, 1]];
        q = forest[[All, 2, All, 2, 1]];
        tt = forest[[All, 1]];
        pt = MapThread[#1[[#2]] &, {tt, p}];
        qt = MapThread[#1[[#2]] &, {tt, q}];
        kt = tt[[All, nleaf + 1 ;;]];
        ep = Exp[-(kt - pt) theta/2];
        eq = Exp[-(kt - qt) theta/2];
        (*Assume each snp site has two possible alleles: 0 and 1;
          need modification for >2 alleles, or n<4*)
        Do[
         y = Extract[like, Transpose[{Range[nsite], p[[All, i]]}]];
         temp = y ep[[All, i]] + (1 - ep[[All, i]]) Transpose[{y[[All, 2]], y[[All, 1]]}];
         y = Extract[like, Transpose[{Range[nsite], q[[All, i]]}]];
         temp *= y eq[[All, i]] + (1 - eq[[All, i]]) Transpose[{y[[All, 2]],y[[All, 1]]}];
         like[[All, nleaf + i]] = temp, {i, nleaf - 1}];
        res = Log[MapThread[Dot, {like[[All, -1]], afreq}]];
        Switch[modelid,
        	"M1"|"M1S", res,
        	"M2", res-Log[PolymorphicProb[forest,theta,epsilon,2]],
        	"M2S", res-Log[PolymorphicProb[forest,theta,epsilon,1]],
        	_, Print["TreeLogLikelihood: Wrong ascertainmodel =", modelid]; Abort[];
        ]        
    ]   
              
(*calculate p(NMA<m_a |T), for m_a=2 by default. NMA is the number of minor alleles at a site*)              
PolymorphicProb[trees : {_CoalescentTree ..}, theta_?NumericQ, epsilon_?NumericQ, nma:1|2] :=
    Module[ {forest=trees,nleaf,pos,set1,set11,set2,lam,lam11,lam2,temp,temp2,i},
    	nleaf=LeafNumber[forest[[1]]];
        (*lam*)
        lam = TotalBranchHeight[forest] theta/2;
        (*lam11*)
        set1 = Map[SortBy[Flatten[#, 1], First] &, forest[[All, 2]]][[All, ;; nleaf]];
        pos = Position[set1, {_, 2 nleaf - 1}, {2}, Heads -> False];
        set11 = If[ pos === {},
                    set1,
                    ReplacePart[set1, MapThread[#1 -> Sequence @@ #2 &, {pos,forest[[pos[[All, 1]], 2, -1]]}]]
                ];
        temp = MapThread[#1[[#2]] &, {forest[[All, 1]], set11[[All, All, 2]]}];
        temp -= MapThread[#1[[#2]] &, {forest[[All, 1]], set11[[All, All, 1]]}];
        lam11 = Total[temp, {2}] theta/2;        
        (*lam2*)
        set2 = Table[
           temp = Select[SplitBy[SortBy[set1[[i]], Last], Last],Length[#] == 2 &][[All, 1, 2]];
           Select[Flatten[forest[[i, 2]], 1], MemberQ[temp, #[[1]]] &], {i,Length[set1]}];
        temp = MapThread[#1[[#2]] &, {forest[[All, 1]], set2[[All, All, 2]]}];
        temp -= MapThread[#1[[#2]] &, {forest[[All, 1]], set2[[All, All, 1]]}];
        lam2 = Total[temp, {2}] theta/2;        
        (**)
        temp = {Exp[-lam], Exp[-lam] lam11,Exp[-lam] lam2};
        temp2 = {(1 - epsilon)^nleaf, nleaf (1 - epsilon)^(nleaf - 1) epsilon, 
                (1 - epsilon)^(nleaf - 2) epsilon^2,(1 - epsilon)^(nleaf - 3) epsilon^3};
        Switch[nma,
        	1,
        	temp = temp[[1]] temp2[[1]]+temp[[2]] temp2[[2]]/nleaf+temp[[3]] temp2[[3]];
            1-temp,
        	2,
        	temp = temp[[1]] Total[temp2[[;;2]]] + 
              temp[[2]] temp2[[1]]+temp[[2]] temp2[[2]]/nleaf+temp[[2]] temp2[[3]] (nleaf-1)+
              temp[[3]] temp2[[2]] 2/nleaf+temp[[3]] temp2[[3]]+temp[[3]] temp2[[4]] (nleaf-2);
            1-temp,
            _,Print["PolymorphicProb: Wrong minimum number of minor alleles=", nma];Abort[];
        ]        
    ]
              
TreeLogLikelihood1[snps_?MatrixQ, trees : {_CoalescentTree ..}, afreq_?MatrixQ, theta_?NumericQ, epsilon_?NumericQ] :=
    Module[ {nsite, nleaf, na, like, p, q, tt, pt, qt, kt, ep, eq, y,i,temp,pos,temp2,forest = trees,
        set1,set11,lam,lam11,set2,lam2},
        {nsite, nleaf} = Dimensions[snps];
        na = Last[Dimensions[afreq]];
        like = ConstantArray[0, {nsite, 2 nleaf - 1}];
        (*Assume each snp site has two possible alleles: 0 and 1;
        like[[All, ;;nleaf]] to be modified for >2 alleles*)
        like[[All, ;;nleaf]] = Replace[snps, {1 -> {epsilon, 1 - epsilon},0 -> {1 - epsilon, epsilon}}, {2}];
        p = forest[[All, 2, All, 1, 1]];
        q = forest[[All, 2, All, 2, 1]];
        tt = forest[[All, 1]];
        pt = MapThread[#1[[#2]] &, {tt, p}];
        qt = MapThread[#1[[#2]] &, {tt, q}];
        kt = tt[[All, nleaf + 1 ;;]];
        ep = Exp[-(kt - pt) theta/2];
        eq = Exp[-(kt - qt) theta/2];
        (*Assume each snp site has two possible alleles: 0 and 1;
          need modification for >2 alleles, or n<4*)
        Do[
         y = Extract[like, Transpose[{Range[nsite], p[[All, i]]}]];
         temp = y ep[[All, i]] + (1 - ep[[All, i]]) Transpose[{y[[All, 2]], y[[All, 1]]}];
         y = Extract[like, Transpose[{Range[nsite], q[[All, i]]}]];
         temp *= y eq[[All, i]] + (1 - eq[[All, i]]) Transpose[{y[[All, 2]],y[[All, 1]]}];
         like[[All, nleaf + i]] = temp, {i, nleaf - 1}];
        (*lam0 = TotalBranchHeight[forest] theta/2;*)
        (*lam*)
        lam = TotalBranchHeight[forest] theta/2;
        (*lam11*)
        set1 = Map[SortBy[Flatten[#, 1], First] &, forest[[All, 2]]][[All, ;; nleaf]];
        pos = Position[set1, {_, 2 nleaf - 1}, {2}, Heads -> False];
        set11 = If[ pos === {},
                    set1,
                    ReplacePart[set1, MapThread[#1 -> Sequence @@ #2 &, {pos,forest[[pos[[All, 1]], 2, -1]]}]]
                ];
        temp = MapThread[#1[[#2]] &, {forest[[All, 1]], set11[[All, All, 2]]}];
        temp -= MapThread[#1[[#2]] &, {forest[[All, 1]], set11[[All, All, 1]]}];
        lam11 = Total[temp, {2}] theta/2;        
        (*lam2*)
        set2 = Table[
           temp = Select[SplitBy[SortBy[set1[[i]], Last], Last],Length[#] == 2 &][[All, 1, 2]];
           Select[Flatten[forest[[i, 2]], 1], MemberQ[temp, #[[1]]] &], {i,Length[set1]}];
        temp = MapThread[#1[[#2]] &, {forest[[All, 1]], set2[[All, All, 2]]}];
        temp -= MapThread[#1[[#2]] &, {forest[[All, 1]], set2[[All, All, 1]]}];
        lam2 = Total[temp, {2}] theta/2;        
        (**)
        temp = {Exp[-lam], Exp[-lam] lam11, Exp[-lam] lam2};
        temp2 = {(1 - epsilon)^nleaf, nleaf (1 - epsilon)^(nleaf - 1) epsilon, 
                (1 - epsilon)^(nleaf - 2) epsilon^2,(1 - epsilon)^(nleaf - 3) epsilon^3};
        temp = temp[[1]] Total[temp2[[;;2]]] + 
              temp[[2]] temp2[[1]]+temp[[2]] temp2[[2]]/nleaf+temp[[2]] temp2[[3]] (nleaf-1)+
              temp[[3]] temp2[[2]] 2/nleaf+temp[[3]] temp2[[3]]+temp[[3]] temp2[[4]] (nleaf-2); 
        (*temp = temp2[[1]] (temp[[1]]+temp[[2]])+
               temp2[[2]] (temp[[1]]+temp[[2]]/nleaf+temp[[3]] 2/nleaf);*)
        (*temp=(1-temp) Exp[-(lam0-lam)];*)
        Log[MapThread[Dot, {like[[All, -1]], afreq}]/(1-temp)]
    ]      
   


TreeLogLikelihood2[snps_?MatrixQ, trees : {_CoalescentTree ..}, afreq_?MatrixQ, theta_?NumericQ, epsilon_?NumericQ] :=
    Module[ {nsite, nleaf, na, like, p, q, tt, pt, qt, kt, ep, eq, y,i,temp,pos,temp2,forest = trees,
        set1,set11,lam,lam11},
        {nsite, nleaf} = Dimensions[snps];
        na = Last[Dimensions[afreq]];
        like = ConstantArray[0, {nsite, 2 nleaf - 1}];
        (*Assume each snp site has two possible alleles: 0 and 1;
        like[[All, ;;nleaf]] to be modified for >2 alleles*)
        like[[All, ;;nleaf]] = Replace[snps, {1 -> {epsilon, 1 - epsilon},0 -> {1 - epsilon, epsilon}}, {2}];
        p = forest[[All, 2, All, 1, 1]];
        q = forest[[All, 2, All, 2, 1]];
        tt = forest[[All, 1]];
        pt = MapThread[#1[[#2]] &, {tt, p}];
        qt = MapThread[#1[[#2]] &, {tt, q}];
        kt = tt[[All, nleaf + 1 ;;]];
        ep = Exp[-(kt - pt) theta/2];
        eq = Exp[-(kt - qt) theta/2];
        (*Assume each snp site has two possible alleles: 0 and 1;
          need modification for >2 alleles, or n<4*)
        Do[
         y = Extract[like, Transpose[{Range[nsite], p[[All, i]]}]];
         temp = y ep[[All, i]] + (1 - ep[[All, i]]) Transpose[{y[[All, 2]], y[[All, 1]]}];
         y = Extract[like, Transpose[{Range[nsite], q[[All, i]]}]];
         temp *= y eq[[All, i]] + (1 - eq[[All, i]]) Transpose[{y[[All, 2]],y[[All, 1]]}];
         like[[All, nleaf + i]] = temp, {i, nleaf - 1}];
        (*lam*)
        lam = TotalBranchHeight[forest] theta/2;
        (*lam11*)
        set1 = Map[SortBy[Flatten[#, 1], First] &, forest[[All, 2]]][[All, ;; nleaf]];
        pos = Position[set1, {_, 2 nleaf - 1}, {2}, Heads -> False];
        set11 = If[ pos === {},
                    set1,
                    ReplacePart[set1, MapThread[#1 -> Sequence @@ #2 &, {pos,forest[[pos[[All, 1]], 2, -1]]}]]
                ];
        temp = MapThread[#1[[#2]] &, {forest[[All, 1]], set11[[All, All, 2]]}];
        temp -= MapThread[#1[[#2]] &, {forest[[All, 1]], set11[[All, All, 1]]}];
        lam11 = Total[temp, {2}] theta/2;        
        (**)
        temp = {Exp[-lam], Exp[-lam] lam11};
        temp2 = {(1 - epsilon)^nleaf, nleaf (1 - epsilon)^(nleaf - 1) epsilon, (1 - epsilon)^(nleaf - 2) epsilon^2};
        (*temp = temp[[1]] Total[temp2[[;;2]]] + temp[[2]] temp2[[1]]+temp[[2]] temp2[[2]]/nleaf+temp[[2]] temp2[[3]] (nleaf-1);*)
        temp = temp[[1]] Total[temp2[[;;2]]] + temp[[2]] temp2[[1]];
        Log[MapThread[Dot, {like[[All, -1]], afreq}]/(1-temp)]
    ]    

(*ParentNode[tree_CoalescentTree, v_Integer?Positive] :=
    If[ v >= 2 LeafNumber[tree] - 1,
        -1,
        Cases[Flatten[tree[[2]], 1], {v, _}][[1, 2]]
    ]*)

(*Clear[compiletransitionN1]
Clear[n, vv, ii, tlst, klst, wlst, kwlst, res, p1, p2, i1, i2, j1, \
j2, kw, pc, pr, p]*)
(*
compiletransitionN1 = Compile[{{tlist, _Real, 1}, {ee, _Integer, 2}},
   n = (Length[tlist] + 1)/2;
   vv = ee /. {-1 -> 2 n, (x_ /; x < n) -> n};
   vv = Reverse[#] & /@ (2 n + 1 - vv);
   If[vv[[1, 2]] < vv[[2, 1]], 0,
    vv[[1, 1]] = Max[vv[[1, 1]], vv[[2, 1]]];
    vv[[2, 2]] = Min[vv[[1, 2]], vv[[2, 2]]];
    ii = vv - Min[vv] + 1;
    ii[[All, 2]] = ii[[All, 2]] - 1;  
    tlst = Append[tlist, Infinity];
    wlst = Reverse[Differences[tlst[[-Max[vv] ;; -Min[vv]]]]];
    klst = Range[Min[vv], Max[vv] - 1];
    kwlst = klst wlst;
    Which[
     ii[[1, 1]] > ii[[2, 2]],
     kw = kwlst[[ii[[2, 1]] ;; ii[[2, 2]]]];
     pc = -Reverse[Most[Prepend[Accumulate[Reverse[kw]], 0]]] +Log[(1 - Exp[-kw])/klst[[ii[[2, 1]] ;; ii[[2, 2]]]]];
     kw = kwlst[[ii[[1, 1]] ;; ii[[1, 2]]]];
     pr = -Most[Prepend[Accumulate[kw], 0]] +Log[(1 - Exp[-kw])/klst[[ii[[1, 1]] ;; ii[[1, 2]]]]];
     p = -Total[kwlst[[ii[[2, 2]] + 1 ;; ii[[1, 1]] - 1]]];
     Total[Exp[Flatten[Outer[Plus, pr, pc] + p]]],
     ii[[1, 1]] <= ii[[2, 2]],
     p1 = Log[(1 - Exp[-kwlst[[ii[[1, 1]] ;; ii[[1, 2]]]]])/klst[[ii[[1, 1]] ;; ii[[1, 2]]]]];
     p2 = Log[(1 - Exp[-kwlst[[ii[[2, 1]] ;; ii[[2, 2]]]]])/klst[[ii[[2, 1]] ;; ii[[2, 2]]]]];
     res = Outer[Plus, p1, p2];
     res -= Transpose[Rest[FoldList[#1 - #2 &, Total[kwlst[[;; ii[[1, 1]] - 1]]] + 
          Prepend[Accumulate[kwlst[[ii[[1, 1]] ;; ii[[1, 2]] - 1]]], 
           0], kwlst[[;; ii[[2, 2]]]]]]];
     i1 = 1;
     i2 = Min[ii[[1, 2]], ii[[2, 2]]] - ii[[1, 1]] + 1;
     j1 = ii[[1, 1]];
     j2 = Min[ii[[1, 2]], ii[[2, 2]]];
     res[[i1 ;; i2, j1 ;; j2]] = LowerTriangularize[res[[i1 ;; i2, j1 ;; j2]], -1] + 
       UpperTriangularize[ConstantArray[-Infinity, {i2 - i1 + 1, j2 - j1 + 1}], 1] + 
       DiagonalMatrix[Log[(Exp[-kwlst[[j1 ;; j2]]] + kwlst[[j1 ;; j2]] - 1)/klst[[j1 ;; j2]]^2]];
     Total[Exp[Flatten[res]]]
     ]
    ], {{n, _Integer}, {vv, _Integer, 2}, {ii, _Integer, 
     2}, {tlst, _Real, 1}, {wlst, _Real, 1}, {klst, _Integer, 
     1}, {kwlst, _Real, 1}, {p1, _Real, 1}, {p2, _Real, 1}, {res, _Real, 
     2}, {i1, _Integer}, {i2, _Integer}, {j1, _Integer}, {j2, _Integer}},
     CompilationTarget -> "C"
   ];

TreeTransitionProb4[tree_CoalescentTree, 
  op : {{_Integer, _Integer}, -1, {_Integer, _Integer}}] := 
 compiletransitionN1[tree[[1]], op[[{1, 3}]]]/TotalBranchHeight[tree]
*)


TreeTransitionProb::dis = "The distance between trees `1` and `2` must be no greater than one transition!"
TreeTransitionProb[tree1_CoalescentTree,tree2_CoalescentTree] :=
    Module[ {type,trans},
        {type,trans} = TreeTransitionTypeI[tree1,tree2];
        Switch[type,
            0,
            TreeTransitionProb[tree1],
            1,
            Total[(TreeTransitionProb[tree1,#]&/@trans[[All,{1,2}]])],            
            _,
            Message[TreeTransitionProb::dis,tree1,tree2];
            Print["{tree1,tree2}=",{tree1,tree2}];
            Return[$Failed]
        ]
    ]


TreeTransitionProb[tree_CoalescentTree, op : {{_Integer, _Integer}, -1, {_Integer, _Integer}}] :=
    Module[ {n, vv, ii, tlist, wlist, klist, kwlist, kw, pc, pr, p, res, 
      i, j, p1, p2},
          (*vv is count from the first node at infinity far way, and the MCRA is the second node, and so on. 
         All the leaf nodes are label as n th node;
         Only list nodes covered by the edges of op;
         klist is the number of lineages starting from the oldest node of op for each interval;
         wlist is the length of time interval for each interval;
         ii is the indices of interval for edges in op;
        *)
        n = LeafNumber[tree];
        vv = op[[{1, 3}]] /. {-1 -> 2 n, (x_ /; x < n) -> n};
        vv = Reverse[#] & /@ (2 n + 1 - vv);
        If[ vv[[1, 2]] < vv[[2, 1]],
            Return[0]
        ];
        vv[[1, 1]] = Max[vv[[1, 1]], vv[[2, 1]]];
        vv[[2, 2]] = Min[vv[[1, 2]], vv[[2, 2]]];
        ii = vv - Min[vv] + 1;
        ii[[All, 2]] -= 1;
        tlist = Append[tree[[1]], Infinity];
        wlist = Reverse[Differences[tlist[[-Max[vv] ;; -Min[vv]]]]];
        klist = Range[Min[vv], Max[vv] - 1];
        kwlist = klist wlist;
        res = Which[
         ii[[1, 1]] > ii[[2, 2]],
         kw = kwlist[[Span @@ ii[[2]]]];
         pc = -Reverse[Most[Prepend[Accumulate[Reverse[kw]], 0]]] +Log[(1 - Exp[-kw])/klist[[Span @@ ii[[2]]]]];
         kw = kwlist[[Span @@ ii[[1]]]];
         pr = -Most[Prepend[Accumulate[kw], 0]] + Log[(1 - Exp[-kw])/klist[[Span @@ ii[[1]]]]];
         p = -Total[kwlist[[ii[[2, 2]] + 1 ;; ii[[1, 1]] - 1]]];
         Total[Exp[Flatten[Outer[Plus, pr, pc] + p]]]/TotalBranchHeight[tree],
         ii[[1, 1]] <= ii[[2, 2]],
         {p1, p2} = Log[(1 - Exp[-kwlist[[Span @@ #]]])/klist[[Span @@ #]] & /@ ii];
         res = Outer[Plus, p1, p2];
         res -= Transpose[Rest[FoldList[#1 - #2 &, Total[kwlist[[;; ii[[1, 1]] - 1]]] + 
              Prepend[Accumulate[kwlist[[ii[[1, 1]] ;; ii[[1, 2]] - 1]]],0], kwlist[[;; ii[[2, 2]]]]]]];
         i = 1 ;; Min[ii[[1, 2]], ii[[2, 2]]] - ii[[1, 1]] + 1;
         j = ii[[1, 1]] ;; Min[ii[[1, 2]], ii[[2, 2]]];
         res[[i, j]] = LowerTriangularize[res[[i, j]], -1] +
             UpperTriangularize[ConstantArray[-Infinity, Dimensions[res[[i, j]]]], 1] + 
               DiagonalMatrix[Log[(Exp[-kwlist[[j]]] + kwlist[[j]] - 1)/klist[[j]]^2]];
         Total[Exp[Flatten[res]]]/TotalBranchHeight[tree],
         True,
         0
         ];
        If[ !NumericQ[res],
            Print["Wrong in TreeTransitionProb! {tree,op}=",{tree,op}];
        ];
        res
    ]
(*
calcaulation the prability that the tree remains after one recombination on it.
each branch of the tree is divived into segments according coalescent times.
p0: transitions between same segmetns
*)
TreeTransitionProb[tree_CoalescentTree] :=
    Module[ {tlist, elist, n, ksq, wsq, kwsq, xysq, p0, ed, cr, a},
        {tlist, elist} = List @@ tree;
        n = LeafNumber[tree];
        ksq = Range[n, 2, -1];
        wsq = Differences[Drop[tlist, n - 1]];
        kwsq = ksq wsq;
        xysq = (Exp[kwsq] - 1)/ksq;
        p0 = Total[(Exp[-kwsq]+kwsq-1)/ksq];
        ed = Flatten[elist, 1];
        ed[[All, 1]] = Replace[ed[[All, 1]], (x_ /; x < n) -> n, {1}];
        ed = 2 n - Pick[ed, ed[[All, 2]] - ed[[All, 1]] - 2, _?NonNegative];
        ed[[All, 2]] += 1;
        ed = n + 1 - ed;
        ed = Tally[ed];
        cr = SparseArray[Thread[ed[[All, 1]] -> ed[[All, 2]]], {n - 1, n - 1}];
        cr = UpperTriangularize[Accumulate[cr], 1];
        cr = Reverse[Transpose[cr]];
        cr = Transpose[Reverse[Accumulate[cr]]];
        cr = UpperTriangularize[cr, 1];
        cr *= UpperTriangularize[KroneckerProduct[xysq, xysq], 1];
        a = UpperTriangularize[Table[kwsq, {Length[cr]}]];
        a = Transpose[Accumulate[Transpose[a]]];
        cr *= Exp[-a];
        (p0 + Total[cr, 2])/TotalBranchHeight[tree]
    ]       
    
TreeTransitionProb[tree_CoalescentTree, transform: {{_Integer, _Integer}, _?NonNegative}] :=
    Module[ {n, tlist, v1, v2, t2, tc, k1, k2, wsq, ksq, res},
        n = LeafNumber[tree];
        tlist = tree[[1]];
        {{v1, v2}, tc} = transform;
        t2 = tlist[[v2]];
        (*k1 is the number of lineages just after node v1*)
        k1 = If[ v1 <= n,
                 n,
                 2 n - v1
             ];
        (*tlist keeps the coalescent times for the remaining k1 lineages,including t1 at v1*)
        tlist = Take[tlist, -k1];
        If[ tc <= t2,
            wsq = Reverse[Differences[Append[Select[tlist, (# < tc &)], tc]]];
            ksq = Range[k1 - Length[wsq] + 1, k1];
            res = Total[Exp[-Most[FoldList[Plus, 0, wsq ksq]]] (1-Exp[-ksq wsq])/ksq];
            res/TotalBranchHeight[tree],
            (*k2 is the number of lineages just after node v2*)
            k2 = 2 n - v2;
            (*wsq is the time interval from t2 at v2 to t1 at v1*)
            wsq = Reverse[Differences[Drop[tlist, -(k2 - 1)]]];
            (*ksq is the number of lineages at each interval corresponding to wsq*)
            ksq = Range[k1 - Length[wsq] + 1, k1];
            res = Total[Exp[-ReplacePart[RotateRight[Accumulate[wsq ksq]], 1 -> 0]] (1-Exp[-ksq wsq])/ksq];
            wsq = Reverse[Differences[Append[Select[Drop[tlist, k1 - k2], (# < tc &)], tc]]];
            ksq = Range[k2 - Length[wsq] + 1, k2];
            Exp[-Total[wsq ksq]] res/TotalBranchHeight[tree]
        ]
    ]    

              
TreeLogPriorProb[tree_CoalescentTree] :=
    First[TreeLogPriorProb[{tree}]];
TreeLogPriorProb::dim = "Trees must have same number of nodes!";
TreeLogPriorProb[trees : {_CoalescentTree ..}] :=
    Module[ {tlists, n},
        tlists = trees[[All, 1]];
        If[ MatrixQ[tlists],
            tlists = Transpose[tlists],
            Message[TreeLogPriorProb::dim];
            Return[$Failed]
        ];
        n = LeafNumber[First[trees]];
        -Total[Differences[Drop[tlists, n - 1]] Range[n, 
            2, -1] Range[n - 1, 1, -1]/2]
    ]
    
TreeTransitionTypeI[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {tlist1, tlist2, elist1, elist2, al1, al2, n, t10, v10, v11, 
      v12, v13, t30, v30, v31, v32, v33, nv12,nv13,share,res},
        If[ SameCoalescentTreeQ[tree1, tree2],
            Return[{0, 0}]
        ];
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        t10 = Complement[tlist1, tlist2];
        If[ Length[t10] == 1,
            t10 = First[t10],
            Return[{-1, -1}]
        ];
        n = LeafNumber[tree1];
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        v10 = Position[tlist1, t10][[1, 1]];
        {v11, v12} = elist1[[v10 - n, All, 1]];
        v13 = al1[[v10]];
        t30 = First[Complement[tlist2, tlist1]];
        v30 = Position[tlist2, t30][[1, 1]];
        {v31, v32} = elist2[[v30 - n, All, 1]];
        v33 = al2[[v30]];
        {v31, v32, v33} = {v31, v32,v33} /. {(x_ /; x > n) :> Position[tlist1, tlist2[[x]]][[1, 1]]};
        share = Intersection[{v11, v12}, {v31, v32}];
        res = Switch[Length[share],
         2,
         If[ v33 == v13,
             {1, If[ t30 < t10,
                     {{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}},
                     {{{v11, v10}, t30, {v10, v13}}, {{v12, v10}, t30, {v10, v13}}}
                 ]},
             {-1, -1}
         ],
         1,
         If[ v31 =!= share[[1]],
             {v31, v32} = {v32, v31}
         ];
         If[ v11 =!= share[[1]],
             {v11, v12} = {v12, v11}
         ];
         {nv12, nv13} = {v12,v13} /. {(x_ /; x > n) :> Position[tlist2, tlist1[[x]]][[1, 1]]};
         If[ al1[[v32]] == v33&&al2[[nv12]]==nv13,
             {1, {{{v31, v10}, t30, {v32, v33}}}},
             {-1, -1}
         ],
         0, {-1, -1}
        ];
        res
    ]  
    
contrasttree[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {tlist1, tlist2, elist1, elist2, n, al1, al2, tnew, newnode, 
      newlabel, newchild, newparent, told, oldnode, oldchild, oldparent, 
      old, new, newlabelrule},
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        tnew = Sort[Complement[tlist2, tlist1]];
        told = Sort[Complement[tlist1, tlist2]];
        If[ Length[tnew] > 2 || Length[told] > 2,
            Return[{-1, {{-1, -1}, {-1, -1}}}]
        ];
        If[ Length[tnew] != Length[told],
            Print["Case in compareTree to be done!", "{tree1,tree2}=", {tree1, tree2}];
            Return[$Failed]
        ];
        n = LeafNumber[tree1];
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        newnode = Flatten[Position[tlist2, #] & /@ tnew];
        newlabelrule = Thread[newnode -> Take[{"A", "B"}, Length[newnode]]];
        newlabel = newnode /. newlabelrule;
        newchild = elist2[[newnode - n, All, 1]] /. newlabelrule;
        newparent = al2[[newnode]] /. newlabelrule;
        {newchild, newparent} = {newchild, newparent} /. {(x_ /; x > n) :> 
            Position[tlist1, tlist2[[x]]][[1, 1]]};
        oldnode = Flatten[Position[tlist1, #] & /@ told];
        oldchild = elist1[[oldnode - n, All, 1]];
        oldparent = al1[[oldnode]];
        If[ Length[oldchild] == 2 && MemberQ[oldchild[[2]], oldnode[[1]]] &&oldchild[[2, 1]] == oldnode[[1]],
            oldchild[[2]] = Reverse[oldchild[[2]]]
        ];
        If[ Length[newchild] == 2 && MemberQ[newchild[[2]], newlabel[[1]]] &&newchild[[2, 1]] == newlabel[[1]],
            newchild[[2]] = Reverse[newchild[[2]]]
        ];
        old = Transpose[{told, oldchild, oldparent, oldnode, oldnode}];
        new = Transpose[{tnew, newchild, newparent, newlabel, newnode}];
        {{old, new}, {al1, al2}}
    ]
    
midtree11[old_, new_, al1_, al2_,tree1_, tree2_] :=
    Module[ {v11, v12, v10, v13, t10, v31, v32, v33, t30, share,tlist1, tlist2, nv11, nv12, nv13,nv33,nvy, vy,n,group,res},
        {t10, {v11, v12}, v13, v10} = Take[Flatten[old, 1], 4];
        {t30, {v31, v32}, v33} = Take[Flatten[new, 1], 3];
        share = Intersection[{v11,v12},{v31,v32}];
        tlist1 = First[tree1];
        tlist2 = First[tree2];
        n = LeafNumber[tree1];
        res = Switch[Length[share],          
            2,
            group = "11";
            If[ v33 == v13,
                {1, If[ t30 < t10,
                        {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}}},
                        {{{{v11, v10}, t30, {v10, v13}}, {{v12, v10},t30, {v10, v13}}}}
                    ]},
                {-1, -1}
            ],
            1,
            If[ v31!=share[[1]],
                {v31,v32} = {v32,v31}
            ];
            If[ v11!=share[[1]],
                {v11,v12} = {v12,v11}
            ];
            vy = al1[[v32]];            
            (*nv: label in tree 2; v: label in tree1*)
            {nv12, nv13,nvy} = {v12,v13,vy} /. {(x_ /; x > n) :>Position[tlist2, tlist1[[x]]][[1, 1]]};
            Which[
                vy==v33&&al2[[nv12]]==nv13,
                group = "11";
                {1,{{{{v11,v10},t30,{v32,v33}}}}},
                vy!=v13&&v13==v33&&al2[[nv12]]==nvy,
                group = "21";
                {2,{{{{v12,v10},-1,{v32,vy}}}}},                
                True,
                group = "-1";
                {-1,-1}
            ],
            0,
            group = "21";
            If[ !MemberQ[al1[[{v31,v32}]],v33],
                Return[{-1,"-1",-1}]
            ];
            If[ al1[[v32]]!=v33,
                {v31,v32} = {v32,v31};
            ];
            vy = al1[[v31]];
            (*nv: label in tree 2; v: label in tree1*)
            {nv11, nv12, nv13,nv33,nvy} = {v11, v12, v13,v33,vy} /. {(x_ /; x > n) :>Position[tlist2, tlist1[[x]]][[1, 1]]};
            If[ !MemberQ[al2[[{nv11,nv12}]],nv13],
                Return[{-1,"-1",-1}]
            ];
            If[ al2[[nv11]]!=nv13,
                {v11,v12} = {v12,v11};
                {nv11,nv12} = {nv12,nv11}
            ];
            Which[
                vy==v33&&al2[[nv12]]==nv33,                   
                {2, {{{{v12, v10}, -1, {v31, v33}}},{{{v12, v10}, -1, {v32, v33}}}}},
                vy==v13&&al2[[nv12]]==nv13,
                {2, {{{{v12, v10}, -1, {v31, v13}}},{{{v11, v10}, -1, {v31, v13}}}}},
                (!MemberQ[{v13,v33},vy])&&al2[[nv12]]==nvy,
                {2, {{{{v12, v10}, -1, {v31, vy}}}}},
                True,
                {-1,-1}                    
            ]                       
        ];
        If[ res===-1,
            group = "-1"
        ];
        Prepend[res,group]
    ]        

(*midtree112 have more operators than midtree11 when there are underdetermined coalescent times (denoted by -1)*)
midtree112[old_, new_, al1_, al2_,tree1_, tree2_] :=
    Module[ {v11, v12, v10, v13, t10, v31, v32, v33, t30, v34,share,tlist1, tlist2, 
        nv11, nv12, nv13,nv33,nvy, vy,n,res,group},
        {t10, {v11, v12}, v13, v10} = Take[Flatten[old, 1], 4];
        {t30, {v31, v32}, v33} = Take[Flatten[new, 1], 3];
        share = Intersection[{v11,v12},{v31,v32}];
        tlist1 = First[tree1];
        tlist2 = First[tree2];
        n = LeafNumber[tree1];
        res = Switch[Length[share],          
            2, 
            group = "11";
            If[ v33 == v13,
                {1, If[ t30 < t10,
                        {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}}},
                        {{{{v11, v10}, t30, {v10, v13}}, {{v12, v10},t30, {v10, v13}}}}
                    ]},
                {-1, -1}
            ],
            1,
            If[ v31!=share[[1]],
                {v31,v32} = {v32,v31}
            ];
            If[ v11!=share[[1]],
                {v11,v12} = {v12,v11}
            ];
            vy = al1[[v32]];            
            (*nv: label in tree 2; v: label in tree1*)
            {nv12, nv13,nvy} = {v12,v13,vy} /. {(x_ /; x > n) :>Position[tlist2, tlist1[[x]]][[1, 1]]};
            Which[
                vy==v33&&al2[[nv12]]==nv13,
                group = "11";
                {1,{{{{v11,v10},t30,{v32,v33}}}}},
                (!MemberQ[{v13,v11},vy])&&{v13,al2[[nv12]]}=={v33,nvy},
                (*t30=!=t10*)
                group = "21";
                {2,If[ t30>t10,
                       {{{{v12,v10},-1,{v32,vy}}},{{{v32,vy},t30,{v10,v13}}}},
                       {{{{v12,v10},-1,{v32,vy}}},{{{v32,vy},t30,{v11,v10}}}}
                   ]
                },
                vy==v11&&{v13,al2[[nv12]]}=={v33,nvy},
                group = "21";
                v34 = First[Complement[tree1[[2, vy-n, All, 1]],{v32}]];
                {2,If[ t30>t10,
                       {{{{v12,v10},-1,{v32,vy}}},{{{v32,vy},t30,{v10,v13}}},{{{v34,vy},t30,{v10,v13}}}},
                       {{{{v12,v10},-1,{v32,vy}}},{{{v32,vy},t30,{v12,v10}}},{{{v32,vy},t30,{v11,v10}},{{v34,vy},t30,{v11,v10}}}}
                   ]
                },
                True,
                group = "-1";
                {-1,-1}
            ],
            0,
            group = "21";
            If[ !MemberQ[al1[[{v31,v32}]],v33],
                Return[{-1,-1}]
            ];
            If[ al1[[v31]]!=v33,
                {v31,v32} = {v32,v31};
            ];
            vy = al1[[v32]];
            (*nv: label in tree 2; v: label in tree1*)
            {nv11, nv12, nv13,nv33,nvy} = {v11, v12, v13,v33,vy} /. {(x_ /; x > n) :>Position[tlist2, tlist1[[x]]][[1, 1]]};
            If[ !MemberQ[al2[[{nv11,nv12}]],nv13],
                Return[{-1,-1}]
            ];
            If[ al2[[nv11]]!=nv13,
                {v11,v12} = {v12,v11};
                {nv11,nv12} = {nv12,nv11}
            ];
            Which[
                vy==v33&&al2[[nv12]]==nv33,                   
                {2, Join[
                    {{{{v12, v10}, -1, {v31, v33}}},
                    {{{v12, v10}, -1, {v32, v33}}},
                    {{{v31,v33},t30,{v32,v33}},{{v32,v33},t30,{v31,v33}}}},
                    If[ MemberQ[{v31,v32},v13],
                        {{{{First[Complement[tree1[[2,v13-n,All,1]],{v10}]],v13},t30,
                        {First[Complement[{v31,v32},{v13}]],v33}}}},
                        {}
                    ]
                    ]
                    },
                vy==v13&&al2[[nv12]]==nv13,                
                {2, Join[
                    {{{{v12, v10}, -1, {v32, v13}}},
                    {{{v11, v10}, -1, {v32, v13}}},
                    If[ v31==v13,
                        {{{v32,v13},t30,{v13,v33}},{{v10,v13},t30,{v13,v33}}},
                        {{{v32,  v13}, t30, {v31, v33}}}
                    ]},
                    If[ MemberQ[{v11,v12},v33],
                        {{{{v31, v33}, t30,{v32,v13}}}},
                        {}
                    ]]
                    },
                (!MemberQ[{v13,v33},vy])&&al2[[nv12]]==nvy,                
                (*If[{v13,al1[[v13]]}=={v33,vy},{{{{v31, v33}, t30,{v32,vy}}}},{}]*)                
                {2, Join[
                    {{{{v12, v10}, -1, {v32, vy}}},{{{v32,  vy}, t30, {v31, v33}}}},
                    If[ v12==v33,
                        {{{{v31, v33}, t30,{v32,vy}}}},
                        {}
                    ], 
                    If[ v33=!=v32&&al1[[v13]]==vy&&v13==v33,
                        {{{{v31, v33}, t30,{v32,vy}}}},
                        {}
                    ],                                         
                    If[ v32==v13,
                        {{{{First[Complement[tree1[[2,v13-n,All,1]],{v10}]],v13},t30,{v31,v33}}}},
                        {}
                    ]
                    ]
                },
                True,
                {-1,-1}                    
            ]                       
        ];
        If[ res==={-1,-1},
            {"-1",-1,-1},
            Prepend[res,group]
        ]
    ]    

midtree2211[old_, new_, al1_] :=
    Module[ {t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, 
        v41, v42, v43,res,group},
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        res = Which[
            v23==v43&&v21==v41&&Complement[{v31, v32}, {v11, v12}] === {},
            group = "22TT1";
            Join[
                Which[
                t30 < t10, {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}}},
                t10 < t30 < t20, {{{{v11, v10}, t30, {v10, v20}}, {{v12, v10}, t30, {v10, v20}}}},
                t30 > t20, {{{{v11, v10}, t30, {v20, v23}}}, {{{v12, v10}, t30, {v20, v23}}}}                   
                ],
                Which[
                t40 < t10, {{{{v21, v20}, t40, {v12, v10}}}, {{{v21, v20}, t40, {v11, v10}}}},
                t10 < t40 < t20, {{{{v21, v20}, t40, {v10, v20}}, {{v10, v20}, t40, {v21, v20}}}},
                t40 > t20, {{{{v21, v20}, t40, {v20, v23}}, {{v10, v20}, t40, {v20, v23}}}, {{{v11, v10}, t40, {v20, v23}}}, {{{v12, v10}, t40, {v20, v23}}}}
                ],
                If[ t40 < t20,
                    {{{{v12, v10}, t40, {v21, v20}}}, {{{v11, v10}, t40, {v21, v20}}}},
                    {}
                ]
            ],          
            v23==v43&&MemberQ[{v31,v32},v21]&&Complement[{v41,v31, v32}, {v21,v11, v12}] === {},
            group = "22TT2";
            If[ v31 != v21,
                {v31, v32} = {v32, v31}
            ];
            Join[
                 Which[
                  t30 > t20, {{{{v32, v10}, t30, {v20, v23}}},{{{v21, v20}, t30, {v20, v23}}, {{v10, v20}, t30, {v20, v23}}}},
                  t10 < t30 < t20, {{{{v21, v20}, t30, {v10, v20}}, {{v10, v20},t30, {v21, v20}}}},
                  t30 < t10, {{{{v21, v20}, t30, {v32, v10}}}}                                    
                  ],
                 If[ t30<t20,
                     {{{{v32, v10}, t30, {v21, v20}}}},
                     {}
                 ],
                 Which[
                  t40 > t20, {{{{v21, v20}, t40, {v20, v23}}, {{v10, v20},t40, {v20, v23}}},{{{v32, v10}, t40, {v20, v23}}},{{{v41, v10}, t40, {v20, v23}}}},
                  t10 < t40 < t20, {{{{v21, v20}, t40, {v10, v20}}, {{v10, v20}, t40, {v21, v20}}},{{{v11, v10}, t40, {v10, v20}}, {{v12, v10},t40, {v10, v20}}}},                          
                  t40 < t10, {{{{v11, v10}, t40, {v12, v10}}, {{v12, v10}, t40, {v11, v10}}},{{{v21,v20},t40,{v41,v10}}}}                          
                  ],
                 If[ t40<t20,
                     {{{{v41, v10}, t40, {v21, v20}}}},
                     {}
                 ]       
                ],
           ((al1[[v31]]==v43&&Complement[{v41, v32}, {v21,v11, v12}] === {})||(al1[[v32]]==v43&&Complement[{v41, v31}, {v21,v11, v12}] === {})),           
           group = "22TT3";
           If[ al1[[v32]] != v43,
               {v31, v32} = {v32, v31}
           ];
           Join[
                {{{{v41, al1[[v41]]}, t40, {v32, v43}}},{{{v31, al1[[v31]]},t30, {v32, v43}}}},                
                If[ Complement[{v11, v12}, {v31, v41}] === {}&& t30>t10,
                    {{{{v10, v20}, t30, {v32, v43}}}},
                    {}
                ],
                If[ Complement[{v11, v12}, {v31, v41}] === {}&& t40>t10,
                    {{{{v10, v20}, t40, {v32, v43}}}},
                    {}
                ]                
            ],
            al1[[v41]] == v43&&Complement[{v31, v32}, {v21,v11, v12}] === {},  
            group = "22TT4";
            If[ MemberQ[{v31, v32}, v21],
                Join[
                 {{{{v31, al1[[v31]]}, t40, {v41, v43}}},{{{v32, al1[[v32]]}, t40, {v41, v43}}}},                 
                 If[ t30 < t10,
                     {{{{v21, v20},  t30, {First[Complement[{v31, v32}, {v21}]], v10}}}},
                     {}
                 ],
                 If[ t30 < t20,
                     {{{{First[Complement[{v31, v32}, {v21}]], v10},t30, {v21, v20}}}},
                     {}
                 ]
                 ],
                Join[
                 {{{{v11, v10}, t40, {v41, v43}}}, {{{v12, v10}, t40, {v41, v43}}}}, 
                 Which[
                     t30 < t10, {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10},t30, {v11, v10}}}},
                     t10 < t30 <t20, {{{{v11, v10}, t30, {v10, v20}}, {{v12, v10},t30, {v10, v20}}}},                     
                     True, {}
                 ],
                 If[ t40 > t10,
                     {{{{v10, v20}, t40, {v41, v43}}}},
                     {}
                 ]                     
                 ]
            ],           
           True,-1
        ];
        If[ res===-1,
            {"-1",-1,-1},
            {group,2,res}
        ]
    ]

midtree2201[old_, new_, al1_] :=
    Module[ {len,t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, v41, v42, v43,group,res},
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        len = Length[#]&/@{Intersection[{v11, v12}, {v41, v31, v32}],Intersection[{v21, v22}, {v41, v31, v32}]};
        If[ len=={1,2},
            {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = {{t20, {v21, v22}, v23, v20}, {t10, {v11, v12}, v13, v10}};
            len = {2,1};
        ];
        res = Which[
            len=={1,1}&&al1[[v41]] == v43,
            group = "22FT1";
            If[ ! MemberQ[{v11, v12}, v31],
                {v31, v32} = {v32, v31}
            ];
            If[ MemberQ[{v11, v12}, v31]&&MemberQ[{v21, v22}, v32],
                Join[
                     {{{{v32, v20}, t40, {v41, v43}}}}, 
                     {{{{v31, v10}, t40, {v41, v43}}}},
                     If[ t30 < t20,
                         {{{{v31, v10}, t30, {v32, v20}}}},
                         {}
                     ],
                     If[ t30 < t10,
                         {{{{v32, v20}, t30, {v31, v10}}}},
                         {}
                     ]
                ],
                -1
            ],
            len=={1,1}&&(al1[[v31]] == v43 || al1[[v32]] == v43),
            group = "22FT2";
            If[ al1[[v32]] != v43,
                {v31, v32} = {v32, v31}
            ];
            If[ ! MemberQ[{v11, v12},  v31],
                {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = {{t20, {v21, v22}, v23, v20}, {t10, {v11, v12}, v13, v10}}
            ];
            If[ MemberQ[{v11, v12},  v31]&&MemberQ[{v21, v22}, v41],
                {{{{v31, v10}, t30, {v32, v43}}},{{{v41, v20}, t40, {v32, v43}}}},
                -1
            ],
            len=={2,1}&&v13==v43&&MemberQ[{v21,v22},v41],
            group = "22FT3";
            Join[
                Which[
                    t30 < t10, {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}}},
                    t30 > t10, {{{{v11, v10}, t30, {v10, v13}}, {{v12, v10}, t30, {v10, v13}}}}
                ], 
                Which[ 
                    t40<t10, {{{{v41, v20}, t40, {v12, v10}}}, {{{v41, v20}, t40, {v11, v10}}}},
                    t40>t10, {{{{v41, v20}, t40, {v10, v13}}}}
                ]
            ],
            len=={2,1}&&v13==v43&&MemberQ[{v11,v12},v41],
            group = "22FT4";
            If[ !MemberQ[{v21,v22},v31],
                {v31, v32} = {v32, v31}
            ];
            Join[
                Which[
                    t40 <t10, {{{{v11, v10}, t40, {v12, v10}}, {{v12, v10}, t40, {v11, v10}}}},
                    t40 >t10, {{{{v11, v10}, t40, {v10, v13}}, {{v12, v10}, t40, {v10, v13}}}}
                ], 
                Which[
                    t40 <t10, {{{{v31, v20},t40, {v41, v10}}}},
                    t40 >t10, {{{{v31, v20}, t40, {v10, v13}}}}
                ],                                           
                Which[ 
                    t30<t10,  {{{{v31, v20}, t30, {v32, v10}}}},
                    t30>t10,  {{{{v31, v20}, t30, {v10, v13}}}}                  
                ],
                If[ t30 < t20,
                    {{{{v32, v10}, t30, {v31, v20}}}},
                    {}
                ]
            ],
            True,-1
        ];
        If[ res===-1,
            {"-1",-1,-1},
            {group,2,res}
        ]
    ]

midtree2210[old_, new_, al1_] :=
    Module[ {len,t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, v41, v42, v43, x,group,res},
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        len = Length[#]&/@{Intersection[{v31,v32},{v21,v11,v12}],Intersection[{v41,v42},{v21,v11,v12}]};
        If[ len=={1,2},
            {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = {{t40, {v41, v42}, v43}, {t30, {v31, v32}, v33}};
            len = {2,1};
        ];
        If[ (len[[1]]==1)&&(! MemberQ[{v21, v11, v12}, v31]),
            {v31, v32} = {v32, v31}
        ];
        If[ (len[[2]]==1)&&(! MemberQ[{v21, v11, v12}, v41]),
            {v41, v42} = {v42, v41}
        ];
        res = Which[
            len=={1,1}&& al1[[{v32, v42}]] == {v33, v43}&&(!MemberQ[{v31, v41},v21]),
            group = "22TF1";
            Join[
                {{{{v31, al1[[v31]]}, t30, {v32, v33}}}}, 
                {{{{v41, al1[[v41]]}, t40, {v42, v43}}}}, 
                If[ t10 < t30,
                    {{{{v10, v20}, t30, {v32, v33}}}},
                    {}
                ], 
                If[ t10 < t40,
                    {{{{v10, v20}, t40, {v42, v43}}}},
                    {}
                ]
            ],
            len=={1,1}&& al1[[{v32, v42}]] == {v33, v43}&&MemberQ[{v31, v41},v21],
            group = "22TF2";
            {{{{v31, al1[[v31]]}, t30, {v32, v33}}},{{{v41, al1[[v41]]}, t40, {v42, v43}}}},
            len=={2,1}&&{v23,al1[[v42]],v41}=={v33,v43,v21},
            group = "22TF3";
            Join[
                 Which[
                     t30 < t10, {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}}},
                     t10 < t30 <  t20, {{{{v11, v10}, t30, {v10, v20}}, {{v12, v10},t30, {v10, v20}}}},
                     t30 > t20, {{{{v11, v10}, t30, {v20, v23}}}, {{{v12, v10}, t30, {v20, v23}}}}
                  ],
                 {{{{v21, v20}, t40, {v42, v43}}}}
            ],
            len=={2,1}&&{v23,al1[[v42]]}=={v33,v43}&&MemberQ[{v11,v12},v41],
            group = "22TF4";
            x = First[Complement[{v11, v12}, {v41}]];
            Join[
                {{{{v41, v10}, t40, {v42, v43}}}},
                Which[
                    t30 > t20, {{{{v21, v20}, t30, {v20, v23}}, {{v10, v20}, t30, {v20, v23}}}, {{{x, v10}, t30, {v20, v23}}}}, 
                    t10 < t30 < t20, {{{{v21, v20}, t30, {v10, v20}}, {{v10, v20}, t30, {v21, v20}}}},
                    t30 < t10, {{{{v21, v20}, t30, {x, v10}}}}
                ],                        
                If[ t30 < t20,
                    {{{{x, v10}, t30, {v21, v20}}}},
                    {}
                ],
                
                If[ t40 > t10,
                    {{{{v10, v20}, t40, {v42, v43}}}},
                    {}
                ]
            ],
            True,-1           
        ];
        If[ res===-1,
            {"-1",-1,-1},
            {group,2,res}
        ]
    ]    
        
midtree2200[inputold_, inputnew_,al1_] :=
    Module[ {old,new,cross,t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, v41, v42, v43,group,res},
        {old,new} = {inputold,inputnew};
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        cross = Outer[Intersection, old[[All, 2]], new[[All, 2]], 1, 1];
        res = Switch[Count[Flatten[cross, 1], {}],
         2, 
         group = "22FF1";
         If[ MatchQ[cross, {{{}, {__}}, {{__}, {}}}],
             old = Reverse[old];
             cross = Reverse[cross];
             {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
             {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
         ];
         If[ MatchQ[cross, {{{__}, {}}, {{}, {__}}}],
             If[ Switch[Length[cross[[1, 1]]], 2, v13 == v33, 1, al1[[Complement[{v31, v32}, cross[[1, 1]]]]] == {v33}, _, False] && 
                 Switch[Length[cross[[2, 2]]], 2, v23 == v43, 1, al1[[Complement[{v41, v42}, cross[[2, 2]]]]] == {v43}, _, False],
                 {firstTransition[{t10, {v11, v12}, v13, v10}, {t30, {v31, v32}, v33}],
                     firstTransition[{t20, {v21, v22}, v23, v20}, {t40, {v41, v42}, v43}]},
                 -1
             ],
             -1
         ],
         1,
         group = "22FF2";
         If[ MatchQ[cross, {{{_}, {}}, {{_}, {_}}}] || MatchQ[cross, {{{}, {_}}, {{_}, {_}}}],
             old = Reverse[old];
             cross = Reverse[cross];
         ];
         If[ MatchQ[cross, {{{_}, {_}}, {{}, {_}}}],
             new = Reverse[new];
             cross = Transpose[Reverse[Transpose[cross]]]
         ];
         If[ MatchQ[cross, {{{_}, {_}}, {{_}, {}}}],
             {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
             {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
             If[ First[cross[[1,1]]]!=v11,
                 {v11,v12} = {v12,v11}
             ];
             If[ First[cross[[2,1]]]!=v21,
                 {v21,v22} = {v22,v21}
             ];
             If[ v13==v33,
                 Join[
                     If[ t30>t10,
                         {{{{v21, v20}, t30, {v10, v13}}}},
                         {{{{v21, v20}, t30, {v11, v10}}}}
                     ],
                    {{{{v12,  v10}, t40, {First[Complement[{v41, v42}, {v12}]], v43}}}}
                 ],
                 -1
             ]
         ],
         0,
         group = "22FF3";
         If[ old[[All, 3]] =!= new[[All, 3]],
             old = Reverse[old];
             cross = Reverse[cross];
         ];
         old[[All, 2]] = Flatten[#] & /@ cross;
         new[[All, 2]] = Flatten[#] & /@ Transpose[cross];
         {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
         {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
         Join[
             {If[ t30 > t10,
                  {{{v21, v20}, t30, {v10, v13}}},
                  {{{v21, v20}, t30, {v11, v10}}}
              ], 
             If[ t40 > t20,
                 {{{v12, v10}, t40, {v20, v23}}},
                 {{{v12, v10}, t40, {v22, v20}}}
             ]},
             If[ v33 ==  v43,
                 {If[ t30 > t20,
                      {{{v11, v10}, t30, {v20, v23}}},
                      {{{v11, v10}, t30, {v21, v20}}}
                  ], 
                 If[ t40 > t10,
                     {{{v22, v20}, t40, {v10, v13}}},
                     {{{v22, v20}, t40, {v12, v10}}}
                 ]},
                 {}
             ]
         ],
         _, -1
         ];
        If[ res===-1,
            {"-1",-1,-1},
            {group,2,res}
        ]
    ]

TreeTransitionTypeII2[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {old, new, al1,al2, oldQ, newQ,temp,res},
        If[ SameCoalescentTreeQ[tree1,tree2],
            {0,0},
            {{old, new}, {al1, al2}} = contrasttree[tree1, tree2];
            res = Switch[{Length[old],Length[new]},         
             {1,1},
             midtree11[old, new, al1, al2,tree1, tree2],
             {2,2},
             oldQ = old[[2, 2, 2]] == old[[1, -2]];
             newQ = new[[2, 2, 2]] == new[[1, -2]];
             temp = Which[
               oldQ && newQ, midtree2211[old, new, al1],
               (! oldQ) && newQ, midtree2201[old, new, al1],
               oldQ && (! newQ), midtree2210[old, new, al1],
               (! oldQ) && (! newQ),midtree2200[old, new,al1]
             ],
             _, {"-1",-1,-1}
            ];
            Rest[res]
        ]
    ]       

(*TreeTransitionTypeII differs from TreeTransitionTypeII2 by
changing midtree11 into midtree112*)        
TreeTransitionTypeII[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {old, new, al1,al2, oldQ, newQ,temp,res},
        If[ SameCoalescentTreeQ[tree1,tree2],
            {"0",0,0},
            {{old, new}, {al1, al2}} = contrasttree[tree1, tree2];
            res = Switch[{Length[old],Length[new]},         
             {1,1},
             midtree112[old, new, al1, al2,tree1, tree2],
             {2,2},
             oldQ = old[[2, 2, 2]] == old[[1, -2]];
             newQ = new[[2, 2, 2]] == new[[1, -2]];
             temp = Which[
               oldQ && newQ, midtree2211[old, new, al1],
               (! oldQ) && newQ, midtree2201[old, new, al1],
               oldQ && (! newQ), midtree2210[old, new, al1],
               (! oldQ) && (! newQ),midtree2200[old, new,al1]
             ],
             _, {"-1",-1,-1}
            ];
            res
        ]
    ]      

midtreecase2211[old_, new_, al1_,tree1_] :=
    Module[ {t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, 
        v41, v42, v43, v44,v24,vx,condx,condy},
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        Which[
            v23==v43&&v21==v41&&Complement[{v31, v32}, {v11, v12}] === {},
            Join[
                Which[
                t30 < t10, {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}}},
                t10 < t30 < t20, {{{{v11, v10}, t30, {v10, v20}}, {{v12, v10}, t30, {v10, v20}}}},
                t30==t20, {{{{v11, v10}, -1, {v21, v20}}}, {{{v12, v10}, -1, {v21, v20}}}},
                t30 > t20, {{{{v11, v10}, t30, {v20, v23}}}, {{{v12, v10}, t30, {v20, v23}}}}                   
                ],
                Which[
                t40 < t10, {{{{v21, v20}, t40, {v12, v10}}}, {{{v21, v20}, t40, {v11, v10}}}},
                t40==t10, {{{{v21, v20}, -1, {v12, v10}}}, {{{v21, v20}, -1, {v11, v10}}}},
                t10 < t40 < t20, {{{{v21, v20}, t40, {v10, v20}}, {{v10, v20}, t40, {v21, v20}}}},
                t40 > t20, {{{{v21, v20}, t40, {v20, v23}}, {{v10, v20}, t40, {v20, v23}}}, {{{v11, v10}, t40, {v20, v23}}}, {{{v12, v10}, t40, {v20, v23}}}}
                ],
                If[ t40 < t20,
                    {{{{v12, v10}, t40, {v21, v20}}}, {{{v11, v10}, t40, {v21, v20}}}},
                    {}
                ]
            ],          
            v23==v43&&MemberQ[{v31,v32},v21]&&Complement[{v41,v31, v32}, {v21,v11, v12}] === {},
            If[ v31 != v21,
                {v31, v32} = {v32, v31}
            ];
            (*vx ==v32*)
            vx = First[Complement[{v11, v12}, {v41}]];
            Join[
                 Which[
                  t30 > t20, {{{{vx, v10}, t30, {v20, v23}}},{{{v21, v20}, t30, {v20, v23}}, {{v10, v20}, t30, {v20, v23}}}},
                  t10 < t30 < t20, {{{{v21, v20}, t30, {v10, v20}}, {{v10, v20},t30, {v21, v20}}}},
                  t30==t10,{{{{v21,v20},-1,{v41,v10}}}},
                  t30 < t10, {{{{v21, v20}, t30, {vx, v10}}}}                                    
                  ],
                 If[ t30<t20,
                     {{{{vx, v10}, t30, {v21, v20}}}},
                     {}
                 ],
                 Which[
                  t40 > t20, {{{{v21, v20}, t40, {v20, v23}}, {{v10, v20},t40, {v20, v23}}},{{{vx, v10}, t40, {v20, v23}}},{{{v41, v10}, t40, {v20, v23}}}},
                  t10 < t40 < t20, {{{{v21, v20}, t40, {v10, v20}}, {{v10, v20}, t40, {v21, v20}}},{{{v11, v10}, t40, {v10, v20}}, {{v12, v10},t40, {v10, v20}}}},                          
                  t40 < t10, {{{{v11, v10}, t40, {v12, v10}}, {{v12, v10}, t40, {v11, v10}}},{{{v21,v20},t40,{v41,v10}}}}                          
                  ],
                 If[ t40<t20,
                     {{{{v41, v10}, t40, {v21, v20}}}},
                     {}
                 ]       
                ],
            al1[[v41]] == v43&&Complement[{v31, v32}, {v21,v11, v12}] === {},                        
            If[ MemberQ[{v31, v32}, v21],
                condx = If[ al1[[v31]]==v10,
                            t40!=t20,
                            t40!=t10
                        ];
                condy = If[ al1[[v32]]==v10,
                            t40!=t20,
                            t40!=t10
                        ];
                Join[
                 If[ condx,
                     {{{{v31, al1[[v31]]}, t40, {v41, v43}}}},
                     {}
                 ],
                 If[ condy,
                     {{{{v32, al1[[v32]]}, t40, {v41, v43}}}},
                     {}
                 ],
                 If[ t30 < t10,
                     {{{{v21, v20},  t30, {First[Complement[{v31, v32}, {v21}]], v10}}}},
                     {}
                 ],
                 If[ t30 < t20,
                     {{{{First[Complement[{v31, v32}, {v21}]], v10},t30, {v21, v20}}}},
                     {}
                 ],
                 Which[
                     v23==v43&&t30==t20&&MemberQ[{v31,v32},v21], 
                     vx = First[Complement[{v11,v12},{v31,v32}]];
                     {{{{vx,v10},-1,{v41,v23}}},{{{vx,v10},-1,{v20,v23}}},{{{v41,v23},t40,{v20,v23}},{{v20,v23},t40,{v41,v23}}}},                     
                     v23=!=v43&&t30==t20&&MemberQ[{v31,v32},v21],                     
                     v24 = First[Complement[tree1[[2,v23-LeafNumber[tree1],All,1]], {v20}]];
                     vx = First[Complement[{v11,v12},{v31,v32}]];
                     Join[
                         {{{{vx,v10},-1,{v20,v23}}}},                    
                         If[ v41==v23,
                             {{{{v20,v23},t40,{v41,v43}},{{v24,v23},t40,{v41,v43}}},{{{vx,v10},t40,{v41,v43}}}},
                             {{{{v20,v23},t40,{v41,v43}}}}
                         ]
                     ],        
                     True,
                     {}
                 ]
                 ],
                Join[
                 If[ t40!=t20,
                     {{{{v11, v10}, t40, {v41, v43}}}, {{{v12, v10}, t40, {v41, v43}}}},
                     {}
                 ], 
                 Which[
                     t30 < t10, {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10},t30, {v11, v10}}}},
                      t10 < t30 <t20, {{{{v11, v10}, t30, {v10, v20}}, {{v12, v10},t30, {v10, v20}}}},
                      True, {}
                      ],
                 If[ t40 > t10,
                     {{{{v10, v20}, t40, {v41, v43}}}},
                     {}
                 ]                     
                 ]
            ],
           ((al1[[v31]]==v43&&Complement[{v41, v32}, {v21,v11, v12}] === {})||(al1[[v32]]==v43&&Complement[{v41, v31}, {v21,v11, v12}] === {})),           
           If[ al1[[v32]] != v43,
               {v31, v32} = {v32, v31}
           ];
           condx = If[ al1[[v41]]==v10,
                       t40!=t20,
                       t40!=t10
                   ];
           condy = If[ al1[[v31]]==v10,
                       t30!=t20,
                       t30!=t10
                   ];
           Join[
                If[ condx,
                    {{{{v41, al1[[v41]]}, t40, {v32, v43}}}},
                    {}
                ],
                If[ condy,
                    {{{{v31, al1[[v31]]},t30, {v32, v43}}}},
                    {}
                ],
                If[ Complement[{v11, v12}, {v31, v41}] === {}&& t30>t10,
                    {{{{v10, v20}, t30, {v32, v43}}}},
                    {}
                ],
                If[ Complement[{v11, v12}, {v31, v41}] === {}&& t40>t10,
                    {{{{v10, v20}, t40, {v32, v43}}}},
                    {}
                ],
                vx = First[Complement[{v21,v11, v12},{v41, v31}]];
                Which[
                    v23==v43&&t40==t20&&v41==v21,
                    Join[{{{{vx,v10},-1,{v32,v23}}}},
                        Which[
                            t30>t10,{{{{v32,v23},t30,{v10,v20}}}},
                            t30<t10,{{{{v32,v23},t30,{v31,v10}}}},
                            True,{}
                        ]
                    ],
                    v23 == v43&&t40==t20&&v21==v31,
                    Join[{{{{vx, v10}, -1, {v32, v23}}}},
                        If[ MemberQ[{t10, t20}, t30],
                            {},
                            {{{{v32, v23}, t30, {v21, v20}}}}
                        ]],                    
                    v23==v43&&t40==t10&&v21==vx,
                    {{{{v21,v20},-1,{v32,v23}}},{{{v10,v20},-1,{v32,v23}}},{{{v32,v23},t30,{v31,v10}}}},
                    v23=!=v43&&t40==t10&&v21==vx,
                    v44 = First[Complement[tree1[[2, v43 - LeafNumber[tree1], All, 1]],{v32}]];
                    Join[{{{{v10,v20},-1,{v32,v43}}},{{{v32,v43},t30,{v31,v10}}}},
                        If[ v21==v43,
                            {{{{v44,v43},tree1[[1,v43]],{v10,v20}}}},
                            {}
                        ]
                        ],
                    True,
                    {}
                ]
            ],
           True,-1
        ]
    ]

midtreecase2201[old_, new_, al1_] :=
    Module[ {len,t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, v41, v42, v43},
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        len = Length[#]&/@{Intersection[{v11, v12}, {v41, v31, v32}],Intersection[{v21, v22}, {v41, v31, v32}]};
        Which[
            MatchQ[len, {1, 2} | {2, 1}],
            If[ len=={1,2},
                {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = {{t20, {v21, v22}, v23, v20}, {t10, {v11, v12}, v13, v10}};
            ];
            If[ v13==v43,
                If[ MemberQ[{v21,v22},v41],
                    Join[
                        If[ t30 < t10,
                            {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}}},
                            {{{{v11, v10}, t30, {v10, v13}}, {{v12, v10}, t30, {v10, v13}}}}
                        ],
                        If[ t40 < t10,
                            {{{{v41, v20}, t40, {v12, v10}}}, {{{v41, v20}, t40, {v11, v10}}}},
                            {{{{v41, v20}, t40, {v10, v13}}}}
                        ]
                    ],
                    If[ !MemberQ[{v21,v22},v31],
                        {v31, v32} = {v32, v31}
                    ];
                    Join[
                        If[ t40 <t10,
                            {{{{v11, v10}, t40, {v12, v10}}, {{v12, v10}, t40, {v11, v10}}}},
                            {{{{v11, v10}, t40, {v10, v13}}, {{v12, v10}, t40, {v10, v13}}}}
                        ], 
                        If[ t40 <t10,
                            {{{{v31, v20},t40, {v41, v10}}}},
                            {{{{v31, v20}, t40, {v10, v13}}}}
                        ],                                           
                        If[ t30<t10,
                            {{{{v31, v20}, t30, {v32, v10}}}},
                            {{{{v31, v20}, t30, {v10, v13}}}}
                        ],
                        If[ t30 < t20,
                            {{{{v32, v10}, t30, {v31, v20}}}},
                            {}
                        ]
                    ]
                ],
                -1
            ],
            len=={1,1},
            Which[
                al1[[v41]] == v43,
                If[ ! MemberQ[{v11, v12}, v31],
                    {v31, v32} = {v32, v31}
                ];
                If[ MemberQ[{v11, v12}, v31]&&MemberQ[{v21, v22}, v32],
                    Join[
                     {{{{v32, v20}, t40, {v41, v43}}}}, 
                     {{{{v31, v10}, t40, {v41, v43}}}},
                     If[ t30 < t20,
                         {{{{v31, v10}, t30, {v32, v20}}}},
                         {}
                     ],
                     If[ t30 < t10,
                         {{{{v32, v20}, t30, {v31, v10}}}},
                         {}
                     ]],
                    -1
                ],
                al1[[v31]] == v43 || al1[[v32]] == v43,
                If[ al1[[v32]] != v43,
                    {v31, v32} = {v32, v31}
                ];
                If[ ! MemberQ[{v11, v12},  v31],
                    {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = {{t20, {v21, v22}, v23, v20}, {t10, {v11, v12}, v13, v10}}
                ];
                If[ MemberQ[{v11, v12},  v31]&&MemberQ[{v21, v22}, v41],
                    {{{{v31, v10}, t30, {v32, v43}}},{{{v41, v20}, t40, {v32, v43}}}},
                    -1
                ],
                True,-1
            ],
            True,-1
        ]
    ]

midtreecase2210[old_, new_, al1_] :=
    Module[ {len,t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, v41, v42, v43, x},
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        len = Length[#]&/@{Intersection[{v31,v32},{v21,v11,v12}],Intersection[{v41,v42},{v21,v11,v12}]};
        Which[
            MatchQ[len, {1, 2} | {2, 1}],
            If[ len=={1,2},
                {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = {{t40, {v41, v42}, v43}, {t30, {v31, v32}, v33}}
            ];
            If[ ! MemberQ[{v21, v11, v12}, v41],
                {v41, v42} = {v42, v41}
            ];
            If[ v23==v33&&al1[[v42]]==v43,
                If[ v41 == v21,
                    Join[
                     Which[
                         t30 < t10, {{{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}}},
                         t10 < t30 <  t20, {{{{v11, v10}, t30, {v10, v20}}, {{v12, v10},t30, {v10, v20}}}},
                         t30 > t20, {{{{v11, v10}, t30, {v20, v23}}}, {{{v12, v10}, t30, {v20, v23}}}}
                      ],
                     {{{{v21, v20}, t40, {v42, v43}}}}
                     ], 
                    (*case MemberQ[{v11,v12},v41]*)
                    x = First[Complement[{v11, v12}, {v41}]];
                    Join[
                            Which[
                            t30 > t20, {{{{v21, v20}, t30, {v20, v23}}, {{v10, v20}, t30, {v20, v23}}}, {{{x, v10}, t30, {v20, v23}}}}, 
                            t10 < t30 < t20, {{{{v21, v20}, t30, {v10, v20}}, {{v10, v20}, t30, {v21, v20}}}},
                            t30 < t10, {{{{v21, v20}, t30, {x, v10}}}}
                            ],                        
                            If[ t30 < t20,
                                {{{{x, v10}, t30, {v21, v20}}}},
                                {}
                            ],
                            {{{{v41, v10}, t40, {v42, v43}}}},
                            If[ t40 > t10,
                                {{{{v10, v20}, t40, {v42, v43}}}},
                                {}
                            ]
                        ]
                ],
                -1
            ],
            len=={1,1},
            If[ ! MemberQ[{v21, v11, v12}, v31],
                {v31, v32} = {v32, v31}
            ];
            If[ ! MemberQ[{v21, v11, v12}, v41],
                {v41, v42} = {v42, v41}
            ];
            If[ al1[[{v32, v42}]] == {v33, v43},
                Join[
                    {{{{v31, al1[[v31]]}, t30, {v32, v33}}}}, 
                    {{{{v41, al1[[v41]]}, t40, {v42, v43}}}}, 
                    If[ Complement[{v11, v12}, {v31, v41}] === {} && t10 < t30,
                        {{{{v10, v20}, t30, {v32, v33}}}},
                        {}
                    ], 
                    If[ Complement[{v11, v12}, {v31, v41}] === {} && t10 < t40,
                        {{{{v10, v20}, t40, {v42, v43}}}},
                        {}
                    ]
                ],
                -1
            ],
            True,-1                        
        ]
    ]

midtreecase2200[inputold_, inputnew_,al1_,tree1_] :=
    Module[ {old,new,cross,t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, v41, v42, v43,v14},
        {old,new} = {inputold,inputnew};
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        cross = Outer[Intersection, old[[All, 2]], new[[All, 2]], 1, 1];
        Switch[Count[Flatten[cross, 1], {}],
         2,         
         If[ MatchQ[cross, {{{}, {__}}, {{__}, {}}}],
             old = Reverse[old];
             cross = Reverse[cross];
             {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
             {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
         ];
         If[ MatchQ[cross, {{{__}, {}}, {{}, {__}}}],
             If[ Switch[Length[cross[[1, 1]]], 2, v13 == v33, 1, al1[[Complement[{v31, v32}, cross[[1, 1]]]]] == {v33}, _, False] && 
                 Switch[Length[cross[[2, 2]]], 2, v23 == v43, 1, al1[[Complement[{v41, v42}, cross[[2, 2]]]]] == {v43}, _, False],
                 Join[
                     If[ t30 != t20,
                         {firstTransition[{t10, {v11, v12}, v13, v10}, {t30, {v31, v32}, v33}]},
                         {}
                     ],
                     If[ t40 !=t10,
                         {firstTransition[{t20, {v21, v22}, v23, v20}, {t40, {v41, v42}, v43}]},
                         {}
                     ]
                 ],
                 -1
             ],
             -1
         ],
         1,
         If[ MatchQ[cross, {{{_}, {}}, {{_}, {_}}}] || MatchQ[cross, {{{}, {_}}, {{_}, {_}}}],
             old = Reverse[old];
             cross = Reverse[cross];
         ];
         If[ MatchQ[cross, {{{_}, {_}}, {{}, {_}}}],
             new = Reverse[new];
             cross = Transpose[Reverse[Transpose[cross]]]
         ];
         If[ MatchQ[cross, {{{_}, {_}}, {{_}, {}}}],
             {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
             {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
             If[ First[cross[[1,1]]]!=v11,
                 {v11,v12} = {v12,v11}
             ];
             If[ First[cross[[2,1]]]!=v21,
                 {v21,v22} = {v22,v21}
             ];
             If[ v13==v33,
                 Join[
                     Which[ 
                     t30 >t10, {{{{v21, v20}, t30, {v10, v13}}}},
                     t30<t10,  {{{{v21, v20}, t30, {v11, v10}}}},
                     t30==t10, {{{{v21, v20}, -1, {v12, v10}}}}
                      ],
                    If[ t40!=t20,
                        {{{{v12,  v10}, t40, {First[Complement[{v41, v42}, {v12}]], v43}}}},
                        {}
                    ]
                 ],
                 -1
             ]
         ],
         0,
         If[ old[[All, 3]] =!= new[[All, 3]],
             old = Reverse[old];
             cross = Reverse[cross];
         ];
         old[[All, 2]] = Flatten[#] & /@ cross;
         new[[All, 2]] = Flatten[#] & /@ Transpose[cross];
         {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
         {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
         Join[
             {Which[ 
                 t30 > t10, {{{v21, v20}, t30, {v10, v13}}},
                 t30<t10,  {{{v21, v20}, t30, {v11, v10}}},
                 t30==t10,{{{v21, v20}, -1, {v12, v10}}}
              ], 
             Which[ 
                 t40 > t20,  {{{v12, v10}, t40, {v20, v23}}},
                t40<t20, {{{v12, v10}, t40, {v22, v20}}},
                t40==t20,{{{v12, v10}, -1, {v21, v20}}}
             ]},
             If[ MemberQ[{v21,v22},v13]&&t40==t20,
                 If[ v22!=v13,
                     {v21,v22} = {v22,v21}
                 ];
                 v14 = First[Complement[tree1[[2,v13-LeafNumber[tree1],All,1]],{v10}]];
                 {{{{v14,v13},tree1[[1,v13]],{v21,v20}}}, {{{First[Intersection[{v11,v12},{v31,v32}]],v10},t30,{v21,v20}}}},
                 {}
             ],
             If[ v33 ==  v43,
                 {Which[ 
                     t30 > t20,  {{{v11, v10}, t30, {v20, v23}}},
                     t30==t20,{{{v11, v10}, -1, {v22, v20}}},
                    t30<t20, {{{v11, v10}, t30, {v21, v20}}}
                  ], 
                 Which[ 
                     t40 > t10, {{{v22, v20}, t40, {v10, v13}}},
                     t40==t10,{{{v22, v20}, -1, {v11, v10}}},
                    t40<t10,{{{v22, v20}, t40, {v12, v10}}}
                 ]},
                 {}
             ]
         ],
         _, -1
         ]
    ]
    
treediff[tree1_, tree2_] :=
    Module[ {n, oldnode, newnode, eleft1, eleft2, tlist1, tlist2, elist1,
       elist2, told, tnew, cc12, cc21, rule12, rule21, diff},
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        AppendTo[tlist1, Infinity];
        AppendTo[tlist2, Infinity];
        elist1 = Sort[#] & /@ elist1;
        elist2 = Sort[#] & /@ elist2;
        n = LeafNumber[tree1];
        elist1 = Append[elist1, {{-2, -1}, {2 n - 1, -1}}];
        elist2 = Append[elist2, {{-2, -1}, {2 n - 1, -1}}];
        tnew = Complement[tlist2, tlist1];
        told = Complement[tlist1, tlist2];
        oldnode = Flatten[Position[tlist1, #] & /@ told];
        newnode = Flatten[Position[tlist2, #] & /@ tnew];
        eleft1 = Delete[Transpose[{Range[2 n - 1], Most[tlist1]}], Partition[oldnode, 1]];
        eleft2 = Delete[Transpose[{Range[2 n - 1], Most[tlist2]}], Partition[newnode, 1]];
        diff =  Select[Transpose[{eleft1[[n + 1 ;;, 1]], eleft2[[n + 1 ;;, 1]]}],Unequal @@ # &];
        rule12 = Join[Thread[oldnode -> Take[CharacterRange["A", "Z"], Length[oldnode]]],Thread[diff[[All, 1]] -> diff[[All, 2]]]];
        rule21 = Join[Thread[newnode -> Take[CharacterRange["A", "Z"], Length[newnode]]], Thread[diff[[All, 2]] -> diff[[All, 1]]]];
        cc12 = SortBy[Complement[elist1,Replace[elist2, rule21, {3}]], tlist1[[#[[1, 2]]]] &];
        cc21 = SortBy[Complement[elist2, Replace[elist1, rule12, {3}]], tlist2[[#[[1, 2]]]] &];
        {tlist1, tlist2, elist1, elist2, told, tnew, cc12, cc21, rule12,rule21}
    ];
       
(*oldtreediff[tree1_, tree2_] :=
    Module[ {tlist1, tlist2, elist1, elist2, n, tnew, newnode, 
      newlabelrule, temp, rule21, cc12},
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        tlist1 = Append[tlist1, Infinity];
        tlist2 = Append[tlist2, Infinity];
        elist1 = Sort[#] & /@ elist1;
        elist2 = Sort[#] & /@ elist2;
        n = LeafNumber[tree1];
        elist1 = Append[elist1, {{-2, -1}, {2 n - 1, -1}}];
        elist2 = Append[elist2, {{-2, -1}, {2 n - 1, -1}}];
        tnew = Complement[tlist2, tlist1];
        newnode = Flatten[Position[tlist2, #] & /@ tnew];
        newlabelrule = 
         Thread[newnode -> Take[CharacterRange["A", "Z"], Length[newnode]]];
        temp = Complement[Range[n + 1, 2 n - 1], newnode];
        rule21 = 
         Thread[temp -> (Position[tlist1, #][[1, 1]] & /@ tlist2[[temp]])];
        rule21 = Join[newlabelrule, rule21];
        cc12 = Complement[elist1, elist2 /. rule21];
        cc12 = SortBy[cc12, tlist1[[#[[1, 2]]]] &];
        {tlist1, elist1, tnew, cc12, rule21}
    ]    
*)

gettreeoperator1[{tlist1_, c12_, rule21_, al1_}, {tlist2_, c21_,rule12_, al2_}] :=
    Module[ {cc12,cc21,cc21to1, v10, v20, v30, v40, v50, v60, t10, v21, v22, v31, v51,
      newt, v11, v12, b1, b2, set, set1, set2, temp1, temp2},
        cc12 = c12;
        cc21 = c21;
        cc21to1 = cc21 /. rule21;
        Switch[{Length[cc12], Length[cc21]},
         {2, 2},
         If[ al1[[cc12[[1, 1, 2]]]] == cc12[[2, 1, 2]] ==cc21to1[[2, 1, 2]]&&al2[[cc21[[1, 1, 2]]]] == cc21[[2, 1, 2]],
             {v10, v20} = cc12[[All, 1, 2]];
             v21 = First[Complement[cc12[[2, All, 1]], {v10}]];
             t10 = tlist1[[v10]];
             newt = tlist2[[cc21[[1, 1, 2]]]];
             Which[
              MemberQ[cc21to1[[2, All, 1]], v21],
              {v11, v12} = cc12[[1, All, 1]];
              Which[
               newt < t10, {{{v11, v10}, newt, {v12, v10}}, {{v12, v10},newt, {v11, v10}}},
               newt > t10, {{{v11, v10}, newt, {v10, v20}}, {{v12, v10},newt, {v10, v20}}},
               _, -1
               ],
              MemberQ[cc21to1[[1, All, 1]], v21],
              v11 = First[Complement[cc21to1[[1, All, 1]], {v21}]];
              {{{v11, v10}, newt, {v21, v20}}},
              True, -1
              ],
             -1
         ],
         {3, 3},
         b1 = {Most[al1[[cc12[[All, 1, 2]]]]] == Rest[cc12[[All, 1, 2]]], 
           cc12[[;; 2, 1, 2]] === cc12[[-1, All, 1]]};
         b2 = {Most[al2[[cc21[[All, 1, 2]]]]] == Rest[cc21[[All, 1, 2]]], 
           cc21[[;; 2, 1, 2]] === cc21[[-1, All, 1]]};
         Switch[{b1, b2},
          {{True, False}, {False, True}},
          If[ cc12[[-1, 1, 2]] == cc21to1[[-1, 1, 2]]&&MemberQ[cc21to1[[;; 2, 1, 2]], cc12[[2, 1, 2]]],
              If[ cc12[[2, 1, 2]]!= cc21to1[[2, 1, 2]],
                  cc21to1[[{1,2}]] = cc21to1[[{2,1}]];
                  cc21[[{1,2}]] = cc21[[{2,1}]];
              ];
              {v10, v20,v30} = cc12[[All, 1, 2]];
              v21 = First[Complement[cc12[[2, All, 1]], {v10}]];
              v31 = First[Complement[cc12[[3, All, 1]], {v20}]];
              If[ MemberQ[cc21to1[[2,All,1]],v21]&&MemberQ[cc21to1[[1,All,1]],v31],
                  {v11, v31} = First[Intersection[cc21to1[[1, All, 1]], #]] & /@cc12[[{1, 3}, All, 1]];
                  newt = tlist2[[cc21[[1, 1, 2]]]];
                  {{{v11, v10}, newt, {v31, v30}}},
                  -1
              ],
              -1
          ],
          {{False, True}, {True, False}},
          If[ cc12[[-1, 1, 2]] == cc21to1[[-1, 1, 2]]&&MemberQ[cc12[[;; 2, 1, 2]], cc21to1[[2, 1, 2]]],
              If[ cc21to1[[2, 1, 2]]!=cc12[[2, 1, 2]],
                  cc12[[{1,2}]] = cc12[[{2,1}]];
              ];
              v51 = First[Complement[cc21to1[[2, All, 1]], cc21to1[[;; 2, 1, 2]]]];
              If[ MemberQ[cc12[[2, All, 1]], v51]&&Length[Intersection[cc12[[1, All, 1]],cc21to1[[1, All, 1]]]]==1,
                  {v11, v22} = First[Intersection[#, cc21to1[[1, All, 1]]]] & /@cc12[[;; 2, All, 1]];
                  {v10, v20} = cc12[[;; 2, 1, 2]];
                  newt = tlist2[[cc21[[1, 1, 2]]]];
                  {{{v11, v10}, newt, {v22, v20}}},
                  -1
              ],
              -1
          ],
          {{True, False}, {True, False}},
          {v10, v20, v30} = cc12[[All, 1, 2]];
          {v40, v50, v60} = cc21to1[[All, 1, 2]];
          v21 = First[Complement[cc12[[2, All, 1]], {v10}]];
          v31 = First[Complement[cc12[[3, All, 1]], {v20}]];
          Which[
           !(v30 == v60&&MemberQ[cc21to1[[1, All, 1]], v21]&&MemberQ[cc21to1[[3, All, 1]], v31]), -1,
           v20==v40,
           v11 = First[Complement[cc21to1[[2, All, 1]], {v40, v50}]];
           newt = tlist2[[cc21[[2, 1, 2]]]];
           {{{v11, v10}, newt, {v20, v30}}},
           v10==v50,
           v12 = First[Intersection[cc21to1[[1, All, 1]], cc12[[1, All, 1]]]];
           newt = tlist2[[cc21[[1, 1, 2]]]];
           {{{v21, v20}, newt, {v12, v10}}},
           True, -1
           ],
          {{False, False}, {False, False}},
          set = Drop[Range[3], {#}] & /@ Reverse[Range[3]];
          set1 = Select[set,al1[[cc12[[First[#], 1, 2]]]] == cc12[[Last[#], 1, 2]] &];
          set1 = {#, Complement[Range[3], #]} & /@ set1;
          set2 = Select[set, al2[[cc21[[First[#], 1, 2]]]] == cc21[[Last[#], 1, 2]] &];
          set2 = {#, Complement[Range[3], #]} & /@ set2;
          If[ Length[set1] == Length[set2] == 1 &&cc12[[set1[[1, All, -1]], 1, 2]] == 
            Reverse[cc21to1[[set2[[1, All, -1]], 1, 2]]],
              v10 = cc12[[set1[[1, 1, 1]], 1, 2]];
              v20 = cc12[[set1[[1, 2, 1]], 1, 2]];
              temp1 = 
               Outer[Intersection, 
                cc12[[{set1[[1, 1, 1]], set1[[1, 2, 1]]}, All, 1]], 
                cc21to1[[{set2[[1, 1, 1]], set2[[1, 2, 1]]}, All, 1]], 1, 1];
              temp2 = Intersection[cc12[[set1[[1, 1, 2]], All, 1]],cc21to1[[set2[[1, 2, 1]], All, 1]]];
              If[ MatchQ[temp1, {{{_}, {_}}, {{_}, {}}}] && Length[temp2] == 1,
                  {v11, v22} = Flatten[temp1[[All, 1]]];
                  newt = tlist2[[cc21[[set2[[1, 1, 1]], 1, 2]]]];
                  {{{v11, v10}, newt, {v22, v20}}},
                  -1
              ],
              -1
          ],
          _,-1
          ],
          _,-1
         ]
    ]
    
nodediffI[{tlist1_,elist1_, cc12_, rule21_, al1_}, {tlist2_,elist2_, cc21_, rule12_, al2_}] :=
    Module[ {cc21to1, oldnode, newnode, cross, b1, b2, set, set1, set2,v31, v61, vvv1, vvv2,n,temp,
        set11,set12,set21,set22,v10,v20,v30,v40,v60,v70,v80},
        cc21to1 = cc21 /. rule21;
        oldnode = newnode = {};
        Switch[{Length[cc12], Length[cc21]},
             {2, 2},             
             If[ tlist1==tlist2&&al1[[cc12[[1, 1, 2]]]]!=cc12[[2,1,2]]&&al2[[cc21[[1, 1, 2]]]]!=cc21[[2,1,2]],
                 cross = Outer[Intersection, cc12[[All, All, 1]], cc21to1[[All, All, 1]], 1,1];
                 If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}}],
                     oldnode = cc12[[All, 1, 2]];
                     newnode = cc21[[All, 1, 2]]
                 ]
             ],
             {3, 3},
             b1 = {Most[al1[[cc12[[All, 1, 2]]]]] == Rest[cc12[[All, 1, 2]]], cc12[[;; 2, 1, 2]] === cc12[[-1, All, 1]]};
             b2 = {Most[al2[[cc21[[All, 1, 2]]]]] == Rest[cc21[[All, 1, 2]]], cc21[[;; 2, 1, 2]] === cc21[[-1, All, 1]]};
             Switch[{b1, b2},
                  {{True, False}, {True, False}},
                  oldnode = cc12[[;; 2, 1, 2]];
                  newnode = cc21[[;; 2, 1, 2]],
                  {{False, True}, {False, True}},
                  cross = Outer[Intersection, cc12[[;; 2, All, 1]], cc21to1[[;; 2, All, 1]],1, 1];
                  If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}}|{{{_,_}, {}}, {{}, {_,_}}}|{{{}, {_,_}}, {{_,_}, {}}}],
                      oldnode = cc12[[;; 2, 1, 2]];
                      newnode = cc21[[;; 2, 1, 2]]
                  ],
                  {{False, False}, {False, False}},
                  set = Drop[Range[3], {#}] & /@ Reverse[Range[3]];
                  set1 = Select[set,al1[[cc12[[First[#], 1, 2]]]] == cc12[[Last[#], 1, 2]] &];
                  set1 = {#, Complement[Range[3], #]} & /@ set1;
                  set2 = Select[set,al2[[cc21[[First[#], 1, 2]]]] == cc21[[Last[#], 1, 2]] &];
                  set2 = {#, Complement[Range[3], #]} & /@ set2;
                  Switch[{Length[set1], Length[set2]},
                     {0,0},
                     If[ tlist1 == tlist2,
                         vvv1 = Intersection[al1[[cc12[[All, 1, 2]]]], 
                           Flatten[cc12[[All, All, 1]]]];
                         vvv2 = Intersection[al2[[cc21[[All, 1, 2]]]], 
                           Flatten[cc21to1[[All, All, 1]]]];
                         cross = Outer[Intersection, cc12[[All, All, 1]], 
                           cc21to1[[All, All, 1]], 1, 1];
                         cross = Map[Length, cross, {2}] /. {x_ /; x > 1 :> 100};
                         If[ vvv1 === vvv2 && Length[vvv1] == 1 && 
                           Total[cross] === Total[cross, {2}] === {2, 2, 2},
                             n = Length[tlist1]/2;
                             temp = elist1[[First[vvv1] - n, All, 1]];
                             If[ Complement[temp, elist2[[First[vvv2] - n, All, 1]]] === {},
                                 If[ Length[Intersection[temp, cc12[[All, 1, 2]]]] == 1,
                                     oldnode = newnode = Union[vvv1, Intersection[temp, cc12[[All, 1, 2]]]]
                                 ]
                             ]
                         ]
                     ],
                     {1,1},
                     If[ cc12[[set1[[1, All, -1]], 1, 2]] ==cc21to1[[set2[[1, All, -1]], 1, 2]],
                         cross = Outer[Intersection, 
                           cc12[[{set1[[1, 1, 1]], set1[[1, 2, 1]]}, All, 1]], 
                           cc21to1[[{set2[[1, 1, 1]], set2[[1, 2, 1]]}, All, 1]], 1, 1];
                         If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}}],
                             oldnode = cc12[[{set1[[1, 1, 1]], set1[[1, 2, 1]]}, 1, 2]];
                             newnode = cc21[[{set2[[1, 1, 1]], set2[[1, 2, 1]]}, 1, 2]];
                         ]
                     ]
                  ]
                         
              ],
             {4, 4},
             set = Drop[Range[4], {#}] & /@ Reverse[Range[4]];
             set1 = Select[set,al1[[cc12[[Most[#], 1, 2]]]] == cc12[[Rest[#], 1, 2]] &];
             set1 = {#, Complement[Range[4], #]} & /@ set1;
             set2 = Select[set,al2[[cc21[[Most[#], 1, 2]]]] == cc21[[Rest[#], 1, 2]] &];
             set2 = {#, Complement[Range[4], #]} & /@ set2;
             Switch[{Length[set1], Length[set2]},
                  {0,0},                  
                  b1 = Select[Transpose[{cc12[[All, 1, 2]], al1[[cc12[[All, 1, 2]]]]}],First[#] != -1 && MemberQ[cc12[[All, 1, 2]], Last[#]] &];
                  b2 = Select[Transpose[{cc21[[All, 1, 2]], al2[[cc21[[All, 1, 2]]]]}],First[#] != -1 && MemberQ[cc21[[All, 1, 2]], Last[#]] &];
                  If[ Length[b1] == Length[b2] == 2 &&Complement[cc12[[All, 1, 2]], Flatten[b1]] === {} && 
                     Complement[cc21[[All, 1, 2]], Flatten[b2]] === {}&&Complement[b1[[All, 2]], b2[[All, 2]] /. rule21] === {},
                      cross = Outer[Intersection,
                        Select[cc12, MemberQ[b1[[All, 1]], #[[1, 2]]] &][[All, All, 1]],
                        Select[cc21, MemberQ[b2[[All, 1]], #[[1, 2]]] &][[All, All, 1]] /. rule21, 1, 1];
                      If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}}],
                          oldnode = b1[[All, 1]];
                          newnode = b2[[All, 1]]
                      ];
                  ],       
                  {1, 1},
                  0,
                  {1, 2},
                  temp = Select[set2, cc21to1[[#[[1, -1]], 1, 2]] == cc12[[set1[[1, 2, 1]], 1, 2]] &];
                  If[ temp==={},
                      set2 = Select[set2,cc21to1[[#[[1, 2]], 1, 2]] == cc12[[set1[[1, 2, 1]], 1, 2]] &],
                      set2 = temp
                  ],
                  {2, 1},
                  temp = Select[set1, cc12[[#[[1, -1]], 1, 2]] == cc21to1[[set2[[1, 2, 1]], 1, 2]] &];
                  If[ temp === {},
                      set1 = Select[set1,cc12[[#[[1, 2]], 1, 2]] == cc21to1[[set2[[1, 2, 1]], 1, 2]] &],
                      set1 = temp
                  ],
                  {2, 2},
                  b1 = {Most[al1[[cc12[[All, 1, 2]]]]] == Rest[cc12[[All, 1, 2]]], cc12[[;; 2, 1, 2]] === cc12[[3, All, 1]]};
                  b2 = {Most[al2[[cc21[[All, 1, 2]]]]] == Rest[cc21[[All, 1, 2]]], cc21[[;; 2, 1, 2]] === cc21[[3, All, 1]]};
                  Switch[{b1, b2},
                       {{True, False}, {True, False}},
                       Which[
                        cc12[[-2, 1, 2]] == cc21to1[[1, 1, 2]],
                        set1 = {{{1, 2, 3}, {4}}};
                        set2 = {{{2, 3, 4}, {1}}},
                        cc12[[1, 1, 2]] == cc21to1[[-2, 1, 2]],
                        set1 = {{{2, 3, 4}, {1}}};
                        set2 = {{{1, 2, 3}, {4}}};
                        ],
                       {{True, False}, {False, True}},
                       If[ MemberQ[cc21to1[[;; 2, 1, 2]], cc12[[-2, 1, 2]]],
                           set1 = {{{1, 2, 3}, {4}}};
                           set2 = If[ cc21to1[[set2[[1, 1, 1]], 1, 2]] ==cc12[[-2, 1, 2]],
                                      {set2[[2]]},
                                      {set2[[1]]}
                                  ]
                       ],
                       {{False, True}, {True, False}},
                       If[ MemberQ[cc12[[;; 2, 1, 2]], cc21to1[[-2, 1, 2]]],
                           set1 = If[ cc12[[set1[[1, 1, 1]], 1, 2]] == cc21[[-2, 1, 2]],
                                      {set1[[2]]},
                                      {set1[[1]]}
                                  ];
                           set2 = {{{1, 2, 3}, {4}}};
                       ]
                   ];
              ];
             If[ Length[set1] == Length[set2] == 1,
                 Which[
                 cc21to1[[set2[[1, 1, -1]], 1, 2]] ==cc12[[set1[[1, 2, 1]], 1, 2]]&&
                 cc12[[set1[[1, 1, -1]], 1, 2]] ==cc21to1[[set2[[1, 2, 1]], 1, 2]],
                 v31 = First[Complement[cc12[[set1[[1, 1, -1]], All, 1]], cc12[[set1[[1, 1]], 1, 2]]]];
                 v61 = First[Complement[cc21to1[[set2[[1, 1, -1]], All, 1]],cc21to1[[set2[[1, 1]], 1, 2]]]];
                 vvv1 = Complement[Flatten[cc12[[set1[[1, 1, ;; 2]], All, 1]]],cc12[[set1[[1, 1]], 1, 2]]];
                 vvv2 = Complement[Flatten[cc21to1[[set2[[1, 1, ;; 2]], All, 1]]], cc21to1[[set2[[1, 1]], 1, 2]]];
                 If[ MemberQ[cc21to1[[set2[[1, 2, 1]], All, 1]], v31]&&MemberQ[cc12[[set1[[1, 2, 1]], All, 1]], v61]&&
                     Length[Intersection[vvv1, vvv2]] == 2,
                     oldnode = cc12[[Most[set1[[1, 1]]], 1, 2]];
                     newnode = cc21[[Most[set2[[1, 1]]], 1, 2]];
                 ],
                 cc12[[Rest[set1[[1, 1]]], 1, 2]] ==cc21to1[[Rest[set2[[1, 1]]], 1, 2]],
                 cross = Outer[Intersection, cc12[[set1[[1, All, 1]], All, 1]], cc21to1[[set2[[1, All, 1]], All, 1]], 1, 1];
                 If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}} | {{{}, {_}}, {{_}, {_}}}],
                     oldnode = cc12[[set1[[1, All, 1]], 1, 2]];
                     newnode = cc21[[set2[[1, All, 1]], 1, 2]];
                 ],
                 tlist1==tlist2&&cc12[[-1, 1, 2]] == cc21to1[[-1, 1, 2]] && 
                 cc12[[set1[[1, 1, 1]], 1, 2]] == cc21to1[[set2[[1, 1, 2]], 1, 2]],
                 cross = Outer[Intersection, cc12[[set1[[1, All, 1]], All, 1]], 
                 cc21to1[[set2[[1, All, 1]], All, 1]], 1, 1];
                 If[ MatchQ[cross, {{{_}, {}}, {{_}, {_}}}],
                     oldnode = Sort[cc12[[{set1[[1, 1, 2]], set1[[1, 2, 1]]}, 1, 2]]];
                     newnode = Sort[cc21to1[[set2[[1, All, 1]], 1, 2]]];
                 ],
                 tlist1==tlist2&&cc12[[-1, 1, 2]] == cc21to1[[-1, 1, 2]] && 
                 cc21to1[[set2[[1, 1, 1]], 1, 2]] == cc12[[set1[[1, 1, 2]], 1, 2]],
                 cross = Outer[Intersection, cc12[[set1[[1, All, 1]], All, 1]], 
                 cc21to1[[set2[[1, All, 1]], All, 1]], 1, 1];
                 If[ MatchQ[cross, {{{_}, {_}}, {{}, {_}}}],
                     oldnode = Sort[cc12[[set1[[1, All, 1]], 1, 2]]];
                     newnode = Sort[cc21to1[[{set2[[1, 1, 2]], set2[[1, 2, 1]]}, 1, 2]]];
                 ],
                 tlist1 == tlist2 && 
                 cc12[[set1[[1, 1, 2]], 1, 2]] == cc21[[set2[[1, 2, 1]], 1, 2]] &&
                 al1[[cc12[[set1[[1, 1, 3]], 1, 2]]]] == 
                 al2[[cc21[[set2[[1, 2, 1]], 1, 2]]]] == 
                 cc12[[set1[[1, 2, 1]], 1, 2]] == cc21[[set2[[1, 1, 3]], 1, 2]],
                 cross = Outer[
                 Intersection, {cc12[[set1[[1, 1, 1]], All, 1]], 
                 Complement[Flatten[cc12[[set1[[1, 1, 2 ;;]], All, 1]]], 
                 cc12[[set1[[1, 1]], 1, 2]]]}, cc21[[set2[[1, All, 1]], All, 1]], 
                 1, 1];
                 If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}} | {{{_}, {_}}, {{}, {_}}}],
                     oldnode = newnode = Sort[cc12[[set1[[1, 1, {1, 3}]], 1, 2]]];
                 ],
                 tlist1 == tlist2 && 
                 cc21[[set2[[1, 1, 2]], 1, 2]] == cc12[[set1[[1, 2, 1]], 1, 2]] &&
                 al2[[cc21[[set2[[1, 1, 3]], 1, 2]]]] == 
                 al1[[cc12[[set1[[1, 2, 1]], 1, 2]]]] == 
                 cc21[[set2[[1, 2, 1]], 1, 2]] == cc12[[set1[[1, 1, 3]], 1, 2]],
                 cross = Outer[Intersection, 
                 cc12[[set1[[1, All, 1]], All, 
                 1]], {cc21[[set2[[1, 1, 1]], All, 1]], 
                 Complement[Flatten[cc21[[set2[[1, 1, 2 ;;]], All, 1]]], 
                 cc21[[set2[[1, 1]], 1, 2]]]}, 1, 1];
                 If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}}| {{{_}, {}}, {{_}, {_}}}],
                     oldnode = newnode = Sort[cc21[[set2[[1, 1, {1, 3}]], 1, 2]]];
                 ]
                 ]
             ],
             {5,5},
             If[ tlist1 == tlist2,
                 n = Length[elist1];
                 set11 = Select[cc12[[All, 1, 2]], MemberQ[cc12[[All, 1, 2]], al1[[#]]] &];
                 set12 = Select[set11, MemberQ[set11, al1[[#]]] &];
                 set21 = Select[cc21[[All, 1, 2]], MemberQ[cc21[[All, 1, 2]], al2[[#]]] &];
                 set22 = Select[set21, MemberQ[set21, al2[[#]]] &];
                 If[ Length[set11]==Length[set21]==3&&Length[set12]==Length[set22] == 1 &&Complement[set11, set21] === {},
                     v10 = First[set12];
                     v20 = al1[[v10]];
                     v30 = al1[[v20]];
                     v40 = First[Complement[set11, {v10, v20}]];
                     v60 = First[set22];
                     v70 = al2[[v60]];
                     v80 = First[Complement[set21, {v60, v70}]];
                     Which[
                      v10 == v80 && v20 == v60 && v40 == v70 &&v30 == al1[[v20]] == al2[[v70]] == al2[[v80]]&&
                       MemberQ[elist1[[v40 - n, All, 1]],First[Complement[elist2[[v70 - n, All, 1]], {v60}]]],
                      cross = Outer[Intersection, {elist1[[v10 - n, All, 1]],
                          Complement[Flatten[elist1[[{v30, v20} - n, All, 1]]], {v10, v20}]}, 
                        elist2[[{v60, v80} - n, All, 1]], 1, 1];
                      If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}}],
                          oldnode = newnode = {v10, v40};
                      ],
                      v10 == v70 && v20 == v80 && v40 == v60 &&v30 == al1[[v20]] == al1[[v40]] == al2[[v70]]&&
                       MemberQ[elist2[[v80 - n, All, 1]],First[Complement[elist1[[v20 - n, All, 1]], {v10}]]],
                      cross = Outer[Intersection,elist1[[{v10, v40} - n, All, 1]], {elist2[[v60 - n, All, 1]], 
                         Complement[Flatten[elist2[[{v30, v70} - n, All, 1]]], {v60, v70}]}, 1, 1];
                      If[ MatchQ[cross, {{{_}, {_}}, {{_}, {_}}}],
                          oldnode = newnode = {v20, v40};
                      ]
                      ]
                 ]
             ]  
             ];
        {oldnode, newnode}
    ] 
  
nodediffII[{tlist1_,elist1_, c12_, rule21_, al1_}, {tlist2_,elist2_, c21_, rule12_, al2_},isprint_:False] :=
    Module[ {n, oldnode, cc12, remain, edges, cond,set, i, j, v10, v13, vv, pos},
        n = Length[tlist1]/2;
        oldnode = Flatten[Position[tlist1, #] & /@ Complement[tlist1, tlist2]];
        cc12 = Select[c12, ! (MemberQ[oldnode, #[[1, 2]]] || #[[1, 2]] == -1) &];
        remain = Flatten[elist2 /. rule21, 1];
        remain = remain //. (Rule @@ # & /@ Select[remain, (StringQ[First[#]] &)]);
        edges = al1[[cc12[[All, 1, 2]]]];
        edges = Thread[#] & /@Transpose[{Transpose[Append[Transpose[cc12[[All, All, 1]]], cc12[[All, 1, 2]]]],edges}];
        If[ isprint,
            Print["edges=",edges]
        ];
        cond = Map[MemberQ[remain, #] &, edges, {2}];
        If[ isprint,
            Print["cond=",cond]
        ];
        set = Flatten[Position[cond, {False, True, True} | {True, False, True}, {1},Heads -> False]];
        Do[
         j = set[[i]];
         {v10, v13} = edges[[j, -1]];
         v13 = NestWhile[al1[[#]] &, v13, (MemberQ[oldnode, #] &)];
         If[ v13!=edges[[j,-1,-1]],
             Print["edges=",edges,";oldnode=",oldnode]
         ];
         If[ v13==-1,
             Continue[]
         ];
         vv = First[Complement[elist1[[v13 - n, All, 1]], {v10}]];
         vv = If[ vv > n,
                  Prepend[elist1[[vv - n, All, 1]], vv],
                  {vv}
              ];
         vv = Select[vv /. rule12, NumericQ];
         vv = (NestList[al2[[#]] &, #, 2] & /@ vv) /. rule21;
         If[ MemberQ[Flatten[vv], v10],
             cond[[j, 3]] = False
         ], {i,Length[set]}];
        If[ isprint,
            Print["cond=",cond]
        ];
        pos = (Xor[#[[1]], #[[2]]]) && (! #[[3]]) & /@ cond;
        pos = Flatten[Position[pos, True, {1}, Heads -> False]];
        If[ pos =!= {},
            oldnode = Union[oldnode, cc12[[pos, 1, 2]]];
        ];
        oldnode
    ]
  
nodediffIII[oldnode_, newnode_, n_, tlist1_,elist1_, al1_,rule12_, tlist2_, elist2_,al2_, rule21_] :=
    Module[ {v30, v31, v32, v33, v10, v11, v12, v13, nv31,nv32,nv33, nv11,nv12,nv13, oldnode2,newnode2},
        {oldnode2, newnode2} = {oldnode, newnode};
        v30 = First[newnode];
        {v31, v32} = elist2[[v30 - n, All, 1]];
        v33 = al2[[v30]];
        {nv31,nv32,nv33} = {v31,v32,v33}/.rule21;
        v10 = First[oldnode];
        {v11, v12} = elist1[[v10 - n, All, 1]];
        v13 = al1[[v10]];
        {nv11,nv12,nv13} = {v11,v12,v13}/.rule12;
        Switch[al1[[#]]==nv33&/@ {nv31, nv32},
         {False, True},
         oldnode2 = Union[oldnode2, {al1[[nv31]]}],
         {True, False},
         oldnode2 = Union[oldnode2, {al1[[nv32]]}]
         ];
        Switch[al2[[#]]==nv13&/@ {nv11, nv12},
         {False, True},
         newnode2 = Union[newnode2, {al2[[nv11]]}],
         {True, False},
         newnode2 = Union[newnode2, {al2[[nv12]]}]
         ];
        If[ (al2[[{nv11, nv12}]] /. rule21) == al1[[{nv31, nv32}]] &&Equal @@ al1[[{nv31, nv32}]],
            oldnode2 = Union[oldnode2, {al1[[nv31]]}];
            newnode2 = Union[newnode2, {al2[[nv11]]}];
        ];
        oldnode2 = Union[oldnode2, Flatten[Position[tlist1, #] & /@Complement[tlist1, Complement[tlist2, tlist2[[newnode2]]]]]];
        newnode2 = Union[newnode2,Flatten[Position[tlist2, #] & /@Complement[tlist2, Complement[tlist1, tlist1[[oldnode2]]]]]];
        {oldnode2, newnode2}
    ]
    
TreeTransitionType1[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {tlist1, tlist2, elist1, elist2, n, t10, v10, v11, v12, v13,t30, v30, v31, v32, 
        v33, share, x, set, al1, al2, cc12, cc21,pos},
        If[ SameCoalescentTreeQ[tree1, tree2],
            Return[{0, 0}]
        ];
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        n = LeafNumber[tree1];
        set = Complement[tlist1, tlist2];
        Switch[Length[set],
             0,
             elist1 = Sort[#] & /@ elist1;
             elist2 = Sort[#] & /@ elist2;
             cc12 = SortBy[Complement[elist1, elist2], #[[1, 2]] &];
             cc21 = SortBy[Complement[elist2, elist1], #[[1, 2]] &];
             Switch[{Length[cc12], Length[cc21]},
              {2, 2},
              If[ al1[[cc12[[1, 1, 2]]]] === cc12[[2, 1, 2]] && 
                al2[[cc21[[1, 1, 2]]]] === cc21[[2, 1, 2]],
                  t10 = t30 = tlist1[[cc12[[1, 1, 2]]]],
                  Return[{-1, -1}]
              ],
              {3, 3},
              pos = Flatten[Position[MapThread[SameQ, {al1[[cc12[[All, 1, 2]]]], al2[[cc21[[All, 1, 2]]]]}],False]];
              If[ Length[pos] == 1,
                  t10 = t30 =  tlist1[[cc12[[First[pos], 1, 2]]]],
                  Return[{-1, -1}]
              ],
              _, Return[{-1, -1}]
              ],
             1,
             t10 = First[set];
             t30 = First[Complement[tlist2, tlist1]],
             _, Return[{-1, -1}]
         ];
        v10 = Position[tlist1, t10][[1, 1]];
        {v11, v12} = elist1[[v10 - n, All, 1]];
        v13 = al1[[v10]];
        v30 = Position[tlist2, t30][[1, 1]];
        {v31, v32} = elist2[[v30 - n, All, 1]];
        v33 = al2[[v30]];
        {v31, v32, v33} = {v31, v32, v33} /. {(x_ /; x > n) :> Position[tlist1, tlist2[[x]]][[1, 1]]};
        share = Intersection[{v11, v12}, {v31, v32}];
        Switch[Length[share],
              2,
              {1, If[ t30 < t10,
                      {{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, t30, {v11, v10}}},
                      {{{v11, v10}, t30, {v10, v13}}, {{v12, v10}, t30, {v10, v13}}}
                  ]
                },
              1,
              x = First[Complement[{v31, v32}, share]];
              If[ {x, v33} == {2 n - 1, -1} || 
                MemberQ[Flatten[elist1, 1], {x, v33}],
                  {1, {{{First[share], v10}, t30, {x, v33}}}},
                  {-1, -1}
              ],
              0, {-1, -1}
          ]
    ]        
         
TreeTransitionType1Robust[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {al1, al2, n,oldnode,newnode,eleft1,eleft2,diff, tlist1, elist1, tnew, cc12, rule21, tlist2, elist2,
       told, cc21, rule12, opp},
        If[ SameCoalescentTreeQ[tree1, tree2],
            Return[{0, 0}]
        ];
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        AppendTo[tlist1, Infinity];
        AppendTo[tlist2, Infinity];
        elist1 = Sort[#] & /@ elist1;
        elist2 = Sort[#] & /@ elist2;
        n = LeafNumber[tree1];
        elist1 = Append[elist1, {{-2, -1}, {2 n - 1, -1}}];
        elist2 = Append[elist2, {{-2, -1}, {2 n - 1, -1}}];
        tnew = Complement[tlist2, tlist1];
        told = Complement[tlist1, tlist2];
        If[ Length[tnew] > 2 || Length[told] > 2,
            {-1, -1},
            oldnode = Flatten[Position[tlist1, #] & /@ told];
            newnode = Flatten[Position[tlist2, #] & /@ tnew];
            eleft1 = Delete[Transpose[{Range[2 n - 1], Most[tlist1]}], Partition[oldnode, 1]];
            eleft2 = Delete[Transpose[{Range[2 n - 1], Most[tlist2]}], Partition[newnode, 1]];
            diff =  Select[Transpose[{eleft1[[n + 1 ;;, 1]], eleft2[[n + 1 ;;, 1]]}],Unequal @@ # &];
            rule12 = Join[Thread[oldnode -> Take[CharacterRange["A", "Z"], Length[oldnode]]],Thread[diff[[All, 1]] -> diff[[All, 2]]]];
            rule21 = Join[Thread[newnode -> Take[CharacterRange["A", "Z"], Length[newnode]]], Thread[diff[[All, 2]] -> diff[[All, 1]]]];
            cc12 = SortBy[Complement[elist1,Replace[elist2, rule21, {3}]], tlist1[[#[[1, 2]]]] &];
            cc21 = SortBy[Complement[elist2, Replace[elist1, rule12, {3}]], tlist2[[#[[1, 2]]]] &];
            opp = gettreeoperator1[{tlist1, cc12, rule21, al1}, {tlist2, cc21, rule12, al2}];
            If[ opp === -1,
                {-1, -1},
                {1, opp}
            ]
        ]
    ]         

TreeTransitionType2[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {al1, al2, eleft1,eleft2,diff, tlist1, elist1, tnew, cc12, rule21, tlist2, elist2,
       told, cc21, rule12, opp, oldnode, newnode, n, newchild, newparent,
       newlabelrule, newnodelabel, oldchild, oldparent, old,
       new, oldQ, newQ},
        If[ SameCoalescentTreeQ[tree1, tree2],
            Return[{0, 0}]
        ];
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        AppendTo[tlist1, Infinity];
        AppendTo[tlist2, Infinity];
        elist1 = Sort[#] & /@ elist1;
        elist2 = Sort[#] & /@ elist2;
        n = LeafNumber[tree1];
        elist1 = Append[elist1, {{-2, -1}, {2 n - 1, -1}}];
        elist2 = Append[elist2, {{-2, -1}, {2 n - 1, -1}}];
        tnew = Complement[tlist2, tlist1];
        told = Complement[tlist1, tlist2];
        If[ Length[tnew] > 2 || Length[told] > 2,
            Return[{-1, -1}]
        ];
        oldnode = Flatten[Position[tlist1, #] & /@ told];
        newnode = Flatten[Position[tlist2, #] & /@ tnew];
        If[ Length[newnode]<2&&Length[oldnode]<2,
            eleft1 = Delete[Transpose[{Range[2 n - 1], Most[tlist1]}], Partition[oldnode, 1]];
            eleft2 = Delete[Transpose[{Range[2 n - 1], Most[tlist2]}], Partition[newnode, 1]];
            diff =  Select[Transpose[{eleft1[[n + 1 ;;, 1]], eleft2[[n + 1 ;;, 1]]}],Unequal @@ # &];
            rule12 = Join[Thread[oldnode -> Take[CharacterRange["A", "Z"], Length[oldnode]]],Thread[diff[[All, 1]] -> diff[[All, 2]]]];
            rule21 = Join[Thread[newnode -> Take[CharacterRange["A", "Z"], Length[newnode]]], Thread[diff[[All, 2]] -> diff[[All, 1]]]];
            cc12 = SortBy[Complement[elist1,Replace[elist2, rule21, {3}]], tlist1[[#[[1, 2]]]] &];
            cc21 = SortBy[Complement[elist2, Replace[elist1, rule12, {3}]], tlist2[[#[[1, 2]]]] &];
            opp = gettreeoperator1[{tlist1, cc12, rule21, al1}, {tlist2, cc21,rule12, al2}];
            If[ opp =!= -1,
                Return[{1, {opp}}]
            ];
            {oldnode, newnode} = nodediffI[{tlist1, elist1,cc12, rule21, al1}, {tlist2,elist2, cc21, rule12, al2}];
            If[ oldnode === newnode === {},
                oldnode = nodediffII[{tlist1, elist1, cc12, rule21, al1}, {tlist2, elist2, cc21, rule12, al2}];
                newnode = nodediffII[{tlist2, elist2, cc21, rule12, al2}, {tlist1, elist1, cc12, rule21, al1}];
                oldnode = Union[oldnode,Flatten[Position[tlist1, #]&/@Complement[tlist1, Complement[tlist2, tlist2[[newnode]]]]]];
                newnode = Union[newnode,Flatten[Position[tlist2, #]&/@Complement[tlist2, Complement[tlist1, tlist1[[oldnode]]]]]];
            ];
            If[ Length[oldnode]==1,
                {oldnode, newnode} = nodediffIII[oldnode, newnode, n, tlist1,elist1, al1,rule12, tlist2, elist2,al2, rule21];
            ];
            If[ Length[newnode] > 2 || Length[oldnode] > 2 ||Length[newnode] != Length[oldnode],
                Return[{-1, -1}]
            ];
            tnew = tlist2[[newnode]];
            told = tlist1[[oldnode]];
        ];
        newchild = elist2[[newnode - n, All, 1]];
        newparent = al2[[newnode]];
        newlabelrule = Thread[newnode -> Take[CharacterRange["A", "Z"], Length[newnode]]];
        newnodelabel = newnode /. newlabelrule;
        {newparent, newchild} = {newparent, newchild} /.newlabelrule /. {(x_ /; x > n) :>Position[tlist1, tlist2[[x]]][[1, 1]]};
        oldchild = elist1[[oldnode - n, All, 1]];
        oldparent = al1[[oldnode]];
        If[ Length[oldchild] == 2 && MemberQ[oldchild[[2]], oldnode[[1]]] &&oldchild[[2, 1]] == oldnode[[1]],
            oldchild[[2]] = Reverse[oldchild[[2]]]
        ];
        If[ Length[newchild] == 2 &&MemberQ[newchild[[2]], newnodelabel[[1]]] &&newchild[[2, 1]] == newnodelabel[[1]],
            newchild[[2]] = Reverse[newchild[[2]]]
        ];
        old = Transpose[{told, oldchild, oldparent, oldnode, oldnode}];
        new = Transpose[{tnew, newchild, newparent, newnodelabel, newnode}];
        If[ Length[old] == Length[new] == 2,
            oldQ = old[[2, 2, 2]] == old[[1, -2]];
            newQ = new[[2, 2, 2]] == new[[1, -2]];
            opp = Which[
             oldQ && newQ, midtreecase2211[old, new, al1,tree1],
             (! oldQ) && newQ, midtreecase2201[old, new, al1],
             oldQ && (! newQ), midtreecase2210[old, new, al1],
             (! oldQ) && (! newQ), midtreecase2200[old, new, al1,tree1]
             ];
            If[ opp===-1,
                {-1,-1},
                {2,opp}
            ],
            Print["Case to be done!"];
            Print["length of cc12 and cc21:",{Length[cc12],Length[cc21]}];
            {-1, -1}
        ]
    ]
    
TreeTransitionType1C[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {x, tlist1, elist1, tlist2, elist2, al1, al2, n, told, t10, 
      v10, v11, v12, v13, t30, v30, v31, v32, v33, share, cc12, cc21, 
      v20, v21, b1, b2, v40, v50, v60, v51, v61, v22, set, set1, 
      set2, temp1, temp2, newt},
        If[ SameCoalescentTreeQ[tree1, tree2],
            Return[{0, 0}]
        ];
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        n = LeafNumber[tree1];
        told = Complement[tlist1, tlist2];
        Switch[Length[told],
         1,
         t10 = First[told];
         v10 = Position[tlist1, t10][[1, 1]];
         {v11, v12} = elist1[[v10 - n, All, 1]];
         v13 = al1[[v10]];
         t30 = First[Complement[tlist2, tlist1]];
         v30 = Position[tlist2, t30][[1, 1]];
         {v31, v32} = elist2[[v30 - n, All, 1]];
         v33 = al2[[v30]];
         {v31, v32, 
           v33} = {v31, v32, 
            v33} /. {(x_ /; x > n) :> 
             Position[tlist1, tlist2[[x]]][[1, 1]]};
         share = Intersection[{v11, v12}, {v31, v32}];
         Switch[Length[share],
          2,
          {1, If[ t30 < t10,
                  {{{v11, v10}, t30, {v12, v10}}, {{v12, v10}, 
                  t30, {v11, v10}}},
                  {{{v11, v10}, 
                  t30, {v10, v13}}, {{v12, v10}, t30, {v10, v13}}}
              ]},
          1,
          x = First[Complement[{v31, v32}, share]];
          If[ al1[[x]] == 
            v33,
              {1, {{{First[share], v10}, t30, {x, v33}}}},
              {-1, -1}
          ],
          0,
          {-1, -1}
          ],
         0,
         AppendTo[tlist1, Infinity];
         AppendTo[tlist2, Infinity];
         elist1 = Sort[#] & /@ elist1;
         elist2 = Sort[#] & /@ elist2;
         n = LeafNumber[tree1];
         elist1 = Append[elist1, {{-2, -1}, {2 n - 1, -1}}];
         elist2 = Append[elist2, {{-2, -1}, {2 n - 1, -1}}];
         cc12 = SortBy[Complement[elist1, elist2], tlist1[[#[[1, 2]]]] &];
         cc21 = SortBy[Complement[elist2, elist1], tlist2[[#[[1, 2]]]] &];
         Switch[{Length[cc12], Length[cc21]},
          {2, 2},
          If[ al1[[cc12[[1, 1, 2]]]] == cc12[[2, 1, 2]] == cc21[[2, 1, 2]] &&
             al2[[cc21[[1, 1, 2]]]] == cc21[[2, 1, 2]],
              {v10, v20} = cc12[[All, 1, 2]];
              v21 = First[Complement[cc12[[2, All, 1]], {v10}]];
              If[ MemberQ[cc21[[1, All, 1]], v21],
                  v11 = First[Complement[cc21[[1, All, 1]], {v21}]];
                  {1, {{{v11, v10}, tlist1[[v10]], {v21, v20}}}},
                  {-1, -1}
              ],
              {-1, -1}
          ],
          {3, 3},
          b1 = {Most[al1[[cc12[[All, 1, 2]]]]] == Rest[cc12[[All, 1, 2]]], 
            cc12[[;; 2, 1, 2]] === cc12[[-1, All, 1]]};
          b2 = {Most[al2[[cc21[[All, 1, 2]]]]] == Rest[cc21[[All, 1, 2]]], 
            cc21[[;; 2, 1, 2]] === cc21[[-1, All, 1]]};
          Switch[{b1, b2},
           {{True, False}, {False, True}},
           If[ cc12[[All, 1, 2]] == cc21[[All, 1, 2]],
               {v10, v20, v30} = cc12[[All, 1, 2]];
               v21 = First[Complement[cc12[[2, All, 1]], {v10}]];
               v31 = First[Complement[cc12[[3, All, 1]], {v20}]];
               If[ MemberQ[cc21[[2, All, 1]], v21] && 
                 MemberQ[cc21[[1, All, 1]], v31],
                   {v11, v31} = 
                    First[Intersection[cc21[[1, All, 1]], #]] & /@ 
                     cc12[[{1, 3}, All, 1]];
                   {{{v11, v10}, tlist1[[cc12[[1, 1, 2]]]], {v31, v30}}},
                   {-1, -1}
               ],
               {-1, -1}
           ],
           {{False, True}, {True, False}},
           If[ cc12[[All, 1, 2]] == cc21[[All, 1, 2]],
               {v40, v50, v60} = cc21[[All, 1, 2]];
               v51 = First[Complement[cc21[[2, All, 1]], {v40}]];
               v61 = First[Complement[cc21[[3, All, 1]], {v50}]];
               If[ MemberQ[cc12[[2, All, 1]], v51] && 
                 MemberQ[cc12[[1, All, 1]], v61],
                   {v11, v22} = 
                    First[Intersection[#, cc21[[1, All, 1]]]] & /@ 
                     cc12[[;; 2, All, 1]];
                   {{{v11, v40}, tlist1[[cc12[[1, 1, 2]]]], {v22, v50}}},
                   {-1, -1}
               ],
               {-1, -1}
           ],
           {{False, False}, {False, False}},
           set = Drop[Range[3], {#}] & /@ Reverse[Range[3]];
           set1 = 
            Select[set, 
             al1[[cc12[[First[#], 1, 2]]]] == cc12[[Last[#], 1, 2]] &];
           set1 = {#, Complement[Range[3], #]} & /@ set1;
           set2 = 
            Select[set, 
             al2[[cc21[[First[#], 1, 2]]]] == cc21[[Last[#], 1, 2]] &];
           set2 = {#, Complement[Range[3], #]} & /@ set2;
           If[ Length[set1] == Length[set2] == 1 && 
             cc12[[set1[[1, All, -1]], 1, 2]] == 
              Reverse[cc21[[set2[[1, All, -1]], 1, 2]]],
               v10 = cc12[[set1[[1, 1, 1]], 1, 2]];
               v20 = cc12[[set1[[1, 2, 1]], 1, 2]];
               temp1 = 
                Outer[Intersection, 
                 cc12[[{set1[[1, 1, 1]], set1[[1, 2, 1]]}, All, 1]], 
                 cc21[[{set2[[1, 1, 1]], set2[[1, 2, 1]]}, All, 1]], 1, 1];
               temp2 = 
                Intersection[cc12[[set1[[1, 1, 2]], All, 1]], 
                 cc21[[set2[[1, 2, 1]], All, 1]]];
               If[ MatchQ[temp1, {{{_}, {_}}, {{_}, {}}}] && Length[temp2] == 1,
                   {v11, v22} = Flatten[temp1[[All, 1]]];
                   newt = tlist2[[cc21[[set2[[1, 1, 1]], 1, 2]]]];
                   {1, {{{v11, v10}, newt, {v22, v20}}}},
                   {-1, -1}
               ],
               {-1, -1}
           ],
           _,
           {-1, -1}
           ]
          ],
         _,
         {-1, -1}
         ]
    ]

midtreetest[{tree1_, tree12_, tree2_}] :=
    Module[ {type, type0, tlist1, midtrees, pos1, pos2, pos3, pos4, pos, 
      include, op,i,j},
        type = type0 = TreeTransitionTypeII2[tree1, tree2];
        Switch[ First[type],
            0|1, Return[True],
            -1,Return[False]            
        ];
        If[ Count[Flatten[type[[2, All, All, 2]]], -1] > 0,
            tlist1 = tree1[[1]];
            type[[-1]] = 
             Map[If[ #[[2]] == -1,
                     ReplacePart[#, 
                      2 -> RandomReal[{Max[tlist1[[{#[[1, 1]], #[[3, 1]]}]]], 
                         tlist1[[#[[3, 2]]]]}]],
                     #
                 ] &, type[[-1]], {2}];
        ];
        midtrees = Map[NextLocalTree[tree1, #] &, Last[type], {2}];
        pos1 = Flatten[
          Position[
           Table[Switch[Length[midtrees[[i]]], 1, True, 2, 
             SameCoalescentTreeQ @@ midtrees[[i]], 3, False], {i, 
             Length[midtrees]}], False]];
        If[ pos1 =!= {},
            Print["Midtrees from ith transformations  (i\[Element]" <> 
              ToString[pos1] <> ") are different!" <> "\n transformations=", 
             type, "{tree1,tree2}=", {tree1, tree2}]
        ];
        midtrees = midtrees[[All, 1]];
        pos2 = Select[
           Flatten[Table[{{i, j}, 
              SameCoalescentTreeQ[midtrees[[i]], midtrees[[j]]]}, {i, 
              Length[midtrees] - 1}, {j, i + 1, Length[midtrees]}], 1], 
           Last[#] &][[All, 1]];
        If[ pos2 =!= {},
            Print["Midtrees from ith and jth transformations ({i,j}\[Element]" \
            <> ToString[pos2] <> ") are same!" <> "\n transformations=", type, 
             "{tree1,tree2}=", {tree1, tree2}]
        ];
        pos3 = Flatten[
          Position[(TreeTransitionTypeII2[tree1, #] & /@ midtrees)[[All, 
             1]], _?(# != 1 &)]];
        If[ pos3 =!= {},
            Print["The distance between tree1 and midtree from ith tansformation (i\[Element]" <> ToString[pos3] <> ") \[NotEqual]1!" <>
               "\n transformations=", type, "{tree1,tree2}=", {tree1, tree2}]
        ];
        pos4 = Flatten[
          Position[(TreeTransitionTypeII2[#, tree2] & /@ midtrees)[[All, 
             1]], _?(# != 1 &)]];
        If[ pos4 =!= {},
            Print["The distance between midtree and tree2 from ith transformation (i\[Element]" <> ToString[pos4] <> ") \[NotEqual]1!" <>
               "\n transformations=", type, "{tree1,tree2}=", {tree1, tree2}]
        ];
        op = Last[TreeTransitionType1Robust[tree1, tree12]];
        include = If[ Complement[op, Flatten[Last[type0], 1]] === {},
                      True,
                      pos = 
                       Flatten[Position[type0[[-1, All, All, {1, -1}]], 
                         op[[All, {1, -1}]]]];
                      pos =!= {} && type0[[-1, First[pos], 1, 2]] === -1
                  ];
        If[ ! include,
            Print["tree12 is not included in the midtrees!\n \
{tree1,tree12,tree2}=", {tree1, tree12, tree2}]
        ];
        pos1 === pos2 === pos3 === pos4==={} && include
    ]
  

findoldnode[tree1_, tree2_, bool_: True, isprint_: False] :=
    Module[ {temp, tlist1, elist1, tlist2, elist2, n, al1, told, 
      tnew, newnode, newlabelrule, rule2to1, elist2to1, res, c12, remain,
       pos, oldnode},
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        n = LeafNumber[tree1];
        elist1 = Append[Sort[#] & /@ elist1, {{-2, -1}, {2 n - 1, -1}}];
        elist2 = Append[Sort[#] & /@ elist2, {{-2, -1}, {2 n - 1, -1}}];
        told = Complement[tlist1, tlist2];
        tnew = Complement[tlist2, tlist1];
        al1 = TreeToAdjacencyList[tree1];
        oldnode = If[ bool,
                      Flatten[Position[tlist1, #] & /@ told],
                      {}
                  ];
        newnode = Flatten[Position[tlist2, #] & /@ tnew];
        newlabelrule = Thread[newnode -> Take[{"E", "F"}, Length[newnode]]];
        temp = Complement[Range[n + 1, 2 n - 1], newnode];
        rule2to1 = 
         Thread[temp -> (Position[tlist1, #][[1, 1]] & /@ tlist2[[temp]])];
        rule2to1 = Join[newlabelrule, rule2to1];
        elist2to1 = elist2 /. rule2to1;
        c12 = SortBy[Complement[elist1, elist2to1], #[[1, 2]] &];
        res = {al1, c12, rule2to1};
        c12 = Select[
          c12, ! (MemberQ[oldnode, #[[1, 2]]] || #[[1, 2]] == -1) &];
        remain = Flatten[elist2to1, 1];
        remain = Append[remain, {remain[[-1, -1]], -1}];
        remain = remain //. (Rule @@ # & /@ 
             Select[remain, (StringQ[First[#]] &)]);
        While[c12=!={},
         If[ isprint,
             Print["oldnode=", oldnode, ";c12=", c12]
         ];
         temp = NestWhile[al1[[#]] &, #, (MemberQ[oldnode, #] &)] & /@al1[[c12[[All, 1, 2]]]];
         temp = Thread[#] & /@ Transpose[{Transpose[Append[Transpose[c12[[All, All, 1]]], c12[[All, 1, 2]]]], temp}];
         If[ isprint,
             Print["branches around node=", temp]
         ];
         temp = Map[MemberQ[remain, #] &, temp, {2}];
         If[ isprint,
             Print["branches around node=", temp]
         ];
         temp = (#[[1]] || #[[2]]) && (! #[[3]]) & /@ temp;
         If[ isprint,
             Print["branches around node=", temp]
         ];
         pos = Flatten[Position[temp, True, {1}, Heads -> False]];
         If[ pos === {},
             Break[],
             oldnode = Union[oldnode, c12[[pos, 1, 2]]];
             c12 = Select[c12, ! MemberQ[oldnode, #[[1, 2]]] &]
         ];
         ];
        Prepend[res, oldnode]
    ]
   
finddiffnode[tree1_, tree2_, al1_, al2_, bool_: True, 
   isprint_: False] :=
    Module[ {temp, tlist1, elist1, tlist2, elist2, n, told, tnew, 
      newnode, newlabelrule, rule2to1, invrule2to1, elist2to1, res, c12,
       remain, pos, oldnode,cond,set},
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        n = LeafNumber[tree1];
        elist1 = Append[Sort[#] & /@ elist1, {{-2, -1}, {2 n - 1, -1}}];
        elist2 = Append[Sort[#] & /@ elist2, {{-2, -1}, {2 n - 1, -1}}];
        told = Complement[tlist1, tlist2];
        tnew = Complement[tlist2, tlist1];
        oldnode = If[ bool,
                      Flatten[Position[tlist1, #] & /@ told],
                      {}
                  ];
        newnode = Flatten[Position[tlist2, #] & /@ tnew];
        newlabelrule = Thread[newnode -> Take[{"E", "F"}, Length[newnode]]];
        temp = Complement[Range[n + 1, 2 n - 1], newnode];
        rule2to1 = 
         Thread[temp -> (Position[tlist1, #][[1, 1]] & /@ tlist2[[temp]])];
        rule2to1 = Join[newlabelrule, rule2to1];
        invrule2to1 = Reverse[#] & /@ rule2to1;
        elist2to1 = elist2 /. rule2to1;
        c12 = SortBy[Complement[elist1, elist2to1], #[[1, 2]] &];
        res = {c12, rule2to1};
        c12 = Select[
          c12, ! (MemberQ[oldnode, #[[1, 2]]] || #[[1, 2]] == -1) &];
        remain = Flatten[elist2to1, 1];
        remain = Append[remain, {remain[[-1, -1]], -1}];
        remain = 
         remain //. (Rule @@ # & /@ Select[remain, (StringQ[First[#]] &)]);
        While[c12 =!= {},
         If[ isprint,
             Print["oldnode=", oldnode, "newnode=", newnode, ";c12=", c12]
         ];
         temp = NestWhile[al1[[#]] &, #, (MemberQ[oldnode, #] &)] & /@ 
           al1[[c12[[All, 1, 2]]]];
         temp = 
          Thread[#] & /@ 
           Transpose[{Transpose[
              Append[Transpose[c12[[All, All, 1]]], c12[[All, 1, 2]]]], 
             temp}];
         If[ isprint,
             Print["branches around node=", temp]
         ];
         cond = Map[MemberQ[remain, #] &, temp, {2}];
         If[ isprint,
             Print["branches around node=", cond]
         ];
         pos = (Xor[#[[1]], #[[2]]]) && (! #[[3]]) & /@ cond;
         If[ isprint,
             Print["branches around node=", pos]
         ];
         pos = Flatten[Position[pos, True, {1}, Heads -> False]];
         If[ pos === {},
             Break[],
             set = 
              Complement[
                Flatten[Pick[c12[[pos, All, 1]], cond[[pos, ;; 2]], False]], 
                oldnode] /. invrule2to1;
             newnode = Union[newnode, al2[[set]]];
             oldnode = Union[oldnode, c12[[pos, 1, 2]]];
             c12 = Select[c12, ! MemberQ[oldnode, #[[1, 2]]] &]
         ];
         ];
        Join[{oldnode, newnode}, res]
    ];
      
comparetree[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {tlist1, tlist2, elist1, elist2, n, told, tnew, oldnode, al1, 
      cc12, rule2to1, newnode, al2, cc21, rule1to2, oldchild, oldparent, 
      old, newchild, newparent, newnodelabel, new,cross,newlabelrule,newnode2,oldnode2},
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        elist1 = Sort[#] & /@ elist1;
        elist2 = Sort[#] & /@ elist2;
        If[ tlist1 == tlist2 && elist1 == elist2,
            Return[0]
        ];
        n = LeafNumber[tree1];
        told = Complement[tlist1, tlist2];
        tnew = Complement[tlist2, tlist1];
        If[ Length[tnew] > 2 || Length[told] > 2,
            (*Print["case to be excluded!"];*)
            Return[-1]
        ];
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        {oldnode, newnode, cc12, rule2to1} = finddiffnode[tree1, tree2, al1,al2];
        {newnode2, oldnode2, cc21, rule1to2} = finddiffnode[tree2, tree1,al2,al1];
        oldnode = Union[oldnode,oldnode2];
        newnode = Union[newnode,newnode2];
        oldnode = Union[oldnode, Flatten[Position[tlist1, #] & /@Complement[tlist1, Complement[tlist2, tlist2[[newnode]]]]]];
        newnode = Union[newnode, Flatten[Position[tlist2, #] & /@Complement[tlist2, Complement[tlist1, tlist1[[oldnode]]]]]];
        If[ Length[newnode] > 2 || Length[oldnode] > 2 ||Length[newnode] != Length[oldnode],
            Return[-1]
        ];
        Which[
             Length[newnode] == 0&&Length[cc12] == 2,
             cross = Outer[Intersection, cc12[[All,All,1]],cc21[[All,All,1]], 1, 1];
             Which[
                 MatchQ[cross, {{{_}, {_}}, {{_}, {_}}}], 
                 oldnode = newnode = 
                     If[ MemberQ[cc12[[2, All, 1]], cc12[[1, 1, 2]]] && 
                         MemberQ[cc21[[2, All, 1]], cc21[[1, 1, 2]]],
                         {cc12[[1,1,2]]},
                         cc12[[All,1,2]]
                     ],
                 True,
                  Print["case00 to be done!"];
                  Return[-1]
             ],    
             (*Length[cc12] == 3,
             Which[
                 al1[[Most[cc12[[All, 1, 2]]]]] === Rest[cc12[[All, 1, 2]]]&&
                 al2[[Most[cc21[[All, 1, 2]]]]] === Rest[cc21[[All, 1, 2]]],
                 Print["cc12=",cc12,";cc21=",cc21];
                 oldnode=Most[cc12[[All, 1, 2]]];
                 newnode=Most[cc21[[All, 1, 2]]],
                 cc12[[;; 2, 1, 2]] === cc12[[-1, All, 1]]&&
                 cc21[[;; 2, 1, 2]] === cc21[[-1, All, 1]],
                 Print["cc12=",cc12,";cc21=",cc21];
                 oldnode=Most[cc12[[All, 1, 2]]];
                 newnode=Most[cc21[[All, 1, 2]]]                 
             ],*)                      
             Length[newnode] == 0,
             Print["case00 to be done!"];
             Return[-1]            
            ];
        If[ Length[newnode] ==0,
            old = new = {},
            tnew = tlist2[[newnode]];
            newchild = elist2[[newnode - n, All, 1]];
            newparent = al2[[newnode]];
            newlabelrule = Thread[newnode -> Take[{"E", "F"}, Length[newnode]]];
            newnodelabel = newnode /. newlabelrule;
            {newparent, newchild} = {newparent, newchild} /.newlabelrule /. {(x_ /; x > n) :>Position[tlist1, tlist2[[x]]][[1, 1]]};
            told = tlist1[[oldnode]];
            oldchild = elist1[[oldnode - n, All, 1]];
            oldparent = al1[[oldnode]];
            If[ Length[oldchild] == 2 && MemberQ[oldchild[[2]], oldnode[[1]]] && oldchild[[2, 1]] == oldnode[[1]],
                oldchild[[2]] = Reverse[oldchild[[2]]]
            ];
            If[ Length[newchild] == 2 && MemberQ[newchild[[2]], newnodelabel[[1]]] &&newchild[[2, 1]] == newnodelabel[[1]],
                newchild[[2]] = Reverse[newchild[[2]]]
            ];
            old = Transpose[{told, oldchild, oldparent, oldnode, oldnode}];
            new = Transpose[{tnew, newchild, newparent, newnodelabel, newnode}];
        ];
        {Length[new], {old, new}, {al1, al2}}
    ]   
   
    
comparetree2[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {tlist1, tlist2, elist1, elist2, n, al1,al2, tnew, newnode, 
      newlabel, newchild, newparent, told, oldnode, oldchild, oldparent, 
      old, new, newlabelrule,cc12,cc21,groupid,pos},
        {tlist1, elist1} = List @@ tree1;
        {tlist2, elist2} = List @@ tree2;
        tnew = Sort[Complement[tlist2, tlist1]];
        told = Sort[Complement[tlist1, tlist2]];
        If[ Length[tnew] > 2 || Length[told] > 2,
            Return[{-1,{{-1,-1},{-1,-1}}}]
        ];
        If[ Length[tnew]==Length[told],
            groupid = Length[tnew],
            Print["Case in compareTree to be done!"];
            Return[$Failed]
        ];
        n = LeafNumber[tree1];
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        If[ told === tnew === {},
            cc12 = SortBy[Complement[elist1, elist2], #[[1, 2]] &];
            cc21 = SortBy[Complement[elist2, elist1], #[[1, 2]] &];
            If[ cc12[[All, 1, 2]] != cc21[[All, 1, 2]],
                Return[{-1, {{-1, -1}, {-1, -1}}}]
            ];
            Which[
             Length[cc12] == Length[cc21] == 2,
             If[ al1[[cc12[[1, 1, 2]]]] === cc12[[2, 1, 2]],
                 tnew = told = tlist1[[cc12[[All, 1, 2]]]],
                 Print["specisal cases d=0."];
                 Return[{-1, {{-1, -1}, {-1, -1}}}]
             ],
             Length[cc12] == Length[cc21] == 3, 
             pos = Flatten[Position[MapThread[SameQ, {al1[[cc12[[All, 1, 2]]]], al2[[cc21[[All, 1, 2]]]]}],False]];
             If[ Length[pos] == 1,
                 groupid = 1;
                 tnew = told = tlist1[[cc12[[pos, 1, 2]]]],
                 Print["specisal cases d=0."];
                 Return[{-1, {{-1, -1}, {-1, -1}}}]
             ],
             True,
             Print["specisal cases d=0."];
             Return[{-1, {{-1, -1}, {-1, -1}}}]             
            ]
        ];
        newnode = Flatten[Position[tlist2, #] & /@ tnew];
        newlabelrule = Thread[newnode -> Take[{"E", "F"}, Length[newnode]]];
        newlabel = newnode /. newlabelrule;
        newchild = elist2[[newnode - n, All, 1]] /. newlabelrule;
        newparent = al2[[newnode]] /. newlabelrule;
        {newchild, newparent} = {newchild, newparent} /. {(x_ /; x > n) :>Position[tlist1, tlist2[[x]]][[1, 1]]};
        oldnode = Flatten[Position[tlist1, #] & /@ told];
        oldchild = elist1[[oldnode - n, All, 1]];
        oldparent = al1[[oldnode]];
        If[ Length[oldchild] == 2 && MemberQ[oldchild[[2]], oldnode[[1]]] &&  oldchild[[2, 1]] == oldnode[[1]],
            oldchild[[2]] = Reverse[oldchild[[2]]]
        ];
        If[ Length[newchild] == 2 && MemberQ[newchild[[2]], newlabel[[1]]] && newchild[[2, 1]] == newlabel[[1]],
            newchild[[2]] = Reverse[newchild[[2]]]
        ];
        old = Transpose[{told, oldchild, oldparent, oldnode, oldnode}];
        new = Transpose[{tnew, newchild, newparent, newlabel, newnode}];
        {groupid,{old, new}, {al1,al2}}
    ]
  

firstTransition::wrongshare = "Wrong share = `1`";
firstTransition[old : {_, {_, _}, _, _}, new : {_, {_, _}, _}] :=
    Module[ {v11, v12, v10, v13, t10, v31, v32, v33, t30, share},
        {t10, {v11, v12}, v13, v10} = old;
        {t30, {v31, v32}, v33} = new;
        share = Intersection[{v11, v12}, {v31, v32}];
        Switch[Length[share],
         1, {{{First[share],v10}, t30, {First[Complement[{v31, v32}, share]], v33}}},
         2,
         If[ t30 < t10,
             {{{v11,v10},t30,{v12,v10}},{{v12,v10},t30,{v11,v10}}},
             {{{v11,v10},t30,{v10,v13}},{{v12,v10},t30,{v10,v13}}}
         ]
        ]
    ]

midtreecase002[old_, new_] :=
    Module[ {t10, v11, v12, v13, v10, t20, v21, v22, v23, v20, t30, v31, v32, v33, t40, v41, v42, v43},
        {{t10, {v11, v12}, v13, v10}, {t20, {v21, v22}, v23, v20}} = old[[All, ;; 4]];
        {{t30, {v31, v32}, v33}, {t40, {v41, v42}, v43}} = new[[All, ;; 3]];
        If[ {t10,t20}=={t30,t40}&&v23==v43&&MemberQ[{v31,v32},v21]&&Complement[{v41,v31, v32}, {v21,v11, v12}] === {},
            If[ v31 != v21,
                {v31, v32} = {v32, v31}
            ];
            {1,{{{{v32,v10},t10,{v21,v20}}}}},
            {-1,-1}
        ]
    ];                

midtreecase11old[old_, new_, al1_, al2_,tree1_, tree2_] :=
    Module[ {v11, v12, v10, v13, t10, v31, v32, v33, t30, share,x,tlist1, tlist2, nv11, nv12, nv13, y,n},
        {t10, {v11, v12}, v13, v10} = Take[Flatten[old, 1], 4];
        {t30, {v31, v32}, v33} = Take[Flatten[new, 1], 3];
        share = Intersection[{v11,v12},{v31,v32}];
        Switch[Length[share],          
            2, 
            If[ v13==v33,
                {1,If[ t30 < t10,
                       {{{{v11,v10},t30,{v12,v10}},{{v12,v10},t30,{v11,v10}}}},
                       {{{{v11,v10},t30,{v10,v13}},{{v12,v10},t30,{v10,v13}}}}
                   ]},
                -1
            ],
            1,
            If[ !MemberQ[{v11,v12},v31],
                {v31,v32} = {v32,v31}
            ];
            If[ v31!=v11,
                {v11,v12} = {v12,v11}
            ];
            x = al1[[v32]];
            If[ x==v33,
                {1,{{{{v11,v10},t30,{v32,v33}}}}},
                If[ v13==v33,
                    tlist1 = First[tree1];
                    {2,Join[
                       If[ tlist1[[v12]]<tlist1[[x]],
                           {{{{v12,v10},-1,{v32,x}}}},
                           {}
                       ],
                       If[ t30>t10,
                           {{{{v32,x},t30,{v10,v13}}}},
                           {{{{v32,x},t30,{v11,v10}}}}
                       ]]},
                    -1
                ]
            ],
            0,
            tlist1 = First[tree1];
            tlist2 = First[tree2];
            n = LeafNumber[tree1];
            (*nv: label in tree 2; v: label in tree1*)
            {nv11, nv12, nv13} = {v11, v12, v13} /. {(x_ /; x > n) :>Position[tlist2, tlist1[[x]]][[1, 1]]};
            If[ MemberQ[al2[[{nv11,nv12}]],nv13],
                (*set v11 connected to v13 in new tree*)
                If[ al2[[nv11]]!=nv13,
                    {v11,v12} = {v12,v11};
                    {nv11,nv12} = {nv12,nv11}
                ];
                (*-1 denotes any time within the branch*)
                Switch[Count[al1[[{v31, v32}]], v33],
                    1,
                    If[ al1[[v31]]!=v33,
                        {v31,v32} = {v32,v31}
                    ];
                    y = al1[[v32]];
                    (*  If[ y ==v31&&t30>t10,
                            {{{{v32,  y}, t30, {v31, v33}},{{First[Complement[tree1[[2,v31 - n, All, 1]], {v32}]], v31}, t30, {v31, v33}}}},
                            {{{{v32,  y}, t30, {v31, v33}}}}
                        ]*)
                    {2, Join[
                        If[ tlist1[[v12]]<tlist1[[y]],
                            {{{{v12, v10}, -1, {v32, y}}}},
                            {}
                        ],
                        If[ y ==v31&&t30>t10,
                            {{{{v32,  y}, t30, {v31, v33}},{{First[Complement[tree1[[2,v31 - n, All, 1]], {v32}]], v31}, t30, {v31, v33}}}},
                            {{{{v32,  y}, t30, {v31, v33}}}}
                        ],
                        If[ y == v13&&al2[[nv12]]==nv13,
                            {{{{v11, v10}, -1, {v32, y}}}},
                            {}
                        ]                        
                    ]},
                    2,
                    {2, Join[
                        If[ tlist1[[v12]]<tlist1[[v33]],
                            {{{{v12, v10}, -1, {v31, v33}}},{{{v12, v10}, -1, {v32, v33}}}},
                            {}
                        ],
                        {{{{v31,v33},t30,{v32,v33}},{{v32,v33},t30,{v31,v33}}}}]},
                    _, -1
                ],
                -1
            ]           
        ]
    ]

midtreecase11[old_, new_, al1_, al2_,tree1_, tree2_] :=
    Module[ {v11, v12, v10, v13, t10, v31, v32, v33, t30, share,tlist1, tlist2, nv11, nv12, nv13,nv33,nvy, vy,n},
        {t10, {v11, v12}, v13, v10} = Take[Flatten[old, 1], 4];
        {t30, {v31, v32}, v33} = Take[Flatten[new, 1], 3];
        share = Intersection[{v11,v12},{v31,v32}];
        Switch[Length[share],          
            2, 
            If[ v13==v33,
                Which[ 
                    t30 < t10,  {1,{{{{v11,v10},t30,{v12,v10}},{{v12,v10},t30,{v11,v10}}}}},
                    t30>t10,    {1,{{{{v11,v10},t30,{v10,v13}},{{v12,v10},t30,{v10,v13}}}}},
                    True,
                    {-1,-1}
                 ],
                {-1,-1}
            ],
            1,
            If[ !MemberQ[{v11,v12},v31],
                {v31,v32} = {v32,v31}
            ];
            If[ v31!=v11,
                {v11,v12} = {v12,v11}
            ];
            vy = al1[[v32]];
            tlist1 = First[tree1];
            tlist2 = First[tree2];
            n = LeafNumber[tree1];
            If[ vy==v33,
                {1,{{{{v11,v10},t30,{v32,v33}}}}},
                (*nv: label in tree 2; v: label in tree1*)
                {nv12, nvy} = {v12,vy} /. {(x_ /; x > n) :>Position[tlist2, tlist1[[x]]][[1, 1]]};
                If[ v13==v33&&al2[[nv12]]==nvy,
                    {2,Join[{{{{v12,v10},-1,{v32,vy}}}},
                       Which[ 
                        t30>t10, {{{{v32,vy},t30,{v10,v13}}}},
                        t30<t10, {{{{v32,vy},t30,{v11,v10}}}},
                        True, {}
                       ]]},
                    {-1,-1}
                ]
            ],
            0,
            tlist1 = First[tree1];
            tlist2 = First[tree2];
            n = LeafNumber[tree1];
            If[ al1[[v32]]!=v33,
                {v31,v32} = {v32,v31};
            ];
            vy = al1[[v31]];
            (*nv: label in tree 2; v: label in tree1*)
            {nv11, nv12, nv13,nv33,nvy} = {v11, v12, v13,v33,vy} /. {(x_ /; x > n) :>Position[tlist2, tlist1[[x]]][[1, 1]]};
            If[ al2[[nv11]]!=nv13,
                {v11,v12} = {v12,v11};
                {nv11,nv12} = {nv12,nv11}
            ];
            Which[
                al1[[v32]]!=v33,{-1,-1},
                al1[[v31]]==v33&&Complement[al2[[{nv11,nv12}]],{nv33,nv13}]==={},                   
                {2, Join[{{{{v12, v10}, -1, {v31, v33}}},{{{v12, v10}, -1, {v32, v33}}}},
                    {{{{v31,v33},t30,{v32,v33}},{{v32,v33},t30,{v31,v33}}}}]
                 },
                al1[[v31]]==v13&&al2[[nv11]]==nv13&&al2[[nv12]]==nv13,
                {2, {{{{v12, v10}, -1, {v31, v13}}},{{{v11, v10}, -1, {v31, v13}}}, {{{v31,  v13}, t30, {v32, v33}}}}},
                (!MemberQ[{v13,v33},vy])&&Complement[al2[[{nv11,nv12}]],{nvy,nv13}]==={},
                {2, {{{{v12, v10}, -1, {v31, vy}}}, {{{v31,  vy}, t30, {v32, v33}}}}},
                True,
                {-1,-1}                    
            ]                       
        ]
    ]


TreeTransitionType22[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {old, new, al1,al2, oldQ, newQ,temp,groupid},
        If[ SameCoalescentTreeQ[tree1,tree2],
            {0,0},
            {groupid,{old, new}, {al1,al2}} = comparetree2[tree1,tree2];
            Switch[groupid,
             0,Print["case00 to be done!"];
               {-1,-1},             
             1,
             midtreecase11old[old, new, al1, al2,tree1, tree2],
             2,
             oldQ = old[[2, 2, 2]] == old[[1, -2]];
             newQ = new[[2, 2, 2]] == new[[1, -2]];
             temp = Which[
               oldQ && newQ, midtreecase2211[old, new, al1],
               (! oldQ) && newQ, midtreecase2201[old, new, al1],
               oldQ && (! newQ), midtreecase2210[old, new, al1],
               (! oldQ) && (! newQ),midtreecase2200[old, new,al1]
               ];
             If[ temp===-1,
                 {-1,-1},
                 {2,temp}
             ],
             True, {-1, -1}
             ]
        ]
    ]
    
TreeTransitionType222[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {diff,old, new, al1,al2, oldQ, newQ,groupid,temp},
        diff = comparetree[tree1, tree2];
        Which[
            diff===0,{0,0},
            diff===-1,{-1,-1},
            True,
            {groupid,{old, new}, {al1,al2}} = diff;
            Switch[groupid,
             0,Print["case00 to be done!"];
               {-1,-1},             
             1,
             midtreecase11[old, new, al1, al2,tree1, tree2],
             2,
             oldQ = old[[2, 2, 2]] == old[[1, -2]];
             newQ = new[[2, 2, 2]] == new[[1, -2]];
             temp = Which[
               oldQ && newQ, midtreecase2211[old, new, al1,tree1],
               (! oldQ) && newQ, midtreecase2201[old, new, al1],
               oldQ && (! newQ), midtreecase2210[old, new, al1],
               (! oldQ) && (! newQ),midtreecase2200[old, new,al1,tree1]
               ];
             If[ temp===-1,
                 {-1,-1},
                 {2,temp}
             ],
             True, {-1, -1}
             ]
        ]
    ]
    
nodediffIIold[{tlist1_,elist1_, c12_, rule21_, al1_}, {tlist2_,elist2_, c21_, rule12_, al2_},isprint_:False] :=
    Module[ {n, oldnode, newnode, cc12, remain, edges, cond,set, i, j, v10, v13, vv, pos},
        n = Length[tlist1]/2;
        oldnode = newnode = {};
        cc12 = Select[c12, ! (MemberQ[oldnode, #[[1, 2]]] || #[[1, 2]] == -1) &];
        remain = Flatten[elist2 /. rule21, 1];
        remain = remain //. (Rule @@ # & /@ Select[remain, (StringQ[First[#]] &)]);
        While[cc12 =!= {},
         edges = NestWhile[al1[[#]] &, #, (MemberQ[oldnode, #] &)] & /@al1[[cc12[[All, 1, 2]]]];
         edges = Thread[#] & /@Transpose[{Transpose[Append[Transpose[cc12[[All, All, 1]]], cc12[[All, 1, 2]]]],edges}];
         If[ isprint,
             Print["edges=",edges]
         ];
         cond = Map[MemberQ[remain, #] &, edges, {2}];
         If[ isprint,
             Print["cond=",cond]
         ];
         set = Flatten[Position[cond, {False, True, True} | {True, False, True}, {1},Heads -> False]];
         Do[
          j = set[[i]];
          {v10, v13} = edges[[j, -1]];
          If[ v13==-1,
              Continue[]
          ];
          vv = First[Complement[elist1[[v13 - n, All, 1]], {v10}]];
          vv = If[ vv > n,
                   Prepend[elist1[[vv - n, All, 1]], vv],
                   {vv}
               ];
          vv = Select[vv /. rule12, NumericQ];
          vv = (NestList[al2[[#]] &, #, 2] & /@ vv) /. rule21;
          If[ MemberQ[Flatten[vv], v10],
              cond[[j, 3]] = False
          ], {i,Length[set]}];
         If[ isprint,
             Print["cond=",cond]
         ];
         pos = (Xor[#[[1]], #[[2]]]) && (! #[[3]]) & /@ cond;
         pos = Flatten[Position[pos, True, {1}, Heads -> False]];
         If[ pos === {},
             Break[],
             set = Complement[Flatten[Pick[cc12[[pos, All, 1]], cond[[pos, ;; 2]], False]],oldnode];
             set = Select[set /. rule12,!StringQ[#]&]; 
             (*set = Select[set, Complement[elist1[[# - n, All, 1]],elist2[[# - n, All, 1]] /. rule21] === {} &];*)
             If[ isprint,
                 Print["set=",set]
             ];
             newnode = Union[newnode, al2[[set]]];
             oldnode = Union[oldnode, cc12[[pos, 1, 2]]];
             If[ isprint,
                 Print["oldnode=",oldnode,";newnode=",newnode]
             ];
             cc12 = Select[cc12, ! MemberQ[oldnode, #[[1, 2]]] &]
         ];
         ];
        {oldnode, newnode}
    ]


TreeTransitionType2222[tree1_CoalescentTree, tree2_CoalescentTree] :=
    Module[ {al1, al2, tlist1, elist1, tnew, cc12, rule21, tlist2, elist2,
       told, cc21, rule12, opp, oldnode, newnode, n, newchild, newparent,
       newlabelrule, newnodelabel, oldchild, oldparent, old,
       new, oldQ, newQ},
        If[ SameCoalescentTreeQ[tree1, tree2],
            Return[{0, 0}]
        ];
        al1 = TreeToAdjacencyList[tree1];
        al2 = TreeToAdjacencyList[tree2];
        {tlist1, tlist2, elist1, elist2, told, tnew, cc12, cc21, rule12,rule21} = treediff[tree1,tree2];
        If[ Length[tnew] > 2 || Length[told] > 2,
            Return[{-1, -1}]
        ];
        opp = gettreeoperator1[{tlist1, cc12, rule21, al1}, {tlist2, cc21,rule12, al2}];
        If[ opp =!= -1,
            Return[{1, {opp}}]
        ];
        {oldnode, newnode} = nodediffI[{tlist1, elist1,cc12, rule21, al1}, {tlist2,elist2, cc21, rule12, al2}];
        If[ oldnode === newnode === {},
            {oldnode, newnode} = MapThread[Union, {nodediffIIold[{tlist1, elist1, cc12, rule21, al1}, {tlist2,elist2, cc21, rule12, al2}],
            Reverse[nodediffIIold[{tlist2, elist2, cc21, rule12, al2}, {tlist1,elist1, cc12, rule21, al1}]]}]
        ];
        oldnode = Union[oldnode, Flatten[Position[tlist1, #] & /@Complement[tlist1, Complement[tlist2, tlist2[[newnode]]]]]];
        newnode = Union[newnode, Flatten[Position[tlist2, #] & /@Complement[tlist2, Complement[tlist1, tlist1[[oldnode]]]]]];
        If[ Length[newnode] > 2 || Length[oldnode] > 2 ||Length[newnode] != Length[oldnode],
            Return[{-1, -1}]
        ];
        n = LeafNumber[tree1];
        tnew = tlist2[[newnode]];
        newchild = elist2[[newnode - n, All, 1]];
        newparent = al2[[newnode]];
        newlabelrule = Thread[newnode -> Take[CharacterRange["A", "Z"], Length[newnode]]];
        newnodelabel = newnode /. newlabelrule;
        {newparent, newchild} = {newparent, newchild} /.newlabelrule /. {(x_ /; x > n) :>Position[tlist1, tlist2[[x]]][[1, 1]]};
        told = tlist1[[oldnode]];
        oldchild = elist1[[oldnode - n, All, 1]];
        oldparent = al1[[oldnode]];
        If[ Length[oldchild] == 2 && MemberQ[oldchild[[2]], oldnode[[1]]] &&oldchild[[2, 1]] == oldnode[[1]],
            oldchild[[2]] = Reverse[oldchild[[2]]]
        ];
        If[ Length[newchild] == 2 &&MemberQ[newchild[[2]], newnodelabel[[1]]] &&newchild[[2, 1]] == newnodelabel[[1]],
            newchild[[2]] = Reverse[newchild[[2]]]
        ];
        old = Transpose[{told, oldchild, oldparent, oldnode, oldnode}];
        new = Transpose[{tnew, newchild, newparent, newnodelabel, newnode}];
        If[ Length[old] == Length[new] == 2,
            oldQ = old[[2, 2, 2]] == old[[1, -2]];
            newQ = new[[2, 2, 2]] == new[[1, -2]];
            opp = Which[
             oldQ && newQ, midtreecase2211[old, new, al1],
             (! oldQ) && newQ, midtreecase2201[old, new, al1],
             oldQ && (! newQ), midtreecase2210[old, new, al1],
             (! oldQ) && (! newQ), midtreecase2200[old, new, al1]
             ];
            If[ opp===-1,
                {-1,-1},
                {2,opp}
            ],
            Print["Case to be done!"];
            {-1, -1}
        ]
    ]   
    


End[] (* End Private Context *)

EndPackage[]