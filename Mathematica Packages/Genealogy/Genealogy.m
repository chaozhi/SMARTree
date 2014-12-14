(* Mathematica Package *)

(* Created by the Wolfram Workbench Jul 28, 2010 *)

(* 	Title: 		Mathematica Package Genealogy
	Context: 	Genealogy`.
	Author:		Chaozhi Zheng	
	Copyright:	Copyright 2010-, Chaozhi Zheng
	Mathematica version 6+
	Keywords:
	Summary: 
*)

BeginPackage["Genealogy`",{"Genealogy`CoalescentTree`","Combinatorica`"}]

RandomARG::usage = "RandomARG  "

ARGToLocalTrees::usage = "ARGToLocalTrees  "

GetLocalTree::usage = "GetLocalTree  "

ResetTreeX::usage = "ResetTreeX  "

ToCoalescentTree::usage = "ToCoalescentTree  "

(*SimpleTree::usage = "SimpleTree  "*)

(* Exported symbols added here with SymbolName::usage *) 

Begin["`Private`"]
(* Implementation of the package *)
 (*
1) I denote by a graph the genealogy of n sequences for a Wright-Fisher dipolid population of size N. 
I scale the time by 2 N. We assume the sequences are continuous from 0 to 1. I set vertex weights and edge weights as follows. 
2) I set the edge weight by the set of form {{{x11_,x12_},k1}..} 
where the ancestral material of the lineage (edge) in interval {xi1,xi2} have ki copied in the sampled sequences.
when ki=n, the interval has reached the MRCA and it will be dropped. 
3) I set the vertex weight for the vertex with recombination. The weight has form {p,type} where p is the recombination
position along the sequence, and type is the recombination type.
type=1: recombination in the ancestral material
type=2: recombination in the non-ancestral trapped material
type=3: recombination in the non-ancestral that has ancestral material only to the right
type=4: recombination in the non-ancestral that has ancestral material only to the right
type=5: recombination in a sequence that carries no ancestral material  
*)   

(*n:sample size;
r=4Ne c,scaled recombination rate,where Ne is effective size of a \
diploid population,and r is the recombination rate per base pair per \
generation.Ne=10^4 and c=10^(-8)(ie 1 cM/Mb) are appropriate for the \
human genome;
r=r lenchro is the normalized scaled recombination rate;lenchro is \
the length of chromsome in unit of base pair.For Set lenchro=10kb,r=4 \
10^4 10^(-8) 10^4=4;
*)   
(***************************Illustration****************************)
 
AncestralMerge[{},n_Integer] :=
    {}; 
AncestralMerge[x : {{{_?NumericQ, _?NumericQ}, _Integer} ..}, 
   n_Integer] :=
    Module[ {p, res, temp,i,j},
        p = Union[Flatten[x[[All, 1]]]];
        res = Flatten[
          Table[temp = Select[p, x[[i, 1, 1]] <= # <= x[[i, 1, 2]] &];
                Table[{{temp[[j]], temp[[j + 1]]}, x[[i, 2]]}, {j, 
                  Length[temp] - 1}], {i, Length[x]}], 1];
        res = Gather[res, #1[[1]] == #2[[1]] &];
        res = {#[[1, 1]], Total[#[[All, 2]]]} & /@ res;
        DeleteCases[res, {{_, _}, n}]
    ];
    
AncestralSplit[{},p_?NumericQ] :=
    {{},{}};    

AncestralSplit[x : {{{_?NumericQ, _?NumericQ}, _Integer} ..}, 
   p_?NumericQ] :=
    Module[ {y,i},
        y = Transpose[{{{Min[#[[1]]], Min[p, Max[#[[1]]]]}, #[[2]]}, {{Max[p, Min[#[[1]]]],Max[#[[1]]]}, #[[2]]}} & /@ x];
        Table[Select[y[[i]], #[[1, 1]] < #[[1, 2]] &], {i, 2}]
    ];

AncestralWeight[{}] :=
    0;    
AncestralWeight[x : {{{_?NumericQ, _?NumericQ}, _Integer} ..}] :=
    Total[Abs[#[[2]] - #[[1]]] & /@ x[[All, 1]]];
    
 
AncestralMemberQ[x : {{{_?NumericQ, _?NumericQ}, _Integer} ..}, 
   p_?NumericQ] :=
    Module[ {inter},
        inter = Interval @@ x[[All, 1]];
        IntervalMemberQ[inter, p]
    ];
    
RandomARG::notGMRCA = "ARG has not reach the Grand Most Recent Common Ancestor.";

(*rectype=1: recombination in ancestral material
  rectype=2: recombination in non-ancestral trapped material*)    
RandomARG[n_Integer, r_?NonNegative] :=
    Module[ {vertices, states, edges = {}, cp = 0,
      count = 1, k, weights,totw, totr, isrec, pos, s, recloc,i,rectype},
        vertices = Table[{{i, 0}}, {i, n}];
        states = Table[{i, {{{0, 1}, 1}}}, {i, n}];
        weights = Table[1,{n}];
        totw = n;
        k = n;
        While[totw > 0,         
         totr = r totw;
         cp += RandomReal[ExponentialDistribution[k (k - 1)/2 + totr/2]];
         isrec = If[ RandomReal[] < totr/(k (k - 1) + totr),
                     True,
                     False
                 ];
         If[ isrec,
             pos = RandomChoice[weights -> Range[Length[states]]];
             s = states[[pos, 2, All, 1]];
             recloc = RandomReal[{Min[s], Max[s]}];
             rectype =If[ !AncestralMemberQ[states[[pos, 2]],recloc],2,1];
             AppendTo[vertices, {{vertices[[states[[pos,1]], 1,1]], cp},VertexWeight->{recloc,rectype}}];
             AppendTo[edges, {{states[[pos, 1]], n + count},EdgeWeight -> states[[pos, 2]]}];
             states = Join[Drop[states, {pos}],Transpose[{Table[n + count, {2}],AncestralSplit[states[[pos, 2]], recloc]}]],
             pos = RandomSample[Range[Length[states]], 2];
             AppendTo[vertices, {{N[Mean[vertices[[states[[pos,1]], 1,1]]]], cp}}];
             edges = Join[edges,Table[{{states[[pos[[i]], 1]], n + count},EdgeWeight -> states[[pos[[i]], 2]]}, {i, 2}]];
             states = Append[Delete[states, Partition[pos, 1]], {n + count, AncestralMerge[Join @@ states[[pos, 2]], n]}];
         ];
         count++;
         k = Length[states];
         weights = AncestralWeight[#] & /@ states[[All, 2]];
         totw = Total[weights];
         ];
        If[ states[[1, 2]]=!={},
            Message[RandomARG::notGMRCA]
        ];
        Graph[edges, vertices, EdgeDirection -> True]
    ];   
    
GetTreeTopology[tr_Graph] := 
  Module[{al, verticesY, ls},
   al = ToAdjacencyMatrix[tr];
   verticesY = VerticesY[tr];
   (*each inner node has two offspring nodes*)
   ls = Partition[SortBy[Position[al, 1], verticesY[Last[#]] &], 2];
   ls = #[[All, 1]] -> #[[1, 2]] & /@ ls;
   ls
   ];

ResetTreeX[tr_Graph] :=
    Module[ {nest, rule, vls, leaves, n},
        rule = GetTreeTopology[tr];
        nest=rule[[-1, -1]] //. (Reverse[#] & /@ rule);
        vls = Vertices[tr];
        leaves = Level[nest, {-1}];
        n = Length[leaves];
        vls[[leaves, 1]] = Range[n];        
        Scan[(vls[[#[[2]], 1]] = N[Mean[vls[[#[[1]], 1]]]]) &, rule];
        ChangeVertices[tr, vls]
    ];    
    
   
GetDegree::direction = "The direction must be \"Out\" or \"In\".";
GetDegree[g_Graph, direction_String] :=
    Module[ {i, res, ed, temp},
        i = Switch[direction, "Out", 1, "In", 2, _, 
          Message[GetDegree::direction];
          Return[$Failed]];
        res = Table[0, {V[g]}];
        ed = Edges[g];
        temp = Transpose[Tally[ed[[All, i]]]];
        res[[temp[[1]]]] = temp[[2]];
        res
    ];   
    
SimpleTree[g_Graph] :=
    Module[ {g2 = g, al, inout, root, leaves, innernodes0, innernodes,
       nodes, edges, rule},
        al = ToAdjacencyLists[g2];
        (*inout=Transpose[{InDegree[g2],OutDegree[g2]}];*)
        inout = Transpose[{GetDegree[g2, "In"], GetDegree[g2, "Out"]}];
        {leaves, innernodes0, root} = 
         Flatten[Position[inout, #]] & /@ {{0, 1}, {2, 1}, {2, 0}};
        innernodes = Join[innernodes0, root];
        nodes = Join[leaves, innernodes];
        edges = 
         Prepend[NestWhile[Flatten[al[[#]], 1] &, al[[#]], 
             Intersection[#, nodes] === {} &], #] & /@ 
          Join[leaves, innernodes0];
        rule = MapIndexed[#1 -> First[#2] + Length[leaves] &, innernodes];
        edges = Partition[edges /. rule, 1];
        nodes = Vertices[g2, All][[nodes]];
        Graph[edges, nodes, EdgeDirection -> True]
    ];
   
ARGToLocalTrees[arg_Graph,indices_: All] :=
    Module[ {rp, p, ew, inter, sel, ed, gs, trls},
        rp = Prepend[SortBy[Cases[GetVertexWeights[arg], {_, x_ /; x <= 2}],First], {0, 0}];
        p = MovingAverage[Append[rp[[All, 1]], 1], 2];
        rp = rp[[indices]];
        p = p[[indices]];
        ew = Transpose[{Edges[arg], GetEdgeWeights[arg]}];
        inter = (Interval @@ #[[2, All, 1]]) & /@ ew;
        sel = Transpose[IntervalMemberQ[#, p] & /@ inter];
        ed = Pick[ew[[All, 1]], #, True] & /@ sel;
        gs = Graph[Partition[#, 1], Partition[Vertices[arg], 1], EdgeDirection -> True] & /@ ed;
        trls = SimpleTree[#] & /@ gs;
        Transpose[{rp, trls}]
    ];
       

GetLocalTree[arg_Graph,p_/;0<=p<=1] :=
    Module[ {rp, ew,i, se, g},
        rp = SortBy[Cases[GetVertexWeights[arg], {_, _}], First];
        rp = Append[Prepend[rp, {0, 0}], {1, 0}];
        ew = Transpose[{Edges[arg], GetEdgeWeights[arg]}];
        i = If[ p < 1,
                Count[rp[[All, 1]], _?(# <= p &)],
                Length[rp] - 1
            ];
        se = Select[ew, AncestralMemberQ[#[[2]], p] &][[All, 1]];
        g = Graph[Partition[se, 1], Partition[Vertices[arg], 1], 
          EdgeDirection -> True];
        SimpleTree[g]
    ];    
    
ToCoalescentTree[g_Graph] := Module[{tlist, elist},
   tlist = Vertices[g][[All, 2]];
   elist = Split[SortBy[Edges[g], Last], Last[#1] == Last[#2] &];
   elist = elist /. Dispatch[Thread[Ordering[tlist] -> Range[Length[tlist]]]];
   tlist = Sort[tlist];
   CoalescentTree[tlist, elist]
   ]; 
        
End[]

EndPackage[]

