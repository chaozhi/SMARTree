(* Mathematica Package *)

(* Created by the Wolfram Workbench Oct 24, 2010 *)

BeginPackage["RandomGenerator`"]
(* Exported symbols added here with SymbolName::usage *) 

DiscreteRatioOfUniform::usage = "DiscreteRatioOfUniform[f,lb,ub,S, xm,ym,cdfm, n] gives n random draws from the distributed f(x), truncted to be in the range [lb, ub]. ym=f(xm) and cdfm=CDF(xm) where ym is the maximum of f(xm). f(x) is assumed to be T_ (-1/2) concave, and 1/S is the normalization constant of f(x). By defult, cdfm is not given (-1) and n=1."

Begin["`Private`"]
(* Implementation of the package *)

(*Homann, W., etal 2004, Automatic nonuniform random variate \
generation. 
Berlin Heidelberg: Springer-Verlag. Chap.10. Page231. 
where truncation[u/v] was revised into round[u/v]
Assumpitions: fun is a T_ (-1/2) concave, or log-concave (stronger) \
function
fun: unnormalized pdf
lb (ub): [lb,ub] is the domain
xm (ym): mode and fun value at mode
SS: sum of fun over all possible values
cdfm: cdf at mode F (xm)
size: sample size
*)
DiscreteRatioOfUniform[
    fun_Function, lb_, ub_, 
    SS_?NumericQ, xm_?NumericQ, ym_?NumericQ, 
    cdfm_: - 1, size_: 1] :=
    Module[ {vm, u1, u2, u, v, xprop,res,i},
        vm = Sqrt[ym];
        {u1, u2} = If[ cdfm >= 0,
                       {-cdfm SS, (1 - cdfm) SS}/vm,
                       {-SS, SS}/vm
                   ];
        res = Table[
          While[True,
           u = RandomReal[{u1, u2}];
           v = RandomReal[{0, vm}];
           xprop = Round[u/v] + xm;
           If[ (lb <= xprop <= ub) && (v^2 <= fun[xprop]),
               Break[]
           ];
           ];
          xprop, {i, size}];
        If[ size==1,
            First[res],
            res
        ]
    ]
    

End[]

EndPackage[]

