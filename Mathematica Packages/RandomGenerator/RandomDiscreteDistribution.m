(* Mathematica Package *)

BeginPackage["RandomGenerator`RandomDiscreteDistribution`"]
(* Exported symbols added here with SymbolName::usage *)  

DiscreteRatioOfUniform::usage = "DiscreteRatioOfUniform[f,lb,ub,S, xm,ym,cdfm, n] gives n random draws from the distributed f(x), truncted to be in the range [lb, ub]. ym=f(xm) and cdfm=CDF(xm) where ym is the maximum of f(xm). f(x) is assumed to be T_ (-1/2) concave, and 1/S is the normalization constant of f(x). By defult, cdfm is not given (-1) and n=1."

RandomNegativeBinomial::usage = "RandomNegativeBinomial[n, p, lb, s] gives s random draws from the negative binomial distribution with parameters n and p, truncted to be not less than lb. The defaults lb=0 and s=1." 

RandomPoisson::usage = "RandomPoisson[lam, lb, s] gives s random draws from the poisson distribution with mean lam, truncted to be not less than lb. The defaults lb=0 and s=1."

RandomBinomial::usage="RandomBinomial[n,p,lb,s] gives s random draws from the binomial distribution with parameters n and p, truncted to be not less than lb. The default lb=0 and s=1.";

Begin["`Private`"] (* Begin Private Context *) 
    

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
    
        
logpdfNegativeBinomial[n_, p_,x_] :=
    x Log[(1 - p)] + n Log[ p] + LogGamma[N[n + x]] - LogGamma[N[n]] - LogGamma[N[x + 1]];
(*CDF[NegativeBinomialDistribution[n,p],lb-1]*)
RandomNegativeBinomial[n_Integer?Positive, p_/;1>p>0, lb_: 0, size_: 1] :=
    Module[ {xm, ym, SS},
        xm = Floor[(n - 1) (1 - p)/p];
        xm = Max[lb, xm];
        ym = logpdfNegativeBinomial[n, p,xm];
        SS = If[ lb == 0,
                 1,
                 1 - BetaRegularized[p, n, Floor[lb]]
             ];
        DiscreteRatioOfUniform[Exp[logpdfNegativeBinomial[n, p,#]]&, lb,
          Infinity, SS, xm, Exp[ym], -1, size]
    ];

logpdfPoisson[lam_,x_] :=
    -lam + x Log[lam] - LogGamma[x + 1.];
(*CDF[PoissonDistribution[lam],lb-1]*)
RandomPoisson[lam_?Positive, lb_: 0, size_: 1] :=
    Module[ {xm, ym, SS},
        xm = If[ IntegerPart[lam] == lam,
                 IntegerPart[lam] - 1,
                 Floor[lam]
             ];
        xm = Max[lb, xm];
        ym = logpdfPoisson[lam,xm];
        SS = If[ lb == 0,
                 1,
                 1 - GammaRegularized[Floor[lb], lam]
             ];
        DiscreteRatioOfUniform[Exp[logpdfPoisson[lam,#]]&, lb, Infinity,
          SS, xm, Exp[ym], -1, size]
    ];
    

RandomBinomial::usage = "RandomBinomial[n, p, lb, s] gives s random draws from the binomial distribution with parameters n and p, truncted to be not less than lb. The defaults lb=0 and s=1."
logpdfBinomial[n_, p_,x_] :=
    (n - x) Log[1 - p] + x Log[p] + LogGamma[n + 1.] - LogGamma[n - x + 1.] - LogGamma[x + 1.];
RandomBinomial[n_Integer, p_/;1>=p>=0, lb_: 0, size_: 1] :=
    Module[ {xm, ym, SS},
        xm = Floor[(n + 1) p];
        xm = Min[n, Max[lb, xm]];
        ym = logpdfBinomial[n, p,xm];
        SS = If[ lb == 0,
                 1,
                 1 - BetaRegularized[1 - p, 1 + n - Floor[lb], Floor[lb]]
             ];
        DiscreteRatioOfUniform[Exp[logpdfBinomial[n, p,#]]&,  lb, 
         Infinity, SS, xm, Exp[ym], -1, size]
    ];
       
       
End[] (* End Private Context *)

EndPackage[]