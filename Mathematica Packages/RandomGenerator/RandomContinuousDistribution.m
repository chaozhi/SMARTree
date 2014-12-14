(* Mathematica Package *)

BeginPackage["RandomGenerator`RandomContinuousDistribution`"]
(* Exported symbols added here with SymbolName::usage *)  

RandomMultiNormVariance::usage = "RandomMultiNormVariance[\[Mu],\[CapitalSigma],n] gives n random draws from the multivariate normal distribution with mean \[Mu] and variance matrix \[CapitalSigma].  The optional n=1."

RandomMultiNormPrecision::usage = "RandomMultiNormPrecision[\[Mu],Q,n] gives n random draws from the multivariate normal distribution with mean \[Mu] and precision matrix Q. The optional n=1."

RandomCanonicalMultiNorm::usage = "RandomCanonicalMultiNorm[b,Q,n] gives n random draws from the multivariate normal distribution with mean Q^(-1)b and precision matrix Q. The optional n=1"

RandomCanonicalMultiNorm2::usage = "RandomCanonicalMultiNorm2[b,Q] gives Length[b] random draws from the multivaratie normal distribution with mean Q^(-1)b[i] and precision matrix Q, for i=1,..., Length[b]."

RandomOneSideTruncNorm::usage = "RandomOneSideTruncNorm[truncside, b, \[Mu]] gives a random draw for n normal distributions with mean \[Mu][[i]] and standard deviation 1, truncated to be right of b if truncside[[i]]=1 and to be left of b[[i]] if truncside[[i]]=-1, for i=1, ..., n. \n RandomOneSideTruncNorm[truncside, b,\[Mu],n] gives n random draws for the normal distriubtion with mean \[Mu] and standard deviation 1, truncted to be right of b if truncside=1 and to be left of b if truncside=-1. RandomOneSideTruncNorm[truncside, b, \[Mu]] for n=1. The optional n=1. RandomOneSideTruncNorm[truncside, b] for a standard normal distribution."

Begin["`Private`"] (* Begin Private Context *) 

SymmetricMatrix[A_?MatrixQ] :=
    (A + Transpose[A])/2;

(*RandomMultiNormVariance[mu_,V_,n_:1] gives n random draw from \
MultiNorm[mu,V], where mu is mean, V is covariance matrix*)
RandomMultiNormVariance[mu_/;VectorQ[mu,NumericQ], V_/;MatrixQ[V,NumericQ], n_: 1] :=
    Module[ {down, z, res},
        down = Transpose[CholeskyDecomposition[SymmetricMatrix[V]]];
        z = RandomReal[NormalDistribution[0, 1], {n, Length[V]}];
        res = mu + down.# & /@ z;
        If[ n == 1,
            Flatten[res],
            res
        ]
    ]

RandomMultiNormPrecision[mu_/;VectorQ[mu,NumericQ], Q_/;MatrixQ[Q,NumericQ], n_: 1] :=
    Module[ {z, res, LS},
        z = RandomReal[NormalDistribution[], {n, Length[Q]}];
        LS = LinearSolve[CholeskyDecomposition[SymmetricMatrix[Q]]];
        res = (mu + LS[#]) & /@ z;
        If[ n == 1,
            Flatten[res],
            res
        ]
    ]

RandomCanonicalMultiNorm[imu_/;VectorQ[imu,NumericQ], iV_/;MatrixQ[iV,NumericQ], size_: 1] :=
    Module[ {imu2, temp},
        imu2 = Table[imu, {size}];
        temp = RandomCanonicalMultiNorm2[imu2, iV];
        If[ size == 1,
            temp[[1]],
            temp
        ]
    ]

RandomCanonicalMultiNorm2[imu_/;MatrixQ[imu,NumericQ], iV_/;MatrixQ[iV,NumericQ]] :=
    Module[ {LS, LS2, iU, mu, rnd},
        LS = LinearSolve[SymmetricMatrix[iV], Method -> "Cholesky"];
        mu = LS[#] & /@ imu;
        rnd = RandomReal[NormalDistribution[0, 1], Dimensions[imu]];
        iU = LS[[2, 3, 1]];
        LS2 = LinearSolve[iU];
        mu + (LS2[#] & /@ rnd)
    ]

(*Geweke J, 1991,Computer Sciences and Statistics:Proc .23d \
Symp.Interface,pp .571\[Dash]577*)
RandomOneSideTruncNorm[truncside:1|-1, b_?NumericQ] :=
    Module[ {lb, z, p, u, res},
        lb = truncside b;
        If[ lb <= 0.45,
            While[True, z = RandomReal[NormalDistribution[0, 1]];
                        If[ z >= lb,
                            res = z;
                            Break[]
                        ]],
            While[True,
             z = -Log[RandomReal[UniformDistribution[{0, 1}]]]/lb + lb;
             p = Exp[-(z - lb)^2/2];
             u = RandomReal[UniformDistribution[{0, 1}]];
             If[ u <= p,
                 res = z;
                 Break[]
             ]
             ]
        ];
        truncside res
    ]
    
RandomOneSideTruncNorm[truncside_ 1 | -1, b_?NumericQ, mu_?NumericQ] :=
    RandomOneSideTruncNorm[truncside, b - mu] + mu   

RandomOneSideTruncNorm[
	truncside:1 | -1, 
	b_?NumericQ, 
	mu_?NumericQ,
	n_Integer?Positive,
	r_:0.05] :=RandomOneSideTruncNorm[Table[truncside,{n}],Table[b,{n}],Table[mu,{n}],r]
   
RandomOneSideTruncNorm[
    truncside_/;VectorQ[truncside, MatchQ[#, _ 1 | -1] &],    
    b_/;VectorQ[b,NumericQ], 
    mu_/;VectorQ[mu,NumericQ],
    r_:0.05] :=
    Module[ {n, res, pos, sample, subpos,ratio,j},
        n = Length[b];
        res = Table[0, {n}];
        pos = Range[n];
        While[True,
         sample = mu[[pos]]+RandomReal[NormalDistribution[0, 1], Length[pos]];
         subpos = Flatten[
            Position[truncside[[pos]]  (sample - b[[pos]]), _?NonNegative]];
         res[[pos[[subpos]]]] = sample[[subpos]];
         ratio = N[Length[subpos]/Length[pos]];
         pos = Complement[pos, pos[[subpos]]];
         If[ ratio < r  || Length[pos] == 0,
             Break[]
         ];
         ];
        If[ Length[pos] > 0,
            res[[pos]] = 
             Table[
             	j = pos[[i]];
                RandomOneSideTruncNorm[truncside[[j]], b[[j]],mu[[j]]], {i,Length[pos]}];
        ];
        res
    ]/;Length[truncside]===Length[b]===Length[mu]
    
End[] (* End Private Context *)

EndPackage[]