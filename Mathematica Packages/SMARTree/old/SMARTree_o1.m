(* Mathematica Package *)

(* Created by the Wolfram Workbench 2013-2-6 *)

BeginPackage["SMARTree`",{"Toolbox`","Metropolis`","RandomGenerator`","Genealogy`CoalescentTree`","Genealogy`LocalTrees`","SMARTree`Likelihood`"}]

initializeSMARTchain::usage = "initializeSMARTchain[inputsnpls,nsq,nbp,initheta,inirho,inieps] generates initial values for one MCMC chain, and returns {{theta, epsilon, rho}, treels, adpaptive parameters, snpls}. In the last column of snpls, the log likelihood given the initialized treels over-writes the values in inputsnpls."

saveSMARTchain::usage = "saveSMARTchain[mcstate,nsq,adpheat,outputconcise,outputcold] saves mcstates and adpheat into output files outputconsise and outputcold."

swapSMARTchain::usage = "swapSMARTchain[mcstate[[i]],nbp,isadp,adprate,adpheat,size] updates mcstate[[i]] states for group i by swapping the temperatures of two random chosen chains within the group in a number of size times. If isadp=True, the difference of inverse temperatures in group i is adpated according the accept ratio of swapping. Returns updated mcstate[[i]] and  adpheat."

permuteSMARTchain::usage = "permuteSMARTchain[mcstate[[i]],nbp,isadp,adprate,adpheat,size] updates mcstate[[i]] states for group i by permuting the temperatures of chains within the group in a number of size times. If isadp=True, the difference of inverse temperatures in group i is adpated according the accept ratio of permuting. Returns updated mcstate[[i]] and  adpset."

updateSMARTchain::usage = "updateSMARTchain[chainstate,nsq,nbp, pritheta,prirho,prieps,isadp,adprate,rjfreq,isprint] updates chainstate for one MCMC chain."

testtreels::usage = "testtreels[treels,msg] gives True if the total branch lengths (column 3 of treels) for trees (column 2) are consistent and the log transtion probabilities between sucessive trees  (or the log prior probability for the first tree; column 4) are consistent, otherwise give False and print the message msg."

testsnpls::usage = "testsnpls[snpls,treels,nbp,theta,epsilon,modelid,msg] gives True if the log likelihoods of SNP data given the local trees (column 2 of treels) are consistent, and otherwise gives False and prints message msg."

(* Exported symbols added here with SymbolName::usage *) 

Begin["`Private`"]
(* Implementation of the package*)

p00 = 0.25;

initializeSMARTchain[inputsnpls_,nsq_, nbp_,initheta_,inirho_,inieps_,modelid_] :=
    Module[ {inithetashape, inithetascale, inithetamin, inithetamax, inirhoshape, inirhoscale, inirhomin, inirhomax, 
    	iniepsalpha, iniepsbeta, iniepsmin, iniepsmax,snpls = inputsnpls,theta, epsilon, rho, treels, adpscale, 
    	adpthetaerr,adphyperscale,adpshift,adpposseglen,adprho, adphypershift,indtr,temp,treels2,treemap,ytree},
        {inithetashape, inithetascale, inithetamin, inithetamax} = initheta;
        {inirhoshape, inirhoscale, inirhomin, inirhomax} = inirho;
        {iniepsalpha, iniepsbeta, iniepsmin, iniepsmax} = inieps;
        epsilon = RandomReal[TruncatedDistribution[{iniepsmin, iniepsmax},BetaDistribution[iniepsalpha,iniepsbeta]]];
        theta =RandomReal[If[inithetascale == Infinity, 
 				UniformDistribution[{inithetamin, inithetamax}],
 				TruncatedDistribution[{inithetamin, inithetamax},GammaDistribution[inithetashape, inithetascale]]]];            
        rho =RandomReal[If[inirhoscale == Infinity, 
 				UniformDistribution[{inirhomin, inirhomax}],
 				TruncatedDistribution[{inirhomin, inirhomax},GammaDistribution[inirhoshape, inirhoscale]]]]; 
        treels =  Most[RandomSequentialTrees[RandomCoalescentTree[nsq], rho, 1, nbp][[All, {1, -1}]]];
        temp = {TotalBranchHeight[treels[[All, 2]]],
              Prepend[Log[TreeTransitionProb @@ # & /@Partition[treels[[All, 2]], 2, 1]], TreeLogPriorProb[treels[[1, 2]]]]};
         (*treels column name: 
           1: location of a change-point;
           2: the tree from the change point (including the point);
           3: totalbranchheight of the tree;
           4: log transition probabiltiy from tree at the prevous chainge-point to tree at current pooint
              or log prior probability for the first tree. The probability is unconditional.          
         *)
        If[ Length[treels]==1,
            treels[[1]] = Join[treels[[1]],temp[[All,1]]],
            treels = Transpose[Join[Transpose[treels], temp]];
        ];
        treels2 = Append[treels, ReplacePart[treels[[-1]], 1 -> nbp + 1]];
        treemap = IndexByInterpolation[treels2[[All, 1]]];
        ytree = treels2[[treemap[snpls[[All, 1]]], 2]];
        snpls[[All, -1]] = TreeLogLikelihood[snpls[[All, 2]], ytree, snpls[[All, 3]], theta, epsilon,modelid];
        adpthetaerr = InitialProposal[0.234, {theta, epsilon}];
        (*adpthetaerr[[-1, 1, 1]]=theta;
        adpthetaerr[[-1, 2, 2]]=epsilon;
        *)
        adpscale = InitialProposalScale[0.44];
        adpshift = InitialProposalScale[0.44];
        adphyperscale = InitialProposal[0.234, {0,theta, epsilon}];
        adphypershift = InitialProposal[0.234, {Mean[treels[[All,3]]]/nsq, theta, epsilon}];
        adpposseglen = InitialProposalScale[0.234];
        adpposseglen[[-1]] = 5;
        adprho = InitialProposal[0.44,rho];
        adprho[[-1]]=rho;
        indtr = Table[0,{11}];
        {{theta, epsilon,rho}, treels, {adprho,adpthetaerr,adpscale,adphyperscale,adpshift,adphypershift,adpposseglen},indtr,{},snpls}
    ]
    
(*
mcstate:{heat,logposterior (with heat = 1),{theta,epsilon,rho},treels,
{adpscale,adpthetaerr,adphyper,adpshift,adphypershift,indtr},snpls};
*)     
saveSMARTchain[mcstate_, nsq_,adpheat_,outputconcise_, outputcold_] :=
    Module[ {mcstate2, mcstate3,i,ind},
    	ind=Flatten[Transpose[Table[Range[Length[mcstate]], {Dimensions[mcstate][[2]]}]]];
        mcstate2 = Flatten[mcstate, 1][[All, ;;-2]];
        Do[
          mcstate2[[i, 5]] = Prepend[mcstate2[[i, 5]], adpheat[[ind[[i]]]]];
          mcstate2[[i, 4]] = Length[mcstate2[[i, 4]]], {i, Length[mcstate2]}];
        mcstate2[[All, 5]] = mcstate2[[All, 5, All, {4, 5}]];
        PutAppend[mcstate2, outputconcise];
        mcstate3 = First[Select[#[[All, ;; -2]], First[#] == 1 &]] & /@ mcstate;
        Do[
          mcstate3[[i, 5]] = Prepend[mcstate3[[i, 5]], adpheat[[ind[[i]]]]];
          mcstate3[[i, 4, All, 2]] = {#[[1, nsq ;;]], TreeToAdjacencyList[#]}&/@mcstate3[[i, 4, All, 2]], {i, Length[mcstate3]}];
        mcstate3[[All, 5]] = mcstate3[[All, 5, All, {4, 5}]];
        PutAppend[mcstate3, outputcold];
    ]        
      
updateSMARTchain[chainstate_, nsq_, nbp_, pritheta_,prirho_,prieps_,modelid_,isadp_,adprate_,rjfreq_,isprint_] :=
    Module[ {prithetashape, prithetascale, prithetamin, prithetamax, prirhoshape, prirhoscale, prirhomin, prirhomax, 
    	priepsalpha, priepsbeta, priepsmin, priepsmax, heat,logl, theta, epsilon,rho, treels, adpscale, adprho,
    	adpthetaerr,adphyperscale,adpshift,adphypershift,adpposseglen, segarlist,indtr, deltd,snpls,direct,istest},
        {heat, logl, {theta, epsilon,rho}, treels,{adprho,adpthetaerr,adpscale,adphyperscale,adpshift,adphypershift,adpposseglen},
        	indtr,segarlist,snpls} = chainstate;            
        {prithetashape, prithetascale, prithetamin, prithetamax} = pritheta;
        {prirhoshape, prirhoscale, prirhomin, prirhomax} = prirho;
        {priepsalpha, priepsbeta, priepsmin, priepsmax} = prieps;
        direct = RandomChoice[{1,-1}];
        If[ direct == -1,
            {snpls, treels} = reversechromosome[snpls, treels, nbp]
        ];
        istest = RandomReal[]<0.1;
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"0:"<>ToString[heat]]
            ];
        ];
        If[ isprint,
            Print["update rj..."]
        ];
        {snpls,treels,adpposseglen,segarlist,indtr} = rjupdatetreels[snpls, treels, theta, epsilon, nsq,nbp, 
            prirhoshape, prirhoscale,prirhomin,prirhomax,modelid,heat, rjfreq,isadp,adprate,adpposseglen];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"1: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];            
        If[ isprint,
            Print["update timescale..."]
        ];
        {snpls,treels, adpscale} = updatetimescale[snpls, treels, theta, epsilon, nsq, nbp,prirhoshape, prirhoscale,prirhomin,prirhomax, 
        	modelid,heat, isadp, adprate, adpscale];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"2: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];
        If[ isprint,
            Print["update timeshift..."]
        ];
        {snpls,treels, adpshift} = updateshift[snpls, treels, theta, epsilon, nsq, nbp,prirhoshape, prirhoscale,prirhomin,prirhomax, 
        	modelid,heat, isadp, adprate, adpshift];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"3: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];
        If[ isprint,
            Print["update ThetaEpsilon..."]
        ];    
        {snpls,{theta, epsilon}, adpthetaerr} = updateThetaEpsilon[snpls, treels, theta, epsilon,rho,nbp, 
            prithetashape,prithetascale, prithetamin,prithetamax,priepsalpha, priepsbeta,priepsmin,priepsmax,prirhoshape, prirhoscale,prirhomin,prirhomax,
            modelid,heat, isadp, adprate,adpthetaerr];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"4: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];
        If[ isprint,
            Print["update hyperscale..."]
        ];
        {snpls, treels, {theta, epsilon}, adphyperscale} = updatehyperscale[snpls, treels, theta, epsilon, nsq, nbp, 
             prirhoshape, prirhoscale,prirhomin,prirhomax, prithetashape,prithetascale,prithetamin, prithetamax,priepsalpha, priepsbeta, priepsmin,priepsmax,
             modelid,heat, isadp, adprate, adphyperscale];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"5: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];
        If[ isprint,
            Print["update hypershift..."]
        ];  
        {snpls, treels, {theta, epsilon}, adphypershift} = updatehypershift[snpls, treels, theta, epsilon, nsq, nbp, 
            prirhoshape, prirhoscale,prirhomin,prirhomax, prithetashape,prithetascale,prithetamin, prithetamax,priepsalpha, priepsbeta, priepsmin,priepsmax,
            modelid,heat, isadp, adprate, adphypershift];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"6: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];
        If[ isprint,
            Print["update rho..."];
        ];  
        (*rho = updateRho[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];*)
        {rho,adprho}=updateRho[treels, rho, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat,isadp, adprate, adprho, heat];
         If[ isprint,
            Print["calculate logl..."]
        ]; 
        logl = Total[snpls[[All,4]]]+Total[treels[[All, -1]]]+Log[Times@@treels[[;;-2,3]]];
        logl+=(priepsalpha-1) Log[epsilon]+(priepsbeta-1) Log[1-epsilon]+
                (prithetashape-1)Log[theta]-theta/prithetascale+
                (prirhoshape-1)Log[rho]-rho/prirhoscale;
        deltd = Append[Differences[treels[[All, 1]]], nbp - treels[[-1, 1]]];
        logl+=(Length[treels] - 1) Log[rho/2]-(deltd.treels[[All, 3]]) rho/2;
        logl+=logpriSNP[snpls[[All, 1]],treels,1, Length[treels], theta,epsilon, nbp,modelid];
        If[ direct == -1,
            {snpls, treels} = reversechromosome[snpls, treels, nbp]
        ];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"7: "<>ToString[heat]]
            ];
        ];
        {heat, logl, {theta, epsilon,rho}, treels,{adprho,adpthetaerr,adpscale,adphyperscale,adpshift,adphypershift,adpposseglen}, indtr,segarlist,snpls}
    ]
    
permuteSMARTchain[inputmcstate_, nbp_,isadp_, adprate_, inputadpset_,size_] :=
    Module[ {mcstate = inputmcstate,adpset = inputadpset, ii, ar,ind = 0},
        Do[
            ii = RandomSample[Range[Length[mcstate]]];
            ar = mcstate[[ii, 1]].mcstate[[All, 2]]-mcstate[[All, 1]].mcstate[[All, 2]];
            ar = Min[1, Exp[ar]];
            ind+=ar;
            If[ RandomReal[] < ar,
                mcstate[[All, 1]] = mcstate[[ii, 1]];
            ],{size}];
        mcstate = SortBy[mcstate,First];
        adpset = AdaptiveProposalScale[isadp, ind/size, adprate, adpset];
        adpset[[-1]] = Min[adpset[[-1]],1./Length[mcstate]];
        If[ isadp,
            mcstate[[All,1]] = 1 - Last[adpset] Range[Length[mcstate] - 1, 0, -1];
        ];
        {mcstate, adpset}
    ]    

swapSMARTchain[inputmcstate_, nbp_, isadp_, adprate_,inputadpset_,size_] :=
    Module[ {mcstate = inputmcstate, adpset = inputadpset, ii,k, ar, ind},
        ind = Table[0,{size}];
        Do[
         ii = RandomSample[Range[Length[mcstate]], 2];
         ar = mcstate[[Reverse[ii], 1]].mcstate[[ii, 2]]-mcstate[[ii, 1]].mcstate[[ii, 2]];
         ar = Min[1, Exp[ar]];
         ind[[k]] = ar;
         If[ RandomReal[] < ar,  
             mcstate[[ii, 1]] = mcstate[[Reverse[ii], 1]];
            ], {k,size}]; 
        mcstate = SortBy[mcstate, First];
        adpset = AdaptiveProposalScale[isadp, Mean[ind], adprate, adpset];
        adpset[[-1]] = Min[adpset[[-1]], 1./Length[mcstate]];
        If[ isadp,
            mcstate[[All, 1]] = 1 - Last[adpset]  Range[Length[mcstate]-1,0, -1];
        ];
        {mcstate, adpset}
    ]  
    
reversechromosome[inputsnpls_, inputtreels_, nbp_] :=
    Module[ {snpls = inputsnpls, treels = inputtreels, temp,istest},
        istest = RandomReal[]<0.1;
        If[ istest,
            testtreels[treels,"It is before reverse chromosome!"];
        ];
        If[ snpls =!= {},
            snpls[[All, 1]] = nbp + 1 - snpls[[All, 1]];
            snpls = Reverse[snpls];
        ];
        treels[[All, 1]] = nbp + 1 - Append[Rest[treels[[All, 1]]] - 1, nbp];
        temp = TreeLogPriorProb[treels[[All, 2]]];
        treels[[All, -1]] = Append[treels[[2 ;;, -1]] - Differences[Log[treels[[All, 3]]]] - 
           Differences[temp], Last[temp]];
        treels = Reverse[treels];
        If[ istest,
            testtreels[treels,"It is after reverse chromosome!"];
        ];
        {snpls, treels}
    ]    

testtreels[treels_, msg_] :=
    Module[ {temp, b1, b2},
        temp = Prepend[Log[TreeTransitionProb @@ # & /@ Partition[treels[[All, 2]], 2, 1]],
        			TreeLogPriorProb[treels[[1, 2]]]];
        b1 = If[ Abs[Total[temp] - Total[treels[[All, -1]]]] > 10^(-6.),
                 Print["{treelstrue,treels(All,-1)}=", Transpose[{temp, treels[[All, -1]]}]];
                 Print["Wrong in treels[[All,-1]]! " <> msg];
                 False,
                 True
             ];
        temp = TotalBranchHeight[treels[[All, 2]]];
        b2 = If[ Abs[Total[temp] - Total[treels[[All, 3]]]] > 10^(-6.),
                 Print["{treelstrue,treels(All,3)}=", Transpose[{temp, treels[[All, 3]]}]];
                 Print["Wrong in treels[[All,3]]! " <> msg];
                 False,
                 True
             ];
        b1 && b2
    ];      

testsnpls[snpls_,treels_,nbp_,theta_,epsilon_,modelid_,msg_] :=
    Module[ {ls,trpos,temp},
        ls = Append[treels[[All, 1]], nbp + 1];
        trpos = IndexByInterpolation[ls][snpls[[All, 1]]];
        temp = TreeLogLikelihood[snpls[[All, 2]], treels[[trpos, 2]],snpls[[All, 3]], theta, epsilon,modelid];
        If[ (Abs[Total[temp]-Total[snpls[[All, -1]]]]>10^(-4.)),
            Print["{snplstrue,snpls (All,-1)}=",Transpose[{temp,snpls[[All,-1]]}]];
            (*Print["{snpls,treels,nbp,theta,epsilon}=",{snpls,treels,nbp,theta,epsilon}];*)
            Print["Wrong in snpls! "<>msg];
            False,
            True
        ]
    ]      
      
updateThetaEpsilon[inputsnpls_,treels_,inputtheta_, inputepsilon_, rho_,nbp_,prithetashape_, 
    prithetascale_,prithetamin_,prithetamax_,priepsalpha_,priepsbeta_,priepsmin_,priepsmax_,
    prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_,isadpt_, adprate_, inputadpset_] :=
    Module[ {snpls = inputsnpls,adpset = inputadpset,theta = inputtheta,epsilon = inputepsilon,
        logprihyper,treels2,treemap,ytree,deltxls,size = 2,res,propx,propsnpls,ar,i},
        logprihyper = Function[{x},
            (prithetashape-1)Log[x[[1]]]-x[[1]]/prithetascale + 
            (priepsalpha - 1) Log[x[[2]]] + (priepsbeta - 1) Log[1 - x[[2]]]
            ];
        treels2 = Append[treels, ReplacePart[treels[[-1]], 1 -> nbp + 1]];
        treemap = IndexByInterpolation[treels2[[All, 1]]];
        ytree = treels2[[treemap[snpls[[All, 1]]], 2]];
        deltxls = StudentTDisplacements[adpset[[-3]], adpset[[-1]]+10^(-10.) IdentityMatrix[2], size,
            SamplingPattern -> {"Eigen", Infinity, True}];
        res = Table[0, {size}];
        Do[
          propx = {theta, epsilon} + deltxls[[i]];
          If[ prithetamin<propx[[1]]<prithetamax && priepsmin < propx[[2]] <priepsmax,
              propsnpls = snpls;
              If[ propsnpls==={},
                  ar = 0,
                  propsnpls[[All, -1]] = TreeLogLikelihood[propsnpls[[All, 2]], ytree, propsnpls[[All, 3]],propx[[1]], propx[[2]],modelid];
                  ar = Total[propsnpls[[All, -1]]] - Total[snpls[[All, -1]]];                  
              ];
              ar += logprihyper[propx] - logprihyper[{theta, epsilon}]; 
              ar+= logpriSNP[snpls[[All, 1]],treels, 1, Length[treels], propx[[1]],propx[[2]], nbp,modelid]-
              		logpriSNP[snpls[[All, 1]],treels, 1, Length[treels], theta,epsilon, nbp,modelid];             
              ar=heat ar+logpritree[treels, nbp, propx[[1]],prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                  		logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
              ar = Min[1, Exp[ar]];
              If[ RandomReal[] < ar,
                  {theta, epsilon} = propx;
                  snpls = propsnpls;
              ],
              ar = 0
          ];
          res[[i]] = {ar, {theta, epsilon}}, {i, size}];
        adpset = AdaptiveProposal[isadpt, Mean[res[[All, 1]]], adprate,adpset, res[[All, -1]]];
        {snpls, {theta,epsilon}, adpset}
    ]
    
(*updateRho[treels_, nbp_, theta_,prilambdamin_,prilambdamax_,heat_] :=
    Module[ {deltd,a,b,rho},
        deltd = Append[Differences[treels[[All, 1]]], nbp - treels[[-1, 1]]];
        a = heat (Length[treels]-1);
        b = heat (Total[deltd treels[[All, 3]]]/2);
        rho = RandomReal[GammaDistribution[a+1,1/b]];  
        While[!(theta prilambdamin<rho<theta prilambdamax),
            rho = RandomReal[GammaDistribution[a+1,1/b]]
        ];
        rho
    ]
*)    
       
updateRho[treels_, rho_, nbp_, theta_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,heat_,isadp_, adprate_, adprho_, heat_] :=
    Module[ {deltd,a,b,logpdf, res},
        deltd = Append[Differences[treels[[All, 1]]], nbp - treels[[-1, 1]]];
        a = heat (Length[treels] + prirhoshape-2);
        b = heat (1/prirhoscale + Total[deltd treels[[All, 3]]]/2);
        If[prirhomin==0&&prirhomax==Infinity,
        	{RandomReal[GammaDistribution[a+1,1/b]],adprho},
	        logpdf = Function[{x}, a Log[x] - x b];
	        res = AdaptiveMetropolis[logpdf, rho, 1, 1, Boole[isadp], adprate, adprho, DomainConstraint -> (prirhomin<#<prirhomax&)];
	        res[[{1, -1}]]
        ]
    ]    
    
logpritree[treels_, nbp_, theta_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,heat_] :=
    Module[ {deltd, a,b,res},
        deltd = Append[Differences[treels[[All, 1]]], nbp-treels[[-1, 1]]];
        a = heat (Length[treels] + prirhoshape-2);
        b = heat (1/prirhoscale + Total[deltd treels[[All, 3]]]/2);
        res=-heat (Length[treels]-1) Log[2.]-(a+1) Log[b];
        res+Which[
        	prirhomin==0&&prirhomax==Infinity,LogGamma[a+1],
        	prirhomax> (a+1)/b,
        	Log[Gamma[a+1]-Gamma[a+1,prirhomax b]-Gamma[a+1,0,prirhomin b]],
        	True,
        	Log[Gamma[a+1,prirhomin b,prirhomax b]]
        ]
    ]    
    
logpritbl[treels_, nbp_, theta_]:=Module[{forest,nleaf,tbl,set},
	forest = treels[[All, 2]];
    nleaf = LeafNumber[forest[[1]]];
    tbl = treels[[All, 3]];
	If[nleaf >=10,	 
	 (*snp acertainment for n=40 depends on # ma among initial 10 chromosomes,asysmetric*)
	 nleaf = 10;
	 forest = SubCoalescentTree[#, Range[nleaf]] & /@ forest;
	 tbl = TotalBranchHeight[forest];
	 ];
	set = Map[SortBy[Flatten[#, 1], First] &, forest[[All, 2]]][[All, ;; nleaf]];
	tbl -= Total[MapThread[#1[[#2]] &, {forest[[All, 1]], set[[All, All, 2]]}],{2}];
	tbl
]    

(*This version of logpriSNP does not depend on epsilon*)
(*logpriSNP1[snpx_,treels_, istart_, iend_, theta_, epsilon_, nbp_] :=
    Module[ {rbp,dls, tbl},
        If[ snpx==={},
            0,
            rbp = If[ iend >= Length[treels],
                      nbp + 1,
                      treels[[iend+1,1]]
                  ];
            dls = Append[treels[[istart;;iend,1]],rbp];
            dls = Differences[dls] - BinCounts[snpx, {dls}];
            tbl = logpritbl[treels[[istart ;; iend]], nbp, theta];
            -dls.tbl theta/2
        ]
    ]   
*)

(*This version of logpriSNP does depend on epsilon*)    
logpriSNP[snpx_,treels_, istart_, iend_, theta_, epsilon_, nbp_, modelid_] :=
    Module[ {dls, forest},
        If[ snpx==={}||modelid=="M2",
            0,            
            dls = If[ iend >= Length[treels],
                      Append[treels[[istart ;;, 1]], nbp + 1],
                      treels[[istart ;; iend + 1,1]]
                  ];
            dls = Differences[dls] - BinCounts[snpx, {dls}];
            forest = treels[[istart ;; iend, 2]];
            dls.Log[1-PolymorphicProb[forest,theta,epsilon]]
        ]
    ]

updatetimescale[snpls_,treels_, theta_, epsilon_,nsq_,nbp_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,
    heat_,isadp_, adprate_, inputadpset_] :=
    Module[ {dim,props,proptreels,propsnpls,ls,map,trees,ar,adpset = inputadpset},
        dim = Length[Union[Flatten[treels[[All, 2, 1, nsq + 1 ;;]]]]];
        props = RandomReal[NormalDistribution[0, Last[adpset]]];
        proptreels = treels;
        proptreels[[All, 2, 1, nsq + 1 ;;]]*= Exp[props];
        proptreels[[All, 3]] *= Exp[props];
        proptreels[[All, -1]] = Prepend[Log[TreeTransitionProb @@ # & /@
            Partition[proptreels[[All, 2]], 2, 1]],TreeLogPriorProb[proptreels[[1, 2]]]];
        propsnpls = snpls;
        If[ propsnpls==={},
            ar = 0,
            ls = Append[proptreels, ReplacePart[proptreels[[-1]], 1 -> nbp + 1]];
            map = IndexByInterpolation[ls[[All, 1]]];
            trees = ls[[map[propsnpls[[All, 1]]], 2]];
            propsnpls[[All, -1]] = TreeLogLikelihood[propsnpls[[All, 2]], trees, propsnpls[[All, 3]],theta, epsilon,modelid];
            ar = Total[propsnpls[[All, -1]]] - Total[snpls[[All, -1]]];
        ];
        ar += (Total[proptreels[[All, -1]]] + Log[Times @@ proptreels[[;; -2, 3]]])-
                (Total[treels[[All, -1]]] + Log[Times @@ treels[[;; -2, 3]]]);
        ar+= logpriSNP[snpls[[All, 1]],proptreels, 1,Length[proptreels], theta, epsilon, nbp,modelid]-
        		logpriSNP[snpls[[All, 1]],treels, 1,Length[treels], theta, epsilon, nbp,modelid];
        ar = heat ar+logpritree[proptreels, nbp, theta, prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                logpritree[treels, nbp, theta, prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
        ar = Min[1, Exp[ar + dim props]];
        adpset = AdaptiveProposalScale[isadp, ar, adprate, adpset];
        If[ RandomReal[] < ar,
            {propsnpls, proptreels, adpset},
            {snpls, treels, adpset}
        ]
    ]    
       
updateshift[snpls_, treels_, theta_, epsilon_, nsq_, nbp_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,
  heat_, isadp_, adprate_, inputadpset_] :=
    Module[ {props, proptreels, propsnpls, ls, map, trees, ar,adpset = inputadpset},
        props = RandomReal[NormalDistribution[0, Last[adpset]]];
        If[ Min[treels[[All, 2, 1, nsq + 1]]] + props <= 0,
            ar = 0,
            proptreels = treels;
            proptreels[[All, 2, 1, nsq + 1 ;;]] += props;
            proptreels[[All, 3]] += nsq props;
            proptreels[[All, -1]] = Prepend[Log[TreeTransitionProb @@ # & /@ 
                Partition[proptreels[[All, 2]], 2, 1]],TreeLogPriorProb[proptreels[[1, 2]]]];
            propsnpls = snpls;
            If[ propsnpls === {},
                ar = 0,
                ls = Append[proptreels, ReplacePart[proptreels[[-1]], 1 -> nbp + 1]];
                map = IndexByInterpolation[ls[[All, 1]]];
                trees = ls[[map[propsnpls[[All, 1]]], 2]];
                propsnpls[[All, -1]] = TreeLogLikelihood[propsnpls[[All, 2]], trees,propsnpls[[All, 3]], theta, epsilon,modelid];
                ar = Total[propsnpls[[All, -1]]] - Total[snpls[[All, -1]]];
            ];
            ar += (Total[proptreels[[All, -1]]] + Log[Times @@ proptreels[[;; -2, 3]]])-
                  (Total[treels[[All, -1]]] + Log[Times @@ treels[[;; -2, 3]]]);
            ar+= logpriSNP[snpls[[All, 1]],proptreels, 1,Length[proptreels], theta, epsilon, nbp,modelid]-
            		logpriSNP[snpls[[All, 1]],treels, 1,Length[treels], theta, epsilon, nbp,modelid];
            ar = heat ar+logpritree[proptreels, nbp,theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                 logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
            ar = Min[1, Exp[ar]];
        ];
        adpset = AdaptiveProposalScale[isadp, ar, adprate, adpset];
        If[ RandomReal[] < ar,
            {propsnpls, proptreels, adpset},
            {snpls, treels, adpset}
        ]
    ]

updatehyperscale[inputsnpls_, inputtreels_, inputtheta_, 
   inputepsilon_, nsq_,nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,prithetashape_,
   prithetascale_,prithetamin_,prithetamax_, priepsalpha_, priepsbeta_,priepsmin_, priepsmax_,modelid_,
   heat_, isadpt_, adprate_, inputadpset_] :=
    Module[ {adpset = inputadpset, snpls = inputsnpls, 
      treels = inputtreels, theta = inputtheta, 
      epsilon = inputepsilon, logprihyper, dim, temp, map, 
      trpos, size = 1, deltxls, res, logmtbh, propx, proptreels, 
      propsnpls, ar,i},
        logprihyper = Function[{x},
            (prithetashape-1) Log[x[[2]]]-x[[2]]/prithetascale + 
            (priepsalpha - 1) Log[x[[3]]] + (priepsbeta - 1) Log[1 - x[[3]]]
          ];
        dim = Length[Union[Flatten[treels[[All, 2, 1, nsq + 1 ;;]]]]];
        temp = Append[treels, ReplacePart[treels[[-1]], 1 -> nbp + 1]];
        map = IndexByInterpolation[temp[[All, 1]]];
        trpos = map[snpls[[All, 1]]];
        deltxls = StudentTDisplacements[adpset[[-3]], adpset[[-1]]+10^(-10.) IdentityMatrix[3], size,
            SamplingPattern -> {"Cholesky", Infinity, False}];
        res = Table[0, {size}];
        logmtbh = Log[Mean[treels[[All, 3]]]];
        Do[
         propx = {logmtbh, theta, epsilon} + deltxls[[i]];
         If[ prithetamin<propx[[2]] <prithetamax &&priepsmin < propx[[-1]] < priepsmax,
             proptreels = treels;
             proptreels[[All, 2, 1, nsq + 1 ;;]] *= Exp[propx[[1]] - logmtbh];
             proptreels[[All, 3]] *= Exp[propx[[1]] - logmtbh];
             proptreels[[All, -1]] = Prepend[Log[TreeTransitionProb @@ # & /@Partition[proptreels[[All, 2]], 2, 1]], 
               TreeLogPriorProb[proptreels[[1, 2]]]];
             propsnpls = snpls;
             If[ propsnpls === {},
                 ar = 0,
                 propsnpls[[All, -1]] = TreeLogLikelihood[propsnpls[[All, 2]], proptreels[[trpos, 2]],
                     propsnpls[[All, 3]], propx[[2]], propx[[3]],modelid];
                 ar = Total[propsnpls[[All, -1]]] - Total[snpls[[All, -1]]];
             ];
             ar += (Total[proptreels[[All, -1]]] +Log[Times @@ proptreels[[;; -2, 3]]])-
                     (Total[treels[[All, -1]]] + Log[Times @@ treels[[;; -2, 3]]]);
             ar += logprihyper[propx] - logprihyper[{logmtbh, theta, epsilon}];
             ar+= logpriSNP[snpls[[All, 1]],proptreels, 1,Length[proptreels], propx[[2]],propx[[3]], nbp,modelid]-
             		logpriSNP[snpls[[All, 1]],treels, 1,Length[treels], theta, epsilon, nbp,modelid];
             ar = heat ar+logpritree[proptreels, nbp, propx[[2]],prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                     logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             ar = Min[1, Exp[ar + dim (propx[[1]] - logmtbh)]];
             If[ RandomReal[] < ar,
                 {theta, epsilon} = Rest[propx];
                 snpls = propsnpls;
                 treels = proptreels;
                 logmtbh = Log[Mean[proptreels[[All, 3]]]];
             ],
             ar = 0
         ];
         res[[i]] = {ar, {logmtbh, theta,epsilon}}, {i, size}];
        adpset = AdaptiveProposal[isadpt, Mean[res[[All, 1]]], adprate, adpset,res[[All, -1]]];
        {snpls, treels, {theta, epsilon}, adpset}
    ]    

updatehypershift[inputsnpls_, inputtreels_, inputtheta_, 
   inputepsilon_, nsq_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,prithetashape_,
   prithetascale_,prithetamin_,prithetamax_, priepsalpha_, priepsbeta_, priepsmin_,priepsmax_, modelid_,heat_, 
   isadpt_, adprate_, inputadpset_] :=
    Module[ {adpset = inputadpset, snpls = inputsnpls, 
      treels = inputtreels, theta = inputtheta, 
      epsilon = inputepsilon, size = 1, logprihyper, temp, map, 
      trpos, deltxls, res, mh, propx, proptreels, propsnpls, ar,i},
        logprihyper = Function[{x}, 
            (prithetashape-1) Log[x[[2]]]-x[[2]]/prithetascale + 
            (priepsalpha - 1) Log[x[[3]]] + (priepsbeta - 1) Log[1 - x[[3]]]
            ];
        temp = Append[treels, ReplacePart[treels[[-1]], 1 -> nbp + 1]];
        map = IndexByInterpolation[temp[[All, 1]]];
        trpos = map[snpls[[All, 1]]];
        deltxls = StudentTDisplacements[adpset[[-3]],adpset[[-1]] + 10^(-10.) IdentityMatrix[3], size, 
          SamplingPattern -> {"Cholesky", Infinity, False}];
        res = Table[0, {size}];
        mh = Mean[treels[[All,2,1,nsq + 1]]];
        Do[         
         propx = {mh, theta, epsilon} + deltxls[[i]];
         If[ propx[[1]] - mh + Min[treels[[All, 2, 1, nsq + 1]]] > 0 && 
             prithetamin<propx[[2]]<prithetamax &&priepsmin <propx[[-1]]<priepsmax,
             proptreels = treels;
             proptreels[[All, 2, 1, nsq + 1 ;;]] += propx[[1]] - mh;
             proptreels[[All, 3]] +=nsq (propx[[1]] - mh);
             proptreels[[All, -1]] = Prepend[Log[TreeTransitionProb @@ # & /@Partition[proptreels[[All, 2]], 2, 1]], 
               TreeLogPriorProb[proptreels[[1, 2]]]];
             propsnpls = snpls;
             If[ propsnpls === {},
                 ar = 0,
                 propsnpls[[All, -1]] = TreeLogLikelihood[propsnpls[[All, 2]], proptreels[[trpos, 2]],
                    propsnpls[[All, 3]], propx[[2]], propx[[3]],modelid];
                 ar = Total[propsnpls[[All, -1]]] - Total[snpls[[All, -1]]];
             ];
             ar += (Total[proptreels[[All, -1]]] + Log[Times @@ proptreels[[;; -2, 3]]])-
                    (Total[treels[[All, -1]]] + Log[Times @@ treels[[;; -2, 3]]]);
             ar += logprihyper[propx] - logprihyper[{mh, theta, epsilon}];
             ar+= logpriSNP[snpls[[All, 1]],proptreels, 1,Length[proptreels], propx[[2]],propx[[3]], nbp,modelid]-
             		logpriSNP[snpls[[All, 1]],treels, 1,Length[treels], theta, epsilon, nbp,modelid];
             ar = heat ar+logpritree[proptreels, nbp, propx[[2]],prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                     logpritree[treels, nbp,theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             ar = Min[1, Exp[ar]];
             If[ RandomReal[] < ar,
                 {theta, epsilon} = Rest[propx];
                 snpls = propsnpls;
                 treels = proptreels;
                 mh = Mean[treels[[All,2,1,nsq + 1]]];
             ],
             ar = 0
         ];
         res[[i]] = {ar, {mh, theta, epsilon}}, {i, size}];
        adpset = AdaptiveProposal[isadpt, Mean[res[[All, 1]]], adprate, adpset,res[[All, -1]]];
        {snpls, treels, {theta, epsilon}, adpset}
    ]    

segloglike[inputsnpls_, tree_, xstart_, xend_, theta_, epsilon_,modelid_] :=
    Module[ {pos, newsnpls = inputsnpls},
        pos = Flatten[Position[inputsnpls[[All, 1]], _?(xstart <= # < xend &), {1}, Heads -> False]];
        If[ pos === {},
            {newsnpls,0},
            newsnpls[[pos, -1]] = TreeLogLikelihood[newsnpls[[pos, 2]], Table[tree, {Length[pos]}], 
              newsnpls[[pos, 3]], theta, epsilon,modelid];
            {newsnpls,Total[newsnpls[[pos, -1]]] - Total[inputsnpls[[pos, -1]]]}
        ]
    ]
    
    
    
(*type=TreeTransitionProbII[tree1,tree3],and 
  type[[1]]=1;
  and tree2=!=tree1,tree2=!=tree3*)
caltype1ff[tree1_,type_] :=
    Module[ {ff, n, v11, v10, v12, v13,v32,v33,e00,opp},
        opp = type[[-1,1]];
        e00 = opp[[All,1]];
        ff = Partition[tree1[[1, Flatten[e00]]], 2];
        ff = (Total[ff[[All, 2]] - ff[[All, 1]]])/TotalBranchHeight[tree1];
        ff -= Total[TreeTransitionProb[tree1, {#, -1, #}] & /@ e00];
        If[ Length[e00] == 2,
            {ff, 0, 0,0},
            n = LeafNumber[tree1];
            {{v11, v10}} = e00;
            {v32,v33} = opp[[1,-1]];
            v12 = Complement[tree1[[2, v10 - n, All, 1]], {v11}][[1]];
            v13 = TreeNodeParent[tree1, v10];
            {ff, 
             TreeTransitionProb[tree1, {{v12, v10}, -1, {v11, v10}}],
             TreeTransitionProb[tree1, {{v12, v10}, -1, {v10, v13}}],
             If[ v33==v13,
                 TreeTransitionProb[tree1, {{v12, v10}, -1, {v32, v13}}],
                 0
             ]
            }
        ]
    ]    
    
(*type=TreeTransitionProbII[tree1,tree3]*)
caljumpprob[tree1_,tree2_,tree3_,type_] :=
    Module[ {ff},
        Switch[First[type],
            0,
            1,
            1,
            If[ SameCoalescentTreeQ[tree1, tree2] || SameCoalescentTreeQ[tree3, tree2],
                Log[p00],
                ff = caltype1ff[tree1,type];
                Log[TreeTransitionProb[tree1,tree2]] +Log[(1 - 2 p00)/Total[ff]]
            ],
            2,
            If[ Union[type[[-1, All, All, 2]]] == {{-1}},
                ff = If[ Length[type[[-1]]]==1,
                         {1},
                         TreeTransitionProb[tree1,#]&/@type[[-1,All,1]]
                     ];
                Log[TreeTransitionProb[tree1,tree2]]-Log[Total[ff]],
                -Log[Length[type[[-1]]]]
            ],
            _,
            Print["Wrong type in caljumprob! {tree1,tree2,tree3,type}=",{tree1,tree2,tree3,type}];
            0
        ]
    ]        
    
(*type = TreeTransitionTypeII2[tree1, tree2]*)    
randommidtree[tree1_, tree2_,type_] :=
    Module[ {ff,op,e00,n,v10,v11,v12,v13,v32,v33},
        Switch[First[type],
            0,
            {tree1,type,1},
            1,
            Switch[RandomChoice[{p00, p00, 1 - 2 p00} -> {1, 2, 3}],
                1, {tree1,type,1},
                2, {tree2,type,1},
                3,
                e00 = type[[-1, 1, All, 1]];
                ff = Partition[tree1[[1, Flatten[e00]]], 2];
                ff = (Total[ff[[All, 2]] - ff[[All, 1]]])/TotalBranchHeight[tree1];
                ff -= Total[TreeTransitionProb[tree1, {#, -1, #}] & /@ e00];
                ff = If[ Length[e00] == 2,
                         {ff, 0, 0,0},
                         n = LeafNumber[tree1];
                         {{v11, v10}} = e00;
                         {v32,v33} = type[[-1,1,1,-1]];
                         v12 = Complement[tree1[[2, v10 - n, All, 1]], {v11}][[1]];
                         v13 = TreeNodeParent[tree1, v10];
                         {ff, 
                          TreeTransitionProb[tree1, {{v12, v10}, -1, {v11, v10}}],
                          TreeTransitionProb[tree1, {{v12, v10}, -1, {v10, v13}}],
                          If[ v33==v13,
                              TreeTransitionProb[tree1, {{v12, v10}, -1, {v32, v13}}],
                              0
                          ]
                         }
                     ];
                {Switch[RandomChoice[ff->Range[4]],
                    1,
                    NestWhile[NextLocalTree[#, mSMCTransition[#, e00]]&,tree1,SameCoalescentTreeQ[tree1, #]&],
                    2,
                    NextLocalTree[tree1, mSMCTransition[tree1, {{v12, v10}, -1, {v11, v10}}]],
                    3,
                    NextLocalTree[tree1, mSMCTransition[tree1, {{v12, v10}, -1, {v10, v13}}]],
                    4,
                    NextLocalTree[tree1, mSMCTransition[tree1, {{v12, v10}, -1, {v32, v13}}]]                                       
                ],type,ff}
            ],
            2,
            Which[
                Union[type[[-1, All, All, 2]]] == {{-1}},
                ff = If[ Length[type[[-1]]]==1,
                         {1},
                         TreeTransitionProb[tree1,#]&/@type[[-1,All,1]]
                     ];
                op = RandomChoice[ff->type[[-1,All,1]]];
                {NextLocalTree[tree1, mSMCTransition[tree1, op]],type,ff},
                !MemberQ[Flatten[Union[type[[-1, All, All, 2]]]],-1],
                {NextLocalTree[tree1, RandomChoice[RandomChoice[type[[-1]]]]],type,1},
                True,
                Print["wrong type in randommidtree!","\n{tree1,tree2,type}=",{tree1,tree2,type}];
                Return[$Failed]
            ],
            True,
            Print["wrong type in randommidtree!","\n{tree1,tree2,type}=",{tree1,tree2,type}];
            Return[$Failed]
          ]
    ]    
        
randommidtree[tree1_, tree2_] :=
    randommidtree[tree1, tree2,TreeTransitionTypeII2[tree1, tree2]]

updatecoaltime[snpls_, treels_, theta_, epsilon_, nsq_, nbp_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_, modelid_,heat_] :=
    Module[ {ct,temp,minpos,maxpos,coalt2,propct,proptreels,propsnpls,ar,ar2},
        ct = RandomChoice[Union[Flatten[treels[[All, 2, 1, nsq + 1 ;;]]]]];
        temp = Flatten[Position[treels[[All, 2, 1,nsq;;]], _?(MemberQ[#, ct] &), {1},Heads -> False]];
        minpos = Min[temp];
        maxpos = Max[temp];
        coalt2 = Append[Union[Flatten[treels[[minpos ;; maxpos, 2, 1,nsq;;]]]], Infinity];
        propct = 
         If[ ct == coalt2[[-2]],
             RandomReal[ExponentialDistribution[1]] + coalt2[[-3]],
             temp = Position[coalt2, ct][[1, 1]];
             RandomReal[coalt2[[{temp - 1, temp + 1}]]]
         ];
        proptreels = treels;
        proptreels[[minpos ;; maxpos, 2]] = proptreels[[minpos ;; maxpos, 2]] /. {ct -> propct};
        proptreels[[minpos ;; maxpos, 3]] = TotalBranchHeight[proptreels[[minpos ;; maxpos, 2]]];
        temp = Log[TreeTransitionProb @@ # & /@Partition[proptreels[[minpos ;; maxpos, 2]], 2, 1]];
        proptreels[[minpos ;; maxpos, -1]] = Prepend[temp, 
          If[ minpos == 1,
              TreeLogPriorProb[proptreels[[1, 2]]],
              Log[TreeTransitionProb @@ proptreels[[minpos - 1 ;; minpos, 2]]]
          ]];
        If[ maxpos < Length[proptreels],
            proptreels[[maxpos + 1, -1]] = Log[TreeTransitionProb @@ proptreels[[maxpos ;; maxpos + 1, 2]]]
        ];
        {propsnpls, ar} = segloglike2[snpls, proptreels, minpos, maxpos, theta, epsilon,nbp,modelid];
        temp = Min[maxpos + 1, Length[treels]];
        ar += (Total[proptreels[[minpos ;; temp, -1]]]+Log[Times @@ proptreels[[Max[minpos - 1, 1] ;; temp - 1, 3]]])-
              (Total[treels[[minpos ;; temp, -1]]]+Log[Times @@ treels[[Max[minpos - 1, 1] ;; temp - 1, 3]]]);
        ar+= logpriSNP[snpls[[All, 1]],proptreels, minpos, maxpos, theta, epsilon, nbp,modelid]-
        		logpriSNP[snpls[[All, 1]],treels, minpos, maxpos, theta, epsilon, nbp,modelid];
        ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
            logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
        ar2 = If[ ct == coalt2[[-2]],
                  (-2 propct) - (-2 ct),
                  0
              ];
        ar = Min[1, Exp[ar - ar2]];
        If[ RandomReal[] < ar,
            {propsnpls, proptreels, ar},
            {snpls, treels, ar}
        ]
    ]
    
updatetreepos[snpls_,treels_, foc_,theta_, epsilon_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_] :=
    Module[ {ar,s1, s2, s3, props, proptreels,propsnpls,tree1,tree2},
        If[ Length[treels] == 1,
            ar = 0,
            {{s1, tree1}, {s2, tree2}} = treels[[foc - 1 ;; foc, {1, 2}]];
            s3 = If[ foc == Length[treels],
                     nbp + 1,
                     treels[[foc + 1, 1]]
                 ];
            props = RandomInteger[{s1 + 1, s3 - 1}];
            proptreels = treels;
            proptreels[[foc, 1]] = props;
            {propsnpls,ar} = If[ s2 < props,
                                 segloglike[snpls, tree1, s2, props, theta, epsilon,modelid],
                                 segloglike[snpls, tree2, props,s2, theta, epsilon,modelid]
                             ];
            ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
               logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
            ar = Min[1, Exp[ar]]
        ];
        If[ RandomReal[] < ar,
            {propsnpls,proptreels,ar},
            {snpls,treels,ar}
        ]
    ]

updatetreepostree[snpls_,treels_,foc_,theta_, epsilon_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_] :=
    Module[ {s1, tree1, proptr, proptreels,propsnpls,ar, ar2, s2, tree2, s3,tree3, type, ff,direct,xleft,xright,
        i1,i2,temp,ls,pos,map,trees,props2,i1m,i2m},
        (*Print["{foc,Length[treels]}=",{foc,Length[treels]}];*)
        Which[
         Length[treels] == 1,
         {s1, tree1} = treels[[1, {1, 2}]];
         proptr = Last[RandomNextTree[tree1]];
         proptreels = treels;
         proptreels[[1, 2 ;;]] = {proptr, TotalBranchHeight[proptr], TreeLogPriorProb[proptr]};
         {propsnpls,ar} = segloglike[snpls, proptr, s1, nbp + 1, theta, epsilon,modelid];
         ar = heat ar+ logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                 logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
         ar2 = Log[TreeTransitionProb[tree1, proptr]/(1-TreeTransitionProb[tree1])]-
                 Log[TreeTransitionProb[proptr, tree1]/(1-TreeTransitionProb[proptr])];
         ar = Min[1, Exp[ar - ar2]],
         foc == 1,
         {{s1, tree1}, {s2, tree2}} = treels[[foc ;; foc + 1, {1, 2}]];
         proptr = Last[RandomNextLocalTree[tree2]];
         If[ SameCoalescentTreeQ[tree1,proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr], TreeLogPriorProb[proptr]};
             proptreels[[foc + 1, -1]] = Log[TreeTransitionProb[proptr, tree2]];
             {propsnpls,ar} = segloglike[snpls, proptr, s1,s2,  theta, epsilon,modelid];
             ar += (Total[proptreels[[foc ;; foc + 1, -1]]] + Log[proptreels[[foc, 3]]])-
                  (Total[treels[[foc ;; foc + 1, -1]]] +Log[treels[[foc, 3]]]);
             ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                     logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             ar2 = Log[TreeTransitionProb[tree2, proptr]] -Log[TreeTransitionProb[tree2, tree1]];
             ar = Min[1, Exp[ar - ar2]]
         ],
         foc == Length[treels], 
         {{s1, tree1}, {s2, tree2}} = treels[[foc - 1 ;; foc, {1, 2}]];
         proptr = Last[RandomNextLocalTree[tree1]];
         props2 = RandomInteger[{s1+1,nbp}];
         proptreels = treels;
         proptreels[[foc]] = {props2,proptr, TotalBranchHeight[proptr], Log[TreeTransitionProb[tree1, proptr]]};
         propsnpls = snpls;
         If[ props2<=s2,
             If[ !SameCoalescentTreeQ[tree2,proptr],
                 {propsnpls,ar} = segloglike[propsnpls,proptr,props2,nbp+1,theta,epsilon,modelid],
                 {propsnpls,ar} = segloglike[propsnpls,proptr,props2,s2,theta,epsilon,modelid]
             ],
             {propsnpls,ar} = segloglike[propsnpls,tree1,s2,props2,theta,epsilon,modelid];
             If[ !SameCoalescentTreeQ[tree2,proptr],
                 {propsnpls,temp} = segloglike[propsnpls,proptr,props2,nbp+1,theta,epsilon,modelid];
                 ar+=temp;
             ];
         ];
         ar+=proptreels[[foc, -1]]-treels[[foc, -1]];
         ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                 logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
         ar2 = proptreels[[foc, -1]] - treels[[foc, -1]];
         ar = Min[1, Exp[ar - ar2]], 
         True, 
         If[ SameCoalescentTreeQ[treels[[foc - 1, 2]], treels[[foc + 1, 2]]],
             If[ RandomReal[] < 0.5,
                 i1 = NestWhile[# - 1 &, foc, # >= 1&&SameCoalescentTreeQ[treels[[foc, 2]], treels[[#, 2]]] &] + 1;
                 i2 = foc + 1,
                 i1 = foc;
                 i2 = NestWhile[# + 1 &,foc, # <= Length[treels]&&SameCoalescentTreeQ[treels[[foc, 2]], treels[[#, 2]]]&]
             ];
             Which[
               i1 == 1,
               temp = NestList[Last[RandomNextLocalTree[#]] &, treels[[i2, 2]], i2 - i1];
               ar2 = (Total[Log[TreeTransitionProb @@ ## & /@ Partition[temp, 2, 1]]]) -
                 (Total[Log[TreeTransitionProb @@ ## & /@Partition[treels[[i2 ;; i1 ;; -1, 2]], 2, 1]]]);
               proptr = Reverse[Rest[temp]];
               props2 = Sort[RandomSample[Range[treels[[i1,1]]+1,treels[[i2,1]]-1],i2-i1-1]];
               proptreels = treels;
               proptreels[[i1+1 ;; i2 - 1, 1]] = props2;
               proptreels[[i1 ;; i2 - 1, 2]] = proptr;
               proptreels[[i1 ;; i2 - 1, 3]] = TotalBranchHeight[proptr];
               proptreels[[i1, -1]] = TreeLogPriorProb[First[proptr]];
               proptreels[[i1 + 1 ;; i2 - 1, -1]] = Log[TreeTransitionProb @@ ## & /@ Partition[proptr, 2, 1]];
               proptreels[[i2, -1]] = Log[TreeTransitionProb[Last[proptr], proptreels[[i2, 2]]]],
               i2 == Length[treels] + 1,
               proptr = Rest[NestList[Last[RandomNextLocalTree[#]] &, treels[[i1 - 1, 2]],i2 - i1]];
               props2 = Sort[RandomSample[Range[treels[[i1-1,1]]+1,nbp],i2-i1]];
               proptreels = treels;
               proptreels[[i1 ;;, 1]] = props2;
               proptreels[[i1 ;;, 2]] = proptr;
               proptreels[[i1 ;;, 3]] = TotalBranchHeight[proptr];
               proptreels[[i1, -1]] = Log[TreeTransitionProb[proptreels[[i1 - 1, 2]], proptreels[[i1, 2]]]];
               proptreels[[i1 + 1 ;;, -1]] = Log[TreeTransitionProb @@ ## & /@ Partition[proptr, 2, 1]];
               ar2 = Total[proptreels[[i1 ;;, -1]]] - Total[treels[[i1 ;;, -1]]],
               True,
               direct = If[ RandomReal[] < 0.5,
                            1,
                            -1
                        ];
               If[ direct == 1,
                   temp = NestList[randommidtree[First[#], treels[[i2, 2]]] &, {treels[[i1 - 1, 2]], 0, 0}, i2 - i1];
                   proptr = temp[[2 ;;, 1]],
                   temp = NestList[randommidtree[First[#], treels[[i1 - 1, 2]]] &, {treels[[i2, 2]],0, 0}, i2 - i1];
                   proptr = Reverse[temp[[2 ;;, 1]]];
               ];
               props2 = Sort[RandomSample[Range[treels[[i1-1,1]]+1,treels[[i2,1]]-1],i2-i1]];
               proptreels = treels;
               proptreels[[i1 ;; i2 - 1, 1]] = props2;
               proptreels[[i1 ;; i2 - 1, 2]] = proptr;
               proptreels[[i1 ;; i2 - 1, 3]] = TotalBranchHeight[proptr];
               proptreels[[i1, -1]] = Log[TreeTransitionProb[proptreels[[i1 - 1, 2]], proptreels[[i1, 2]]]];
               proptreels[[i1 + 1 ;; i2 - 1, -1]] = Log[TreeTransitionProb @@ ## & /@ Partition[proptr, 2, 1]];
               proptreels[[i2, -1]] = Log[TreeTransitionProb[Last[proptr], proptreels[[i2, 2]]]];
               If[ direct == 1,
                   ls = Transpose[Join[Transpose[Partition[temp[[All, 1]], 2, 1]], {temp[[2 ;;, 3]], 
                       proptreels[[i1 ;; i2 - 1, -1]]}]];
                   ar2 = Total[
                      Which[
                        SameCoalescentTreeQ[#[[1]], #[[2]]]&&SameCoalescentTreeQ[treels[[i2, 2]], #[[2]]],
                        0,
                        SameCoalescentTreeQ[#[[1]], #[[2]]]||SameCoalescentTreeQ[treels[[i2, 2]], #[[2]]],
                        Log[p00],
                        True,
                        #[[-1]] + Log[(1 - 2 p00)/Total[#[[3]]]]
                      ] & /@ls];
                   ar2-= If[ SameCoalescentTreeQ[treels[[i1 - 1, 2]],treels[[i1, 2]]],
                             (i2 - i1) Log[p00],
                             Log[p00]
                         ],
                   ls = Partition[temp[[All, 1]], 2, 1];
                   ls = Transpose[Join[Transpose[ls], {temp[[2 ;;, 3]],Log[TreeTransitionProb @@ # & /@ ls]}]];
                   ar2 = Total[
                       Which[
                          SameCoalescentTreeQ[#[[1]], #[[2]]]&&SameCoalescentTreeQ[treels[[i1 - 1, 2]], #[[2]]],
                          0,
                          SameCoalescentTreeQ[#[[1]], #[[2]]]||SameCoalescentTreeQ[treels[[i1 - 1, 2]], #[[2]]],
                        Log[p00],
                        True,
                        #[[-1]] + Log[(1 - 2 p00)/Total[#[[3]]]]
                      ] & /@ ls];
                   ar2-= If[ SameCoalescentTreeQ[treels[[i2, 2]],treels[[i2-1, 2]]],
                             (i2 - i1) Log[p00],
                             Log[p00]
                         ];
               ]
             ];
             propsnpls = snpls;
             i1m = Max[1,i1-1];
             i2m = Min[i2, Length[proptreels]];
             xleft = treels[[i1m, 1]];
             xright = If[ i2 == Length[treels] + 1,
                          nbp + 1,
                          treels[[i2, 1]]
                      ];
             pos = Flatten[Position[propsnpls[[All, 1]], _?(xleft <= # < xright &), {1},Heads -> False]];
             ls = If[ i2m==Length[proptreels],
                      Append[proptreels[[i1m ;; i2m]],ReplacePart[proptreels[[-1]],1->nbp+1]],
                      proptreels[[i1m ;; i2m]]
                  ];
             (*Print[{i1m,i2m,i1,i2,xleft,xright,ls[[All,1]]}];*)
             map = IndexByInterpolation[ls[[All, 1]]];
             trees = ls[[map[propsnpls[[pos, 1]]], 2]];
             propsnpls[[pos, -1]] = TreeLogLikelihood[propsnpls[[pos, 2]], trees, propsnpls[[pos, 3]],theta, epsilon,modelid];
             ar = (Total[propsnpls[[pos, -1]]] - Total[snpls[[pos, -1]]]);
             ar += (Total[proptreels[[i1 ;; i2m, -1]]] + Log[Times @@proptreels[[i1m ;; i2m - 1, 3]]]) - 
                 (Total[treels[[i1 ;; i2m, -1]]] +Log[Times @@ treels[[i1m ;; i2m - 1, 3]]]);
             ar = heat ar + logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
                     logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             ar = Min[1, Exp[ar - ar2]],
             {{s1, tree1},{s2,tree2}, {s3, tree3}} = treels[[foc - 1;;foc + 1, {1, 2}]];
             props2 = RandomInteger[{s1+1,s3-1}];
             {proptr, type, ff} = randommidtree[tree1, tree3];
             proptreels = treels;
             proptreels[[foc]] = {props2,proptr, TotalBranchHeight[proptr],Log[TreeTransitionProb[tree1, proptr]]};
             proptreels[[foc + 1, -1]] = Log[TreeTransitionProb[proptr, tree3]];
             propsnpls = snpls;
             If[ props2<=s2,
                 If[ !SameCoalescentTreeQ[tree2,proptr],
                     {propsnpls,ar} = segloglike[propsnpls,proptr,props2,s3,theta,epsilon,modelid],
                     {propsnpls,ar} = segloglike[propsnpls,proptr,props2,s2,theta,epsilon,modelid]
                 ],
                 {propsnpls,ar} = segloglike[propsnpls,tree1,s2,props2,theta,epsilon,modelid];
                 If[ !SameCoalescentTreeQ[tree2,proptr],
                     {propsnpls,temp} = segloglike[propsnpls,proptr,props2,s3,theta,epsilon,modelid];
                     ar+=temp;
                 ];
             ];
             ar+=(Total[proptreels[[foc ;; foc + 1, -1]]] +Log[proptreels[[foc, 3]]])-
                 (Total[treels[[foc ;; foc + 1, -1]]] +Log[treels[[foc, 3]]]);
             ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                     logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             Switch[ First[type],                                   
                 0,
                 ar2 = 0,
                 2,
                 ar2 = If[ Union[type[[-1, All, All, 2]]] == {{-1}},
                           proptreels[[foc, -1]]-treels[[foc, -1]],
                           0
                       ],
                 1,
                 ar2 = If[ SameCoalescentTreeQ[tree1,proptr]||SameCoalescentTreeQ[tree3,proptr],
                           Log[p00],
                           proptreels[[foc, -1]] +Log[(1-2 p00)/Total[ff]]
                       ];
                 ar2-=If[ SameCoalescentTreeQ[tree1,tree2]||SameCoalescentTreeQ[tree3,tree2],
                          Log[p00],
                          If[ (SameCoalescentTreeQ[tree1,proptr]||SameCoalescentTreeQ[tree3,proptr]),
                              ff = caltype1ff[tree1,type]
                          ];
                          treels[[foc,-1]]+Log[(1- 2 p00)/Total[ff]]
                      ],
                 _,
                 Print["wrong type in updatetreetree!","\n{tree1,tree2,tree3,type13}=",{tree1,tree2,tree3,type}];
                 ar2 = 0
             ];
             ar = Min[1, Exp[ar - ar2]];
         ];
        ];
        If[ RandomReal[] < ar,
            {propsnpls,proptreels,ar},
            {snpls,treels,ar}
        ]
    ]
         
updatetreetree[snpls_,treels_,foc_,theta_, epsilon_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_] :=
    Module[ {s1, tree1, proptr, proptreels,propsnpls,ar, ar2, s2, tree2, s3,tree3, type, ff,direct,xleft,xright,i2min,
        i1,i2,temp,ls,pos,map,trees},
        Which[
         Length[treels] == 1,         
         tree1 = treels[[1,2]];
         {type,proptr} = randomjumptree[tree1];
         If[ SameCoalescentTreeQ[treels[[1, 2]], proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[1, 2 ;;]] = {proptr, TotalBranchHeight[proptr], TreeLogPriorProb[proptr]};
             {propsnpls, ar} = segloglike[snpls, proptr, treels[[1, 1]], nbp + 1, theta,epsilon,modelid];
             ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat] - 
               logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             ar2 = logpdfjumptree[tree1,type]-logpdfjumptree[proptr,TreeTransitionTypeI[proptr,tree1]];
             ar = Min[1, Exp[ar - ar2]]
         ],
         foc == 1,
         tree1 = treels[[foc+1,2]];
         {type,proptr} = randomjumptree[tree1];
         If[ SameCoalescentTreeQ[treels[[foc,2]],proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr], TreeLogPriorProb[proptr]};
             proptreels[[foc + 1, -1]] = Log[TreeTransitionProb[proptr, treels[[foc+1,2]]]];
             {propsnpls,ar} = segloglike[snpls, proptr, treels[[foc,1]],treels[[foc+1,1]],  theta, epsilon,modelid];
             ar += (Total[proptreels[[foc ;; foc + 1, -1]]] + Log[proptreels[[foc, 3]]])-
                  (Total[treels[[foc ;; foc + 1, -1]]] +Log[treels[[foc, 3]]]);
             ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                     logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             ar2 = logpdfjumptree[tree1,type]-logpdfjumptree[tree1,TreeTransitionTypeI[tree1,treels[[foc,2]]]];
             ar = Min[1, Exp[ar - ar2]]
         ], 
         foc == Length[treels], 
         tree1 = treels[[foc-1,2]];
         {type,proptr} = randomjumptree[tree1];
         If[ SameCoalescentTreeQ[treels[[foc, 2]], proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr], Log[TreeTransitionProb[treels[[foc - 1, 2]], proptr]]};
             {propsnpls, ar} = segloglike[snpls, proptr, treels[[foc, 1]], nbp + 1, theta, epsilon,modelid];
             ar += proptreels[[foc, -1]] - treels[[foc, -1]];
             ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
               logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
             ar2 = logpdfjumptree[tree1,type]-logpdfjumptree[tree1,TreeTransitionTypeI[tree1,treels[[foc,2]]]];
             ar = Min[1, Exp[ar - ar2]]
         ], 
         True, 
         If[ SameCoalescentTreeQ[treels[[foc - 1, 2]], treels[[foc + 1, 2]]],
             If[ RandomReal[] < 0.5,
                 i1 = NestWhile[# - 1 &, foc, # >= 1&&SameCoalescentTreeQ[treels[[foc, 2]], treels[[#, 2]]] &] + 1;
                 i2 = foc + 1,
                 i1 = foc;
                 i2 = NestWhile[# + 1 &,foc, # <= Length[treels]&&SameCoalescentTreeQ[treels[[foc, 2]], treels[[#, 2]]]&]
             ];
             Which[
               i1 == 1,
               temp = NestList[Last[RandomNextLocalTree[#]] &, treels[[i2, 2]], i2 - i1];
               ar2 = (Total[Log[TreeTransitionProb @@ ## & /@ Partition[temp, 2, 1]]]) -
                 (Total[Log[TreeTransitionProb @@ ## & /@Partition[treels[[i2 ;; i1 ;; -1, 2]], 2, 1]]]);
               proptr = Reverse[Rest[temp]];
               proptreels = treels;
               proptreels[[i1 ;; i2 - 1, 2]] = proptr;
               proptreels[[i1 ;; i2 - 1, 3]] = TotalBranchHeight[proptr];
               proptreels[[i1, -1]] = TreeLogPriorProb[First[proptr]];
               proptreels[[i1 + 1 ;; i2 - 1, -1]] = Log[TreeTransitionProb @@ ## & /@ Partition[proptr, 2, 1]];
               proptreels[[i2, -1]] = Log[TreeTransitionProb[Last[proptr], proptreels[[i2, 2]]]],
               i2 == Length[treels] + 1,
               proptr = Rest[NestList[Last[RandomNextLocalTree[#]] &, treels[[i1 - 1, 2]],i2 - i1]];
               proptreels = treels;
               proptreels[[i1 ;;, 2]] = proptr;
               proptreels[[i1 ;;, 3]] = TotalBranchHeight[proptr];
               proptreels[[i1, -1]] = Log[TreeTransitionProb[proptreels[[i1 - 1, 2]], proptreels[[i1, 2]]]];
               proptreels[[i1 + 1 ;;, -1]] = Log[TreeTransitionProb @@ ## & /@ Partition[proptr, 2, 1]];
               ar2 = Total[proptreels[[i1 ;;, -1]]] - Total[treels[[i1 ;;, -1]]],
               True,
               direct = If[ RandomReal[] < 0.5,
                            1,
                            -1
                        ];
               If[ direct == 1,
                   temp = NestList[randommidtree[First[#], treels[[i2, 2]]] &, {treels[[i1 - 1, 2]], 0, 0}, i2 - i1];
                   proptr = temp[[2 ;;, 1]],
                   temp = NestList[randommidtree[First[#], treels[[i1 - 1, 2]]] &, {treels[[i2, 2]],0, 0}, i2 - i1];
                   proptr = Reverse[temp[[2 ;;, 1]]];
               ];
               proptreels = treels;
               proptreels[[i1 ;; i2 - 1, 2]] = proptr;
               proptreels[[i1 ;; i2 - 1, 3]] = TotalBranchHeight[proptr];
               proptreels[[i1, -1]] = Log[TreeTransitionProb[proptreels[[i1 - 1, 2]], proptreels[[i1, 2]]]];
               proptreels[[i1 + 1 ;; i2 - 1, -1]] = Log[TreeTransitionProb @@ ## & /@ Partition[proptr, 2, 1]];
               proptreels[[i2, -1]] = Log[TreeTransitionProb[Last[proptr], proptreels[[i2, 2]]]];
               If[ direct == 1,
                   ls = Transpose[Join[Transpose[Partition[temp[[All, 1]], 2, 1]], {temp[[2 ;;, 3]], 
                       proptreels[[i1 ;; i2 - 1, -1]]}]];
                   ar2 = Total[
                      Which[
                          SameCoalescentTreeQ[#[[1]], #[[2]]]&&SameCoalescentTreeQ[treels[[i2, 2]], #[[2]]],
                          0,
                          SameCoalescentTreeQ[#[[1]], #[[2]]]||SameCoalescentTreeQ[treels[[i2, 2]], #[[2]]],
                        Log[p00],
                        True,
                        #[[-1]] + Log[(1 - 2 p00)/Total[#[[3]]]]
                      ] & /@ls];
                   ar2-= If[ SameCoalescentTreeQ[treels[[i1 - 1, 2]],treels[[i1, 2]]],
                             (i2 - i1) Log[p00],
                             Log[p00]
                         ],
                   ls = Partition[temp[[All, 1]], 2, 1];
                   ls = Transpose[Join[Transpose[ls], {temp[[2 ;;, 3]],Log[TreeTransitionProb @@ # & /@ ls]}]];
                   ar2 = Total[
                       Which[
                          SameCoalescentTreeQ[#[[1]], #[[2]]]&&SameCoalescentTreeQ[treels[[i1 - 1, 2]], #[[2]]],
                          0,
                          SameCoalescentTreeQ[#[[1]], #[[2]]]||SameCoalescentTreeQ[treels[[i1 - 1, 2]], #[[2]]],
                        Log[p00],
                        True,
                        #[[-1]] + Log[(1 - 2 p00)/Total[#[[3]]]]
                      ] & /@ ls];
                   ar2-= If[ SameCoalescentTreeQ[treels[[i2, 2]],treels[[i2-1, 2]]],
                             (i2 - i1) Log[p00],
                             Log[p00]
                         ];
               ]
             ];
             propsnpls = snpls;
             xleft = treels[[i1, 1]];
             xright = If[ i2 == Length[treels] + 1,
                          nbp + 1,
                          treels[[i2, 1]]
                      ];
             pos = Flatten[Position[propsnpls[[All, 1]], _?(xleft <= # < xright &), {1},Heads -> False]];
             i2min = Min[i2, Length[proptreels]];
             ls = If[ i2min==Length[proptreels],
                      Append[proptreels[[i1 ;; i2min]],ReplacePart[proptreels[[-1]],1->nbp+1]],
                      proptreels[[i1 ;; i2min]]
                  ];
             map = IndexByInterpolation[ls[[All, 1]]];
             trees = ls[[map[propsnpls[[pos, 1]]], 2]];
             propsnpls[[pos, -1]] = TreeLogLikelihood[propsnpls[[pos, 2]], trees, propsnpls[[pos, 3]],theta, epsilon,modelid];
             ar = (Total[propsnpls[[pos, -1]]] - Total[snpls[[pos, -1]]]);
             ar += (Total[proptreels[[i1 ;; i2min, -1]]] + Log[Times @@proptreels[[Max[i1 - 1, 1] ;; i2min - 1, 3]]]) - 
                 (Total[treels[[i1 ;; i2min, -1]]] +Log[Times @@ treels[[Max[i1 - 1, 1] ;; i2min - 1, 3]]]);
             ar = heat ar + logpritree[proptreels, nbp,theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
                     logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             ar = Min[1, Exp[ar - ar2]],
             {{s1, tree1},{s2,tree2}, {s3, tree3}} = treels[[foc - 1;;foc + 1, {1, 2}]];
             {proptr, type, ff} = randommidtree[tree1, tree3];
             If[ SameCoalescentTreeQ[tree2,proptr],
                 proptreels = treels;
                 propsnpls = snpls;
                 ar = 0,
                 proptreels = treels;
                 proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr],Log[TreeTransitionProb[tree1, proptr]]};
                 proptreels[[foc + 1, -1]] = Log[TreeTransitionProb[proptr, tree3]];
                 {propsnpls,ar} = segloglike[snpls, proptr, s2,s3,  theta, epsilon,modelid];
                 ar+=(Total[proptreels[[foc ;; foc + 1, -1]]] +Log[proptreels[[foc, 3]]])-
                     (Total[treels[[foc ;; foc + 1, -1]]] +Log[treels[[foc, 3]]]);
                 ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                         logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
             ];
             Switch[ First[type],                                   
                 0,
                 ar2 = 0,
                 2,
                 ar2 = If[ Union[type[[-1, All, All, 2]]] == {{-1}},
                           proptreels[[foc, -1]]-treels[[foc, -1]],
                           0
                       ],
                 1,
                 ar2 = If[ SameCoalescentTreeQ[tree1,proptr]||SameCoalescentTreeQ[tree3,proptr],
                           Log[p00],
                           proptreels[[foc, -1]] +Log[(1-2 p00)/Total[ff]]
                       ];
                 ar2-=If[ SameCoalescentTreeQ[tree1,tree2]||SameCoalescentTreeQ[tree3,tree2],
                          Log[p00],
                          If[ (SameCoalescentTreeQ[tree1,proptr]||SameCoalescentTreeQ[tree3,proptr]),
                              ff = caltype1ff[tree1,type]
                          ];
                          treels[[foc,-1]]+Log[(1- 2 p00)/Total[ff]]
                      ],
                 _,
                 Print["wrong type in updatetreetree!","\n{tree1,tree2,tree3,type13}=",{tree1,tree2,tree3,type}];
                 ar2 = 0
             ];
             ar = Min[1, Exp[ar - ar2]];
         ];
        ];
        If[ RandomReal[] < ar,
            {propsnpls,proptreels,ar},
            {snpls,treels,ar}
        ]
    ]
    

updatetreeinsert[snpls_,treels_,theta_, epsilon_,nbp_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_] :=
    Module[ {foc, s1, tree1, s2, tree2, ar, ar2, props, proptr, type, ff, proptreels,propsnpls},
        foc = RandomInteger[{1, Length[treels]}];
        {s1, tree1} = treels[[foc, {1, 2}]];
        If[ foc == Length[treels],
            If[ nbp - s1 < 1,
                ar = 0,
                props = RandomInteger[{s1 + 1, nbp}];
                proptr = Last[RandomNextLocalTree[tree1]];
                proptreels = Insert[treels, {props, proptr, TotalBranchHeight[proptr],Log[TreeTransitionProb[tree1, proptr]]}, foc + 1];
                {propsnpls,ar} = segloglike[snpls, proptr, props,nbp+1,  theta, epsilon,modelid];
                ar+=proptreels[[foc + 1, -1]] + Log[proptreels[[foc, 3]]];
                ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                    logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
                ar2 = proptreels[[foc + 1, -1]] +Log[1./(nbp - s1)];
                (*Print[{ar,ar2}];*)
                ar = Min[1, Exp[ar - ar2]]
            ],
            {s2, tree2} = treels[[foc + 1, {1, 2}]];
            If[ s2 - s1 < 2,
                ar = 0,
                props = RandomInteger[{s1 + 1, s2 - 1}];
                {proptr, type, ff} = randommidtree[tree1, tree2];
                Switch[First[type],
                    0,
                    (*proptr==tree1==tree2*)
                    proptreels = Insert[treels, {props, proptr, TotalBranchHeight[proptr],treels[[foc+1,-1]]}, foc + 1];
                    propsnpls = snpls;
                    ar = (Total[proptreels[[foc + 1 ;; foc + 2, -1]]] +Total[Log[proptreels[[foc ;; foc + 1, 3]]]])-
                        (treels[[foc + 1, -1]] + Log[treels[[foc, 3]]]);
                    ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                        logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
                    ar2 = Log[1./(s2 - s1 - 1)];
                    ar = Min[1, Exp[ar - ar2]],
                    1,
                    proptreels = Insert[treels, {props, proptr, TotalBranchHeight[proptr],Log[TreeTransitionProb[tree1, proptr]]}, foc + 1];
                    proptreels[[foc + 2, -1]] = Log[TreeTransitionProb[proptr, tree2]];
                    {propsnpls,ar} = segloglike[snpls, proptr, props,s2,  theta, epsilon,modelid];
                    ar+=(Total[proptreels[[foc + 1 ;; foc + 2, -1]]] +Total[Log[proptreels[[foc ;; foc + 1, 3]]]])-
                        (treels[[foc + 1, -1]] + Log[treels[[foc, 3]]]);
                    ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                        logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
                    ar2 = If[ SameCoalescentTreeQ[tree1,proptr]||SameCoalescentTreeQ[tree2,proptr],
                              Log[p00],
                              proptreels[[foc+1, -1]] +Log[(1-2 p00)/Total[ff]]
                          ]+Log[1./(s2 - s1 - 1)];
                    ar = Min[1, Exp[ar - ar2]],
                    2,
                    ar = 0,
                    _,
                    Print["wrong type in updatetreeinsert!","\n{tree1,tree2,type12}=",{tree1,tree2,type}];
                    ar = 0
                ]
            ];
        ];
        (*Print[{ar,type}];*)
        If[ RandomReal[] < ar,
            {propsnpls,proptreels,ar},
            {snpls,treels,ar}
        ]
    ]
   
updatetreedelete[snpls_,treels_, theta_, epsilon_, nbp_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_] :=
    Module[ {ar,ar2,foc,s1,s2,s3,tree1,tree2,tree3,type,proptreels,propsnpls,ff},
        If[ Length[treels] == 1,
            ar = 0,
            foc = RandomInteger[{2, Length[treels]}]; 
            (*Print["{foc,len}=",{foc,Length[treels]}];*)
            If[ foc == Length[treels],
                {s1, s2} = treels[[foc - 1 ;; foc, 1]];
                tree1 = treels[[foc-1,2]];
                proptreels = Delete[treels, foc];
                {propsnpls,ar} = segloglike[snpls, tree1, s2,nbp+1,  theta, epsilon,modelid];
                ar+=0-(treels[[foc, -1]] +Log[treels[[foc - 1, 3]]]);
                ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                    logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
                ar2 = 0 - (treels[[foc, -1]] +Log[1./(nbp - s1)]);
                ar = Min[1, Exp[ar - ar2]],
                {{s1, tree1}, {s2, tree2}, {s3, tree3}} = treels[[foc - 1 ;; foc + 1, {1, 2}]];
                type = TreeTransitionTypeII2[tree1, tree3];
                (*Print["type=",type];*)
                Switch[ First[type],
                    0,
                    If[ !SameCoalescentTreeQ[tree2,tree1],
                        Print["Wrong in delete trees!","\n{tree1,tree2,tree3}=",{tree1,tree2,tree3}];
                        Return[$Failed],
                        proptreels = Delete[treels, foc];
                        propsnpls = snpls;
                        ar = (proptreels[[foc, -1]] +Log[proptreels[[foc - 1, 3]]])-
                        (Total[treels[[foc ;; foc + 1, -1]]] +Total[Log[treels[[foc - 1 ;; foc, 3]]]]);
                        ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                             logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
                        ar2 = 0-Log[1./(s3 - s1 - 1)];
                        ar = Min[1, Exp[ar - ar2]]
                    ],
                    1,
                    proptreels = Delete[treels, foc];
                    proptreels[[foc, -1]] = Log[TreeTransitionProb[tree1, tree3]];
                    {propsnpls,ar} = segloglike[snpls, tree1, s2,s3, theta, epsilon,modelid];
                    ar+=(proptreels[[foc, -1]] +Log[proptreels[[foc - 1, 3]]])-
                        (Total[treels[[foc ;; foc + 1, -1]]] +Total[Log[treels[[foc - 1 ;; foc, 3]]]]);
                    ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                    logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
                    ar2 = 0-(If[ SameCoalescentTreeQ[tree2,tree1]||SameCoalescentTreeQ[tree2,tree3],
                                 Log[p00],
                                 ff = caltype1ff[tree1, type];
                                 treels[[foc, -1]] +Log[(1-2 p00)/Total[ff]]
                             ]+Log[1./(s3 - s1 - 1)]);
                    ar = Min[1, Exp[ar - ar2]],
                    2,
                    ar = 0,
                    _,
                    Print["wrong type in updatetreeinsert!","\n{tree1,tree2,tree3,type13}=",{tree1,tree2,tree3,type}];
                    ar = 0
                ];
            ];
        ];
        (*Print[{ar,foc,type}];*)
        If[ RandomReal[] < ar,
            {propsnpls,proptreels,ar},
            {snpls,treels,ar}
        ]
    ]
    

updatetreeswapch[inputsnpls_, inputtreels_, theta_, epsilon_, nsq_, nbp_, modelid_,heat_] :=
    Module[ {ind = {0, 10^(-6.)},snpls = inputsnpls,treels = inputtreels,size, ls, trpos,pair,segpos,i1,i2,seg,
        proptreels,propsnpls,ar,xleft,xright,pos,i},
        size = Ceiling[nsq/2.];
        ls = Append[treels[[All, 1]], nbp + 1];
        trpos = IndexByInterpolation[ls][snpls[[All, 1]]];
        Do[
         pair = Sort[RandomSample[Range[nsq], 2]];
         ls = Select[#, Max[#] <= nsq &] & /@ treels[[All, 2, 2, All, All, 1]];
         ls = Map[Sort, ls, 2];
         segpos = Flatten[Position[ls, _?(! MemberQ[#, pair] &), {1},Heads -> False]];
         segpos = Split[segpos, (#2 - #1 == 1 &)];
         segpos = Transpose[{segpos[[All, 1]], segpos[[All, -1]]}];         
         If[ Length[segpos]==0,
             ar = 1;
             ind+={ar,1};
             {propsnpls, proptreels} = {snpls, treels}
         ];
         Do[
          {i1, i2} = seg;
          proptreels = treels;
          proptreels[[i1 ;; i2, 2, 2, All, All, 1]] = proptreels[[i1 ;; i2, 2, 2, All, All, 1]] /.Thread[pair -> Reverse[pair]];
          propsnpls = snpls;
          If[ propsnpls === {},
              ar = 1,
              (*Print["pair=",pair,";snpls=",snpls];*)
              xleft = treels[[i1, 1]];
              xright = If[ i2 == Length[treels],
                           nbp + 1,
                           treels[[i2 + 1, 1]]
                       ];
              pos = Flatten[Position[propsnpls[[All, 1]], _?(xleft <= # < xright &), {1}, Heads -> False]];
              If[(!(pair[[1]]<=10&&pair[[2]]>10)),
              	(*snp acertainment for n=40 depends on #ma among initial 10 chromosomes, asysmetric*)
              	pos = Pick[pos, Unequal @@ # & /@ propsnpls[[pos, 2, pair]]]
              	];
              If[ pos === {},
                  ar = 1,
                  propsnpls[[pos, -1]] = TreeLogLikelihood[propsnpls[[pos, 2]],proptreels[[trpos[[pos]], 2]], 
                      propsnpls[[pos, 3]], theta,epsilon,modelid];                  
                  ar = Total[propsnpls[[pos, -1]]] - Total[snpls[[pos, -1]]];
                  ar = Min[1, Exp[heat ar]];
              ];              
          ];         
          ind+={ar, 1};
          If[ RandomReal[] < ar,
              {snpls, treels} = {propsnpls, proptreels}
          ], {seg,segpos}], {i,size}];
        If[ (Divide @@ ind)==0,
            Print["{size,ind,ar}3=",{size,ind,ar}]
        ];
        {snpls, treels, Divide @@ ind}
    ]
    
updatetreechecker[inputsnpls_, inputtreels_, theta_, epsilon_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_, modelid_,heat_] :=
    Module[ {snpls = inputsnpls, treels = inputtreels, ind = {0, 10^(-6.)}, ar, ar2, rho, rg, props, proptr, proptreels,
        trees, pos, logl, tree1, tree2, tree3,tree2p, temp, ap,start,i,types,ff,ff0},
        {snpls, treels, ar} = updatetreepostree[snpls, treels, 1, theta, epsilon, nbp, prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
        ind += {ar, 1};
        {snpls, treels, ar} = updatetreepostree[snpls, treels, Length[treels], theta, epsilon, nbp,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
        ind += {ar, 1};
        Do[
            (*rho = updateRho[treels, nbp, theta,prilambdamin,prilambdamax heat];*)
            (*{rho,adprho}=updateRho[treels, rho, nbp, theta,prilambdamin,prilambdamax,heat,isadp, adprate, adprho, heat];*)
            (*update changepoint at rg, condinional on changepoints at rg-1 and rg+1*)
            rg = Range[start, Length[treels] - 1, 2];
            props = RandomInteger[#]&/@Transpose[{treels[[rg - 1, 1]] + 1, treels[[rg + 1, 1]] - 1}];
            {proptr, types, ff} = Transpose[MapThread[randommidtree, {treels[[rg - 1, 2]], treels[[rg + 1, 2]]}]];
            proptreels = treels;
            proptreels[[rg, 1]] = props;
            proptreels[[rg, 2]] = proptr;
            proptreels[[rg, 3]] = TotalBranchHeight[proptr];
            proptreels[[rg, -1]] = Log[MapThread[TreeTransitionProb, {proptreels[[rg - 1, 2]],proptreels[[rg, 2]]}]];
            proptreels[[rg + 1, -1]] = Log[MapThread[TreeTransitionProb, {proptreels[[rg, 2]],proptreels[[rg + 1, 2]]}]];
            trees = proptreels[[IndexByInterpolation[Append[proptreels[[All, 1]], nbp + 1]][snpls[[All, 1]]],2]];
            logl = TreeLogLikelihood[snpls[[All, 2]], trees,snpls[[All, 3]], theta, epsilon,modelid];
            pos = Table[Flatten[Position[snpls[[All,1]], _?(treels[[i - 1, 1]] <= # <treels[[i + 1, 1]] &), {1}, Heads -> False]], {i, rg}];
            ar = (Total[logl[[#]]] - Total[snpls[[#, -1]]]) & /@ pos;
            ar += proptreels[[rg, -1]] + proptreels[[rg + 1, -1]] +Log[proptreels[[rg,3]]];
            ar -= treels[[rg, -1]] + treels[[rg + 1, -1]] +Log[treels[[rg, 3]]];
            ar += -(proptreels[[rg, 1]] - proptreels[[rg - 1, 1]]) proptreels[[rg - 1, 3]] rho/2 - 
                    (proptreels[[rg + 1, 1]] -proptreels[[rg, 1]]) proptreels[[rg, 3]] rho/2;
            ar -= -(treels[[rg, 1]] - treels[[rg - 1, 1]]) treels[[rg - 1,3]] rho/2 - 
                    (treels[[rg + 1, 1]] -treels[[rg, 1]]) treels[[rg, 3]] rho/2;
            ar2 = Table[
              Switch[types[[i, 1]],
               0, 0,
               1,
               {tree1, tree2, tree3} = treels[[rg[[i]] - 1 ;; rg[[i]] + 1, 2]];
               tree2p = proptreels[[rg[[i]], 2]];
               temp = If[ SameCoalescentTreeQ[tree1, tree2p]||SameCoalescentTreeQ[tree3, tree2p],
                          Log[p00],
                          proptreels[[rg[[i]], -1]] + Log[(1 - 2 p00)/Total[ff[[i]]]]
                      ];
               temp -If[ SameCoalescentTreeQ[tree1, tree2]||SameCoalescentTreeQ[tree3, tree2],
                         Log[p00],
                         ff0 = If[ (SameCoalescentTreeQ[tree1,tree2p]||SameCoalescentTreeQ[tree3,tree2p]),
                                   caltype1ff[tree1,types[[i]]],
                                   ff[[i]]
                               ];
                         treels[[rg[[i]], -1]] +Log[(1 - 2 p00)/Total[ff0]]
                     ],
               2,
               If[ Union[types[[i, -1, All, All, 2]]] == {{-1}},
                   proptreels[[rg[[i]], -1]] - treels[[rg[[i]], -1]],
                   0
               ]
              ], {i, Length[rg]}];
            ar = Min[1, #] & /@ Exp[heat ar - ar2];
            ap = Flatten[Position[Thread[RandomReal[{0, 1}, Length[ar]] < ar], True]];
            If[ ap=!={},
                Do[treels[[rg[[ap]], i]] = proptreels[[rg[[ap]], i]], {i,Dimensions[treels][[2]]}];
                treels[[rg[[ap]] + 1, -1]] = proptreels[[rg[[ap]] + 1, -1]];
                Do[snpls[[pos[[i]], -1]] = logl[[pos[[i]]]], {i, ap}];
            ];
            ind += {Total[ar], Length[ar]}, {start, 2,Min[3, Length[treels]-1]}];
        {snpls, treels, Divide @@ ind}
    ]    

    
blockupdatetreepos[inputsnpls_, inputtreels_, theta_, epsilon_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_, modelid_,heat_] :=
    Module[ {accept = {0, 10^(-6.)}, snpls = inputsnpls,treels = inputtreels, 
        rg, props, proptreels, propsnpls, ar, ind,i, start},
        Do[
         rg = Range[start, Length[treels], 2];
         props = RandomInteger[# + {1, -1}] & /@Partition[Append[treels[[All, 1]], nbp + 1][[start - 1 ;; ;; 2]], 2, 1];
         proptreels = treels;
         proptreels[[rg, 1]] = props;
         propsnpls = snpls;
         ar = 0;
         Do[
          {propsnpls, ind} = If[ treels[[i, 1]] < proptreels[[i, 1]],
                                 segloglike[propsnpls, proptreels[[i - 1, 2]], treels[[i, 1]], proptreels[[i, 1]], theta, epsilon,modelid],
                                 segloglike[propsnpls, proptreels[[i, 2]], proptreels[[i, 1]], treels[[i, 1]], theta, epsilon,modelid]
                             ];
          ar+=ind, {i, rg}];
         ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
               logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
         ar = Min[1, Exp[ar]];
         accept += {ar, 1};
         If[ RandomReal[] < ar,
             snpls = propsnpls;
             treels = proptreels;
         ], {start, 2, Min[3, Length[treels]]}];
        {snpls, treels, Divide @@ accept}
    ]      
    
rjupdatetreels[inputsnpls_, inputtreels_, theta_, epsilon_, nsq_,nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,
    modelid_,heat_,size_,isadp_,adprate_, inputadpposseglen_] :=
    Module[ {adpposseglen = inputadpposseglen,snpls = inputsnpls,treels = inputtreels, ind = Table[{0, 10^(-6.)}, {10}],
        maxseglen,act,actls,ar,segar,segarlist,posseglen,cc,istest=False},
        cc = 0.35;
        actls = RandomChoice[{(1-2 cc)/3,(1-2 cc)/3,(1-2 cc)/3,cc,cc}->Range[5],size];        
        (*cc = 0.25;
        actls = RandomChoice[{(1-3 cc)/2,(1-3 cc)/2,cc,cc,cc}->Range[5],size];
        *)
        posseglen = Max[1,Round[Last[adpposseglen]]];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"rjstart: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];
        segarlist = Reap[Do[
            Check[
            (*maxseglen = Length[treels]+10;*)
            maxseglen =100;
            (*Print["bef; act=",act];*)
            Switch[act,          
               1,{snpls,treels, ar} = updatecoaltime[snpls, treels, theta, epsilon, nsq,nbp, prirhoshape, prirhoscale,prirhomin,prirhomax,modelid,heat],
               2,{snpls, treels, ar} = updatesegtreepos[snpls, treels, posseglen, theta, epsilon, nbp, prirhoshape, prirhoscale,prirhomin,prirhomax,modelid,heat],
               3,{snpls,treels, segar} = updatesegtreetree[snpls, treels, theta, epsilon, nbp, prirhoshape, prirhoscale,prirhomin,prirhomax,modelid,heat,maxseglen],
               4,{snpls,treels, segar} = updatesegtreeinsert[snpls,treels, theta, epsilon,nbp,prirhoshape, prirhoscale,prirhomin,prirhomax,modelid,heat,maxseglen],
               5,{snpls,treels, segar} = updatesegtreedelete[snpls,treels, theta, epsilon, nbp,prirhoshape, prirhoscale,prirhomin,prirhomax,modelid,heat,maxseglen]
             ];
            (*Print["after; act=",act];*)
            If[ MemberQ[{3,4,5},act],  
                (*Print["{act,segar}=",{act,segar}];*)
                Sow[segar];
                ar = segar[[2]]
            ];
            If[ istest,
                If[ snpls=!={},
                    testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"rj: wrong in snpls[[All,-1]]!! act="<>ToString[act]]
                ];
            ];
            ind[[act]]+={ar,1},Print["act=",act]],{act,actls}];
        ];
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"rjend: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];
        If[ ind[[2,2]]>=1,
            adpposseglen = AdaptiveProposalScale[isadp, Divide@@ind[[2]], adprate, adpposseglen];
            adpposseglen[[-1]] = Min[Max[1,Last[adpposseglen]],100];
        ];
        (*Print["swaping tree leaves.."];*)
        {snpls, treels, ar} = updatetreeswapch[snpls, treels, theta, epsilon, nsq, nbp, modelid,heat];
        ind[[10]]+={ar,1};
        If[ istest,
            If[ snpls=!={},
                testsnpls[snpls,treels,nbp,theta,epsilon,modelid,"swapch: wrong in snpls[[All,-1]]!!"<>ToString[heat]]
            ];
        ];
        segarlist = If[ segarlist[[2]]=!={},
                        Mean[#] & /@ SplitBy[SortBy[segarlist[[2, 1]], First], First],
                        {}
                    ];
        {snpls,treels,adpposseglen,segarlist,ind[[All, 1]]/ind[[All, 2]]}
    ]
(****************************************************************************)
(*tree---tree'=randomjumptree[tree]
type=TreeTransitionTypeI[tree,tree']
logpdfjumptree[tree, type] gives the log probability of proposing tree' from tree.
*)  

randomjumptree[tree_, e0 : {_Integer, _Integer}] :=
    Module[ {tlist, elist, n, k, t0, t1, e1, newtree, newtype},
        {tlist, elist} = List @@ tree;
        n = Length[elist] + 1;
        k = RandomInteger[{Max[n, e0[[1]]], e0[[2]] - 1}];
        t0 = RandomReal[tlist[[{k, k + 1}]]];
        {t1, e1} = mSMCCoalescent[tree, t0];
        newtree = NextLocalTree[tree, {e0, t1, e1}];
        newtype = Which[
            e0 == e1,{0, 0}, 
            e0[[2]] == e1[[2]], {1, {{e0, t1, e1}, {e1, t1, e0}}},
            e0[[2]] == e1[[1]], {1, {{e0, t1, e1}, {First[DeleteCases[tree[[2, e0[[2]] - n]], e0]], t1, e1}}},
            True, {1, {{e0, t1, e1}}}
        ];
        {newtype, newtree}
    ]
    
randomjumptree[tree_] :=
    Module[ {n,els, e0},
        n = LeafNumber[tree];
        els = Flatten[tree[[2]], 1];
        e0 = RandomChoice[els[[All, 2]] - Replace[els[[All, 1]], (x_ /; x < n) -> n, {1}] ->els];
        randomjumptree[tree, e0]
    ]

jumptreeprob[tree_] :=
    Module[ {tlist, elist, n, ksq, wsq, kwsq, xysq, p0, ed, cr, p, a},
        {tlist, elist} = List @@ tree;
        n = LeafNumber[tree];
        ksq = Range[n, 2, -1];
        wsq = Differences[Drop[tlist, n - 1]];
        kwsq = ksq wsq;
        xysq = (Exp[kwsq] - 1)/ksq;
        p0 = Total[(Exp[-kwsq] + kwsq - 1)/kwsq];
        ed = Flatten[elist, 1];
        ed[[All, 1]] = Replace[ed[[All, 1]], (x_ /; x < n) -> n, {1}];
        ed = 2 n - Pick[ed, ed[[All, 2]] - ed[[All, 1]] - 2, _?NonNegative];
        ed[[All, 2]] += 1;
        ed = Tally[n + 1 - ed];
        cr = SparseArray[Thread[ed[[All, 1]] -> ed[[All, 2]]], {n - 1, n - 1}];
        cr = UpperTriangularize[Accumulate[cr], 1];
        cr = Reverse[Transpose[cr]];
        cr = Transpose[Reverse[Accumulate[cr]]];
        cr = UpperTriangularize[cr, 1];
        p = cr UpperTriangularize[KroneckerProduct[xysq/wsq, xysq], 1];
        a = UpperTriangularize[Table[kwsq, {n - 1}]];
        a = Transpose[Accumulate[Transpose[a]]];
        p *= Exp[-a];
        (p0 + Total[p, 2])/((n + 2) (n - 1)/2)
    ]

jumptreeprob[tree_, transform : {{_Integer, _Integer}, _?NonNegative}] :=
    Module[ {n, tlist, v1, v2, tc, t2, k1, ii, wsq, ksq, kwsq, res, k2},
        n = LeafNumber[tree];
        tlist = tree[[1]];
        {{v1, v2}, tc} = transform;
        t2 = tlist[[v2]];
        If[ tc < tlist[[v1]],
            0,
               (*k1 is the number of lineages just after node v1*)
            k1 = If[ v1 <= n,
                     n,
                     2 n - v1
                 ];
            (*tlist keeps the coalescent times for the remaining k1 lineages,
            including t1 at v1*)
            tlist = Take[tlist, -k1];
            If[ tc <= t2,
                ii = Position[tlist, _?(# > tc &), {1}, 1, Heads -> False][[1, 1]];
                wsq = Reverse[Differences[Append[tlist[[;; ii - 1]], tc]]];
                ksq = Range[k1 - Length[wsq] + 1, k1];
                res = 
                 Exp[-Most[FoldList[Plus, 0, wsq ksq]]] (1 - Exp[-ksq wsq])/ksq;
                res = 
                 Total[res/ReplacePart[wsq, 1 -> tlist[[ii]] - tlist[[ii - 1]]]],
                (*k2 is the number of lineages just after node v2*)
                k2 = 2 n - v2;
                (*wsq is the time interval from t2 at v2 to t1 at v1*)
                wsq = Reverse[Differences[Drop[tlist, -(k2 - 1)]]];
                (*ksq is the number of lineages at each interval corresponding to  wsq*)
                ksq = Range[k1 - Length[wsq] + 1, k1];
                kwsq = ksq wsq;
                res =  Total[Exp[-ReplacePart[RotateRight[Accumulate[kwsq]], 1 -> 0]] (1 - Exp[-kwsq])/kwsq];
                wsq =  Reverse[Differences[Append[Select[Drop[tlist, k1 - k2], (# < tc &)], tc]]];
                ksq = Range[k2 - Length[wsq] + 1, k2];
                res *= Exp[-Total[wsq ksq]]
            ];
            res/((n + 2) (n - 1)/2)
        ]
    ]

jumptreeprob[tree_, e0 : {_Integer, _Integer}] :=
    Module[ {n, ksq, wsq, kwsq, xysq, p0, p, a},
        n = LeafNumber[tree];
        ksq = Range[2 n - Max[n, e0[[1]]], 2 n - e0[[2]] + 1, -1];
        wsq = Differences[tree[[1, Max[n, e0[[1]]] ;; e0[[2]]]]];
        kwsq = ksq wsq;
        xysq = (Exp[kwsq] - 1)/ksq;
        p0 = Total[(Exp[-kwsq] + kwsq - 1)/(ksq^2 wsq)];
        p = UpperTriangularize[KroneckerProduct[xysq/wsq, xysq], 1];
        a = UpperTriangularize[Table[kwsq, {Length[kwsq]}]];
        a = Transpose[Accumulate[Transpose[a]]];
        p *= Exp[-a];
        (p0 + Total[p, 2])/((n + 2) (n - 1)/2)
    ]
  

logpdfjumptree[tree_, type_] :=
    If[ First[type] == 0,
        Log[jumptreeprob[tree]],
        Log[Total[jumptreeprob[tree,#]&/@type[[-1,All,{1,2}]]]]
    ]        
    
segnumber[n_, els_] :=
    Total[els[[All, 2]]-Replace[els[[All, 1]], (x_ /; x < n) -> n, {1}]];          
    
(*tree1---tree2 along chromosomes;
  |       |
tree1---(to sample)*)
randomjumpcondtree[tree1_, tree2_] :=
    Module[ {n,type, e0ls, newtree, newtype, ar2,sharee0},
        n = LeafNumber[tree1];
        type = TreeTransitionTypeI[tree1, tree2];
        Check[
        e0ls = 
         If[ First[type] == 0,
             Flatten[tree1[[2]], 1],
             type[[-1, All, 1]]
         ];
        {newtype,newtree} = randomjumptree[tree1, RandomChoice[e0ls]];
        Switch[ {First[type],First[newtype]},
          {0,0},
          ar2 = 0,
          {0,1},
          sharee0 = newtype[[-1,All,1]];
          ar2 = logpdfjumptree[tree1,newtype];
          ar2-=Log[Total[jumptreeprob[tree1,#]&/@sharee0]]-Log[segnumber[n,sharee0]/((n + 2) (n - 1)/2)],
          {1,0},
          sharee0 = type[[-1,All,1]];
          ar2 = Log[Total[jumptreeprob[tree1,#]&/@sharee0]]-Log[segnumber[n,sharee0]/((n + 2) (n - 1)/2)];
          ar2-=logpdfjumptree[tree1,type],
          {1,1},          
          e0ls = Select[newtype[[-1, All, ;; 2]],MemberQ[type[[-1, All, 1]], First[#]]&];
          ar2 = Log[Total[jumptreeprob[tree1,#]&/@e0ls]]-Log[segnumber[n,type[[-1, All, 1]]]/((n + 2) (n - 1)/2)];
          e0ls = Select[type[[-1, All, ;; 2]],MemberQ[newtype[[-1, All, 1]], First[#]]&];
          ar2-=Log[Total[jumptreeprob[tree1,#]&/@e0ls]]-Log[segnumber[n,newtype[[-1, All, 1]]]/((n + 2) (n - 1)/2)];          
        ];
        ,Print["{tree1,tree2,type}=",{tree1,tree2,type}];{ar2, newtree}={0,tree2}];
        {ar2, newtree}
    ]
    

(****************************************************************************)    
(*
randomjumptree[tree_, e0 : {_Integer, _Integer}] := 
 Module[{tlist, elist, n, t0, t1, e1, newtree, 
   newtype}, {tlist, elist} = List @@ tree;
  n = Length[elist] + 1;
  t0 = tlist[[e0[[1]]]];
  {t1, e1} = mSMCCoalescent[tree, t0];
  newtree = NextLocalTree[tree, {e0, t1, e1}];
  newtype = 
   Which[e0 == e1, {0, 0}, 
    e0[[2]] == e1[[2]], {1, {{e0, t1, e1}, {e1, t1, e0}}}, 
    e0[[2]] == 
     e1[[1]], {1, {{e0, t1, 
       e1}, {First[DeleteCases[tree[[2, e0[[2]] - n]], e0]], t1, 
       e1}}}, True, {1, {{e0, t1, e1}}}];
  {newtype, newtree}]

randomjumptree[tree_] := 
 randomjumptree[tree, RandomChoice[RandomChoice[tree[[2]]]]]

jumptreeprob[tree_] := 
 Module[{n, tlist, elist, ksq, wsq, xysq, ed, cr, a}, 
  n = LeafNumber[tree];
  {tlist, elist} = List @@ tree;
  ksq = Range[n, 2, -1];
  wsq = Differences[Drop[tlist, n - 1]];
  xysq = (1 - Exp[-ksq wsq])/ksq;
  ed = Flatten[elist, 1];
  ed[[All, 1]] = Replace[ed[[All, 1]], (x_ /; x < n) -> n, {1}];
  ed = 2 n - ed;
  ed[[All, 2]] += 1;
  ed = n + 1 - ed;
  ed = Tally[ed];
  cr = SparseArray[Thread[ed[[All, 1]] -> ed[[All, 2]]], {n - 1, n - 1}];
  cr = Reverse[Transpose[cr]];
  cr = UpperTriangularize[Transpose[Reverse[Accumulate[cr]]]];
  cr = cr Table[xysq, {n - 1}];
  a = UpperTriangularize[Table[PadLeft[Most[ksq wsq], n - 1], {n - 1}], 1];
  a = Transpose[Accumulate[Transpose[a]]];
  cr *= Exp[-a];
  Total[cr, 2]/(2 n - 2)
  ]

(*Probabaiblity for op=[e0,-1,e0]*)
jumptreeprob[tree_, e0 : {_Integer, _Integer}] := 
 Module[{n, k0, wls, kls, xysq},
  n = LeafNumber[tree];
  k0 = 2 n - Max[e0[[1]], n];
  wls = Differences[tree[[1, Max[e0[[1]], n] ;; e0[[2]]]]];
  kls = Range[k0, k0 - Length[wls] + 1, -1];
  xysq = (1 - Exp[-kls wls])/kls;
  Total[Exp[-Accumulate[Prepend[Most[kls wls], 0]]] xysq]/(2 n - 2)
  ]

jumptreeprob[tree_,transform : {{_Integer, _Integer}, _?NonNegative}] := 
 Module[{n, t1, v, tls, kls},
  t1 = transform[[2]];
  If[t1 < tree[[1, transform[[1, 1]]]], 0,
   n = LeafNumber[tree];
   v = If[transform[[1, 1]] <= n, n, transform[[1, 1]]];
   tls = Differences[Append[Select[tree[[1, v ;;]], # < t1 &], t1]];
   kls = Range[2 n - v, 2 n - v - Length[tls] + 1, -1];
   Exp[-Total[tls kls]]/(2 n - 2)
   ]
  ]       

logpdfjumptree[tree_, type_] :=
    If[ First[type] == 0,
        Log[jumptreeprob[tree]],
        Log[Total[jumptreeprob[tree,#]&/@type[[-1,All,{1,2}]]]]
    ]        
    
(*tree1---tree2 along chromosomes;
  |       |
tree1---(to sample)*)
randomjumpcondtree[tree1_, tree2_] :=
    Module[ {n,type, e0ls, newtree, newtype, ar2,sharee0},
        n = LeafNumber[tree1];
        type = TreeTransitionTypeI[tree1, tree2];
        e0ls = 
         If[ First[type] == 0,
             Flatten[tree1[[2]], 1],
             type[[-1, All, 1]]
         ];
        {newtype,newtree} = randomjumptree[tree1, RandomChoice[e0ls]];
        Switch[ {First[type],First[newtype]},
          {0,0},
          ar2 = 0,
          {0,1},
          sharee0 = newtype[[-1,All,1]];
          ar2 = logpdfjumptree[tree1,newtype];
          ar2-=Log[Total[jumptreeprob[tree1,#]&/@sharee0]]-Log[Length[sharee0]/(2 n-2)],
          {1,0},
          sharee0 = type[[-1,All,1]];
          ar2 = Log[Total[jumptreeprob[tree1,#]&/@sharee0]]-Log[Length[sharee0]/(2 n-2)];
          ar2-=logpdfjumptree[tree1,type],
          {1,1},          
          e0ls = Select[newtype[[-1, All, ;; 2]],MemberQ[type[[-1, All, 1]], First[#]]&];
          ar2 = Log[Total[jumptreeprob[tree1,#]&/@e0ls]]-Log[Length[e0ls]/(2 n-2)];
          e0ls = Select[type[[-1, All, ;; 2]],MemberQ[newtype[[-1, All, 1]], First[#]]&];
          ar2-=Log[Total[jumptreeprob[tree1,#]&/@e0ls]]-Log[Length[e0ls]/(2 n-2)];          
        ];
        {ar2, newtree}
    ]
*)    
(****************************************************************************)    
   
(*tree1---tree2 along chromosomes;
   |       |
tree11---(to sample tree22)*)
(*find nodes in tree11 corresponding to v in tree1. n is number of leaves*)
mapnode[tree1_,v_, tree11_] :=
    Module[ {n,pos, vv},
        n = LeafNumber[tree1];
        If[ v <= n,
            v,
            pos = Position[tree11[[1]], tree1[[1, v]], {1}, 1, Heads -> False];
            If[ pos === {},
                vv = tree1[[2, v - n, All, 1]];
                If[ # <= n,
                    #,
                    Position[tree11[[1]], tree1[[1, #]], {1}, 1,Heads -> False]
                ] & /@ vv,
                pos
            ]
        ]
    ]
    
(*
tree11---tree1---tree2
type=TreeTransitionTypeI[tree1,tree2]*)
(*type11to1=TreeTransitionTypeI[tree11,tree1]
tree11set return list of {e0,e1}, so that |({e0,_,e1}tree11)-tree1|=1*)
gete0e1 = Function[{tree11,type11to1},
            Module[ {v11,v10,v13,v32,v33,v12,set},
                {{v11, v10}} = type11to1[[-1, All, 1]];
                {v32, v33} = type11to1[[-1, 1, -1]];
                v12 = Complement[tree11[[2, v10 - LeafNumber[tree11], All, 1]], {v11}][[1]];
                v13 = TreeNodeParent[tree11, v10];
                set = {{{v12, v10}, {v11, v10}}, {{v12, v10}, {v10, v13}}};
                If[ v13 == v33,
                    AppendTo[set,{{v12, v10}, {v32, v13}}]
                ];
                set
            ]
        ];
    
condTreeTransitionTypeII[tree11_, tree1_, tree2_, type_] :=
    Module[ {type2,n, t1,temp,set},
        type2 = TreeTransitionTypeII[tree11, tree2];
        n = LeafNumber[tree11];
        If[ type2[[2]] == 2,
            t1 = type[[-1, 1, 2]];
            temp = TreeTransitionTypeI[tree11, tree1];
            type2[[-1]] = Select[type2[[-1]], (#[[1, 2]]==t1)&&Intersection[#[[All, 1]], temp[[-1, All, 1]]] === {}&];
            If[ Length[temp[[-1]]] == 1,
                set = gete0e1[tree11,temp];
                type2[[-1]] = Select[type2[[-1]], Intersection[#[[All, {1, 3}]], set] === {}&];
            ]
        ];
        type2
    ]    
      
randomjumpcondtree[tree1_, tree2_, tree11_] :=
    Module[ {msg = "Special proposal trees in randomjumpcondtree!", type, type2,t1, 
      e0, e1, tree22, newtype, type22,logratio},
        type = TreeTransitionTypeI[tree1, tree2];
        Catch[
         Switch[First[type],
          0, {0, tree11},
          1,
          type2 = condTreeTransitionTypeII[tree11, tree1, tree2, type];
          Switch[type2[[2]],
           1, {0, tree2},
           2,
           If[ type2[[-1]] === {},
                  (*Print["1: {tree1,tree2,tree11}=", {tree1,tree2, tree11}];*)
               Throw["emptytype2"],
               {e0, t1, e1} = RandomChoice[RandomChoice[type2[[-1]]]];
               tree22 = NextLocalTree[tree11, {e0, t1, e1}];
               newtype = Which[
                 e0[[2]] == e1[[2]], {1, {{e0, t1, e1}, {e1, t1, e0}}},
                 e0[[2]] == e1[[1]], {1, {{e0, t1, e1}, {First[DeleteCases[tree11[[2, e0[[2]] - LeafNumber[tree1]]], e0]], t1, e1}}},
                 True, {1, {{e0, t1, e1}}}
               ];
               type22 = condTreeTransitionTypeII[tree1, tree11, tree22, newtype];
               If[ type22[[-1]] === {},
                   Print["2: {tree1,tree2,tree11,tree22}=", {tree1,tree2, tree11, tree22}];
                   Throw["emptytype22"],
                   If[ type22[[2]]!=2,
                       Print["{newtype,type22}=",{newtype,type22}];
                       Print["22: {tree1,tree2,tree11,tree22}=", {tree1,tree2, tree11, tree22}];
                   ];
                   If[ Intersection[type[[-1]], Flatten[type22[[-1]], 1]] === {},
                       Print["{newtype,type22}=",{newtype,type22}];
                       Print["222: {tree1,tree2,tree11,tree22}=", {tree1,tree2, tree11, tree22}];
                   ];
                   (*If[First[TreeTransitionTypeI[tree2, tree22]]!=1,                   
                           Print["{newtype,type22}=",{newtype,type22}];
                        Print["2222: {tree1,tree2,tree11,tree22}=", {tree1,tree2, tree11, tree22}];
                   ];*)
                   logratio=Log[1./Length[type2[[-1]]]] - Log[1./Length[type22[[-1]]]];
                   If[!NumericQ[logratio],
                   	Print["{logratio,tree1,tree2,tree11,tree22}=",{logratio,tree1,tree2,tree11,tree22}];
                   	Print["cond{type2,type22}=",{type2,type22}];
                   ];
                   {logratio, tree22}
               ]
           ],
           _,
           Print["3: "<>msg];
           Print["3: {tree1,tree2,tree11}=", {tree1,tree2, tree11}];
           Throw["errortype1"]
           ],
          _, 
          Print["4: "<>msg];
          Print["4: {tree1,tree2,tree11}=", {tree1,tree2, tree11}];
          Throw["errortype2"]
          ]
         ]
    ]

(*T1---T2--T3--... along chromosomes;
  |    |
T11---T22--T33--..*)
(*return {ar2,{T11,T22,T33,..}} if exists, otherwise {None,None} if failed*)
randomjumpsegcondtree[trees_, tree11_,maxseglen_] :=
    Module[ {proptrees, newtree,ar2 = 0,i,temp,logratio},
        proptrees = Table[0, {Length[trees]}];
        proptrees[[1]] = tree11;
        Catch[
         Do[          
          If[ SameCoalescentTreeQ[proptrees[[i,1]], trees[[i, 1]]],
              proptrees = Take[proptrees, i];
              Break[]
          ];
          If[ i>=maxseglen+1,
              Throw[{None,None}]
          ];
          temp = randomjumpcondtree[trees[[i, 1]], trees[[i + 1, 1]],proptrees[[i, 1]]];
          (*Print["proptree=",temp];*)
          If[ Head[temp]===String,
              Switch[temp,
                  "errortype1",
                  Print["i=",i];
                  Print["trees=",trees];
                  Print["proptrees=",proptrees],     
                  "errortype2",
                  Print["i=",i];
                  Print["trees=",trees];
                  Print["proptrees=",proptrees];                    
              ];
              Throw[{None,None}],
              {logratio, newtree} = temp;
              ar2+=logratio;
              proptrees[[i + 1]] = {newtree, TotalBranchHeight[newtree], Log[TreeTransitionProb[proptrees[[i, 1]],newtree]]};
          ], {i, Length[trees] - 1}];
         {ar2, proptrees}
         ]
    ]


(****************************************************************************)
    
(*update inputsnpls[[All,-1]] when treels[[istart;;iend]] were revised*)
segloglike2[inputsnpls_, treels_, istart_, iend_, theta_, epsilon_, nbp_,modelid_] :=
    Module[ {xstart, xend, pos, newsnpls = inputsnpls, trees,ls},
        xstart = treels[[istart, 1]];
        xend = If[ iend == Length[treels],
                   nbp + 1,
                   treels[[iend + 1, 1]]
               ];
        pos = Flatten[Position[inputsnpls[[All, 1]], _?(xstart <= # < xend &), {1},Heads -> False]];
        If[ pos === {},
            {newsnpls, 0},
            ls = If[ iend == Length[treels],
                     Append[treels[[istart ;;]],ReplacePart[treels[[-1]], 1 -> nbp + 1]],
                     treels[[istart ;; iend + 1]]
                 ];
            trees = ls[[IndexByInterpolation[ls[[All, 1]]][inputsnpls[[pos, 1]]], 2]];
            newsnpls[[pos, -1]] = TreeLogLikelihood[newsnpls[[pos, 2]], trees, newsnpls[[pos, 3]],theta, epsilon,modelid];
            {newsnpls, Total[newsnpls[[pos, -1]]] - Total[inputsnpls[[pos, -1]]]}
        ]
    ]    

randomsegtreeprop[1, snpls_, treels_, theta_, epsilon_, nbp_,maxseglen_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_] :=
    Module[ {tree1,type,proptr,ar2,proptreels,propsnpls,propsub,seglen,ar,ar2sub},
        tree1 = treels[[1,2]];
        {type,proptr} = randomjumptree[tree1];
        If[ SameCoalescentTreeQ[tree1, proptr],
            {1,1,snpls,treels},
            ar2 = logpdfjumptree[tree1,type]-logpdfjumptree[proptr,TreeTransitionTypeI[proptr,tree1]];
            proptreels = treels;
            proptreels[[1, 2 ;;]] = {proptr, TotalBranchHeight[proptr],TreeLogPriorProb[proptr]};
            {ar2sub,propsub} = randomjumpsegcondtree[treels[[All,2 ;;]], proptreels[[1, 2 ;;]],maxseglen];
            If[ ar2sub === None,
                {0,-1,snpls,treels},
                seglen = Length[propsub];
                proptreels[[;; seglen, 2 ;;]] = propsub;
                ar2 += ar2sub;
                {propsnpls, ar} = segloglike2[snpls, proptreels, 1, seglen, theta, epsilon,nbp,modelid];
                ar += (Total[proptreels[[;; seglen, -1]]] + Log[Times @@ proptreels[[;; seglen - 1, 3]]]) - 
                        (Total[treels[[;; seglen, -1]]] + Log[Times @@ treels[[;; seglen - 1, 3]]]);
                ar+= logpriSNP[snpls[[All,1]],proptreels, 1,seglen, theta, epsilon, nbp,modelid]-
                        logpriSNP[snpls[[All,1]],treels, 1,seglen, theta, epsilon, nbp,modelid];
                ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat]- 
                          logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
                ar = Min[1, Exp[ar - ar2]];
                {ar, seglen, propsnpls,proptreels}
            ]
        ]
    ]

randomsegtreeprop[foc_/;foc>=2, snpls_, treels_, theta_, epsilon_, nbp_,maxseglen_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_] :=
    Module[ {ii,proptr,ar2sub,propsub,proptreels,ar2sub2,propsub2,seglen,ar2,ar,propsnpls,maxpos},
        ii = foc-1;
        proptr = Last[randomjumptree[treels[[ii,2]]]];
        If[ SameCoalescentTreeQ[treels[[ii,2]], proptr],
            {1,1,snpls,treels},
            {ar2sub,propsub} = randomjumpsegcondtree[treels[[ii ;;,2 ;;]], {proptr,-1,-1},maxseglen];
            If[ ar2sub === None,
                {0,-1,snpls,treels},
                proptreels = Insert[treels, {-1, proptr, -1,-1}, ii + 1];
                proptreels[[ii + 1 ;; ii + Length[propsub], 2;;]] = propsub;
                {ar2sub2,propsub2} = randomjumpsegcondtree[proptreels[[ii+1 ;;,2 ;;]], treels[[ii,2;;]],maxseglen];
                If[ ar2sub2 === None,
                    {0,-1,snpls,treels},
                    ar2 = ar2sub + ar2sub2;
                    proptreels = Delete[proptreels,ii+1];
                    proptreels[[ii ;; ii + Length[propsub2] - 1, 2;;]] = propsub2;
                    seglen = Max[Length[#] & /@ {propsub, propsub2}];
                    maxpos = Min[ii+seglen,Length[treels]];
                    {propsnpls, ar} = segloglike2[snpls, proptreels, ii+1, maxpos, theta, epsilon,nbp,modelid];
                    ar += (Total[proptreels[[ii+1;; maxpos, -1]]] + Log[Times @@ proptreels[[ii+1 ;; maxpos - 1, 3]]]) - 
                            (Total[treels[[ii+1;; maxpos, -1]]] + Log[Times @@ treels[[ii+1 ;; maxpos - 1, 3]]]);
                    ar+= logpriSNP[snpls[[All,1]],proptreels, ii+1,maxpos, theta, epsilon, nbp,modelid]-
                            logpriSNP[snpls[[All,1]],treels, ii+1,maxpos, theta, epsilon, nbp,modelid];
                    ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
                      logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
                    ar = Min[1, Exp[ar - ar2]];
                    {ar, seglen, propsnpls,proptreels}
                ]
            ]
        ]
    ]
             
    
updatesegtreetree2[snpls_, treels_, theta_, epsilon_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_,maxseglen_] :=
    Module[ {foc,proptr, proptreels, propsnpls, seglen = 1, ar, ar2,type,tree1},
        foc = RandomChoice[ReplacePart[Table[1, {Length[treels]}], 1 -> LeafNumber[treels[[1,2]]]+1]->Range[Length[treels]]];
        (*Print["segtreetree;foc=",foc];*)
        Which[
         Length[treels] ==foc==1,         
         tree1 = treels[[foc,2]];
         {type,proptr} = randomjumptree[tree1];
         If[ SameCoalescentTreeQ[treels[[foc, 2]], proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr], TreeLogPriorProb[proptr]};
             {propsnpls, ar} = segloglike[snpls, proptr, treels[[foc, 1]], nbp + 1, theta,epsilon,modelid];
             ar+= proptreels[[foc, -1]] - treels[[foc, -1]];
             ar+= logpriSNP[snpls[[All,1]],proptreels, foc,foc, theta, epsilon, nbp,modelid]-
             		logpriSNP[snpls[[All,1]],treels, foc,foc, theta, epsilon, nbp,modelid];
             ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
               logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
             ar2 = logpdfjumptree[tree1,type]-logpdfjumptree[proptr,TreeTransitionTypeI[proptr,tree1]];
             If[ !(NumericQ[ar]&&NumericQ[ar2]),
                 Print["1: {ar,ar2}=",{ar,ar2}]
             ];
             ar = Min[1, Exp[ar - ar2]]
         ],
         foc == Length[treels],
         tree1 = treels[[foc-1,2]];
         {type,proptr} = randomjumptree[tree1];
         If[ SameCoalescentTreeQ[treels[[foc, 2]], proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr], Log[TreeTransitionProb[tree1, proptr]]};
             {propsnpls, ar} = segloglike[snpls, proptr, treels[[foc, 1]], nbp + 1, theta, epsilon,modelid];
             ar += proptreels[[foc, -1]] - treels[[foc, -1]];
             ar+= logpriSNP[snpls[[All,1]],proptreels, foc,foc, theta, epsilon, nbp,modelid]-
             		logpriSNP[snpls[[All,1]],treels, foc,foc, theta, epsilon, nbp,modelid];
             ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
               logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
             ar2 = logpdfjumptree[tree1,type]-logpdfjumptree[tree1,TreeTransitionTypeI[tree1,treels[[foc,2]]]];
             If[ !(NumericQ[ar]&&NumericQ[ar2]),
                 Print["2: {ar,ar2}=",{ar,ar2}]
             ];
             ar = Min[1, Exp[ar - ar2]]
         ],
         True,
         {ar,seglen,propsnpls,proptreels} = randomsegtreeprop[foc, snpls, treels, theta, epsilon,nbp,maxseglen,prirhoshape, prirhoscale,prirhomin,prirhomax,modelid,heat]
        ];
        If[ RandomReal[] < ar,
            {propsnpls,proptreels, {seglen, ar}},
            {snpls, treels, {seglen, ar}}
        ]
    ]


updatesegtreetree[snpls_, treels_, theta_, epsilon_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_,maxseglen_] :=
    Module[ {foc,proptr, proptreels, propsnpls, seglen = 1, ar, ar2,type,tree1,ar2sub,propsub, maxpos},
        foc = RandomChoice[ReplacePart[Table[1, {Length[treels]}], 1 -> LeafNumber[treels[[1,2]]]+1]->Range[Length[treels]]];
        Which[
         Length[treels] ==foc==1,         
         tree1 = treels[[foc,2]];
         {type,proptr} = randomjumptree[tree1];
         If[ SameCoalescentTreeQ[treels[[foc, 2]], proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr], TreeLogPriorProb[proptr]};
             {propsnpls, ar} = segloglike[snpls, proptr, treels[[foc, 1]], nbp + 1, theta,epsilon,modelid];
             ar+= proptreels[[foc, -1]] - treels[[foc, -1]];
             ar+= logpriSNP[snpls[[All,1]],proptreels, foc,foc, theta, epsilon, nbp,modelid]-
                   logpriSNP[snpls[[All,1]],treels, foc,foc, theta, epsilon, nbp,modelid];
             ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
               logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
             ar2 = logpdfjumptree[tree1,type]-logpdfjumptree[proptr,TreeTransitionTypeI[proptr,tree1]];
             If[ !(NumericQ[ar]&&NumericQ[ar2]),
                 Print["1: {ar,ar2}=",{ar,ar2}]
             ];
             ar = Min[1, Exp[ar - ar2]]
         ],
         foc == Length[treels],
         tree1 = treels[[foc-1,2]];
         {type,proptr} = randomjumptree[tree1];
         If[ SameCoalescentTreeQ[treels[[foc, 2]], proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr], Log[TreeTransitionProb[tree1, proptr]]};
             {propsnpls, ar} = segloglike[snpls, proptr, treels[[foc, 1]], nbp + 1, theta, epsilon,modelid];
             ar += proptreels[[foc, -1]] - treels[[foc, -1]];
             ar+= logpriSNP[snpls[[All,1]],proptreels, foc,foc, theta, epsilon, nbp,modelid]-
             		logpriSNP[snpls[[All,1]],treels, foc,foc, theta, epsilon, nbp,modelid];
             ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
               logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
             ar2 = logpdfjumptree[tree1,type]-logpdfjumptree[tree1,TreeTransitionTypeI[tree1,treels[[foc,2]]]];
             If[ !(NumericQ[ar]&&NumericQ[ar2]),
                 Print["2: {ar,ar2}=",{ar,ar2}]
             ];
             ar = Min[1, Exp[ar - ar2]]
         ],
         True,         
         If[ foc==1,
             tree1 = treels[[foc,2]];
             {type,proptr} = randomjumptree[tree1];
             ar2 = If[ SameCoalescentTreeQ[tree1,proptr],
                       0,
                       logpdfjumptree[tree1,type]-logpdfjumptree[proptr,TreeTransitionTypeI[proptr,tree1]]
                   ],
             {ar2,proptr} = randomjumpcondtree@@treels[[foc-1;;foc,2]];
         ];
         If[ !(NumericQ[ar2]),
             Print["2.5: {ar2}=",{ar2}]
         ];
         If[ SameCoalescentTreeQ[treels[[foc, 2]], proptr],
             proptreels = treels;
             propsnpls = snpls;
             ar = 1,
             proptreels = treels;
             proptreels[[foc, 2 ;;]] = {proptr, TotalBranchHeight[proptr],
                 If[ foc == 1,
                     TreeLogPriorProb[proptr],
                     Log[TreeTransitionProb[treels[[foc - 1, 2]], proptr]]
                 ]};
             {ar2sub,propsub} = randomjumpsegcondtree[treels[[foc ;;,2 ;;]], proptreels[[foc, 2 ;;]],maxseglen];
             If[ ar2sub === None,
                 seglen = -1;
                 ar = 0,
                 seglen = Length[propsub];
                 maxpos = foc + seglen - 1;
                 proptreels[[foc ;; maxpos, 2 ;;]] = propsub;
                 ar2 += ar2sub;
                 {propsnpls, ar} = segloglike2[snpls, proptreels, foc, maxpos, theta, epsilon,nbp,modelid];
                 ar += (Total[proptreels[[foc ;; maxpos, -1]]] + Log[Times @@ proptreels[[foc ;; maxpos - 1, 3]]]) - 
                         (Total[treels[[foc ;; maxpos, -1]]] + Log[Times @@ treels[[foc ;; maxpos - 1, 3]]]);
                 ar+= logpriSNP[snpls[[All,1]],proptreels, foc,maxpos, theta, epsilon, nbp,modelid]-
                 		logpriSNP[snpls[[All,1]],treels, foc,maxpos, theta, epsilon, nbp,modelid];
                 ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
                   logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
                 If[ !(NumericQ[ar]&&NumericQ[ar2]),
                     Print["3: {ar,ar2}=",{ar,ar2}]
                 ];
                 ar = Min[1, Exp[ar - ar2]]
             ];
         ];
         ];
        If[ RandomReal[] < ar,
            {propsnpls,proptreels, {seglen, ar}},
            {snpls, treels, {seglen, ar}}
        ]
    ]


updatesegtreeinsert[inputsnpls_, inputtreels_, theta_, epsilon_, nbp_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_, modelid_,heat_,maxseglen_] :=
    Module[ {snpls = inputsnpls, treels = inputtreels,seglen = 1,foc, proptr, seg0,type, 
      ar, props, proptreels, ar2, propsnpls,ar2sub,propsub, maxpos},
        foc = RandomInteger[{1, Length[treels]}];
        (*Print["foc=",foc];*)
        seg0 = {treels[[foc, 1]] + 1, If[ foc == Length[treels],
                                          nbp,
                                          treels[[foc + 1, 1]] - 1
                                      ]};
        props = RandomInteger[seg0];
        Which[
            Greater @@ seg0,
            ar = 0,            
            foc == Length[treels],
            {type,proptr} = randomjumptree[treels[[foc,2]]];
            proptreels = Insert[treels, {props, proptr, TotalBranchHeight[proptr],
                    Log[TreeTransitionProb[treels[[foc,2]], proptr]]}, foc + 1];
            {propsnpls, ar} = segloglike[snpls, proptr, props, nbp + 1, theta, epsilon,modelid];
            ar += proptreels[[foc + 1, -1]] + Log[proptreels[[foc, 3]]];
            ar+= logpriSNP[snpls[[All,1]],proptreels, foc,foc+1, theta, epsilon, nbp,modelid]-
            		logpriSNP[snpls[[All,1]],treels, foc,foc, theta, epsilon, nbp,modelid];
            ar = heat ar + logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
              logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
            ar2 = logpdfjumptree[treels[[foc,2]],type]+Log[1./(seg0[[2]]-seg0[[1]]+1)];
            If[ !(NumericQ[ar]&&NumericQ[ar2]),
                Print["4: {ar,ar2}=",{ar,ar2}]
            ];
            ar = Min[1, Exp[ar - ar2]],
            True,
            {type,proptr} = randomjumptree[treels[[foc,2]]];
            ar2 = logpdfjumptree[treels[[foc,2]],type]+Log[1./(seg0[[2]]-seg0[[1]]+1)];
            proptreels = Insert[treels, {props, proptr, TotalBranchHeight[proptr],
                Log[TreeTransitionProb[treels[[foc, 2]], proptr]]}, foc + 1];
            {ar2sub,propsub} = randomjumpsegcondtree[treels[[foc ;;,2 ;;]],    proptreels[[foc + 1, 2 ;;]],maxseglen];
            If[ ar2sub === None,
                seglen = -1;
                ar = 0,
                seglen = Length[propsub];
                ar2+=ar2sub;
                maxpos = foc + seglen;
                proptreels[[foc + 1 ;;maxpos, 2 ;;]] = propsub;
                {propsnpls, ar} = segloglike2[snpls, proptreels, foc+1,  maxpos, theta, epsilon,nbp,modelid];
                ar += (Total[proptreels[[foc+1;; maxpos, -1]]] + Log[Times @@ proptreels[[foc;; maxpos - 1, 3]]]) - 
                    (Total[treels[[foc+1;; maxpos - 1, -1]]] +Log[Times @@ treels[[foc;; maxpos - 2, 3]]]);
                ar+= logpriSNP[snpls[[All,1]],proptreels, foc,maxpos, theta, epsilon, nbp,modelid]-
                		logpriSNP[snpls[[All,1]],treels, foc,maxpos-1, theta, epsilon, nbp,modelid];
                ar = heat ar +logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
                  logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
                If[ !(NumericQ[ar]&&NumericQ[ar2]),
                    Print["5: {ar,ar2}=",{ar,ar2}]
                ];
                ar = Min[1, Exp[ar - ar2]];
            ];
        ];
        If[ RandomReal[] < ar,
            {propsnpls,proptreels, {seglen, ar}},
            {snpls, treels, {seglen, ar}}
        ]
    ]    

updatesegtreedelete[inputsnpls_, inputtreels_, theta_, epsilon_, nbp_,prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_, heat_,maxseglen_] :=
    Module[ {snpls = inputsnpls, treels = inputtreels, seglen = 1,ar,type,
      foc, proptreels, ar2, propsub,  maxpos, propsnpls, seg0,ar2sub},
        If[ Length[treels]==1,
            ar = 0,
            foc = RandomInteger[{2, Length[treels]}];
            proptreels = Delete[treels, foc];
            seg0 = {treels[[foc - 1, 1]]+1,If[ foc == Length[treels],
                                               nbp,
                                               treels[[foc + 1, 1]] - 1
                                           ]};
            If[ foc==Length[treels],
                {propsnpls,ar} = segloglike[snpls, treels[[foc-1,2]], treels[[foc,1]],nbp+1,  theta, epsilon,modelid];
                ar+=0-(treels[[foc, -1]] +Log[treels[[foc - 1, 3]]]);
                ar+= logpriSNP[snpls[[All,1]],proptreels, foc-1,foc-1, theta, epsilon, nbp,modelid]-
                		logpriSNP[snpls[[All,1]],treels, foc-1,foc, theta, epsilon, nbp,modelid];
                ar = heat ar+logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat]-
                    logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax,heat];
                type = TreeTransitionTypeI@@treels[[{foc-1,foc},2]];
                ar2 = 0-(logpdfjumptree[treels[[foc-1,2]],type]+Log[1/(seg0[[2]]-seg0[[1]]+1)]);
                If[ !(NumericQ[ar]&&NumericQ[ar2]),
                    Print["6: {ar,ar2}=",{ar,ar2}]
                ];
                ar = Min[1, Exp[ar - ar2]],
                {ar2sub,propsub} = randomjumpsegcondtree[treels[[foc ;;, 2 ;;]],proptreels[[foc - 1, 2 ;;]],maxseglen];
                If[ ar2sub ===None,
                    seglen = -1;
                    ar = 0,
                    seglen = Length[propsub];
                    type = TreeTransitionTypeI@@treels[[{foc-1,foc},2]];
                    ar2 = 0-(logpdfjumptree[treels[[foc-1,2]],type]+Log[1/(seg0[[2]]-seg0[[1]]+1)]);
                    ar2+=ar2sub;
                    maxpos = foc + seglen - 2;
                    proptreels[[foc - 1 ;; maxpos, 2 ;;]] = propsub;
                    {propsnpls, ar} = segloglike2[snpls, proptreels, foc-1, maxpos, theta, epsilon,nbp,modelid];
                    ar+=(Total[proptreels[[foc;; maxpos, -1]]] +Log[Times @@ proptreels[[foc-1;; maxpos-1,3]]]) -
                         (Total[treels[[foc;; maxpos+1, -1]]]+Log[Times @@ treels[[foc-1;; maxpos, 3]]]);
                    ar+= logpriSNP[snpls[[All,1]],proptreels, foc-1,maxpos, theta, epsilon, nbp,modelid]-
                    		logpriSNP[snpls[[All,1]],treels, foc-1,maxpos+1, theta, epsilon, nbp,modelid];
                    ar = heat ar + logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat] - 
                      logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
                    If[ !(NumericQ[ar]&&NumericQ[ar2]),
                        Print["7: {ar,ar2}=",{ar,ar2}]
                    ];
                    ar = Min[1, Exp[ar - ar2]];
                ];
            ];
        ];
        If[ RandomReal[] < ar,
            {propsnpls,proptreels, {seglen, ar}},
            {snpls, treels, {seglen, ar}}
        ]
    ]   
    
updatesegtreepos[snpls_, treels_, segl_, theta_, epsilon_, nbp_, prirhoshape_, prirhoscale_,prirhomin_,prirhomax_,modelid_,heat_] :=
    Module[ {i1, i2, xleft, xright, proptreels, propsnpls, ar},
        If[ Length[treels]==1,
            {snpls, treels, 0.234},
            (*i1 = RandomInteger[{2, Max[2,Length[treels]-segl+1]}];*)
            i1 = RandomInteger[{2, Length[treels]}];
            i2 = Min[i1 + segl - 1, Length[treels]];
            xleft = treels[[i1 - 1, 1]];
            xright = If[ i2 == Length[treels],
                         nbp + 1,
                         treels[[i2 + 1, 1]]
                     ];
            proptreels = treels;
            proptreels[[i1 ;; i2, 1]] = Sort[RandomSample[Range[xleft + 1, xright - 1], i2 - i1 + 1]];
            {propsnpls, ar} = segloglike2[snpls, proptreels, i1-1,i2, theta, epsilon,nbp,modelid];
            ar+= logpriSNP[snpls[[All,1]],proptreels, i1-1,i2, theta, epsilon, nbp,modelid]-
            		logpriSNP[snpls[[All,1]],treels, i1-1,i2, theta, epsilon, nbp,modelid];
            ar = heat ar + logpritree[proptreels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat]-
                logpritree[treels, nbp, theta,prirhoshape, prirhoscale,prirhomin,prirhomax, heat];
            If[ !(NumericQ[ar]),
                Print["1: {ar}=",{ar}]
            ];
            ar = Min[1, Exp[ar]];
            If[ RandomReal[] < ar,
                {propsnpls, proptreels, ar},
                {snpls, treels, ar}
            ]
        ]
    ]    
           
(****************************************************************************)    


End[]

EndPackage[]

