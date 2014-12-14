Steps for infering local trees from SNP data:
1) There are two folders "InferSMARTree" and "Mathematica Packages", and one text file "README_SMARTree.txt".

2) Refer to "SMARTree_MathematicaCode\Mathematica Packages\HTML Documents\SMARTree-HTML\html\guide\SMARTree.html" for the package "SMARTree". Similarly "Genealogy.html" for the package "Genealogy", and for other packages.

3) Run the mathematica notebook "SMARTree_MathematicaCode\InferSMARTree\InferSMARTree.nb" by double clicking it.

4) In the paper, 
	model M1 = primodel "P1"+ ascertainmodel "M1"
	model M2 = primodel "P1"+ ascertainmodel "M2"
	model M1*=  primodel "P2"+ ascertainmodel "M1"
   The default is model M1 for the data set: Data_Standard.txt

4.1) In "InferSMARTree.nb", the default setup for running conditions:

	isrestart = False;
	isparallel = True;
	incorpdata = True;
	primodel="P1";
	ascertainmodel = "M1";
	ngroupofchain = 2;
	nchainofgroup = 16;
	heatingAR = 0.5;
	datafile = "Standard_Data.txt";
	outcompact =  "SMARTree_OutCompact_" <> ascertainmodel <> primodel <> "_" <> datafile
	outcold = "SMARTree_OutCold_" <> ascertainmodel <> primodel <> "_" <> datafile
	outset = "SMARTree_OutSetting_" <> ascertainmodel <> primodel <> "_" <> datafile
	outstate = "SMARTree_MCMCState_" <> ascertainmodel <> primodel <> "_" <> datafile

	rjfreq = Ceiling[nbp/5000];
	adpmax = 2000;
	adprate = 3/adpmax;
	itmin = 1;
	itmax = 1000000
	printfreq = 100;
	savefreq = 2;
	swapfreq = nchainofgroup;
	isprint = False;

4.2) 
        *) If isrestart=True, continue running by using the saved mcmc states from the previous running.
	*) If isparallel=True, parallel compuation for multiple MCMC chains
	*) If iscorpdata=False, the SNP data is not accounted for. It may be used for model check whether the posteriors are same as the priors
	*) There are only two choice for prior model: "P1" and "P2"
	*) There are only two choices for ascertainmodel: "M1" and "M2".
	*) There are in total ngroupofchain*nchainofgroup mcmc chains. They are divided into "ngroupofchain" independent groups. 
	*) Within each group, parallel temperating algorithm is used, and the target accept ratio is heatingAR.
	*) The running conditions are saved in the outset file. And the MCMC state is saved in the outstate file, which may be used to restart/continue the mcmc chains. 
	*) The results are saved in the two files. Only cold chains (temperature=1) are saved in the e.g. "SMARTree_OutCold_M1_Standard_Data.txt". All the chains are saved un the outcompact where estimated local trees are excluded.


4.3) 
	*) rjfreq is the number of reversible jump sampling per iteration
	*) adpmax is the number of iniital iterations, during which the proposal distributions are adpatived.
	*) adprate is the rate of adaptation if isadp=True.
	*) itmax-itmin: is the max number of iteations.
	*) printfreq is the interval number of iterations at which some parameter values are printed on monitor.
	*) savefreq is the interval number of iteration at which results are saved.
	*) swapfreq is the number of swapping of mcmc chains with each group
	*) isprint=True, messages for the detailed steps are printed on monitor.

5) The resoults saved in outcold are analyzed by ShowResults_SMARTree.nb
